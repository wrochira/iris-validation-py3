import os
from math import sin, cos, pi

import svgwrite
from svgwrite.animate import Animate

from iris_validation._defs import COLORS, CHAIN_VIEW_RINGS, CHAIN_VIEW_GAP_ANGLE


class ChainView():
    def __init__(self, data, chain_index, canvas_size=(1000, 1000), hidden=False):
        self.data = data
        self.chain_index = chain_index
        self.canvas_size = canvas_size
        self.hidden = hidden

        self.dwg = None
        self.cfa_cache = { }
        self.num_rings = len(CHAIN_VIEW_RINGS)
        self.num_versions = self.data['num_versions']
        self.num_segments = self.data['aligned_length']
        self.center = (self.canvas_size[0] // 2, self.canvas_size[1] // 2)
        self.full_radius = round(min(self.canvas_size) / 2 - 10, 2)
        self.division_size = round(self.full_radius / (self.num_rings + 2), 2)
        self.angle_delta = (2 * pi - CHAIN_VIEW_GAP_ANGLE) / self.num_segments
        self.svg_id = f'iris-chain-view-{self.chain_index}'

        self._draw()

    def _coords_from_angle(self, angle, radius, gap=True):
        gap_angle = CHAIN_VIEW_GAP_ANGLE if gap else 0.0
        arg_string = str([ self.center, angle, radius, gap_angle ])
        if arg_string in self.cfa_cache:
            coords = self.cfa_cache[arg_string]
        else:
            result_x = self.center[0] + radius * sin(angle + gap_angle/2)
            result_y = self.center[1] - radius * cos(angle + gap_angle/2)
            coords = (round(result_x, 1), round(result_y, 1))
            self.cfa_cache[arg_string] = coords
        return coords

    def _draw(self):
        # Initialise drawing
        self.dwg = svgwrite.Drawing(profile='full')

        # Set HTML attributes
        self.dwg.attribs['viewBox'] = '0 0 ' + ' '.join([ str(x) for x in self.canvas_size ])
        self.dwg.attribs['id'] = self.svg_id
        if self.hidden:
            self.dwg.attribs['style'] = 'display: none;'

        # Draw background
        self.dwg.add(self.dwg.circle(r=self.full_radius,
                                     center=self.center,
                                     fill=COLORS['WHITE'],
                                     fill_opacity=1,
                                     stroke_opacity=0))


        # Draw data rings
        for ring_id, ring_metric in enumerate(CHAIN_VIEW_RINGS):
            self._add_ring(ring_id, ring_metric)

        # Draw missing-data shade
        for version_id, residue_validities in enumerate(self.data['residue_validities']):
            group_opacity = 1 if version_id == self.num_versions-1 else 0
            shade_group = self.dwg.g(id=f'{self.svg_id}-shade-{version_id}', opacity=group_opacity)
            for segment_id, residue_valid in enumerate(residue_validities):
                if not residue_valid:
                    shade_group.add(self.dwg.polygon([ self.center,
                                                       self._coords_from_angle(self.angle_delta * segment_id, self.full_radius+5),
                                                       self._coords_from_angle(self.angle_delta * (segment_id+1), self.full_radius+5) ],
                                                       stroke_opacity=0,
                                                       fill=COLORS['L_PINK'],
                                                       fill_opacity=1))
            self.dwg.add(shade_group)

        # Draw outer rings
        self.dwg.add(self.dwg.circle(r=self.full_radius-24,
                                     center=self.center,
                                     fill_opacity=0,
                                     stroke=COLORS['BLACK'],
                                     stroke_width=1,
                                     stroke_opacity=0.5))
        self.dwg.add(self.dwg.circle(r=self.full_radius-8,
                                     center=self.center,
                                     fill_opacity=0,
                                     stroke=COLORS['BLACK'],
                                     stroke_width=1,
                                     stroke_opacity=0.5))
        for i in range(self.num_segments+1):
            self.dwg.add(self.dwg.line(self._coords_from_angle(self.angle_delta*i, self.full_radius-24),
                         self._coords_from_angle(self.angle_delta*i, self.full_radius-8),
                         stroke=COLORS['BLACK'],
                         stroke_width=1,
                         stroke_opacity=0.5))

        # Draw segment selector
        center_point = self.angle_delta*0.5
        selector_points = (self._coords_from_angle(center_point, self.full_radius-16),
                           self._coords_from_angle(center_point-0.02, self.full_radius-8),
                           self._coords_from_angle(center_point-0.02, self.full_radius+8),
                           self._coords_from_angle(center_point+0.02, self.full_radius+8),
                           self._coords_from_angle(center_point+0.02, self.full_radius-8))
        self.dwg.add(self.dwg.polygon(selector_points,
                                      stroke=COLORS['BLACK'],
                                      stroke_width=2,
                                      stroke_opacity=1,
                                      fill=COLORS['GREY'],
                                      fill_opacity=0.2,
                                      id=f'{self.svg_id}-residue-selector'))

        # Draw interaction segments
        for segment_id in range(self.num_segments):
            self.dwg.add(self.dwg.polygon([ self.center,
                                            self._coords_from_angle(self.angle_delta * segment_id, self.full_radius+5),
                                            self._coords_from_angle(self.angle_delta * (segment_id+1), self.full_radius+5) ],
                                          stroke=COLORS['BLACK'],
                                          stroke_width=1,
                                          stroke_opacity=0,
                                          fill=COLORS['L_GREY'],
                                          fill_opacity=0,
                                          onmousedown=f'handleSegment(1, {segment_id});',
                                          onmouseover=f'handleSegment(2, {segment_id});',
                                          onmouseup=f'handleSegment(3, {segment_id});',
                                          id=f'{self.svg_id}-interaction-segment-{segment_id}'))
        self.dwg.add(self.dwg.circle(r=1.5*self.division_size,
                                     center=self.center,
                                     fill=COLORS['WHITE'],
                                     fill_opacity=1,
                                     stroke_opacity=0))

        # Draw center text
        self.dwg.add(self.dwg.text(text='Iris',
                                   insert=(self.center[0], self.center[1]-24),
                                   font_size=1.5*16,
                                   font_family='Arial',
                                   font_weight='bold',
                                   text_anchor='middle',
                                   alignment_baseline='central'))
        self.dwg.add(self.dwg.text(text='Chain ' + self.data['chain_id'],
                                   insert=(self.center[0], self.center[1]+16),
                                   font_size=16,
                                   font_family='Arial',
                                   text_anchor='middle',
                                   alignment_baseline='central'))
        if self.data['has_molprobity']:
            self.dwg.add(self.dwg.text(text='MolProbity',
                                       insert=(self.center[0], self.center[1]+48),
                                       font_size=16,
                                       font_family='Arial',
                                       text_anchor='middle',
                                       alignment_baseline='central',
                                       fill=COLORS['L_GREY']))

    def _add_ring(self, ring_id, metric):
        datapoints = self.data[metric['type'] + '_values'][metric['id']]

        # Draw axes
        ring_base_radius = (ring_id + 2) * self.division_size
        self.dwg.add(self.dwg.circle(r=ring_base_radius,
                                     center=self.center,
                                     fill_opacity=0,
                                     stroke=metric['ring_color'],
                                     stroke_width=1,
                                     stroke_opacity=1))
        self.dwg.add(self.dwg.polyline([ self._coords_from_angle((CHAIN_VIEW_GAP_ANGLE/25)*(i-(20-1)/2), ring_base_radius, gap=False) for i in range(20) ],
                                       stroke=metric['ring_color'],
                                       stroke_width=3,
                                       stroke_opacity=1,
                                       fill_opacity=0))
        self.dwg.add(self.dwg.text(text=metric['short_name'],
                                   insert=self._coords_from_angle(0, ring_base_radius+12, gap=False),
                                   font_size=16,
                                   font_family='Arial',
                                   text_anchor='middle',
                                   alignment_baseline='central'))

        if metric['type'] == 'discrete':
            for version_id, version_datapoints in enumerate(datapoints):
                version_ring_segments = [ ]
                for segment_id, datapoint in enumerate(version_datapoints):
                    segment_length = 10
                    segment_color = metric['seq_colors'][-1]
                    segment_opacity = 1
                    if datapoint is not None and 0 <= datapoint < len(metric['seq_colors']):
                        segment_color = metric['seq_colors'][datapoint]
                    if segment_color == metric['seq_colors'][-1]:
                        segment_opacity = 0.5
                    segment_points = (self._coords_from_angle(self.angle_delta * (segment_id),
                                                              ring_base_radius - segment_length),
                                      self._coords_from_angle(self.angle_delta * (segment_id),
                                                              ring_base_radius + segment_length),
                                      self._coords_from_angle(self.angle_delta * (segment_id+1),
                                                              ring_base_radius + segment_length),
                                      self._coords_from_angle(self.angle_delta * (segment_id+1),
                                                              ring_base_radius - segment_length))
                    version_ring_segments.append((segment_points, segment_color, segment_opacity))
                group_opacity = 1 if version_id == self.num_versions-1 else 0
                segment_group = self.dwg.g(id=f'{self.svg_id}-discrete-{version_id}-{ring_id}', opacity=group_opacity)
                for segment_points, segment_color, segment_opacity in version_ring_segments:
                    segment_group.add(self.dwg.polyline(segment_points,
                                                        stroke_width=0,
                                                        stroke_opacity=0,
                                                        fill=segment_color,
                                                        fill_opacity=segment_opacity))
                self.dwg.add(segment_group)

        elif metric['type'] == 'continuous':
            # Get mean metric value
            all_valid_values = [ ]
            for version_datapoints in datapoints:
                for datapoint in version_datapoints:
                    if datapoint is None:
                        continue
                    value = datapoint * metric['polarity']
                    all_valid_values.append(value)
            ring_avg = 0
            if len(all_valid_values) == 0:
                return
            ring_avg = sum(all_valid_values) / len(all_valid_values)

            # Calculate deltas from the ring average
            deltas = [ ]
            for version_datapoints in datapoints:
                version_deltas = [ ]
                for datapoint in version_datapoints:
                    delta = None
                    if datapoint is not None:
                        value = datapoint * metric['polarity']
                        delta = value - ring_avg
                    version_deltas.append(delta)
                deltas.append(version_deltas)

            # Calculate average negative delta in the latest dataset
            latest_negative_deltas = [ x for x in deltas[-1] if x is not None and x < 0 ]
            avg_negative_delta = 0
            if len(latest_negative_deltas) > 0:
                avg_negative_delta = sum(latest_negative_deltas) / len(latest_negative_deltas)

            # Subtract the average negative delta from all deltas to calculate 'magnitudes'
            magnitudes = [ ]
            all_valid_magnitudes = [ ]
            for version_deltas in deltas:
                version_magnitudes = [ x - avg_negative_delta if x is not None else None for x in version_deltas ]
                all_valid_magnitudes += [ x for x in version_magnitudes if x is not None ]
                magnitudes.append(version_magnitudes)
            magnitude_min, magnitude_max = (min(all_valid_magnitudes), max(all_valid_magnitudes))

            # Calculate plot magnitudes
            plot_magnitudes = [ ]
            for version_magnitudes in magnitudes:
                version_plot_magnitudes = [ ]
                for magnitude in version_magnitudes:
                    if magnitude is None:
                        version_plot_magnitudes.append(None)
                        continue
                    plot_magnitude = 0
                    if magnitude > 0 and magnitude_max != 0:
                        plot_magnitude = magnitude / magnitude_max * +0.25
                    elif magnitude < 0 and magnitude_min != 0:
                        plot_magnitude = magnitude / magnitude_min * -0.7
                    version_plot_magnitudes.append(plot_magnitude)
                plot_magnitudes.append(version_plot_magnitudes)

            # Calculate plot point coordinates
            line_points = [ ]
            for version_plot_magnitudes in plot_magnitudes:
                version_line_points = [ ]
                zero_point = self._coords_from_angle(self.angle_delta*0.5, ring_base_radius)
                version_line_points.append(zero_point)
                for segment_id, plot_magnitude in enumerate(version_plot_magnitudes):
                    angle = self.angle_delta * (segment_id + 0.5)
                    plot_radius = ring_base_radius
                    if plot_magnitude is not None:
                        plot_radius += self.division_size * plot_magnitude
                    point = self._coords_from_angle(angle, plot_radius)
                    version_line_points.append(point)
                line_points.append(version_line_points)

            # Draw line
            baseline_circle_points = [ ]
            baseline_point_resolution = 200
            for point_id in range(baseline_point_resolution + 1):
                point_angle = (baseline_point_resolution - point_id) * (2*pi - CHAIN_VIEW_GAP_ANGLE) / baseline_point_resolution
                baseline_circle_points.append(self._coords_from_angle(point_angle, ring_base_radius))
            plot_points = line_points[-1] + baseline_circle_points
            ring_line = self.dwg.polyline(plot_points,
                                          stroke=metric['ring_color'],
                                          stroke_width=2,
                                          stroke_opacity=1,
                                          fill=metric['ring_color'],
                                          fill_opacity=0.2)
            for version_id, version_line_points in enumerate(line_points):
                plot_points = version_line_points + baseline_circle_points
                points_string = ' '.join([ ','.join([ str(x) for x in point ]) for point in plot_points ])
                animation = Animate(values=None,
                                    dur='250ms',
                                    begin='indefinite',
                                    fill='freeze',
                                    attributeName='points',
                                    to=points_string,
                                    id=f'{self.svg_id}-animation-{version_id}-{ring_id}')
                ring_line.add(animation)
            self.dwg.add(ring_line)
