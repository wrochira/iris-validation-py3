import os
import json
import math

import svgwrite
from svgwrite.animate import Animate

from iris_validation.graphics.chain import ChainView
from iris_validation.graphics.residue import ResidueView
from iris_validation._defs import COLORS, CHAIN_VIEW_RINGS, RESIDUE_VIEW_BOXES, RESIDUE_VIEW_BARS, CHAIN_VIEW_GAP_ANGLE


JS_PATH = os.path.join(os.path.dirname(__file__), 'js')
JS_CONSTANTS_PATH = os.path.join(JS_PATH, 'constants.js')
JS_INTERACTION_PATH = os.path.join(JS_PATH, 'interaction.js')


class Panel():
    def __init__(self, data, canvas_size=(1500, 1000)):
        self.data = data
        self.canvas_size = canvas_size

        self.dwg = None
        self.javascript = None
        self.chain_views = None
        self.residue_view = None
        self.num_models = self.data[0]['num_versions']
        self.chain_ids = [ chain_data['chain_id'] for chain_data in self.data ]
        self.swtich_colors = [ COLORS['VL_GREY'], COLORS['CYAN'] ]
        self.svg_id = 'iris-panel'

        self._verify_chosen_metrics()
        self._generate_javascript()
        self._generate_subviews()
        self._draw()

    # TODO: Make this nicer
    def _verify_chosen_metrics(self):
        global CHAIN_VIEW_RINGS, RESIDUE_VIEW_BOXES, RESIDUE_VIEW_BARS
        for metric_list in (CHAIN_VIEW_RINGS, RESIDUE_VIEW_BOXES, RESIDUE_VIEW_BARS):
            if not isinstance(metric_list, list):
                raise ValueError('Chosen metrics in the _defs.py file must be lists')
            for metric_index in reversed(range(len(metric_list))):
                if (metric_list[metric_index]['is_covariance'] and not self.data[0]['has_covariance']):
                    del metric_list[metric_index]
                elif (metric_list[metric_index]['is_molprobity'] and not self.data[0]['has_molprobity']):
                    del metric_list[metric_index]
                elif (metric_list[metric_index]['is_reflections'] and not self.data[0]['has_reflections']):
                    del metric_list[metric_index]

    def _generate_javascript(self):
        json_data = json.dumps(self.data)
        num_chains = len(self.chain_ids)
        bar_metric_ids = [ metric['id'] for metric in RESIDUE_VIEW_BARS ]
        box_metric_ids = [ metric['id'] for metric in RESIDUE_VIEW_BOXES ]
        box_colors = json.dumps([ metric['seq_colors'] for metric in RESIDUE_VIEW_BOXES ])
        box_labels = json.dumps([ metric['seq_labels'] for metric in RESIDUE_VIEW_BOXES ])
        gap_degrees = CHAIN_VIEW_GAP_ANGLE * 180 / math.pi

        with open(JS_CONSTANTS_PATH, 'r', encoding='utf8') as infile:
            js_constants = infile.read()

        with open(JS_INTERACTION_PATH, 'r', encoding='utf8') as infile:
            js_interation = infile.read()

        js_constants = js_constants.format(model_data=json_data,
                                           num_chains=num_chains,
                                           bar_metric_ids=bar_metric_ids,
                                           box_metric_ids=box_metric_ids,
                                           box_colors=box_colors,
                                           box_labels=box_labels,
                                           gap_degrees=gap_degrees,
                                           chain_selector_colors=self.swtich_colors)

        self.javascript = js_constants + js_interation

    def _generate_subviews(self):
        self.chain_views = [ ]
        for chain_index, chain_data in enumerate(self.data):
            chain_view = ChainView(chain_data, chain_index, hidden=chain_index>0).dwg
            self.chain_views.append(chain_view)
        self.residue_view = ResidueView().dwg

    def _draw(self):
        middle_gap = 30
        view_border = 10
        view_title_font = 24
        button_width = 38
        button_height = 32
        view_width, view_height = [ dim - view_border for dim in self.canvas_size ]
        view_divider_x = round(2/3 * view_width, 2)
        chain_view_bounds = (view_border,
                             view_border,
                             view_divider_x - round(middle_gap/2, 2),
                             view_height)
        residue_view_bounds = (view_divider_x + round(middle_gap/2, 2),
                               view_border,
                               view_width,
                               view_height)

        # Initialise drawing
        self.dwg = svgwrite.Drawing(profile='full')

        # Disable text selection
        self.dwg.attribs['style'] = 'user-select: none;'

        # Draw background
        self.dwg.add(self.dwg.polygon(points=[ (0, 0),
                                               (0, self.canvas_size[1]),
                                               (self.canvas_size[0], self.canvas_size[1]),
                                               (self.canvas_size[0], 0) ],
                                      fill=COLORS['WHITE'],
                                      fill_opacity=1,
                                      stroke_opacity=0))

        # Set HTML attributes
        self.dwg.attribs['viewBox'] = '0 0 ' + ' '.join([ str(x) for x in self.canvas_size ])
        self.dwg.attribs['id'] = self.svg_id

        # Add JavaScript
        self.dwg.defs.add(self.dwg.script(content=self.javascript))

        # View titles and divider lines
        self.dwg.add(self.dwg.text(text='Chain',
                                   insert=(chain_view_bounds[0], chain_view_bounds[1]+view_title_font),
                                   font_size=view_title_font,
                                   font_family='Arial'))

        self.dwg.add(self.dwg.text(text='Residue',
                                   insert=(residue_view_bounds[0], residue_view_bounds[1]+view_title_font),
                                   font_size=view_title_font,
                                   font_family='Arial',
                                   id=f'{self.svg_id}-residue-summary'))

        self.dwg.add(self.dwg.line((chain_view_bounds[0], chain_view_bounds[1]+40),
                                   (chain_view_bounds[2], chain_view_bounds[1]+40),
                                   stroke=COLORS['BLACK'],
                                   stroke_width=2))

        self.dwg.add(self.dwg.line((residue_view_bounds[0], residue_view_bounds[1]+40),
                                   (residue_view_bounds[2], residue_view_bounds[1]+40),
                                   stroke=COLORS['BLACK'],
                                   stroke_width=2))

        # Chain selector buttons
        for chain_index, chain_id in enumerate(self.chain_ids[:12]):
            selector_color = self.swtich_colors[1] if chain_index == 0 else self.swtich_colors[0]
            self.dwg.add(self.dwg.rect(insert=(chain_view_bounds[0] + 75 + 50*chain_index, chain_view_bounds[1]),
                                       size=(button_width, button_height),
                                       rx=5,
                                       stroke_opacity=0,
                                       fill_opacity=0.5,
                                       fill=selector_color,
                                       id=f'{self.svg_id}-chain-selector-{chain_index}'))

            self.dwg.add(self.dwg.text(text=chain_id,
                                       insert=(chain_view_bounds[0] + 75 + button_width/2 + 50*chain_index, chain_view_bounds[1] + button_height/2),
                                       font_size=view_title_font,
                                       font_family='Arial',
                                       text_anchor='middle',
                                       alignment_baseline='central'))

            self.dwg.add(self.dwg.rect(insert=(chain_view_bounds[0] + 75 + 50*chain_index, chain_view_bounds[1]),
                                       size=(button_width, button_height),
                                       rx=5,
                                       stroke_opacity=0,
                                       fill_opacity=0,
                                       onmouseover='setPointer();',
                                       onmouseout='unsetPointer();',
                                       onclick=f'setChain({chain_index});'))

        # Extra chains dropdown
        # TODO: finish this
        if len(self.chain_ids) > 12:
            chain_index = 12
            selector_color = self.swtich_colors[0]
            self.dwg.add(self.dwg.rect(insert=(chain_view_bounds[0] + 75 + 50*chain_index, chain_view_bounds[1]),
                                       size=(38, 32),
                                       rx=5,
                                       stroke_opacity=0,
                                       fill_opacity=0.5,
                                       fill=selector_color,
                                       id=f'{self.svg_id}-chain-selector-dropdown'))

            self.dwg.add(self.dwg.text(text='...',
                                       insert=(chain_view_bounds[0] + 85 + 50*chain_index, chain_view_bounds[1]+view_title_font),
                                       font_size=view_title_font,
                                       font_family='Arial'))

            self.dwg.add(self.dwg.rect(insert=(chain_view_bounds[0] + 75 + 50*chain_index, chain_view_bounds[1]),
                                       size=(38, 32),
                                       rx=5,
                                       stroke_opacity=0,
                                       fill_opacity=0,
                                       onmouseover='setPointer();',
                                       onmouseout='unsetPointer();',
                                       onclick=f'toggleDropdown();'))

        # Version toggle switch
        self.dwg.add(self.dwg.text(text='Previous',
                                   insert=(chain_view_bounds[2]-215, chain_view_bounds[1]+20),
                                   font_size=16,
                                   font_family='Arial'))

        self.dwg.add(self.dwg.text(text='Latest',
                                   insert=(chain_view_bounds[2]-55, chain_view_bounds[1]+20),
                                   font_size=16,
                                   font_family='Arial'))

        if self.num_models > 1:
            switch_group = self.dwg.g(id=f'{self.svg_id}-switch',
                                      onmouseover='setPointer();',
                                      onmouseout='unsetPointer();',
                                      onclick='toggleVersion();')
        else:
            switch_group = self.dwg.g(id=f'{self.svg_id}-switch')

        switch_rectangle = self.dwg.rect(insert=(chain_view_bounds[2]-140, chain_view_bounds[1]),
                                         size=(70, 30),
                                         rx=15,
                                         stroke_opacity=0,
                                         fill_opacity=1,
                                         fill=self.swtich_colors[1])

        for version_id in range(2):
            animation = Animate(values=None,
                                dur='250ms',
                                begin='indefinite',
                                fill='freeze',
                                attributeName='fill',
                                to=self.swtich_colors[version_id],
                                id=f'{self.svg_id}-switch-color-animation-{version_id}')
            switch_rectangle.add(animation)

        switch_group.add(switch_rectangle)

        switch_circle = self.dwg.circle(r=10,
                                        center=(chain_view_bounds[2]-85, chain_view_bounds[1]+15),
                                        stroke_opacity=0,
                                        fill_opacity=1,
                                        fill=COLORS['WHITE'])

        for version_id in range(2):
            animation = Animate(values=None,
                                dur='250ms',
                                begin='indefinite',
                                fill='freeze',
                                attributeName='cx',
                                to=(chain_view_bounds[2]-125, chain_view_bounds[2]-85)[version_id],
                                id=f'{self.svg_id}-switch-move-animation-{version_id}')
            switch_circle.add(animation)

        switch_group.add(switch_circle)
        self.dwg.add(switch_group)

        # Place sub-views
        canvas_mid_x = self.canvas_size[0] / 2
        # *** Chain view
        view_mid_x = (chain_view_bounds[0] + chain_view_bounds[2]) / 2
        view_adj_x = int(round(view_mid_x - canvas_mid_x))
        view_space_width = chain_view_bounds[2] - chain_view_bounds[0]
        view_space_height = chain_view_bounds[3] - (chain_view_bounds[1] + 50)
        viewbox_width = int(round(1000**2 / view_space_width))
        viewbox_height = int(round(1000**2 / view_space_height))
        width_buffer = -((viewbox_width - 1000) / 2)
        height_buffer = -((viewbox_height - 1000 - 50) / 2 + 50)
        for chain_view in self.chain_views:
            chain_view.attribs['x'] = str(view_adj_x)
            chain_view.attribs['viewBox'] = f'{width_buffer} {height_buffer} {viewbox_width} {viewbox_height}'
            self.dwg.add(chain_view)
        # *** Residue view
        view_mid_x = (residue_view_bounds[0] + residue_view_bounds[2]) / 2
        view_adj_x = int(round(view_mid_x - canvas_mid_x))
        view_space_width = residue_view_bounds[2] - residue_view_bounds[0]
        view_space_height = residue_view_bounds[3] - (residue_view_bounds[1] + 50)
        viewbox_width = int(round(400**2 / view_space_width))
        viewbox_height = int(round(1000**2 / view_space_height))
        width_buffer = -((viewbox_width - 400) / 2)
        height_buffer = -((viewbox_height - 1000 - 50) / 2 + 50)
        self.residue_view.attribs['x'] = str(view_adj_x)
        self.residue_view.attribs['viewBox'] = f'{width_buffer} {height_buffer} {viewbox_width} {viewbox_height}'
        self.dwg.add(self.residue_view)
