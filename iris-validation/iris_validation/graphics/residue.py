import svgwrite
from svgwrite.gradients import LinearGradient

from iris_validation._defs import COLORS, RESIDUE_VIEW_BOXES, RESIDUE_VIEW_BARS


class ResidueView():
    def __init__(self, canvas_size=(400, 1000)):
        self.canvas_size = canvas_size

        self.dwg = None
        self.svg_id = 'iris-residue-view'

        self.box_names = [ metric['short_name'] for metric in RESIDUE_VIEW_BOXES ]
        self.bar_names = [ metric['long_name'] for metric in RESIDUE_VIEW_BARS ]

        # TODO: allow any number of bars
        self.bar_names = self.bar_names[:2]

        self._draw()

    def _draw(self):
        top_margin = 10
        left_indent = 35

        # Initialise drawing
        self.dwg = svgwrite.Drawing(profile='full')

        # Set HTML attributes
        self.dwg.attribs['viewBox'] = '0 0 ' + ' '.join([ str(x) for x in self.canvas_size ])
        self.dwg.attribs['id'] = self.svg_id

        # Draw background
        self.dwg.add(self.dwg.polygon(points=[ (0, 0),
                                               (0, self.canvas_size[1]),
                                               (self.canvas_size[0], self.canvas_size[1]),
                                               (self.canvas_size[0], 0) ],
                                      fill=COLORS['WHITE'],
                                      fill_opacity=1,
                                      stroke_opacity=0))

        # Boxes
        for box_id, box_title in enumerate(self.box_names):
            box_bounds = [ 0.25*self.canvas_size[0]+left_indent,
                               top_margin + 70*box_id,
                               self.canvas_size[0],
                               top_margin + 70*box_id + 50 ]

            self.dwg.add(self.dwg.polygon(points=[ (left_indent, box_bounds[1]),
                                                   (box_bounds[0], box_bounds[1]),
                                                   (box_bounds[0], box_bounds[3]),
                                                   (left_indent, box_bounds[3]) ],
                                          fill=COLORS['VL_GREY'],
                                          fill_opacity=0.5,
                                          stroke=COLORS['BLACK'],
                                          stroke_width=2,
                                          stroke_opacity=1))
            self.dwg.add(self.dwg.polygon(points=[ (box_bounds[0], box_bounds[1]),
                                                   (box_bounds[2], box_bounds[1]),
                                                   (box_bounds[2], box_bounds[3]),
                                                   (box_bounds[0], box_bounds[3]) ],
                                          fill=COLORS['VL_GREY'],
                                          fill_opacity=0.75,
                                          stroke=COLORS['BLACK'],
                                          stroke_width=2,
                                          stroke_opacity=1,
                                          id=f'{self.svg_id}-box-{box_id}'))
            self.dwg.add(self.dwg.text('',
                                       insert=((box_bounds[0]+box_bounds[2])/2, (box_bounds[1]+box_bounds[3])/2),
                                       font_size=20,
                                       font_family='Arial',
                                       font_weight='bold',
                                       fill=COLORS['BLACK'],
                                       fill_opacity=1,
                                       text_anchor='middle',
                                       alignment_baseline='central',
                                       id=f'{self.svg_id}-box-{box_id}-text'))
            self.dwg.add(self.dwg.text(box_title,
                                       insert=(left_indent + 0.125*self.canvas_size[0], (box_bounds[1]+box_bounds[3])/2),
                                       font_size=18,
                                       font_family='Arial',
                                       fill=COLORS['BLACK'],
                                       fill_opacity=1,
                                       text_anchor='middle',
                                       alignment_baseline='central'))

        # Bars
        bar_width = 120
        bar_charts_bounds = (left_indent,
                             70*len(self.box_names)+30,
                             self.canvas_size[0],
                             self.canvas_size[1]-60)

        # Bar chart container
        self.dwg.add(self.dwg.polygon(points=[ (bar_charts_bounds[0], bar_charts_bounds[1]),
                                               (bar_charts_bounds[2], bar_charts_bounds[1]),
                                               (bar_charts_bounds[2], bar_charts_bounds[3]),
                                               (bar_charts_bounds[0], bar_charts_bounds[3]) ],
                                      fill=COLORS['WHITE'],
                                      fill_opacity=0,
                                      stroke=COLORS['BLACK'],
                                      stroke_width=2,
                                      stroke_opacity=1,
                                      id=f'{self.svg_id}-bar-charts-container'))

        # Bar chart axis
        for label_id in range(10+1):
            height = bar_charts_bounds[1] + label_id*(bar_charts_bounds[3]-bar_charts_bounds[1])/10
            self.dwg.add(self.dwg.line((bar_charts_bounds[0]-5, height), (bar_charts_bounds[0]+5, height),
                                       stroke=COLORS['BLACK'],
                                       stroke_width=2,
                                       stroke_opacity=1))
            self.dwg.add(self.dwg.text(str(100-label_id*10),
                                       insert=(bar_charts_bounds[0]-8, height+5),
                                       font_size=18,
                                       font_family='Arial',
                                       fill=COLORS['BLACK'],
                                       fill_opacity=1,
                                       text_anchor='end',
                                       alignment_baseline='central'))
        # Bar chart bottom label
        self.dwg.add(self.dwg.text('Percentiles',
                                   insert=(self.canvas_size[0]/2, bar_charts_bounds[3]+50),
                                   font_size=18,
                                   font_family='Arial',
                                   fill=COLORS['BLACK'],
                                   fill_opacity=1,
                                   text_anchor='middle',
                                   alignment_baseline='central'))

        bar_chart_width = bar_charts_bounds[2] - bar_charts_bounds[0]
        for bar_id, bar_name in enumerate(self.bar_names):
            bar_x = bar_charts_bounds[0] + (bar_chart_width * (2*bar_id+1)/4)

            # Bar label
            self.dwg.add(self.dwg.text(bar_name,
                                       insert=(bar_x, bar_charts_bounds[3]+25),
                                       font_size=18,
                                       font_family='Arial',
                                       fill=COLORS['BLACK'],
                                       fill_opacity=1,
                                       text_anchor='middle',
                                       alignment_baseline='central'))

            # Bar
            self.dwg.add(self.dwg.polygon(points=[ (bar_x-bar_width//2, bar_charts_bounds[3]),
                                                   (bar_x-bar_width//2, bar_charts_bounds[1]),
                                                   (bar_x+bar_width//2, bar_charts_bounds[1]),
                                                   (bar_x+bar_width//2, bar_charts_bounds[3]) ],
                                          fill=COLORS['VL_GREY'],
                                          fill_opacity=0.5,
                                          stroke=COLORS['BLACK'],
                                          stroke_width=2,
                                          stroke_opacity=1))
            # Box plot
            box_plot_group = self.dwg.g(id=f'{self.svg_id}-boxplot-{bar_id}', opacity=0)
            box_plot_group.add(self.dwg.polygon(points=[ (bar_x-bar_width//2, bar_charts_bounds[3]),
                                                         (bar_x-bar_width//2, bar_charts_bounds[1]),
                                                         (bar_x+bar_width//2, bar_charts_bounds[1]),
                                                         (bar_x+bar_width//2, bar_charts_bounds[3]) ],
                                                fill=COLORS['WHITE'],
                                                fill_opacity=1,
                                                stroke=COLORS['BLACK'],
                                                stroke_width=2,
                                                stroke_opacity=1))
            box_plot_group.add(self.dwg.polygon(points=[ (bar_x-bar_width//2, bar_charts_bounds[1]+80),
                                                         (bar_x-bar_width//2, bar_charts_bounds[3]-80),
                                                         (bar_x+bar_width//2, bar_charts_bounds[3]-80),
                                                         (bar_x+bar_width//2, bar_charts_bounds[1]+80) ],
                                                fill=f'url(#{self.svg_id}-gradient-{bar_id})',
                                                fill_opacity=0.8,
                                                stroke=COLORS['BLACK'],
                                                stroke_width=2,
                                                stroke_opacity=0.5,
                                                id=f'{self.svg_id}-boxplot-{bar_id}-box'))

            gradient = LinearGradient(start=(0, 0), end=(0,1), id=f'{self.svg_id}-gradient-{bar_id}')
            gradient.add_stop_color(offset='0%', color=COLORS['BAR_GREEN'])
            gradient.add_stop_color(offset='50%', color=COLORS['BAR_ORANGE'])
            gradient.add_stop_color(offset='100%', color=COLORS['BAR_RED'])
            self.dwg.defs.add(gradient)

            box_plot_group.add(self.dwg.line((bar_x-bar_width//2, bar_charts_bounds[1]+200),
                                             (bar_x+bar_width//2, bar_charts_bounds[1]+200),
                                             stroke=COLORS['BLACK'],
                                             stroke_width=2,
                                             stroke_opacity=0.5,
                                             stroke_dasharray=2,
                                             id=f'{self.svg_id}-boxplot-{bar_id}-line-high'))
            box_plot_group.add(self.dwg.line((bar_x-bar_width//2, (bar_charts_bounds[1]+bar_charts_bounds[3])//2),
                                             (bar_x+bar_width//2, (bar_charts_bounds[1]+bar_charts_bounds[3])//2),
                                             stroke=COLORS['BLACK'],
                                             stroke_width=3,
                                             stroke_opacity=0.8,
                                             stroke_dasharray=5,
                                             id=f'{self.svg_id}-boxplot-{bar_id}-line-mid'))
            box_plot_group.add(self.dwg.line((bar_x-bar_width//2, bar_charts_bounds[3]-200),
                                             (bar_x+bar_width//2, bar_charts_bounds[3]-200),
                                             stroke=COLORS['BLACK'],
                                             stroke_width=2,
                                             stroke_opacity=0.5,
                                             stroke_dasharray=2,
                                             id=f'{self.svg_id}-boxplot-{bar_id}-line-low'))
            box_plot_group.add(self.dwg.line((bar_x-bar_width//2, bar_charts_bounds[3]),
                                             (bar_x+bar_width//2, bar_charts_bounds[3]),
                                             fill_opacity=0,
                                             stroke=COLORS['BLACK'],
                                             stroke_width=4,
                                             stroke_opacity=1,
                                             id=f'{self.svg_id}-bar-{bar_id}-mainline'))
            box_plot_group.add(self.dwg.text('',
                                             insert=(bar_x, bar_charts_bounds[3]),
                                             font_size=20,
                                             font_family='Arial',
                                             font_weight='bold',
                                             fill=COLORS['BLACK'],
                                             fill_opacity=1,
                                             text_anchor='middle',
                                             alignment_baseline='central',
                                             id=f'{self.svg_id}-bar-{bar_id}-label'))
            self.dwg.add(box_plot_group)
