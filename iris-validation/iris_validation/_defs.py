COLORS = { 'BLACK'      : 'rgb(000, 000, 000)',
           'WHITE'      : 'rgb(255, 255, 255)',

           'GREY'       : 'rgb(050, 050, 050)',
           'L_GREY'     : 'rgb(150, 150, 150)',
           'VL_GREY'    : 'rgb(200, 200, 200)',

           'RED'        : 'rgb(200, 080, 080)',
           'ORANGE'     : 'rgb(250, 200, 050)',
           'GREEN'      : 'rgb(050, 200, 050)',

           'BLUE'       : 'rgb(050, 050, 200)',
           'CYAN'       : 'rgb(050, 200, 200)',
           'TEAL'       : 'rgb(000, 120, 120)',
           'SLATE'      : 'rgb(120, 160, 200)',
           'MAGENTA'    : 'rgb(200, 050, 200)',
           'INDIGO'     : 'rgb(080, 000, 120)',
           'L_PINK'     : 'rgb(255, 235, 235)',

           'BAR_GREEN'  : 'rgb(090, 237, 141)',
           'BAR_ORANGE' : 'rgb(247, 212, 134)',
           'BAR_RED'    : 'rgb(240, 106, 111)'
         }

CONTINUOUS_METRICS = ( { 'id'            : 0,
                         'type'          : 'continuous',
                         'long_name'     : 'Avg. B-factor',
                         'short_name'    : 'Avg. B',
                         'ring_color'    : COLORS['CYAN'],
                         'polarity'      : -1,
                         'is_molprobity' : False
                       },
                       { 'id'            : 1,
                         'type'          : 'continuous',
                         'long_name'     : 'Max. B-factor',
                         'short_name'    : 'Max. B',
                         'ring_color'    : COLORS['TEAL'],
                         'polarity'      : -1,
                         'is_molprobity' : False
                       },
                       { 'id'            : 2,
                         'type'          : 'continuous',
                         'long_name'     : 'Std. B-factor',
                         'short_name'    : 'Std. B',
                         'ring_color'    : COLORS['SLATE'],
                         'polarity'      : -1,
                         'is_molprobity' : False
                       },
                       { 'id'            : 3,
                         'type'          : 'continuous',
                         'long_name'     : 'Residue Fit',
                         'short_name'    : 'Res. Fit',
                         'ring_color'    : COLORS['MAGENTA'],
                         'polarity'      : -1,
                         'is_molprobity' : False
                       },
                       { 'id'            : 4,
                         'type'          : 'continuous',
                         'long_name'     : 'Main Chain Fit',
                         'short_name'    : 'M.C. Fit',
                         'ring_color'    : COLORS['BLUE'],
                         'polarity'      : -1,
                         'is_molprobity' : False
                       },
                       { 'id'            : 5,
                         'type'          : 'continuous',
                         'long_name'     : 'Side Chain Fit',
                         'short_name'    : 'S.C. Fit',
                         'ring_color'    : COLORS['INDIGO'],
                         'polarity'      : -1,
                         'is_molprobity' : False
                       }
                     )

DISCRETE_METRICS =   ( { 'id'            : 0,
                         'type'          : 'discrete',
                         'long_name'     : 'Rotamer Classification',
                         'short_name'    : 'Rota.',
                         'ring_color'    : COLORS['L_GREY'],
                         'seq_colors'    : (COLORS['RED'],
                                            COLORS['ORANGE'],
                                            COLORS['GREEN']),
                         'seq_labels'    : ('Outlier',
                                            'Allowed',
                                            'Favoured'),
                         'is_molprobity' : False
                       },
                       { 'id'            : 1,
                         'type'          : 'discrete',
                         'long_name'     : 'Ramachandran Classification',
                         'short_name'    : 'Rama.',
                         'ring_color'    : COLORS['L_GREY'],
                         'seq_colors'    : (COLORS['RED'],
                                            COLORS['ORANGE'],
                                            COLORS['GREEN']),
                         'seq_labels'    : ('Outlier',
                                            'Allowed',
                                            'Favoured'),
                         'is_molprobity' : False
                       },
                       { 'id'            : 2,
                         'type'          : 'discrete',
                         'long_name'     : 'Clash Indicator',
                         'short_name'    : 'Clashes',
                         'ring_color'    : COLORS['L_GREY'],
                         'seq_colors'    : (COLORS['RED'],
                                            COLORS['ORANGE'],
                                            COLORS['GREEN']),
                         'seq_labels'    : ('Multiple Clashes',
                                            'One Clash',
                                            'No Clashes'),
                         'is_molprobity' : True
                        }
                     )

CHAIN_VIEW_GAP_ANGLE = 0.35
RAMACHANDRAN_THRESHOLDS = (0.02, 0.002)

CHAIN_VIEW_RINGS   = [ DISCRETE_METRICS[0],
                       DISCRETE_METRICS[1],
                       DISCRETE_METRICS[2],
                       CONTINUOUS_METRICS[0],
                       CONTINUOUS_METRICS[1],
                       CONTINUOUS_METRICS[4],
                       CONTINUOUS_METRICS[5] ]

RESIDUE_VIEW_BOXES = [ DISCRETE_METRICS[0],
                       DISCRETE_METRICS[1],
                       DISCRETE_METRICS[2] ]

RESIDUE_VIEW_BARS  = [ CONTINUOUS_METRICS[0],
                       CONTINUOUS_METRICS[5] ]
