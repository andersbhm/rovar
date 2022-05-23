'''
Anders Lindanger
2020

Make matplotlib look nicer.
'''

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import copy
import matplotlib as mpl
from matplotlib import rcParams
import matplotlib.lines as mlines

#################### SIZES #####################################################
font_size_title = 17
font_size_labels = 16
font_size_numbers_axis = 16
font_size_legend = 13
###############################################################################
rcParams['axes.autolimit_mode'] = 'data'
rcParams['axes.axisbelow'] = 'line'
rcParams['axes.edgecolor'] = 'black'
rcParams['axes.facecolor'] = 'white'
rcParams['axes.formatter.limits'] = [-5, 6]
rcParams['axes.formatter.min_exponent'] = 0
rcParams['axes.formatter.offset_threshold'] = 4
rcParams['axes.formatter.use_locale'] = False
rcParams['axes.formatter.use_mathtext'] = False
rcParams['axes.formatter.useoffset'] = True
rcParams['axes.grid'] = False
rcParams['axes.grid.axis'] = 'both'
rcParams['axes.grid.which'] = 'major'
rcParams['axes.labelcolor'] = 'black'
rcParams['axes.labelpad'] = 4.0
rcParams['axes.labelsize'] = font_size_labels
rcParams['axes.labelweight'] = 'normal'
rcParams['axes.linewidth'] = 0.8
rcParams['axes.spines.bottom'] = True
rcParams['axes.spines.left'] = True
rcParams['axes.spines.right'] = True
rcParams['axes.spines.top'] = True


try:
    rcParams['axes.titlecolor'] = 'auto'
    rcParams['axes.titlelocation'] = 'center'
except KeyError:
    print("Could not set rcParams 'axes.titlecolor' and 'axes.titlelocation', but that does not matter much.")


rcParams['axes.titlepad'] = 6.0
rcParams['axes.titlesize'] = font_size_title
rcParams['axes.titleweight'] = 'normal'
rcParams['axes.unicode_minus'] = True
rcParams['axes.xmargin'] = 0.0
rcParams['axes.ymargin'] = 0.1
rcParams['axes3d.grid'] = True
rcParams['backend'] = 'TkAgg'
rcParams['backend_fallback'] = True

rcParams['date.autoformatter.day'] = '%Y-%m-%d'
rcParams['date.autoformatter.hour'] = '%m-%d %H'
rcParams['date.autoformatter.microsecond'] = '%M:%S.%f'
rcParams['date.autoformatter.minute'] = '%d %H:%M'
rcParams['date.autoformatter.month'] = '%Y-%m'
rcParams['date.autoformatter.second'] = '%H:%M:%S'
rcParams['date.autoformatter.year'] = '%Y'
rcParams['errorbar.capsize'] = 0.0

rcParams['figure.constrained_layout.h_pad'] = 0.04167
rcParams['figure.constrained_layout.hspace'] = 0.02
rcParams['figure.constrained_layout.use'] = False
rcParams['figure.constrained_layout.w_pad'] = 0.04167
rcParams['figure.constrained_layout.wspace'] = 0.02
rcParams['figure.dpi'] = 100.0
rcParams['figure.edgecolor'] = 'white'
rcParams['figure.facecolor'] = 'white'
rcParams['figure.figsize'] = [8, 5]
rcParams['figure.frameon'] = False
rcParams['figure.subplot.bottom'] = 0.11
rcParams['figure.subplot.hspace'] = 0.2
rcParams['figure.subplot.left'] = 0.125
rcParams['figure.subplot.right'] = 0.9
rcParams['figure.subplot.top'] = 0.88
rcParams['figure.subplot.wspace'] = 0.2
rcParams['figure.titlesize'] = 'large'
rcParams['figure.titleweight'] = 'normal'

rcParams['font.family'] = ['Times New Roman']
rcParams['font.serif'] = ['Times New Roman']
rcParams['font.size'] = 10.0
rcParams['grid.alpha'] = 0.7
rcParams['grid.linestyle'] = '-'
rcParams['grid.linewidth'] = 0.6
rcParams['interactive'] = False

rcParams['legend.borderaxespad'] = 0.5
rcParams['legend.borderpad'] = 0.4
rcParams['legend.columnspacing'] = 2.0
rcParams['legend.facecolor'] = 'inherit'
rcParams['legend.fancybox'] = False
rcParams['legend.fontsize'] = font_size_legend
rcParams['legend.framealpha'] = 0.5
rcParams['legend.frameon'] = True
rcParams['legend.handleheight'] = 0.7
rcParams['legend.handlelength'] = 2.0
rcParams['legend.handletextpad'] = 0.3
rcParams['legend.labelspacing'] = 0.1
rcParams['legend.loc'] = 'best'
rcParams['legend.markerscale'] = 1.0
rcParams['legend.numpoints'] = 1
rcParams['legend.scatterpoints'] = 1
rcParams['legend.shadow'] = False

try:
    rcParams['legend.title_fontsize'] = None
except KeyError:
    print("Could not set rcParams 'legend.title_fontsize', but that does not matter much.")


rcParams['lines.antialiased'] = True
rcParams['lines.color'] = 'C0'
rcParams['lines.dash_capstyle'] = 'butt'
rcParams['lines.dash_joinstyle'] = 'round'
rcParams['lines.dashdot_pattern'] = [6.4, 1.6, 1.0, 1.6]
rcParams['lines.dashed_pattern'] = [3.7, 1.6]
rcParams['lines.dotted_pattern'] = [1.0, 1.65]
rcParams['lines.linestyle'] = '-'
rcParams['lines.linewidth'] = 1
rcParams['lines.marker'] = None
rcParams['lines.markersize'] = 6.0
rcParams['markers.fillstyle'] = 'full'
rcParams['mathtext.bf'] = 'Times New Roman:bold'
rcParams['mathtext.cal'] = 'cursive'
rcParams['mathtext.default'] = 'rm'
#rcParams['mathtext.fallback_to_cm'] = True
rcParams['mathtext.fontset'] = 'custom'
rcParams['mathtext.it'] = 'Times New Roman:italic'
rcParams['mathtext.rm'] = 'Times New Roman:roman'
rcParams['mathtext.sf'] = 'Times New Roman'
rcParams['mathtext.tt'] = 'monospace'


rcParams['savefig.dpi'] = 300.0
rcParams['savefig.edgecolor'] = 'white'
rcParams['savefig.facecolor'] = 'white'
rcParams['savefig.format'] = 'pdf'
rcParams['savefig.orientation'] = 'portrait'
rcParams['savefig.pad_inches'] = 0.1
rcParams['savefig.transparent'] = False

try:
    rcParams['scatter.edgecolors'] = 'face'
except KeyError:
    print("Could not set rcParams 'scatter.edgecolors', but that does not matter much.")

rcParams['scatter.marker'] = 'x'


#rcParams['text.latex.preview'] = False
rcParams['text.usetex'] = False

rcParams['xtick.alignment'] = 'center'
rcParams['xtick.bottom'] = True
rcParams['xtick.color'] = 'black'
rcParams['xtick.direction'] = 'in'
rcParams['xtick.labelbottom'] = True
rcParams['xtick.labelsize'] = font_size_numbers_axis
rcParams['xtick.labeltop'] = False
rcParams['xtick.major.bottom'] = True
rcParams['xtick.major.pad'] = 6.
rcParams['xtick.major.size'] = 6.0
rcParams['xtick.major.top'] = True
rcParams['xtick.major.width'] = 1.0
rcParams['xtick.minor.bottom'] = True
rcParams['xtick.minor.pad'] = 3.4
rcParams['xtick.minor.size'] = 3.0
rcParams['xtick.minor.top'] = True
rcParams['xtick.minor.visible'] = True
rcParams['xtick.minor.width'] = 1.0
rcParams['xtick.top'] = True
rcParams['ytick.alignment'] = 'center_baseline'
rcParams['ytick.color'] = 'black'
rcParams['ytick.direction'] = 'in'
rcParams['ytick.labelleft'] = True
rcParams['ytick.labelright'] = False
rcParams['ytick.labelsize'] = font_size_numbers_axis
rcParams['ytick.left'] = True
rcParams['ytick.major.left'] = True
rcParams['ytick.major.pad'] = 3.5
rcParams['ytick.major.right'] = True
rcParams['ytick.major.size'] = 6.0
rcParams['ytick.major.width'] = 1.0
rcParams['ytick.minor.left'] = True
rcParams['ytick.minor.pad'] = 3.4
rcParams['ytick.minor.right'] = True
rcParams['ytick.minor.size'] = 3.0
rcParams['ytick.minor.visible'] = True
rcParams['ytick.minor.width'] = 1.0
rcParams['ytick.right'] = True

#for i,j in zip(rcParams.keys(), rcParams.values()):
#    print("rcParams['%s'] = %s "%(i,j))



if __name__ == '__main__':
    print("")
