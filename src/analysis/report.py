from __future__ import print_function
import matplotlib.pyplot as plt
from biotools.analysis.plot import smoothed
from os import mkdir


# wraps biotools.analysis.plot.plot()
def plot(plotdata, directory, bottom=True, side=True, legend=True,
         save=True, filename='untitled.pdf', upperbound=0.05, factor=21,
         fig=plt.figure(None, facecolor='w', edgecolor='w'), **kwargs):
    try:
        mkdir(directory)
    except OSError:
        pass

    lowerbound = -upperbound / 6

    # smooth the data
    snt = smoothed(plotdata['nt']['var'], factor)
    lnt = len(plotdata['nt']['var'])
    saa = smoothed(plotdata['aa']['var'], factor)
    laa = len(plotdata['aa']['var'])

    # generate x-ranges so that amino acids
    # and nucleotides align
    xnt = np.arange(lnt) * (0.0 + laa) / lnt + 1
    xaa = np.arange(laa) + 1

    # create the proper sized frame, depending on
    # how we draw the plot
    x = 0.09 if side else 0.02
    y = 0.09 if bottom else 0.04

    ax = fig.add_axes([x, y, 0.98 - x, 0.98 - y], xlim=[0, laa * 1.06],
                      ylim=[lowerbound, upperbound])
    ax.minorticks_on()
    ax.tick_params(axis='x', which='minor', length=3)

    for loc, spine in ax.spines.iteritems():
        if loc in ['right', 'top']:
            spine.set_color('none')
            continue
        if loc == 'bottom':
            ax.xaxis.set_ticks_position('bottom')
            if bottom:
                ax.set_xlabel("Amino acids")
            continue
        if loc == 'left':
            if side:
                ax.set_ylabel("Sequence variance")
                ax.yaxis.set_ticks_position('left')
            else:
                spine.set_color('none')
                ax.tick_params('y', which='both', color='none',
                               labelcolor='none')

    ax.hlines(ax.get_yticks(), 0, laa * 1.06, color='0.75', linestyle='dashed')
    ax.hlines(0, 0, laa * 1.06, color='k', linestyle='solid')

    nt_lines = ax.plot(xnt, snt, color='#0000ff', linestyle='solid')
    aa_lines = ax.plot(xaa, saa, color='#00ff00', linestyle='solid')

    starts = plotdata['aa']['starts']
    ends = plotdata['aa']['ends']
    counts = plotdata['aa']['count']
    scale = laa / max(ends)
    ys = (np.arange(len(starts)) + 1) * lowerbound / 3

    ax.hlines(ys, starts, ends, colors='k', lw=4, linestyle='solid')
    for i, c in enumerate(counts):
        ax.text(laa + 10, lowerbound / 3 * (i + 0.25), c)
    if legend:
        fig.legend((nt_lines, aa_lines), ('Nucleotide', 'Amino acid'),
                   'upper right')

    if save:
        fig.savefig(directory + filename)

    print('=============', filename, '=============')
    print('Average variance: ')
    print('\t', sum(plotdata['nt']['var']) / lnt, 'per base pair')
    print('\t', sum(plotdata['aa']['var']) / laa, 'per amino acid')
