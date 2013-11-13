import numpy as np
import matplotlib.pyplot as plt
import math
from os import sep, mkdir


def smoothed(unsmoothed, factor):
    l = len(unsmoothed)
    r = [math.exp(-(i - l) ** 2 / float(factor)) for i in xrange(1, l)]
    w = np.array(r + [1] + r[::-1])
    c = (w[l - i - 1:2 * l - i - 1] for i in xrange(l))
    wm = np.array([x / sum(x) for x in c])
    return np.dot(unsmoothed, wm.transpose())


def plot(plotdata, directory, bottom=True, side=True, legend=True,
         save=True, filename='untitled.pdf', upperbound=0.05, factor=21,
         fig=None, **kwargs):
    if fig is None:
        fig = plt.figure(None, facecolor='w', edgecolor='w')

    if not directory.endswith(sep):
        directory += sep

    try:
        mkdir(directory)
    except OSError:
        pass

    # plotting data
    ntvar = plotdata['nt']['var']
    aavar = plotdata['aa']['var']

    # gene models
    starts = plotdata['aa']['starts']
    ends = plotdata['aa']['ends']
    counts = plotdata['aa']['count']

    # smooth the data
    snt = smoothed(ntvar, factor)
    lnt = len(ntvar)
    saa = smoothed(aavar, factor)
    laa = len(aavar)

    # bounding rectangle
    bound = [0, laa, -upperbound / 6.0, upperbound]

    # x-values to align nucleotide & amino acids
    xnt = np.arange(lnt) / 3.0 + 1
    xaa = np.arange(laa) + 1

    ax = axes(bottom, side, bound, fig, **kwargs)
    ntl = draw(xnt, snt, ax, '#0000ff', **kwargs)
    aal = draw(xaa, saa, ax, '#00ff00', **kwargs)
    models(starts, ends, counts, bound, ax, **kwargs)
    report(filename, ntvar, aavar, lnt, laa)

    if legend:
        fig.legend((ntl, aal), ('Nucleotide', 'Amino acid'), 'upper right')
    if save:
        fig.savefig(directory + filename)


def axes(bottom, side, bound, fig, **kwargs):
    # create the proper sized frame, depending on
    # how we draw the plot
    x = 0.09 if side else 0.02
    y = 0.09 if bottom else 0.04

    xs = [bound[0], bound[1] * 1.06]
    # construct the axes
    ax = fig.add_axes([x, y, 0.98 - x, 0.98 - y], xlim=xs, ylim=bound[2:])
    ax.minorticks_on()
    ax.tick_params(axis='x', which='minor', length=3)

    # hide the unwanted axis lines (typically the top & right)
    # label the wanted axes and draw them
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

    ax.hlines(ax.get_yticks(), xs[0], xs[1], color='0.75',
              linestyle='dashed')
    ax.hlines(0, xs[0], xs[1], color='k', linestyle='solid')
    return ax


def draw(x, y, ax, color, **kwargs):
    return ax.plot(x, y, color=color, linestyle='solid')


def models(starts, ends, counts, bound, ax, **kwargs):
    lb, l = bound[2], len(starts)
    scale = bound[1] / max(ends)
    ys, i = np.arange(1, l + 1) * lb / 3.0, 0

    # draw the gene models
    ax.hlines(ys, starts, ends, colors='k', lw=4, linestyle='solid')
    for c in counts:
        ax.text(bound[1] + 10, lb / 3.0 * (i + 1.25), int(c))
        i += 1


def report(filename, ntvar, aavar, lnt, laa):
    print(filename.center(80, '='))
    print('Average variance: ')
    print('\t%f per base pair' % (sum(ntvar) / lnt))
    print('\t%f per amino acid' % (sum(aavar) / laa))
