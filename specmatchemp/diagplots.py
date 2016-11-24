"""
@filename diagplots.py

Functions to create library diagnostic plots.
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

from specmatchemp import plots
from specmatchemp import detrend
from specmatchemp.library import Library


def library_comparison(params, param_x, param_y, suffix='_sm', ptlabels=False,
                       legend=True, rescale=True):
    """Create a library comparison plot showing the library values of the
    parameters as points, with lines point to derived parameters.

    Args:
        params (pd.DataFrame): Parameter table with library and derived values
        param_x (str): Parameter to plot on the x-axis.
        param_y (str): Parameter to plot on the y-axis.
        suffix (str): Suffix on columns in parameter table with derived values
        ptlabels (str, optional): Column to annotate points with.
            Defaults to False, denoting no labels.
    """
    if param_x not in Library.STAR_PROPS:
        raise ValueError("{0} is not a valid parameter.".format(param_x))
    if param_y not in Library.STAR_PROPS:
        raise ValueError("{0} is not a valid parameter.".format(param_y))

    plt.plot(params[param_x], params[param_y], 'ko', label='Library value')

    x = params[[param_x+suffix, param_x]]
    y = params[[param_y+suffix, param_y]]
    plt.plot(x.T, y.T, 'r')
    plt.plot(x.iloc[0], y.iloc[0], 'r', label='SpecMatch-Emp value')
    plots.label_axes(param_x, param_y, rescale)
    if legend:
        plt.legend(loc='best')

    if ptlabels is not False and ptlabels in params.columns:
        params.apply(lambda row: plots.annotate_point(
            row[param_x], row[param_y], row[ptlabels]), axis=1)


def library_difference(params, prop, suffix='_sm', ptlabels=False,
                       rescale=True, plt_kw={'color': 'blue'}):
    """Plot the residuals (library-derived) for each star in the library.

    Args:
        params (pd.DataFrame): Parameter table with library and derived values
        prop (str): Parameter to plot on the x-axis.
        suffix (str): Suffix on columns in parameter table with derived values
    """
    if prop == 'radius':
        resid = (params[prop+suffix] - params[prop])/params[prop]
        plt.semilogx(params[prop], resid, 'o', **plt_kw)
        plt.xlim(0.1, 20)
    elif prop == 'mass':
        resid = (params[prop+suffix] - params[prop])/params[prop]
        plt.plot(params[prop], resid, 'o', **plt_kw)
    else:
        resid = params[prop+suffix] - params[prop]
        plt.plot(params[prop], resid, 'o', **plt_kw)

    if ptlabels is not False and ptlabels in params.columns:
        params['resid'] = resid
        params.apply(lambda row: plots.annotate_point(
            row[prop], row['resid'], row[ptlabels]), axis=1)

    mean = np.mean(resid)
    rms = np.sqrt(np.mean(resid**2))

    ax = plt.gca()
    bbox = dict(facecolor='white', edgecolor='none', alpha=0.8)
    plt.text(0.05, 0.1, "Mean Diff: {0:.3g}\nRMS Diff: {1:.3g}"
             .format(mean, rms), transform=ax.transAxes, bbox=bbox)
    plt.axhline(y=0, color='k', linestyle='dashed')

    plots.label_axes(param_x=prop, rescale=rescale)


def five_pane(params, suffix, trend=False, ptlabels=False, rescale=True):
    """Five panel diagnostic plot
    """
    if trend:
        d = detrend.Detrend()

    gs = gridspec.GridSpec(6, 2)

    plt.subplot(gs[0:3, 0])
    library_comparison(params, 'Teff', 'radius', suffix, ptlabels=ptlabels,
                       rescale=rescale)

    plt.subplot(gs[3:6, 0])
    library_comparison(params, 'feh', 'radius', suffix, ptlabels=ptlabels,
                       legend=False, rescale=rescale)

    plt.subplot(gs[0:2, 1])
    library_difference(params, 'Teff', suffix=suffix, ptlabels=ptlabels,
                       rescale=rescale)
    if trend:
        xlim = plt.xlim()
        ylim = plt.ylim()
        d.plot('Teff')
        plt.xlim(xlim)
        plt.ylim(ylim)
    plt.ylabel(r'$\Delta\ T_{\mathrm{eff}}$ (K)')

    plt.subplot(gs[2:4, 1])
    library_difference(params, 'radius', suffix=suffix, ptlabels=ptlabels,
                       rescale=rescale)
    if trend:
        xlim = plt.xlim()
        ylim = plt.ylim()
        d.plot('radius')
        plt.xlim(xlim)
        plt.ylim(ylim)
    plt.ylabel(r'$\Delta R_{\star}/R_{\star}$')

    plt.subplot(gs[4:6, 1])
    library_difference(params, 'feh', suffix=suffix, ptlabels=ptlabels,
                       rescale=rescale)
    if trend:
        xlim = plt.xlim()
        ylim = plt.ylim()
        d.plot('feh')
        plt.xlim(xlim)
        plt.ylim(ylim)
    plt.ylabel(r'$\Delta \mathrm{[Fe/H]}$ (dex)')
