"""
@filename plots.py

Helper functions to plot various data from SpecMatch-Emp
"""

import matplotlib.pyplot as plt
import matplotlib.transforms as transforms


def reverse_x():
    """Reverses the x-axis of the current figure"""
    plt.xlim(plt.xlim()[::-1])


def reverse_y():
    """Reverses the y-axis of the current figure"""
    plt.ylim(plt.ylim()[::-1])


def hide_x_ticks():
    """Hide x label ticks"""
    ax = plt.gca()
    ax.axes.get_xaxis().set_ticks([])


def hide_y_ticks():
    """Hide y label ticks"""
    ax = plt.gca()
    ax.axes.get_yaxis().set_ticks([])


def annotate_point(x, y, text, offset=5, offset_x=None, offset_y=None,
                   text_kw={}):
    """Annotates the point at a given x, y position (in data coordinates),
    at a given pixel offset.

    Args:
        x: x-coordinate of point
        y: y-coordinate of point
        text (str): String to annotate
        offset: (optional) pixel offset to use
        offset_x, offset_y: (optional) pixel offset to use in x, y directions
        text_kw (dict): (optional) any additional keywords to pass to plt.text
    """
    if offset_x is None or offset_y is None:
        offset_x = offset
        offset_y = offset
    ax = plt.gca()
    trans_offset = transforms.offset_copy(ax.transData, units='dots',
                                          x=offset_x, y=offset_y)

    plt.text(x, y, text, transform=trans_offset, **text_kw)


def annotate_spectrum(text, spec_offset=0, offset_x=10, offset_y=5,
                      align='left', text_kw={}):
    """Annotates a spectrum.

    Args:
        text (str): String to annotate
        spec_offset: (optional) Vertical offset of spectrum
        offset_x: (optional) Pixel offset from left/right boundary
        offset_y: (optional) Vertical pixel offset from spectrum
        align: (optional) 'left' or 'right' alignment for text
        text_kw (dict): (optional) any additional keywords to pass to plt.text
    """
    ax = plt.gca()
    xlim = ax.get_xlim()
    if align == 'left':
        xpos = xlim[0]
        offset_x = abs(offset_x)
    elif align == 'right':
        xpos = xlim[1]
        offset_x = -abs(offset_x)
    else:
        return

    # transform to pixel coords
    disp_coords = ax.transData.transform((xpos, spec_offset + 1))
    disp_coords = (disp_coords[0] + offset_x, disp_coords[1] + offset_y)
    # invert transform to go back to data coords
    data_coords = ax.transData.inverted().transform(disp_coords)
    ax_coords = ax.transAxes.inverted().transform(disp_coords)
    # fix y position in data coordinates (fixed offset from spectrum)
    # but allow x position to float so we can pan horizontally
    trans = transforms.blended_transform_factory(ax.transAxes, ax.transData)

    bbox = dict(facecolor='white', edgecolor='none', alpha=0.8)
    plt.text(ax_coords[0], data_coords[1], text, bbox=bbox, transform=trans,
             horizontalalignment=align, **text_kw)


def label_axes(param_x=None, param_y=None, rescale=True):
    """Convenience function for tweaking axes to make plots

    Args:
        param_x (str): Parameter to plot on x-axis
        param_y (str): Parameter to plot on y-axis
        rescale (bool): Whether to rescale
    """
    if param_x is 'Teff':
        reverse_x()
        plt.xlabel('Effective Temperature (K)')
        if rescale:
            plt.xticks([3000, 4000, 5000, 6000, 7000])

    if param_x is 'feh':
        plt.xlabel('[Fe/H] (dex)')

    if param_x is 'radius':
        plt.xlabel(r'$R\ (R_\odot)$')
        if rescale:
            ax = plt.gca()
            ax.set_xscale('log')

    if param_y is 'radius':
        plt.ylabel(r'Stellar Radius (Solar-radii)')
        if rescale:
            ax = plt.gca()
            ax.set_yscale('log')
            yt = [0.1, 0.2, 0.3, 0.4, 0.5, 0.7, 1, 2, 3, 4, 5, 7, 10, 20]
            ax.set_yticks(yt, yt)
            ax.set_ylim(0.1, 20)


def set_tight_lims(data_x, data_y, center_x=None, center_y=None,
                   mode='symmetric', buf=0.3):
    """Sets plot limits around a target subset of data, centered at
    a given point.

    Args:
        data_x (np.ndarray): x-coordinates of data
        data_y (np.ndarray): y-coordinates of data
        center_x (optional [float]): x-coordinate of center point
        center_y (optional [float]): y-coordinate of center point
        mode: (optional) 'symmetric': Make limits symmetric about target
                         'tight': Use asymmetric limits
        buf (float): Buffer radius
    """
    ax = plt.gca()

    if center_x is None:
        maxx = max(data_x)
        minx = min(data_x)
        sepx = maxx - minx
        maxx = maxx + buf * sepx
        minx = minx - buf * sepx
        ax.set_xlim((minx, maxx))
    else:
        distx = data_x - center_x
        maxx = max(max(distx), 0)
        minx = min(min(distx), 0)
        if mode == 'symmetric':
            limx = max(abs(maxx), abs(minx))
            limx = limx + buf * limx
            ax.set_xlim((center_x - limx, center_x + limx))
        elif mode == 'tight':
            maxx = maxx + buf * maxx if maxx != 0 else -buf * minx
            minx = minx + buf * minx if minx != 0 else -buf * maxx
            ax.set_xlim((center_x + minx, center_x + maxx))

    if center_y is None:
        maxy = max(data_y)
        miny = min(data_y)
        sepy = maxy - miny
        maxy = maxy + buf * sepy
        miny = miny - buf * sepy
        ax.set_ylim((miny, maxy))
    else:
        disty = data_y - center_y
        maxy = max(max(disty), 0)
        miny = min(min(disty), 0)
        if mode == 'symmetric':
            limy = max(abs(maxy), abs(miny))
            limy = limy + buf * limy
            ax.set_ylim((center_y - limy, center_y + limy))
        elif mode == 'tight':
            maxy = maxy + buf * maxy if maxy != 0 else -buf * miny
            miny = miny + buf * miny if miny != 0 else -buf * maxy
            ax.set_ylim((center_y + miny, center_y + maxy))
