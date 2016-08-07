"""
@filename plots.py

Helper functions to plot various data from SpecMatch-Emp
"""

import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import matplotlib.transforms as transforms

def reverse_x():
    """Reverses the x-axis of the current figure"""
    plt.xlim(plt.xlim()[::-1])

def reverse_y():
    """Reverses the y-axis of the current figure"""
    plt.ylim(plt.ylim()[::-1])

def annotate_point(x, y, text, offset=5, offset_x=None, offset_y=None, text_kw={}):
    """Annotates the point at a given x, y position (in data coordinates),
    at a given pixel offset.

    Args:
        x: x-coordinate in data space
        y: y-coordinate in data space
        text (str): String to annotate
        offset: (optional) pixel offset to use
        offset_x, offset_y: (optional) pixel offset to use in x, y directions
        text_kw (dict): (optional) any additional keywords to pass to pyplot.text
    """
    ax = plt.gca()
    # transform to pixel coords
    disp_coords = ax.transData.transform((x, y))
    if offset_x is None or offset_y is None:
        offset_x = offset
        offset_y = offset
    disp_coords = (disp_coords[0]+offset_x, disp_coords[1]+offset_y)
    # invert transform to go back to data coords
    data_coords = ax.transData.inverted().transform(disp_coords)
    plt.text(data_coords[0], data_coords[1], text, **text_kw)

def annotate_spectrum(text, spec_offset=0, offset_x=10, offset_y=5, align='left', text_kw={'fontsize':'small'}):
    """Annotates a spectrum.

    Args:
        text (str): String to annotate
        spec_offset: (optional) Vertical offset of spectrum
        offset_x: (optional) Pixel offset from left/right boundary
        offset_y: (optional) Vertical pixel offset from spectrum
        align: (optional) 'left' or 'right' alignment for text
        text_kw (dict): (optional) any additional keywords to pass to pyplot.text
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
    disp_coords = ax.transData.transform((xpos, spec_offset+1))
    disp_coords = (disp_coords[0]+offset_x, disp_coords[1]+offset_y)
    # invert transform to go back to data coords
    data_coords = ax.transData.inverted().transform(disp_coords)
    ax_coords = ax.transAxes.inverted().transform(disp_coords)
    # fix y position in data coordinates (fixed offset from spectrum)
    # but allow x position to float
    trans = transforms.blended_transform_factory(ax.transAxes, ax.transData)

    bbox=dict(facecolor='white', edgecolor='none',alpha=0.8)
    plt.text(ax_coords[0], data_coords[1], text, bbox=bbox, transform=trans, horizontalalignment=align, **text_kw)
    