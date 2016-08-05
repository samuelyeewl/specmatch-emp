.. _quickstart:

Quickstart
==========

Here's how to get up and running with ``specmatch-emp`` 

Library
-------

``specmatch-emp`` comes with a large library high-resolution optical
spectra shifted onto the rest wavelength scale. We'll import the
library, along with some other useful modules.

::

    import specmatchemp.library
    import pandas as pd
  
    from matplotlib.transforms import blended_transform_factory
    from pylab import *
    rc('savefig',dpi=160)
  
    # Set up HR diagram figure
    def hr_diagram():
        semilogy()
        xlim(xlim()[::-1])
        autoscale(tight='y')
        xlabel('Effective Temperature (K)')
        ylabel('Stellar Radius (Rsun)')
  
Now we'll load in the library around the Mgb triplet

::

    lib = specmatchemp.library.read_hdf(
        '/Users/petigura/Dropbox/SpecMatch-Emp/library.h5',wavlim=[5140,5200]
    )  

    
Here's how the library spans the HR diagram.

::

    fig = figure()
    hr_diagram()
    plot(lib.library_params.Teff, lib.library_params.radius,'.')

.. image:: quickstart-library.png

And here's the library with the sources labeled.

::

    fig = figure()
    hr_diagram()
    g = lib.library_params.groupby('source')
    colors = ['Red','Orange','LimeGreen','Cyan','RoyalBlue','Magenta']
    i = 0
    for source, idx in g.groups.iteritems():
        cut = lib.library_params.ix[idx]
        color = colors[i]
        plot(cut.Teff, cut.radius,'.',label=source,color=color,alpha=0.8,ms=5) 
        i +=1
    legend()
    
.. image:: quickstart-library-labeled.png

The parameters are stored in a pandas DataFrame which makes querying
easy. Let's grab some representative dwarf star spectra.

::

    cut = lib.library_params.query('radius < 1.5 and -0.25 < feh < 0.25')
    g = cut.groupby(pd.cut(cut.Teff,bins=arange(3000,7000,500)))
    cut = g.first()
    
    fig = figure()
    hr_diagram()
    plot(lib.library_params.Teff, lib.library_params.radius,'.')
    plot(cut.Teff, cut.radius,'.')


.. image:: quickstart-library-selected-stars.png

Plot the Mgb region for the spectra sorted by effective temperature

::

    fig,ax = subplots(figsize=(8,4))
    trans = blended_transform_factory(ax.transAxes,ax.transData)
    bbox=dict(facecolor='white', edgecolor='none',alpha=0.8)
    step = 1
    shift = 0
    for _,row in cut.iterrows():
        spec = lib.library_spectra[cut.lib_index,0,:]
        plot(lib.wav,spec.T + shift,color='RoyalBlue',lw=0.5)
        s = "{cps_name:s}, Teff={Teff:.0f}".format(**row)    
        text(0.01,i+1,s,bbox=bbox,transform=trans)
        shift+=step
    
    grid()
    xlabel('Wavelength (Angstroms)')
    ylabel('Normalized Flux (Arbitrary Offset)')
    

.. image:: quickstart-spectra-selected-stars.png
