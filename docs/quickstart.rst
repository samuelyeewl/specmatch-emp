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
  
Now we'll load in the library around the Mgb triplet. By default,
SpecMatch create the following directory ``${HOME}/.specmatchemp/`` and
download the library into it.

::

   # Only load spectra from 5140 to 5200 angstroms.
   lib = specmatchemp.library.read_hdf(wavlim=[5140,5200])

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


Matching
--------

To see how to perform a match with SpecMatch, we will look at how to
perform the matching on two example stars - HD 190406, a G0V star, as
well as Barnard's Star (GL 699), an M dwarf.

Once again, we import the library module and read in the library. We
obtain our two target spectra from the library and remove them from
the matching process.

::

    import specmatchemp.library
    lib = specmatchemp.library.read_hdf()
    
    idx1 = lib.get_index('190406')
    G_star = lib.pop(idx1)
    idx2 = lib.get_index('GL699')
    M_star = lib.pop(idx2)


To perform SpecMatch, we import the specmatch module and create a 
SpecMatch object. As an example, we will use the wavelength region
from 5300 - 5400 Angstroms.

The match method first compares the target spectrum against each of
the library spectra. It then synthesizes linear combinations of the
best matching spectra.

::
    
    from specmatchemp.specmatch import SpecMatch
    match_G = SpecMatch(G_star[1], lib, (5300,5400))
    match_G.match()

The final derived parameters can be found in the results attribute.

::

    print('Derived Parameters: ')
    print('Teff: {0:.0f}, Radius: {1:.2f}, [Fe/H]: {2:.2f}'.format(\
            match_G.results['Teff'], match_G.results['radius'], match_G.results['feh']))
    print('Library Parameters: ')
    print('Teff: {0:.0f}, Radius: {1:.2f}, [Fe/H]: {2:.2f}'.format(\
            G_star[0]['Teff'], G_star[0]['radius'], G_star[0]['feh']))

::

    Derived Parameters: 
    Teff: 5855, Radius: 1.36, [Fe/H]: 0.06
    Library Parameters: 
    Teff: 5763, Radius: 1.12, [Fe/H]: 0.03

We can take a closer look at the workings of the matching process. First,
examine the chi-squared surfaces of the match with the library spectra, as
a function of stellar parameters.

::

    fig = figure(figsize=(12,8))
    match_G.plot_chi_squared_surface()
    # Indicate library parameters for target star.
    fig.add_subplot(131)
    axvline(x=G_star[0].Teff, color='k')
    fig.add_subplot(132)
    axvline(x=G_star[0].radius, color='k')
    fig.add_subplot(133)
    axvline(x=G_star[0].feh, color='k')


.. image:: quickstart-Gstar-chisquared-surface.png

The 5 closest matches, which were used to synthesize the linear combinations,
have been highlighted. The position of these stars in the HR-diagram can be
seen below.

::

    match_G.plot_references()