# code-start-imports
from matplotlib.transforms import blended_transform_factory
import pandas as pd
from pylab import *
import specmatchemp.library
rc('savefig',dpi=160)

def hr_diagram():
    semilogy()
    xlim(xlim()[::-1])
    autoscale(tight='y')
    xlabel('Effective Temperature (K)')
    ylabel('Stellar Radius (Rsun)')
# code-stop-imports




# code-start-loadlibrary: load in the library around the Mgb triplet
lib = specmatchemp.library.read_hdf(wavlim=[5140,5200])
# code-stop-loadlibrary




# code-start-library: Here's how the library spans the HR diagram.
fig = figure()
hr_diagram()
plot(lib.library_params.Teff, lib.library_params.radius,'.')
# code-stop-library
fig.savefig('quickstart-library.png')




# code-start-library-labeled
fig = figure()
hr_diagram()
g = lib.library_params.groupby('source')
colors = ['Red','Orange','LimeGreen','Cyan','RoyalBlue','Magenta','ForestGreen']
i = 0
for source, idx in g.groups.items():
    cut = lib.library_params.ix[idx]
    color = colors[i]
    plot(cut.Teff, cut.radius,'.',label=source,color=color,alpha=0.8,ms=5) 
    i +=1
legend()
# code-stop-library-labeled
fig.savefig('quickstart-library-labeled.png')





# code-start-library-selected-stars
cut = lib.library_params.query('radius < 1.5 and -0.25 < feh < 0.25')
g = cut.groupby(pd.cut(cut.Teff,bins=arange(3000,7000,500)))
cut = g.first()

fig = figure()
hr_diagram()
plot(lib.library_params.Teff, lib.library_params.radius,'.')
plot(cut.Teff, cut.radius,'o')
fig.savefig('quickstart-library-selected-stars.png')
# code-stop-library-selected-stars





# code-start-spectra-selected-stars
fig,ax = subplots(figsize=(8,4))
trans = blended_transform_factory(ax.transAxes,ax.transData)
bbox = dict(facecolor='white', edgecolor='none',alpha=0.8)
step = 1
shift = 0
for _,row in cut.iterrows():
    spec = lib.library_spectra[row.lib_index,0,:]
    plot(lib.wav,spec.T + shift,color='RoyalBlue',lw=0.5)
    s = "{cps_name:s}, Teff={Teff:.0f}".format(**row)    
    text(0.01, 1+shift, s, bbox=bbox, transform=trans)
    shift+=step

grid()
xlabel('Wavelength (Angstroms)')
ylabel('Normalized Flux (Arbitrary Offset)')
# code-stop-spectra-selected-stars
fig.set_tight_layout(True)
fig.savefig('quickstart-spectra-selected-stars.png')





# code-start-specmatch-load
lib = specmatchemp.library.read_hdf(wavlim=[5300,5400])

idx1 = lib.get_index('190406')
G_star = lib.pop(idx1)
idx2 = lib.get_index('GL699')
M_star = lib.pop(idx2)
# code-stop-specmatch-load




# code-start-specmatch-match
from specmatchemp.specmatch import SpecMatch
match_G = SpecMatch(G_star[1], lib, (5300,5400))
match_G.match()
# code-stop-specmatch-match



# code-start-specmatch-print
print('Derived Parameters: ')
print('Teff: {0:.0f}, Radius: {1:.2f}, [Fe/H]: {2:.2f}'.format(\
        match_G.results['Teff'], match_G.results['radius'], match_G.results['feh']))
print('Library Parameters: ')
print('Teff: {0:.0f}, Radius: {1:.2f}, [Fe/H]: {2:.2f}'.format(\
        G_star[0]['Teff'], G_star[0]['radius'], G_star[0]['feh']))
# code-stop-specmatch-print


# code-start-plot-chisquared
fig = figure(figsize=(12,8))
match_G.plot_chi_squared_surface()
# Indicate library parameters for target star.
axes = fig.axes
axes[0].axvline(G_star[0]['Teff'], color='k')
axes[1].axvline(G_star[0]['radius'], color='k')
axes[2].axvline(G_star[0]['feh'], color='k')
# code-stop-plot-chisquared

fig.set_tight_layout(True)
fig.savefig('quickstart-Gstar-chisquared-surface.png')


# code-start-plot-references
fig = figure(figsize=(12,10))
match_G.plot_references(verbose=True)
# plot target onto HR diagram
axes = fig.axes
axes[0].plot(G_star[0]['Teff'], G_star[0]['radius'], '*', ms=15, color='red', label='Target')
axes[1].plot(G_star[0]['Teff'], G_star[0]['radius'], '*', ms=15, color='red')
axes[2].plot(G_star[0]['feh'], G_star[0]['radius'], '*', ms=15, color='red')
axes[3].plot(G_star[0]['feh'], G_star[0]['radius'], '*', ms=15, color='red')
axes[0].legend(numpoints=1, fontsize='small', loc='best')
# code-stop-plot-references

fig.savefig('quickstart-Gstar-lincomb-references.png')