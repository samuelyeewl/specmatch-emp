# code-start-imports
import pandas as pd
from pylab import *
import specmatchemp.library
import specmatchemp.plotting.plots as smplot

# code-stop-imports
rc('savefig',dpi=160)




# code-start-loadlibrary: load in the library around the Mgb triplet
lib = specmatchemp.library.read_hdf(wavlim=[5140,5200])
# code-stop-loadlibrary



# code-start-library: Here's how the library spans the HR diagram.
fig = figure()
plot(lib.library_params.Teff, lib.library_params.radius,'k.',)
smplot.label_axes('Teff','radius')
# code-stop-library
fig.savefig('quickstart-library.png')




# code-start-library-labeled
fig = figure()
g = lib.library_params.groupby('source')
colors = ['Red','Orange','LimeGreen','Cyan','RoyalBlue','Magenta','ForestGreen']
i = 0
for source, idx in g.groups.items():
    cut = lib.library_params.ix[idx]
    color = colors[i]
    plot(
        cut.Teff, cut.radius,'+', label=source, color=color, alpha=1, ms=5, 
        mew=1.5
    ) 
    i+=1
legend()
smplot.label_axes('Teff','radius')
# code-stop-library-labeled
fig.savefig('quickstart-library-labeled.png')





# code-start-library-selected-stars
cut = lib.library_params.query('radius < 1.5 and -0.25 < feh < 0.25')
g = cut.groupby(pd.cut(cut.Teff,bins=arange(3000,7000,500)))
cut = g.first()

fig = figure()
plot(lib.library_params.Teff, lib.library_params.radius,'.', label='_nolegend_')
plot(cut.Teff, cut.radius,'o', label='Selected Stars')
legend()
smplot.label_axes('Teff','radius')
fig.savefig('quickstart-library-selected-stars.png')
# code-stop-library-selected-stars





# code-start-spectra-selected-stars
from matplotlib.transforms import blended_transform_factory
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

# code-start-plot-lincomb-spectra
fig = plt.figure(figsize=(12,6))
match_G.plot_lincomb()
# code-end-plot-lincomb-spectra

fig.set_tight_layout(True)
fig.savefig('quickstart-Gstar-lincomb-spectra.png')

# code-start-mstar
match_M = SpecMatch(M_star[1], lib, (5300,5400))
match_M.match()

print('Derived Parameters: ')
print('Teff: {0:.0f}, Radius: {1:.2f}, [Fe/H]: {2:.2f}'.format(\
        match_M.results['Teff'], match_M.results['radius'], match_M.results['feh']))
print('Library Parameters: ')
print('Teff: {0:.0f}, Radius: {1:.2f}, [Fe/H]: {2:.2f}'.format(\
        M_star[0]['Teff'], M_star[0]['radius'], M_star[0]['feh']))

fig1 = figure(figsize=(12,8))
match_M.plot_chi_squared_surface()
# Indicate library parameters for target star.
axes = fig1.axes
axes[0].axvline(M_star[0]['Teff'], color='k')
axes[1].axvline(M_star[0]['radius'], color='k')
axes[2].axvline(M_star[0]['feh'], color='k')

fig2 = figure(figsize=(12,10))
match_M.plot_references(verbose=True)
# plot target onto HR diagram
axes = fig2.axes
axes[0].plot(M_star[0]['Teff'], M_star[0]['radius'], '*', ms=15, color='red', label='Target')
axes[1].plot(M_star[0]['Teff'], M_star[0]['radius'], '*', ms=15, color='red')
axes[2].plot(M_star[0]['feh'], M_star[0]['radius'], '*', ms=15, color='red')
axes[3].plot(M_star[0]['feh'], M_star[0]['radius'], '*', ms=15, color='red')
axes[0].legend(numpoints=1, fontsize='small', loc='best')

fig3 = plt.figure(figsize=(12,6))
match_M.plot_lincomb()
# code-end-mstar

fig1.set_tight_layout(True)
fig1.savefig('quickstart-Mstar-chisquared-surface.png')
fig2.savefig('quickstart-Mstar-lincomb-references.png')
fig3.set_tight_layout(True)
fig3.savefig('quickstart-Mstar-lincomb-spectra.png')