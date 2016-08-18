# code-start-imports
import pandas as pd
from pylab import *
import specmatchemp.library
import specmatchemp.plots as smplot
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


# code-start-pop-library
idx1 = lib.get_index('190406')
G_star = lib.pop(idx1)
idx2 = lib.get_index('GL699')
M_star = lib.pop(idx2)
# code-stop-pop-library


# code-start-read-spectrum-G
from specmatchemp import spectrum
G_spectrum = spectrum.read_hires_fits('../samples/rj59.1923.fits').cut(5130,5210)
G_spectrum.name = 'HD190406'
# code-stop-read-spectrum-G


# code-start-shift-spectrum-G
from specmatchemp.specmatch import SpecMatch
sm_G = SpecMatch(G_spectrum, lib)
sm_G.shift()
# code-stop-shift-spectrum-G


# code-start-plot-shifts-G
fig = plt.figure(figsize=(10,5))
sm_G.target_unshifted.plot(normalize=True, plt_kw={'color':'forestgreen'}, text='Target (unshifted)')
sm_G.target.plot(offset=0.5, plt_kw={'color':'royalblue'}, text='Target (shifted): HD190406')
sm_G.shift_ref.plot(offset=1, plt_kw={'color':'firebrick'}, text='Reference: '+sm_G.shift_ref.name)
plt.xlim(5160,5200)
plt.ylim(0,2.2)
# code-stop-plot-shifts-G
fig.set_tight_layout(True)
fig.savefig('quickstart-Gstar-shifts.png')


# code-start-shift-spectrum-M
# Load spectrum
M_spectrum = spectrum.read_hires_fits('../samples/rj130.2075.fits').cut(5130,5210)
M_spectrum.name = 'GL699'

# Shift spectrum
sm_M = SpecMatch(M_spectrum, lib)
sm_M.shift()

# Plot shifts
fig = plt.figure(figsize=(10,5))
sm_M.plot_shifted_spectrum(wavlim=(5160,5200))
# code-stop-shift-spectrum-M
fig.set_tight_layout(True)
fig.savefig('quickstart-Mstar-shifts.png')


# code-start-match-G
sm_G.match()

# Plot chi-squared surfaces
fig = figure(figsize=(12,8))
sm_G.plot_chi_squared_surface()
# Indicate library parameters for target star.
axes = fig.axes
axes[0].axvline(G_star[0]['Teff'], color='k')
axes[1].axvline(G_star[0]['radius'], color='k')
axes[2].axvline(G_star[0]['feh'], color='k')
# code-stop-match-G

fig.set_tight_layout(True)
fig.savefig('quickstart-Gstar-chisquared-surface.png')


# code-start-lincomb-G
sm_G.lincomb()

print('Derived Parameters: ')
print('Teff: {0:.0f}, Radius: {1:.2f}, [Fe/H]: {2:.2f}'.format(
    sm_G.results['Teff'], sm_G.results['radius'], sm_G.results['feh']))
print('Library Parameters: ')
print('Teff: {0:.0f}, Radius: {1:.2f}, [Fe/H]: {2:.2f}'.format(
    G_star[0]['Teff'], G_star[0]['radius'], G_star[0]['feh']))
# code-stop-lincomb-G


# code-start-plot-lincomb-G
# Plot HR diagram
fig1 = figure(figsize=(12,10))
sm_G.plot_references(verbose=True)
# plot target onto HR diagram
axes = fig1.axes
axes[0].plot(G_star[0]['Teff'], G_star[0]['radius'], '*', ms=15, color='red', label='Target')
axes[1].plot(G_star[0]['Teff'], G_star[0]['radius'], '*', ms=15, color='red')
axes[2].plot(G_star[0]['feh'], G_star[0]['radius'], '*', ms=15, color='red')
axes[3].plot(G_star[0]['feh'], G_star[0]['radius'], '*', ms=15, color='red')
axes[0].legend(numpoints=1, fontsize='small', loc='best')


# Plot reference spectra and linear combinations
fig2 = plt.figure(figsize=(12,6))
sm_G.plot_lincomb()
# code-stop-plot-lincomb-G  
fig1.set_tight_layout(True)
fig1.savefig('quickstart-Gstar-lincomb-references.png')
fig2.set_tight_layout(True)
fig2.savefig('quickstart-Gstar-lincomb-spectra.png')



# code-start-mstar
# Perform match
sm_M.match()

# Plot chi-squared surfaces
fig1 = figure(figsize=(12,8))
sm_M.plot_chi_squared_surface()
# Indicate library parameters for target star.
axes = fig1.axes
axes[0].axvline(M_star[0]['Teff'], color='k')
axes[1].axvline(M_star[0]['radius'], color='k')
axes[2].axvline(M_star[0]['feh'], color='k')

# Perform lincomb
sm_M.lincomb()

print('Derived Parameters: ')
print('Teff: {0:.0f}, Radius: {1:.2f}, [Fe/H]: {2:.2f}'.format(
    sm_M.results['Teff'], sm_M.results['radius'], sm_M.results['feh']))
print('Library Parameters: ')
print('Teff: {0:.0f}, Radius: {1:.2f}, [Fe/H]: {2:.2f}'.format(
    M_star[0]['Teff'], M_star[0]['radius'], M_star[0]['feh']))

# Plot HR diagram
fig2 = figure(figsize=(12,10))
sm_M.plot_references(verbose=True)
# plot target onto HR diagram
axes = fig2.axes
axes[0].plot(M_star[0]['Teff'], M_star[0]['radius'], '*', ms=15, color='red', label='Target')
axes[1].plot(M_star[0]['Teff'], M_star[0]['radius'], '*', ms=15, color='red')
axes[2].plot(M_star[0]['feh'], M_star[0]['radius'], '*', ms=15, color='red')
axes[3].plot(M_star[0]['feh'], M_star[0]['radius'], '*', ms=15, color='red')
axes[0].legend(numpoints=1, fontsize='small', loc='best')

# Plot reference spectra and linear combinations
fig3 = plt.figure(figsize=(12,6))
sm_M.plot_lincomb()
# code-stop-mstar

fig1.set_tight_layout(True)
fig1.savefig('quickstart-Mstar-chisquared-surface.png')
fig2.set_tight_layout(True)
fig2.savefig('quickstart-Mstar-lincomb-references.png')
fig3.set_tight_layout(True)
fig3.savefig('quickstart-Mstar-lincomb-spectra.png')
