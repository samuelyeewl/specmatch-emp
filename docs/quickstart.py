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
plot(lib.library_params.Teff, lib.library_params.radius,'.')
smplot.label_axes('Teff','radius')
# code-stop-library
fig.savefig('quickstart-library.png')




# code-start-library-labeled
fig = figure()
g = lib.library_params.groupby('source')
colors = ['Red','Orange','LimeGreen','Cyan','RoyalBlue','Magenta']
i = 0
for source, idx in g.groups.iteritems():
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
plot(lib.library_params.Teff, lib.library_params.radius,'.')
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
