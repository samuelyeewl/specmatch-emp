# SpecMatch-Emp

SpecMatch-Emp is a tool to extract the fundamental properties (effective temperature, radius, and metallicity) by comparing a target star's spectrum to a library of spectra from stars with known properties. 

## Installation

Run `python setup.py install` in the root directory. This downloads all the necessary data files and installs the command-line tools.

## Documentation

Further documentation and a quickstart walkthrough can be found at [http://specmatch-emp.readthedocs.io/](http://specmatch-emp.readthedocs.io/)

## Technical Details and Attribution

The spectral library and matching algorithm are described in [Yee et al. (2017)](http://adsabs.harvard.edu/cgi-bin/bib_query?arXiv:1701.00922). If you use `SpecMatch-Emp` in your research, please cite this paper.


## Note (added 2021-08-23)

The library star GL570B has been removed from the library -- the library
properties were erroneously for GL570A, while the spectrum was likely that of
the binary pair GL570BC. This star was _not_ included in the error analysis
described in the paper, but unfortunately crept back into the uploaded library.

To remove the star from an existing library, simply run
```
lib = specmatchemp.library.read_hdf()
idx = lib.get_index('GL570B')       # Probably 395
lib.remove(idx)
lib.to_hdf(specmatchemp.library.LIBPATH)
```
or download the library again from the dropbox [link](https://www.dropbox.com/s/po0kzgjn1j9ha2v/library.h5?dl=0).

Thanks to Eric Mamajek for catching this mistake.

