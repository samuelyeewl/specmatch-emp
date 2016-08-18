import ez_setup
ez_setup.use_setuptools()

import os, sys
from setuptools import setup, find_packages

on_rtd = os.environ.get('READTHEDOCS') == 'True'

if on_rtd:
    setup(
        name="SpecMatch-Emp",
        version="0.2",
        packages=find_packages(),
        include_package_data=True
    )
else:
    setup(
        name="SpecMatch-Emp",
        version="0.2",
        packages=find_packages(),
        install_requires=[
            "numpy",
            "scipy",
            "matplotlib",
            "h5py",
            "pandas",
            "lmfit",
            "astropy",
            "astroquery",
            "isochrones",
        ],
        include_package_data=True
    )

    # download library
    HOMEDIR = os.environ['HOME']
    LIBPATH = "{0}/.specmatchemp/library.h5".format(HOMEDIR)

    def reporthook(blocknum, blocksize, totalsize):
        readsofar = blocknum * blocksize
        if totalsize > 0:
            percent = readsofar * 1e2 / totalsize
            s = "\r%5.1f%% %*d / %d" % (
                percent, len(str(totalsize)), readsofar, totalsize)
            sys.stderr.write(s)
            if readsofar >= totalsize:  # near the end
                sys.stderr.write("\n")
        else:  # total size is unknown
            sys.stderr.write("read %d\n" % (readsofar,))

    # liburl = "https://zenodo.org/record/60225/files/library.h5"
    liburl = "https://www.dropbox.com/s/po0kzgjn1j9ha2v/library.h5?dl=0"
    if not os.path.exists(os.path.dirname(LIBPATH)):
        os.mkdir(os.path.dirname(LIBPATH))
    if not os.path.exists(LIBPATH):
        from six.moves import urllib
        print("Downloading library.h5")
        urllib.request.urlretrieve(liburl, LIBPATH, reporthook)

