import ez_setup
ez_setup.use_setuptools()

import os
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