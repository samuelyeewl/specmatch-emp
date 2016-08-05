import ez_setup
ez_setup.use_setuptools()

from setuptools import setup, find_packages

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