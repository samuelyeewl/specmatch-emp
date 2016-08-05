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
        "pandas",
        "json",
        "lmfit",
        "isochrones",
        "h5py",
    ],
    include_package_data=True
)