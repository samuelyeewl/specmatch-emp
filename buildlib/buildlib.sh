#!/bin/bash

# read the catalogs
python buildlib/read_catalogs.py catalogs/ "~/Dropbox/SpecMatch-Emp/spectra/cps_templates.csv" lib/

# fill in parameters with isochrone models
python buildlib/get_isochrones.py lib/libstars.csv lib/isochrone_models/ -m lib/libstars_mask.csv

# build references for bootstrapping
python buildlib/shift_spectrum.py CK01241 "/users/samuel/Dropbox/SpecMatch-Emp/spectra/iodfitsdb/rj125.89.fits" "/Users/samuel/Dropbox/SpecMatch-Emp/nso/nso_std.fits" ./results/
python buildlib/shift_spectrum.py 222368 "/users/samuel/Dropbox/SpecMatch-Emp/spectra/iodfitsdb/rj26.531.fits" "/Users/samuel/Dropbox/SpecMatch-Emp/nso/nso_std.fits" ./results/
# extract_spectrum.py
python buildlib/shift_spectrum.py 216899 "/users/samuel/Dropbox/SpecMatch-Emp/spectra/iodfitsdb/rj59.1926.fits" "/Users/samuel/Dropbox/SpecMatch-Emp/nso/nso_std.fits" ./results/
# extract_spectrum.py

# generate script
python buildlib/generate_shift_script.py

# run script
parallel --slf $PBS_NODEFILE -a shift_script.txt

# combine library
