#!/bin/bash

# read the catalogs
python buildlib/read_catalogs.py catalogs/ "~/Dropbox/SpecMatch-Emp/spectra/cps_templates.csv" lib/

# fill in parameters with isochrone models
python buildlib/get_isochrones.py lib/libstars.csv lib/isochrone_models/ -m lib/libstars_mask.csv

