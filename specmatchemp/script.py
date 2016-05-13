#!/usr/bin/env python

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import shift_spectra
import specmatch_io

lib = pd.read_csv('../starswithspectra.csv',index_col=0)
lib = lib.convert_objects(convert_numeric=True)

nso_path = '/Users/samuel/Dropbox/SpecMatch-Emp/nso/nso_std.fits'
s_nso, w_nso, serr_nso, h_nso = specmatch_io.read_standard_spectrum(nso_path)

# # 4500 < Teff < 6500
# lib_restr = lib.query('4500 < Teff < 6500')
# ref_path = nso_path

# # 3700 < Teff < 4500
# lib_restr = lib.query('3700 < Teff <= 4500')
# ref_path = '../lib/rj55.906_adj.fits'

# # Teff <= 3700
# lib_restr = lib.query('Teff <= 3700')
# ref_path = '../lib/rj59.1926_adj.fits'

# Teff >= 6500
lib_restr = lib.query('Teff >= 6300')
ref_path = '../lib/rj187.479_adj.fits'

for index, row in lib_restr.iterrows():
    spec_path = '/Users/samuel/Dropbox/SpecMatch-Emp/spectra/iodfitsdb/'+row['obs']+'.fits'
    out_path = '../lib/'+row['obs']+'_adj.fits'
    img_path = '../lib/Images/'+str(row['Teff'])+'_'+row['obs']+'.png'
    try:
        s, w, serr = shift_spectra.main(spec_path, 'hires', ref_path, out_path)
    except Exception as e:
        print(e)
        continue

    plt.clf()
    plt.plot(w_nso, s_nso)
    plt.plot(w, s)
    plt.xlim(5180, 5190)
    # plt.show()
    plt.savefig(img_path)