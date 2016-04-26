#!/usr/bin/env python
from adjust_spectra import *


if __name__ == '__main__':
    target_path = '/Users/samuel/Dropbox/SpecMatch-Emp/spectra/iodfitsdb/rj122.761.fits'
    ref_path = '/Users/samuel/Dropbox/SpecMatch-Emp/spectra/iodfitsdb/rj76.283.fits'
    # Read in target and reference spectra
    s, serr, w, hdu = open_spectrum(target_path)
    s_ref, serr_ref, w_ref, hdu_ref = open_spectrum(ref_path)

    s_adj, serr_adj, w_adj = adjust_spectra(s, serr, w, s_ref, serr_ref, w_ref)

    plt.plot(w_adj[2], s_adj[2])
    plt.plot(w_ref[2], s_ref[2])
    plt.show()

    # # save file
    # outfile = os.path.splitext(target_path)[0] + '_adj.fits'
    # hdu[0].data = adj[0]
    # hdu[1].data = adj[1]
    # hdu[2].data = adj[2]
    # hdu.writeto(outfile)
    sys.exit()
