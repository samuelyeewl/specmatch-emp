#!/usr/bin/env python
from adjust_spectra import *

def read_nso_spectrum(path, w_min=3000, w_max=13000, num_pts=64000):
    """
    Reads in the NSO spectrum, located at the specified path.
    Gives the spectrum within the given wavelength bounds.
    Also rescales the wavelengths to have a constant difference in log lambda,
    with the specified number of points

    Args:
        path:
            The path to the NSO spectrum
        w_min, w_max:
            The lower and upper wavelength bounds.
        num_pts:
            The number of sample points between those bounds

    Returns:
        s, w:
            The spectrum and wavelength scale
    """
    nso_hdu = fits.open(path)
    w = nso_hdu[0].data
    s = nso_hdu[1].data

    w_range = np.asarray([True if ww > w_min and ww < w_max else False for ww in w])
    w_clipped = w[w_range]
    s_clipped = s[w_range]

    # Rescale the reference spectrum
    log_w_min = np.log10(w_min)
    log_w_max = np.log10(w_max)
    w_log = np.logspace(log_w_min, log_w_max, num_pts, base=10.0)

    s_log = np.interp(w_log, w_clipped, s_clipped)

    return s_log, w_log


if __name__ == '__main__':
    target_path = '/Users/samuel/Dropbox/SpecMatch-Emp/spectra/iodfitsdb/rj122.761.fits'
    nso_path = '/Users/samuel/Dropbox/SpecMatch-Emp/nso/nso.fits'
    w_min = 4800
    w_max = 6500
    dw = 0.02260 # Approximate median dw
    num_pts = int((w_max-w_min)/dw)

    # Read in target and reference spectra
    s, serr, w, hdu = read_spectrum(target_path)
    s_ref, w_ref = read_nso_spectrum(nso_path, w_min, w_max, num_pts)
    serr_ref = np.zeros_like(s_ref)

    s_adj, serr_adj, w_adj = adjust_spectra(s, serr, w, s_ref, serr_ref, w_ref)
    # adjust_spectra(s, serr, w, s_ref, serr_ref, w_ref)

    plt.plot(w_adj, s_adj)
    plt.plot(w_ref, s_ref)
    plt.xlim(5385,5400)
    plt.ylim(0,1.1)
    plt.show()

    # plt.plot(w_adj[2], s_adj[2])
    # plt.plot(w_ref[2], s_ref[2])
    # plt.show()

    # # # save file
    # # outfile = os.path.splitext(target_path)[0] + '_adj.fits'
    # # hdu[0].data = adj[0]
    # # hdu[1].data = adj[1]
    # # hdu[2].data = adj[2]
    # # hdu.writeto(outfile)
    # sys.exit()
