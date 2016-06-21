"""
@filename specmatchemp/io/spectrum.py

Defines a Spectrum class for file input and output
"""

from astropy.io import fits
import numpy as np
FITS_COLUMNS = [
    ('wav', 'D', 'Wavelength', 'Angstroms'),
    ('flux', 'D', 'Normalized flux', 'relative intensity'),
    ('uflux', 'D', 'Flux uncertainty', 'relative intensity')
]
HEADER_REQUIRED_KEYS = ['name','obs']

class Spectrum(np.recarray):
    """Spectrum
    
    A light superclass on top of numpy record array that stores header
    information and and read and write to fits objects.

    Args:
        wav (array): wavelengths corresponding to each pixel in the
            flux array
        flux (array): continuum-normalized flux as a function of
            rest wavelength
        uflux (array): relative flux uncertainty
        header (dict): dictionary containing metadata associated with the
            observed spectrum. Similar to a header from a fits file.
            Required keys: object, observation
    """

    def __new__(cls, wav, flux, uflux, header):
        # Input array is an already formed ndarray instance. We first
        # cast to be our class type.
        obj = np.rec.fromarrays([wav, flux, uflux], names='wav,flux,uflux')
        for key in HEADER_REQUIRED_KEYS:
            headerkeys = [k.lower() for k in header.keys()]
            haskey = headerkeys.count(key.lower()) > 0
            assert haskey, "{} is required key".format(key)
        obj.header = header
        obj = obj.view(cls)
        return obj

    def __array_finalize__(self, obj):
        # see InfoArray.__array_finalize__ for comments
        self.header = getattr(obj, 'header', None)

    def __repr__(self):
        wavlim = [np.nanmin(self.wav),np.nanmax(self.wav)]
        out = "<Spectrum: {}, {}, ({:.1f}-{:.1f} Ang)>".format(
            self.header['name'],self.header['obs'],*wavlim)
        return  out

    def to_fits(self, outfile, clobber=True):
        """Save to FITS

        Save a Spectrum object as a mutli-extension fits file.

        Args:
            outfile (string): name of output file name
            clobber (bool): if true, will overwrite existing file
            
        """
        columns = []
        for i,col in enumerate(FITS_COLUMNS):
            colinfo = FITS_COLUMNS[i]
            coldata = getattr(self,colinfo[0])
            fitscol = fits.Column(
                array=coldata, format=colinfo[1], name=colinfo[0], 
                unit=colinfo[3]
            )
            columns.append(fitscol)

        table_hdu = fits.BinTableHDU.from_columns(columns)
        fitsheader = fits.Header()
        fitsheader.update(self.header)
        primary_hdu = fits.PrimaryHDU(header=fitsheader)
        hdu_list = fits.HDUList([primary_hdu, table_hdu])
        hdu_list.writeto(outfile, clobber=clobber)

def read_fits(filename):
    """Read spectrum from fits file

    Read in a spectrum as saved by the Spectrum.to_fits method into
    a Spectrum object

    Args:
        filename (string): path to fits file
        
    Returns:
        Spectrum object
        
    """
    hdu = fits.open(filename)
    header = hdu[0].header
    table = hdu[1].data
    record_names = table.dtype.names
    required_cols = [k[0] for k in FITS_COLUMNS]
    for k in required_cols:
        assert k in record_names, "Column {0} not found. {1} are all \
        requried columns in the fits table.".format(k, required_cols)
    spec = Spectrum(table['wav'], table['flux'], table['uflux'], header)
    
    return spec