"""
@filename library.py

Defines the library class which will be used for matching
"""

LIB_COLS = ['cps_name', 'obs', 'Teff', 'u_Teff', 'radius', 'u_radius', 'logg', 'u_logg', 'feh', 'u_feh',
           'mass', 'u_mass', 'age', 'u_age', 'vsini', 'source', 'source_name', 'library_index']

class Library():
    """Library class

    This object is a container for the library spectrum and stellar 
    parameters for the library stars.

    Args: 
        library_params (pd.DataFrame): Pandas DataFrame containing
            the parameters for the library stars. It should have
            the columns specified in LIB_COLS, although values can
            be np.nan. The library_index column specifies the position
            of the star's spectrum in the library_spectra array.

        wav (np.ndarray): Wavelength scale for the library spectra.

        library_spectra (np.ndarray): 2D array containing the library
            spectra ordered according to the library_index column.

        header (dict): (optional) Any additional metadata to store
            with the library.

        wavlim (2-element iterable): (optional) The upper and lower
            wavelength limits to be read.
    """
    target_chunk_bytes = 100e3  # Target number of bytes 

    def __init__(self, library_params, wav, library_spectra, header={}, wavlim = None):
        return

    def __init__(self, wav):
        """
        Creates empty library - use class methods to add library spectra
        """
        return

    def __str__(self):
        """
        String representation of library
        """
        return

    ## Container methods
    def __iter__(self):
        """
        Allow library to be an iterable
        """
        return self

    def __next__(self):
        """
        Next item
        """
        return

    def __len__(self):
        """
        Length
        """
        return

    def __getitem__(self, key):
        """
        Get item at specified library_index
        """
        return

    def __contains__(self, key):
        """
        Check if specified library_index is filled
        """
        return

    def insert(self, params, spectrum):
        """
        Insert spectrum and associated parameters
        """
        return

    def to_hdf(self, outfile):
        """
        Saves library as a HDF file
        """
        return

def read_hdf(infile):
    """
    Reads in a library from a HDF file
    """
    return

