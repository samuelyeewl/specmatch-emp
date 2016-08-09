"""
@filename library.py

Defines the library class which will be used for matching
"""

from __future__ import print_function

import os, sys
import datetime

import numpy as np
import pandas as pd
import h5py

LIB_COLS = ['lib_index','cps_name', 'lib_obs', 'Teff', 'u_Teff', 'radius', 'u_radius', 
            'logg', 'u_logg', 'feh', 'u_feh', 'mass', 'u_mass', 'age', 'u_age', 
            'vsini', 'source', 'source_name', 'snr']
STAR_PROPS = ['Teff', 'radius', 'logg', 'feh', 'mass', 'age']
FLOAT_TOL = 1e-3
HOMEDIR = os.environ['HOME']
LIBPATH = "{0}/.specmatchemp/library.h5".format(HOMEDIR)

def reporthook(blocknum, blocksize, totalsize):
    readsofar = blocknum * blocksize
    if totalsize > 0:
        percent = readsofar * 1e2 / totalsize
        s = "\r%5.1f%% %*d / %d" % (
            percent, len(str(totalsize)), readsofar, totalsize)
        sys.stderr.write(s)
        if readsofar >= totalsize: # near the end
            sys.stderr.write("\n")
    else: # total size is unknown
        sys.stderr.write("read %d\n" % (readsofar,))

liburl = "https://zenodo.org/record/59743/files/library.h5"
if not os.path.exists(os.path.dirname(LIBPATH)):
    os.mkdir(os.path.dirname(LIBPATH))
if not os.path.exists(LIBPATH):
    from six.moves import urllib
    print("Downloading library.h5")
    urllib.request.urlretrieve(liburl, LIBPATH, reporthook)

class Library():
    """A container for a library of spectrum and corresponding stellar parameters.

    The Library class is a container for the library spectrum and 
    stellar parameters for the library stars. The library is indexed 
    using the library_index column.

    Attributes:
        library_params (pd.DataFrame): The parameter table, which can
            be queried directly.
        library_spectra (np.ndarray): 3D array containing the library
            spectra ordered according to the index column.

    Args: 
        library_params (pd.DataFrame): Pandas DataFrame containing
            the parameters for the library stars. It should have
            the columns specified in LIB_COLS, although values can
            be np.nan. The library_index of each row should
            be the index of the specturm in library_spectra.

        wav (np.ndarray): Wavelength scale for the library spectra.

        library_spectra (np.ndarray): 3D array containing the library
            spectra ordered according to the index column.
            Each entry contains the spectrum and its uncertainty.

        header (dict): (optional) Any additional metadata to store
            with the library.

        wavlim (2-element iterable): (optional) The upper and lower
            wavelength limits to be read.
    """
    target_chunk_bytes = 100e3  # Target number of bytes per chunk

    def __init__(self, wav=None, library_spectra=None, library_params=None, header={}, wavlim=None, param_mask=None):
        """
        Creates a fully-formed library from a given set of spectra.
        """
        # If no params included, create empty library
        if library_params is None:
            self.library_params = pd.DataFrame(columns=LIB_COLS)
            self.wav = wav
            self.library_spectra = np.empty((0, 2, len(wav)))
            self.header = {'date_created': str(datetime.date.today())}
            self.wavlim = wavlim
            return

        # otherwise we need to include the provided tables
        # ensure that parameter table has the right columns
        for col in library_params:
            assert col in LIB_COLS, \
                "{0} is not an allowed column".format(col)

        # If no spectra included but params are included, we have a spectrumless library
        if library_spectra is None:
            self.library_params = library_params
            self.wav = wav
            self.library_spectra = np.empty((0, 2, len(wav)))
            self.header = {'date_created': str(datetime.date.today())}
            self.wavlim = wavlim
            return

        # ensure that parameter table, library spectra have same length.
        num_spec = len(library_spectra)
        assert len(library_params) == num_spec,    \
            "Error: Length of parameter table and library spectra are not equal."
        # ensure that all indices in library_params can be found in library_spectra
        for i, row in library_params.iterrows():
            assert row.lib_index < num_spec,     \
            "Error: Index {0:d} is out of bounds in library_spectra".format(i)

        # ensure library_spectra is of right shape
        assert np.shape(library_spectra)[1] == 2 and np.shape(library_spectra)[2] == len(wav), \
            "Error: library_spectra should have shape ({0:d}, 2, {1:d}".format(num_spec, len(wav))

        # set index to be equal to lib_index
        library_params.lib_index = library_params.lib_index.astype(int)
        self.library_params = library_params
        self.library_params.set_index('lib_index', inplace=True, drop=False)

        self.wav = wav
        self.library_spectra = library_spectra
        self.header = header
        header['date_created'] = str(datetime.date.today())
        self.wavlim = wavlim
        self.param_mask = param_mask

    def append(self, params, spectrum, u_spectrum):
        """Adds spectrum and associated stellar parameters into library.

        Args:
            params (pd.Series): A row to be added to the library array. It
                should have the fields specified in LIB_COLS.
            spectrum (np.ndarray): Array containing spectrum to be added.
                The spectrum should have been shifted and interpolated onto
                the same wavelength scale as the library.
            u_spectrum (np.ndarray): Array containing uncertainty in spectrum.
        """
        assert self.library_spectra is not None, \
            "Error: Cannot append to library with no spectra"
        # ensure that parameter table, library spectra have same length.
        assert len(self.library_params) == len(self.library_spectra),    \
            "Error: Length of parameter table and library spectra are not equal."

        # ensure that parameter row has the right columns
        for col in params.columns:
            assert col in LIB_COLS, \
                "{0} is not an allowed column".format(col)

        # ensure that the provided spectrum has the same number of elements
        # as the wavelength array
        assert len(spectrum) == len(self.wav), \
            "Error: spectrum is not the same length as library wavelength array"

        # add new star to library
        params.lib_index = len(self.library_spectra)
        self.library_params = pd.concat((self.library_params, params), ignore_index=True)
        self.library_spectra = np.vstack((self.library_spectra, [[spectrum, u_spectrum]]))
        self.library_params.set_index('lib_index', inplace=True, drop=False)

    def remove(self, index):
        """Removes the spectrum and parameters with the given index from the library.

        Args:
            index (int): Index of spectrum to remove.
        """
        if not self.__contains__(index):
            raise KeyError

        if self.library_spectra is not None:
            self.library_spectra = np.delete(self.library_spectra, [index], axis=0)
        self.library_params = self.library_params[self.library_params.lib_index != index]

        # reset index
        self.library_params.lib_index = self.library_params.lib_index.apply(\
            lambda i: i-1 if i > index else i)
        self.library_params.set_index('lib_index', inplace=True, drop=False)


    def get_index(self, searchstr):
        """Searches the library for the given search string. Checks columns
        lib_obs, cps_name, source_name in order.

        Args:
            searchstr (str): String to search for
        Returns:
            lib_index (int): Library index of the found star. Returns None if no
                object found. 
        """
        pattern='^'+searchstr+'$'
        res = self.library_params[self.library_params.lib_obs.str.match(pattern)]
        if len(res)==1:
            return res.iloc[0].lib_index
        elif len(res)>1:
            return np.array(res.iloc[:].lib_index)
        res = self.library_params[self.library_params.cps_name.str.match(pattern)]
        if len(res)==1:
            return res.iloc[0].lib_index
        elif len(res)>1:
            return np.array(res.iloc[:].lib_index)
        res = self.library_params[self.library_params.source_name.str.match(pattern)]
        if len(res)==1:
            return res.iloc[0].lib_index
        elif len(res)>1:
            return np.array(res.iloc[:].lib_index)
        return None


    def to_hdf(self, path):
        """
        Saves library as a HDF file

        Args:
            path (str): Path to store library
        """
        self.library_params.lib_index = self.library_params.lib_index.astype(int)
        self.library_params.set_index('lib_index', inplace=True, drop=False)
        self.library_params.index.rename('idx', inplace=True)

        # store spectrum
        with h5py.File(path, 'w') as f:
            for key in self.header.keys():
                f.attrs[key] = self.header[key]
            f['wav'] = self.wav

            # convert library_params to record array
            r = self.library_params.to_records()
            dt = r.dtype.descr
            for i in range(len(dt)):
                if dt[i][1] == "|O":
                    # max string length = 100
                    dt[i] = (dt[i][0], 'S100')
            r = np.array(r, dtype=dt)
            f['params'] = r

            f['param_mask'] = self.param_mask

            # Compute chunk size - group wavelenth regions together
            chunk_row = len(self.library_spectra)
            chunk_depth = 2
            chunk_col = int(self.target_chunk_bytes / self.library_spectra[:,:,0].nbytes)
            chunk_size = (chunk_row, chunk_depth, chunk_col)

            print("Storing model spectra with chunks of size {0}".format(chunk_size))
            dset = f.create_dataset('library_spectra', data=self.library_spectra,
                compression='gzip', compression_opts=1, shuffle=True, chunks=chunk_size)

    def __str__(self):
        """
        String representation of library
        """
        outstr = "<specmatchemp.library.Library>\n"
        for key, val in self.header.items():
            outstr += "{0}: {1}\n".format(key, val)

        return outstr

    ## Container methods
    def __iter__(self):
        """
        Allow library to be an iterable
        """
        # iterate over the spectrum table
        # iteration over np.ndarray much faster than pd.DataFrame
        self.__it_counter = 0

        return self

    def __next__(self):
        """
        Next item

        Returns:
            params (pd.Series): Stellar parameters for the next star.
            spectrum (np.ndarray): Spectrum for the next star
        """
        if self.__it_counter >= len(self.library_params):
            raise StopIteration

        idx = self.__it_counter
        self.__it_counter += 1

        if self.library_spectra is None:
            return self.library_params.loc[idx]
        else:
            return self.library_params.loc[idx], self.library_spectra[idx]

    def __len__(self):
        """
        Number of spectra in library.

        Returns:
            Number of spectra stored in library.
        """
        return len(self.library_spectra)

    def __getitem__(self, index):
        """
        Get item at specified library_index

        Args:
            index (int): Library index of desired spectrum.
        """
        # Check if library_index specified is in the container
        if not self.__contains__(index):
            raise KeyError

        if self.library_spectra is None:
            return self.library_params.loc[index]
        else:
            return self.library_params.loc[index], self.library_spectra[index]

    def __delitem__(self, index):
        """
        Deletes the item at the given index. Alias for Library.remove

        Args:
            index (int): Library index of spectrum to remove.
        """
        self.remove(index)

    def __contains__(self, index):
        """
        Check if specified library_index is filled

        Args:
            index (int): Library index to check
        """
        return index in self.library_params.lib_index

def read_hdf(path=None, wavlim='all'):
    """
    Reads in a library from a HDF file

    Args:
        paramfile (str): path to h5 file containing star parameters.
        specfile (str): path to h5 file containing spectra.
        wavlim (2-element iterable): (optional) The upper and lower wavelength
            limits to be read. If 'none', reads in library without spectra.

    Returns:
        lib (library.Library) object
    """
    if path is None:
        return read_hdf(LIBPATH, wavlim)

    with h5py.File(path, 'r') as f:
        print("Reading library from {0}".format(path))
        header = dict(f.attrs)
        wav = f['wav'][:]
        param_mask = f['param_mask'][:]
        library_params = pd.DataFrame.from_records(f['params'][:], index='idx')
        # decode strings
        for (col_name, dt) in library_params.dtypes.iteritems():
            if dt == 'object':
                library_params[col_name] = library_params[col_name].str.decode('utf-8')

        if wavlim == 'all':
            library_spectra = f['library_spectra'][:]
        elif wavlim == 'none':
            library_spectra = None
        else:
            idxwav, = np.where( (wav > wavlim[0]) & (wav < wavlim[1]))
            idxmin = idxwav[0]
            idxmax = idxwav[-1] + 1 # add 1 to include last index when slicing
            library_spectra = f['library_spectra'][:,:,idxmin:idxmax]
            wav = wav[idxmin:idxmax]

    lib = Library(wav, library_spectra, library_params, header=header, wavlim=wavlim, param_mask=param_mask)
    return lib

