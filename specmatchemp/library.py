"""
@filename library.py

Defines the library class which will be used for matching
"""

from __future__ import print_function

import datetime

import numpy as np
import pandas as pd
import h5py
import json
import matplotlib.pyplot as plt

from specmatchemp import LIBPATH
from specmatchemp import plots
from specmatchemp import spectrum
from specmatchemp.spectrum import Spectrum


class Library(object):
    """A container for a library of spectrum and corresponding stellar
    parameters.

    The Library class is a container for the library spectrum and
    stellar parameters for the library stars. The library is indexed
    using the library_index column.

    Attributes:
        library_params (pd.DataFrame): The parameter table, which can
            be queried directly as a pandas DataFrame.

        library_spectra (np.ndarray): 3D array containing the library
            spectra ordered according to the index column.

        wav (np.ndarray): Common wavelength scale for the library spectra.

        param_mask (pd.DataFrame): Boolean dataframe with same rows and
            columns as library_params, marked True for model-independent
            parameters.

        nso (Spectrum): NSO spectrum.

    Args:
        wav (np.ndarray, optional): Wavelength scale for the library spectra.

        library_params (pd.DataFrame, optional): Pandas DataFrame
            containing the parameters for the library stars. It should
            have the columns specified in LIB_COLS, although values can
            be np.nan. The library_index of each row should be the index
            of the specturm in library_spectra.

        library_spectra (np.ndarray, optional): 3D array containing the
            library spectra ordered according to the index column.
            Each entry contains the spectrum, its uncertainty, and a mask
            array.

        header (dict, optional): Any additional metadata to store with
            the library.

        wavlim (tuple, optional): The upper and lower
            wavelength limits to be read.

        param_mask (pd.DataFrame, optional): Boolean dataframe with
            same rows and columns as library_params, marked True for
            model-independent parameters.

        nso (Spectrum): NSO spectrum
    """
    #: Columns required in library
    LIB_COLS = ['lib_index', 'cps_name', 'lib_obs', 'Teff', 'u_Teff',
                'radius', 'u_radius', 'logg', 'u_logg', 'feh', 'u_feh',
                'mass', 'u_mass', 'age', 'u_age', 'vsini', 'source',
                'source_name', 'snr']
    #: Numeric star properties
    STAR_PROPS = ['Teff', 'radius', 'logg', 'feh', 'mass', 'age']

    _target_chunk_bytes = 100e3  # Target number of bytes per chunk

    def __init__(self, wav=None, library_spectra=None, library_params=None,
                 header={}, wavlim=None, param_mask=None, nso=None):
        # If no params included, create empty library
        if library_params is None:
            self.library_params = pd.DataFrame(columns=Library.LIB_COLS)
            self.wav = wav
            self.library_spectra = np.empty((0, 3, len(wav)))
            self.header = {'date_created': str(datetime.date.today())}
            self.wavlim = wavlim if wavlim is not None else (wav[0], wav[-1])
            self.nso = nso
            return

        # Check if the necessary columns have been included.
        for col in Library.LIB_COLS:
            if col not in library_params:
                raise ValueError("{0} is a required column.".format(col))

        # We can create a library with no spectra for increased i/o speed.
        if library_spectra is None:
            self.library_params = library_params
            self.wav = None
            self.library_spectra = np.empty((0, 3, len(wav)))
            self.header = {'date_created': str(datetime.date.today())}
            self.wavlim = None
            self.nso = nso
            return

        # ensure that parameter table, library spectra have same length.
        if len(library_params) != len(library_spectra):
            raise ValueError("Library params and spectra are of " +
                             "unequal length")

        # ensure that all indices in library_params can be found
        # are smaller than the number of library_spectra
        num_spec = len(library_spectra)
        for i, row in library_params.iterrows():
            if row.lib_index >= num_spec:
                raise IndexError("Index {0:d} ".format(row.lib_index) +
                                 "is out of bounds in library_spectra")

        # ensure library_spectra is of right shape
        if np.shape(library_spectra)[1] != 3 or \
                np.shape(library_spectra)[2] != len(wav):
            raise IndexError("library_spectra should have shape " +
                             "({0:d}, 3, {1:d})".format(num_spec, len(wav)))

        # set index to be equal to lib_index
        library_params.lib_index = library_params.lib_index.astype(int)
        library_params.set_index('lib_index', inplace=True, drop=False)

        # assign attributes
        self.wav = wav
        self.library_params = library_params
        self.library_spectra = library_spectra
        self.header = header
        header['date_created'] = str(datetime.date.today())
        self.wavlim = wavlim if wavlim is not None else (wav[0], wav[-1])
        self.param_mask = param_mask
        self.nso = nso

    def copy(self):
        """Return a deep copy of the library object.
        """
        wav = np.copy(self.wav) if self.wav is not None else None
        library_spectra = np.copy(self.library_spectra) \
            if self.library_spectra is not None else None
        library_params = self.library_params.copy() \
            if self.library_params is not None else None
        header = self.header.copy() if self.header is not None else None
        wavlim = self.wavlim
        param_mask = self.param_mask.copy() \
            if self.param_mask is not None else None
        nso = self.nso.copy() if self.nso is not None else None

        return Library(wav, library_spectra, library_params, header,
                       wavlim, param_mask, nso)

    def append(self, params, spectrum=None):
        """Adds spectrum and associated stellar parameters into library.

        Args:
            params (pd.Series): A row to be added to the library array. It
                should have the fields specified in Library.LIB_COLS.
            spectrum (Spectrum): Spectrum object. Spectrum should have been
                shifted and interpolated onto the same wavelength scale as
                the library.
        """
        if self.library_spectra is None and spectrum is not None:
            raise ValueError("Error: Cannot append to library with no spectra")

        # ensure that parameter row has the right columns
        for col in params.columns:
            if col not in Library.LIB_COLS:
                raise KeyError("{0} is not an allowed column".format(col))

        # ensure that the provided spectrum has the same number of elements
        # as the wavelength array
        if len(spectrum.w) != len(self.wav) \
                or not np.allclose(spectrum.w, self.wav):
            raise ValueError("Spectrum should be shifted and interpolated " +
                             "onto library wavelength scale.")

        # add new star to library
        params.lib_index = len(self.library_spectra)
        self.library_params = pd.concat((self.library_params, params),
                                        ignore_index=True)
        if spectrum is not None:
            new_spec = [spectrum.s, spectrum.serr, spectrum.mask]
            self.library_spectra = np.vstack((self.library_spectra, new_spec))
        self.library_params.set_index('lib_index', inplace=True, drop=False)

    def remove(self, index):
        """Removes the spectrum and parameters with the given index from
        the library.

        Args:
            index (int of array of ints): Indices of spectra to remove.
        """
        if not self.__contains__(index):
            raise KeyError

        if self.library_spectra is not None:
            self.library_spectra = np.delete(self.library_spectra, index,
                                             axis=0)
        self.library_params = self.library_params[self.library_params.
                                                  lib_index != index]

        # reset index
        self.library_params.lib_index = self.library_params.lib_index. \
            apply(lambda i: i - 1 if i > index else i)
        self.library_params.set_index('lib_index', inplace=True, drop=False)

    def pop(self, index):
        """Removes the spectrum and parameters with the given index from the library.

        Args:
            index (int): Index of spectrum to remove.
        Returns:
            (pd.Series, spectrum.Spectrum):
                - params (pd.Series): Parameters of star
                - spectrum (Spectrum): Spectrum
        """
        if not self.__contains__(index):
            raise KeyError

        params = self.library_params.loc[index]
        if self.library_spectra.size > 0:
            spectrum = self.get_spectrum(index)

        self.remove(index)

        if self.library_spectra.size > 0:
            return params, spectrum
        else:
            return params

    def get_index(self, searchstr):
        """Searches the library for the given search string. Checks columns
        lib_obs, cps_name, source_name in order.

        Args:
            searchstr (str): String to search for.

        Returns:
            int: Library index of the found star.
            Returns None if no object found.
        """
        pattern = '^' + searchstr + '$'
        lib_params = self.library_params
        # Search columns in order
        res = lib_params[lib_params.lib_obs.str.match(pattern)]
        if len(res) == 1:
            return res.iloc[0].lib_index
        elif len(res) > 1:
            return np.array(res.iloc[:].lib_index)

        res = lib_params[lib_params.cps_name.str.match(pattern)]
        if len(res) == 1:
            return res.iloc[0].lib_index
        elif len(res) > 1:
            return np.array(res.iloc[:].lib_index)

        res = lib_params[lib_params.source_name.str.match(pattern)]
        if len(res) == 1:
            return res.iloc[0].lib_index
        elif len(res) > 1:
            return np.array(res.iloc[:].lib_index)

        return None

    def to_hdf(self, path):
        """
        Saves library as a HDF file

        Args:
            path (str): Path to store library
        """
        self.library_params.lib_index = self.library_params.lib_index \
            .astype(int)
        self.library_params.set_index('lib_index', inplace=True, drop=False)
        self.library_params.index.rename('idx', inplace=True)

        # store spectrum
        with h5py.File(path, 'w') as f:
            # saving attributes
            for key in self.header.keys():
                f.attrs[key] = json.dumps(self.header[key])

            # saving wavlength array
            f['wav'] = self.wav

            if self.nso is not None:
                f.create_group('nso')
                self.nso.to_hdf(f['nso'])

            # convert library_params to record array
            r = self.library_params.to_records()
            dt = r.dtype.descr
            for i in range(len(dt)):
                if dt[i][1] == "|O":
                    # max string length = 100
                    dt[i] = (dt[i][0], 'S100')
            r = np.array(r, dtype=dt)
            f['params'] = r

            if self.param_mask is not None:
                f['param_mask'] = self.param_mask

            # Compute chunk size - group wavelenth regions together
            chunk_row = len(self.library_spectra)
            chunk_depth = 3
            chunk_col = int(self._target_chunk_bytes /
                            self.library_spectra[:, :, 0].nbytes)
            chunk_size = (chunk_row, chunk_depth, chunk_col)

            print("Storing library spectra with chunks of size {0}"
                  .format(chunk_size))
            f.create_dataset('library_spectra', data=self.library_spectra,
                             compression='gzip', compression_opts=1,
                             shuffle=True, chunks=chunk_size)

    def get_spectrum(self, indices):
        """Returns the spectrum at the given index.

        Args:
            indices (int or list): Library indices of spectrum.
        Returns:
            Spectrum or list of Spectrum: Spectra at given indices.
        """
        if isinstance(indices, int) or isinstance(indices, np.int_):
            indices = [indices]

        spectra = []
        for idx in indices:
            s = self.library_spectra[idx, 0]
            serr = self.library_spectra[idx, 1]
            w = self.wav
            name = self.library_params.loc[idx, 'cps_name']
            attrs = {}
            for p in Library.STAR_PROPS:
                attrs[p] = self.library_params.loc[idx, p]
            attrs['obs'] = self.library_params.loc[idx, 'lib_obs']

            spectra.append(Spectrum(w, s, serr, name=name, attrs=attrs))

        if len(spectra) == 1:
            return spectra[0]
        else:
            return spectra

    def query_params(self, querystr):
        """Query the library_params table and return the indices of matching
        stars.

        Args:
            querystr (str): Query string

        Returns:
            int or list of int:
            Index or indices of library entries matching the query string.
        """
        return np.array(self.library_params.query(querystr).index)

    def wav_cut(self, minw, maxw, deepcopy=False):
        """Returns a library with spectra within the given limits.

        Args:
            minw (float): Minimum wavelength
            maxw (float): Maximum wavelength
            deepcopy (bool, optional): True to return a library with all
                attributes copied. Default is False.

        Returns:
            Library: New library object.
        """
        idxwav, = np.where((self.wav > minw) & (self.wav < maxw))
        idxmin = idxwav[0]
        idxmax = idxwav[-1] + 1  # add 1 to include last index when slicing

        if deepcopy:
            return Library(self.wav[idxmin:idxmax],
                           self.library_spectra[:, :, idxmin:idxmax],
                           self.library_params.copy(), self.header.copy(),
                           (minw, maxw), self.param_mask.copy(),
                           self.nso.copy())
        else:
            return Library(self.wav[idxmin:idxmax],
                           self.library_spectra[:, :, idxmin:idxmax],
                           self.library_params, self.header, (minw, maxw),
                           self.param_mask, self.nso)

    def plot(self, paramx, paramy, grouped=False, marker='.',
             ptlabels=False, plt_kw={}):
        """Create a plot of the library in parameter space

        Args:
            paramx (str): Parameter to plot on the x-axis
            paramy (str): Parameter to plot on the y-axis
            grouped (bool, optional): Whether to group the stars by source
            ptlabels (str, optional): Library column to annotate points with
            plt_kw (dict, optional): Additional keyword arguments to pass to
                pyplot.plot
        """
        if grouped:
            g = self.library_params.groupby('source')

            for source in g.groups:
                cut = self.library_params.ix[g.groups[source]]
                plt.plot(cut[paramx], cut[paramy], marker, label=source)
        else:
            plt.plot(self.library_params[paramx], self.library_params[paramy],
                     marker, **plt_kw)

        if ptlabels is not False:
            self.library_params.apply(lambda x: plots.annotate_point(x[paramx],
                                      x[paramy], x[ptlabels]), axis=1)

        plots.label_axes(paramx, paramy)

    def __str__(self):
        """String representation of library
        """
        outstr = "<specmatchemp.library.Library>\n"
        outstr += "Contains {0:d} stars\n".format(len(self.library_params))
        if self.wavlim is not None:
            outstr += "Wavelength region: {0:.2f} - {1:.2f}\n" \
                .format(**self.wavlim)

        # print other attributes
        for key, val in self.header.items():
            outstr += "{0}: {1}\n".format(key, val)

        return outstr

    # Container methods
    def __iter__(self):
        """
        Allow library to be an iterable
        """
        # iterate over the spectrum table
        # iteration over np.ndarray much faster than pd.DataFrame
        self.__it_counter = 0

        return self

    def next(self):
        return self.__next__()

    def __next__(self):
        """
        Next item

        Returns:
            params (pd.Series): Stellar parameters for the next star.
            spectrum (Spectrum): Spectrum for the next star
        """
        if self.__it_counter >= len(self.library_params):
            raise StopIteration

        idx = self.__it_counter
        self.__it_counter += 1

        if self.library_spectra is None:
            return self.library_params.loc[idx]
        else:
            return self.library_params.loc[idx], self.get_spectrum(idx)

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

        if self.library_spectra.size > 0:
            return self.library_params.loc[index], self.get_spectrum(index)
        else:
            return self.library_params.loc[index]

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
    Reads in the library from a HDF file.

    Args:
        path (str, optional): path to h5 file containing library.
            If no path is given, reads in from default directory.
        wavlim (2-element iterable or str, optional):
            The upper and lower wavelength limits to be read.
            If 'none', reads in library without spectra.
            If 'all', reads in complete stored spectra.

    Returns:
        Library: Library object
    """
    if path is None:
        return read_hdf(LIBPATH, wavlim)

    with h5py.File(path, 'r') as f:
        print("Reading library from {0}".format(path))
        header = dict(f.attrs)
        for k in header.keys():
            header[k] = json.loads(header[k])

        wav = f['wav'][:]

        if 'nso' in f:
            nso = spectrum.read_hdf(f['nso'])
            nso.name = 'NSO'

        if 'param_mask' in f:
            param_mask = f['param_mask'][:]
        else:
            param_mask = None

        library_params = pd.DataFrame.from_records(f['params'][:], index='idx')
        # decode strings
        for (col_name, dt) in library_params.dtypes.iteritems():
            if dt == 'object':
                library_params[col_name] = library_params[col_name]\
                    .str.decode('utf-8')

        if wavlim == 'all':
            library_spectra = f['library_spectra'][:]
            wavlim = (wav[0], wav[-1])
        elif wavlim == 'none':
            library_spectra = None
        else:
            idxwav, = np.where((wav > wavlim[0]) & (wav < wavlim[1]))
            idxmin = idxwav[0]
            idxmax = idxwav[-1] + 1  # add 1 to include last index when slicing
            library_spectra = f['library_spectra'][:, :, idxmin:idxmax]
            wav = wav[idxmin:idxmax]

    lib = Library(wav, library_spectra, library_params, header=header,
                  wavlim=wavlim, param_mask=param_mask, nso=nso)
    return lib
