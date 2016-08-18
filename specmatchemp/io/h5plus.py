"""
@filename h5plus.py

Module to augment h5py utility.
"""
import h5py


def read_dict(f, recursive=True):
    """Read out a dictionary from a h5py file object.

    Args:
        f (h5py.File): File object to read from
        recursive (bool, optional): Whether to search subgroups in f.
    """
    d = {}
    def _add_to_dict(k):
        if isinstance(f[k], h5py.Dataset):
            d[k] = f[k].value

    if recursive:
        f.visit(_add_to_dict)
    else:
        for k in f:
            _add_to_dict(k)

    return d
