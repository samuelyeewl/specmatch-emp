"""
@filename detrend.py

Class to carry out detrending of parameters
"""

import os
import csv
import matplotlib.pyplot as plt
from specmatchemp import SPECMATCHDIR


class Detrend(object):
    """Class to carry out detrending of parameters.

    Reads in a CSV file containing the trends, provided as tuples of the form
    (param, uncal_1, cal_1, uncal_2, cal_2), where uncal_i are the
    uncalibrated values, and cal_i are the desired calibrated values.
    uncal_1 and uncal_2 provide the range for which the particular trendline
    can be applied.

    Suppose a given parameter is given as u, such that \
    :math:`uncal_1 \leq u \lt uncal_2.`
    Then the detrended parameter is given by

    .. math::

        c = (u - uncal_1) \cdot \frac{cal_2 - cal_1}{uncal_2 - uncal_1} + cal_1

    In this way, we can fit out a piecewise linear relationship by providing
    multiple entries for the same parameters. If overlapping regions are
    provided, the region with the smallest leftmost boundary is used. Outside
    of the provided regions, the class does nothing.

    Args:
        filename (str, optional): Path to csv file containing detrend
            parameters. If not provided, looks in default working location.
    """
    def __init__(self, filename=""):
        if len(filename) == 0:
            filename = os.path.join(SPECMATCHDIR, 'detrend.csv')

        self._detrendtable = {}
        with open(filename, mode='r') as csvfile:
            has_header = csv.Sniffer().has_header(csvfile.read(1024))
            # Rewind
            csvfile.seek(0)
            # Create reader object
            reader = csv.reader(csvfile)
            # Skip header row
            if has_header:
                next(reader)

            for row in reader:
                param = row[0]
                if param in self._detrendtable:
                    self._detrendtable[param].append((float(row[1]),
                        float(row[2]), float(row[3]), float(row[4])))
                else:
                    self._detrendtable[param] = [(float(row[1]),
                        float(row[2]), float(row[3]), float(row[4]))]

        for p in self._detrendtable:
            self._detrendtable[param].sort()

    def detrend(self, value, param):
        """Calculate the detrended parameter.

        Args:
            param (str): Name of parameter to detrend
            value (float): Value to detrend.
        """
        # Do nothing if no trend provided for the given region.
        if param not in self._detrendtable:
            return value
        else:
            # find appropriate interval
            interval = None
            for row in self._detrendtable[param]:
                if value >= row[0] and value < row[2]:
                    interval = row
                    break
            # no interval found, do nothing
            if interval is None:
                return value

            uncal_1 = row[0]
            cal_1 = row[1]
            uncal_2 = row[2]
            cal_2 = row[3]

            # detrend in radius is for delta r / r
            if param == 'radius':
                dr_r_1 = (uncal_1 - cal_1) / cal_1
                dr_r_2 = (uncal_2 - cal_2) / cal_2

                dr_r_target = ((value - uncal_1) *
                               ((dr_r_2 - dr_r_1) / (uncal_2 - uncal_1)) +
                               dr_r_1)

                return value - dr_r_target * value
            else:
                return ((value - uncal_1) *
                        ((cal_2 - cal_1) / (uncal_2 - uncal_1)) +
                        cal_1)
            # return ((value - uncal_1) *
            #         ((cal_2 - cal_1) / (uncal_2 - uncal_1)) +
            #         cal_1)

    def plot(self, param):
        """Plot the trendline(s) for the given parameter.

        Args:
            param (str): Name of parameter to plot
        """
        if param in self._detrendtable:
            for row in self._detrendtable[param]:
                if param == 'radius':
                    plt.plot([row[1], row[3]],
                        [(row[0] - row[1])/row[1], (row[2] - row[3])/row[3]],
                        'r-', linewidth=1.0)
                else:
                    plt.plot([row[1], row[3]],
                        [row[0] - row[1], row[2] - row[3]], 'r-', linewidth=1)
