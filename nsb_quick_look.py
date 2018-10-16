#!/usr/bin/env python3

from astropy.io.fits import getdata, getval
import numpy as np
import os, sys
from astropy.stats import sigma_clipped_stats
import argparse

def _parse_arguments():
    """
    Parser to get arguments from the command line.
    """
    parser = argparse.ArgumentParser()

    parser.add_argument(
        '-f',
        dest='files',
        help='files to work on',
        nargs='+'
    )
    parser.add_argument(
        '-b',
        dest='master_bias',
        help='master bias to use'
    )
    parser.add_argument(
        '-c',
        nargs='+',
        dest='bias_files',
        help='bias files to create master bias',
        default=None
    )
    parser.add_argument(
        '-r',
        dest='reduced_already',
        help='flag if data is already bias subtracted',
        action='store_true'
    )

    args = parser.parse_args()

    files = args.files
    master_bias = args.master_bias
    bias_files = args.bias_files
    reduced_already = args.reduced_already

    return files, master_bias, bias_files, reduced_already

def _make_master_bias(bias_files):
    """
    Make master bias.
    """
    print('\n#--------- Making master bias')

    bias_amount = len(bias_data)
    print('Using {} bias images for master bias'.format(bias_amount))

    master_bias = np.median(bias_data, axis=0)

    std = np.std(master_bias)
    mean = np.mean(master_bias)
    median = np.median(master_bias)
    variance = np.var(master_bias)

    biashead = fits.Header()
    biashead['IMAGETYP'] = 'Master Bias'
    biashead['MEAN'] = (mean, 'mean')
    biashead['MEDIAN'] = (median, 'median')
    biashead['STD'] = (std, 'standard dev')
    biashead['VAR'] = (variance, 'variance')
    biasfits = fits.PrimaryHDU(master_bias, header=biashead)

    print('Saving master bias ----> masterbias.fits')
    biasfits.writeto('./masterbias.fits', overwrite=True)

    print('\nDone...')

    return None

def _subtract_master_bias(data, master_bias):
    """
    Subtract master bias from the data.
    """
    reduced_data = np.subtract(data, master_bias)

    return reduced_data

def _get_background(data, sigma_clipping=True):
    """
    Measure the background stats for the reduced data.
    """
    if sigma_clipping:
        mean, median, std = sigma_clipped_stats(
            data,
            sigma=3.0,
            iters=5
        )

    else:
        mean = np.mean(reduced_data)
        median = np.median(reduced_data)
        std = np.std(reduced_data)

    return mean, median, std

def _print_results(file, mean, median, std):
    """
    Print the results of the analysis to the termainl.
    """

    exposure_time = getval(file, 'EXPTIME')

    print('\nFile:\t\t{}'.format(file))
    print('\nExposure:\t{}'.format(exposure_time))
    print('Mean:\t\t{:.2f}'.format(mean))
    print('Median:\t\t{:.2f}'.format(median))
    print('Std:\t\t{:.2f}'.format(std))

    return None

def main():
    """
    Program that checks to make sure the counts above the noise is good
    enough for nsb measurements.Has the option to create a master bias
    if needed as well.
    """

    files, master_bias, bias_files, reduced_already = _parse_arguments()

    if bias_files:
        _make_master_bias(bias_files)

    else:
        for file in files:
            data = getdata(file).astype(np.float)
            if not reduced_already:
                master_bias_data = getdata(master_bias).astype(np.float)
                data = _subtract_master_bias(data, master_bias_data)
            mean_bkg, median_bkg, std_bkg = _get_background(data)
            _print_results(file, mean_bkg, median_bkg, std_bkg)

    print('\nFINISHED\n')
    sys.exit(0)

if __name__ == '__main__':
    main()
