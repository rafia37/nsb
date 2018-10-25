#!/usr/bin/env python3
import sys, os
import argparse
from astropy.io import fits
from astropy.io.fits import getheader, getval, getdata
import numpy as np
import pdb

def _get_biases(files):
    """
    Grabs the biases by looking at the IMAGETYP header in the fits
    files.
    """
    bias_data = []
    bias_files = []

    for filename in files:
        try:
            type = getval(filename, 'IMAGETYP')
        except OSError:
            continue

        if type == 'Bias Frame':
            bias_data.append(getdata(filename))
            #temp_data = getdata(filename)
            #print(np.mean(temp_data))
            #print(np.std(temp_data))
            #print(np.median(temp_data))
            bias_files.append(filename)


    with open('bias.txt', 'w') as file:
        for bias in bias_files:
            line = '{}\n'.format(bias)
            file.write(line)

    return bias_data

def _get_flats(files):
    """
    Gets the flats and puts them into dictionaries depending on the
    filter type.
    """
    flats = {}
    for filename in files:
        try:
            file_type = getval(filename, 'IMAGETYP')
        except OSError:
            continue

        if file_type == 'Flat Field':
            filter_used = getval(filename, 'FILTER')

            if filter_used not in flats:
                flats[filter_used] = []

            flats[filter_used].append((filename, getdata(filename)))

    return flats

def _make_master_bias(bias_data):
    """
    Creates the master bias using bias data.
    """
    print('\n#--------- Making master bias')

    bias_amount = len(bias_data)
    print('Using {} bias images for master bias'.format(bias_amount))

    master_bias = np.median(bias_data, axis=0)
    # get bias stats and inject into header
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

    return master_bias

def _make_master_flat(flats, master_bias, only_flat=False):
    """
    Makes the master flat. It normalizes the flat data then takes the median
    of the normalized flat to remove stars. Writes a file named
    "masterflat.fits".
    """
    for key, values in flats.items():
        print('\n#--------- Creating {} master flat'.format(key))
        normalized_flats = []

        print('Subtracting bias and normalizing flat data')
        for flat in values:
            filename = flat[0]
            data = flat[1]

            if not only_flat:
                data = data - master_bias
            normalized_flats.append(
                np.divide(data, np.median(data[500:2556,500:2556]))
            )
        print('Using {} flats'.format(len(normalized_flats)))
        print('Taking the median of the normalized flats')
        master_flat_data = np.median(normalized_flats, axis=0)
        master_flat_header = fits.Header()
        master_flat_header['IMAGETYP'] = 'Master Flat'
        master_flat_header['FILTER'] = key
        new_fits_file = fits.PrimaryHDU(
            master_flat_data,
            header=master_flat_header
        )

        master_flat_file = 'masterflat-{}.fits'.format(key)
        print('Saving master flat ----> {}'.format(master_flat_file))
        new_fits_file.writeto(
            master_flat_file,
            overwrite=True
        )

        print('\nDone...')

    return None

def _parse_arguments():
    """
    Parses the command line arguments.
    """
    parser = argparse.ArgumentParser()
    parser.add_argument(
        'files',
        nargs='+',
        help='files to reduce'
    )
    parser.add_argument(
        '-f',
        '--flat',
        help='do flats only',
        action='store_true'
    )
    args = parser.parse_args()
    files = args.files
    only_flat = args.flat

    return files, only_flat

def main():
    """
    Main loop that creates the master image reduction files.
    """
    print('\n########## nsb_masters.py')
    print('Creates master bias and master flats for each filter')

    files, only_flat = _parse_arguments()
    master_bias = None

    if not only_flat:
        print('doing bias stuff')
        bias_data = _get_biases(files)
        master_bias = _make_master_bias(bias_data)


    flats = _get_flats(files)
    _make_master_flat(flats, master_bias, only_flat)

    print('\nFINISHED\n')

    sys.exit(0)

if __name__ == '__main__':
    main()

