#!/usr/env/python
import sys, os
import argparse
from astropy.io import fits
from astropy.io.fits import getheader, getval, getdata
import numpy as np
import shutil
from tqdm import tqdm


def _get_biases(files):

    if not os.path.exists('./bias'):
        os.mkdir('./bias')

    bias_data = []
    bias_files = []
    for filename in files:
        type = getval(filename, 'IMAGETYP')

        if type == 'Bias Frame':
            bias_data.append(getdata(filename))
            bias_files.append(filename)


    with open('bias.txt', 'w') as file:
        for bias in bias_files:
            line = '{}\n'.format(bias)
            file.write(line)

    return bias_data

# show the telescope side
def _get_pierside_info(files):
    """
    Get the pierside info. This will be moved to the setup section when bias
    and darks are introduced. Flats and biases are seperated.
    """
    # create lists
    west_pierside = []
    east_pierside = []
    flat_data = []
    flat_files = []
    bias_data = []
    bias_files = []

    print('\n### Checking images for rotation and if can be used for flats')
    for filename in tqdm(files):
        pierside = getval(filename, 'PIERSIDE')
        altitude = getval(filename, 'ALTITUDE')


        if pierside == 'WEST' and altitude >= 75:
            data = np.rot90(getdata(filename), 2)
            flat_data.append(data)
            flat_files.append(filename)
            west_pierside.append(filename)

        elif pierside =='WEST':
            west_pierside.append(filename)

        elif pierside == 'EAST' and altitude >= 75:
            data = getdata(filename)
            flat_files.append(filename)
            flat_data.append(data)
            east_pierside.append(filename)

        elif pierside == 'EAST':
            east_pierside.append(filename)

    with open('west_pierside.txt', 'w') as file:
        for pierside in west_pierside:
            line = '{}\n'.format(pierside)
            file.write(line)

    with open('east_pierside.txt', 'w') as file:
        for pierside in east_pierside:
            line = '{}\n'.format(pierside)
            file.write(line)

    with open('flats.txt', 'w') as file:
        for flat in flat_files:
            line = '{}\n'.format(flat)
            file.write(line)


    print('Done...')
    return west_pierside, east_pierside, flat_data

def _make_master_bias(bias_data):
    print('##### Making master bias')

    bias_amount = len(bias_data)
    print('Using {} bias images for master bias'.format(bias_amount))

    master_bias = np.median(bias_data, axis=0)
    biashead = fits.Header()
    biashead['IMAGETYP'] = 'Master Bias'
    biasfits = fits.PrimaryHDU(master_bias, header=biashead)

    print('Saving master bias ----> masterbias.fits')
    biasfits.writeto('./masterbias.fits', overwrite=True)

    return master_bias

def _make_master_flat(flat_data, master_bias):
    """
    Makes the master flat. It normalizes the flat data then takes the median
    of the normalized flat to remove stars. Writes a file named
    "masterflat.fits".
    """
    print('\n### Creating master flat')

    flat_amount = len(flat_data)
    #print('Using {} images to make master flat'.format(flat_amount))

    if flat_amount <= 4:
        print('Not enough data to make master flat!')
        sys.exit('Exiting...')

    norm_flat = []
    for flat in tqdm(flat_data):
        flat_bias_sub = flat - master_bias

        norm_flat.append(np.divide(flat_bias_sub, np.median(flat_bias_sub)))

    print('Saving master flat ----> masterflat.fits')
    master_flat = np.median(norm_flat, axis=0)
    flathead = fits.Header()
    flathead['IMAGETYP'] = 'Master Flat'
    newflat = fits.PrimaryHDU(master_flat, header=flathead)

    newflat.writeto('./masterflat.fits', overwrite=True)

    print('Done...')
    return master_flat

def _image_reduce(west_pierside, east_pierside, master_bias, master_flat):
    """
    Goes through the east and west pierside images and divides them by the
    master flat. If a west pierside image, it must be rotated by 180 degrees.
    """
    print('\n### Subtracting bias, dividing flat for west images')
    for filename in tqdm(west_pierside):
        #print('Subtracting master bias, dividing master_flat for {}'.format(filename))
        with fits.open(filename, mode='update') as file:
            science_data = np.rot90(file[0].data, 2)
            science_header = file[0].header

            science_data = science_data - master_bias
            science_data = np.divide(science_data, master_flat)
            science_header['NSB_BIAS'] = 'True'
            science_header['NSB_FLAT'] = 'True'

            file[0].data = science_data

            file.flush()
    print('Done...')

    print('\n### Subtracting bias, dividing flat for east images')
    for filename in tqdm(east_pierside):
        #print('Subtracting master bias, dividing master_flat for {}'.format(filename))
        with fits.open(filename, mode='update') as file:
            science_data = file[0].data
            science_header = file[0].header

            science_data = science_data - master_bias
            science_data = np.divide(science_data, master_flat)

            science_header['NSB_BIAS'] = 'True'
            science_header['NSB_FLAT'] = 'True'

            file[0].data = science_data

            file.flush()
    print('Done...\n')

if __name__ == '__main__':
    """
    Main loop that goes through, seperates files by pierside, and creates a
    master flat with images taken above 75 degrees altitude. This is because
    uniformity is perfect at the zenith and degrades to about two percent per
    degree at a zenith angle near 70 degrees.
    """
    print('\n####### nsb_reduce.py')
    parser = argparse.ArgumentParser()
    parser.add_argument('files', nargs='+', help='files to reduce')
    parser.add_argument('-b', '--bias', help='apply master bias to frames')

    args = parser.parse_args()
    files = args.files
    bias = args.bias

    if bias:
        master_bias = getdata(bias)


    else:
        bias_data = _get_biases(files)
        master_bias = _make_master_bias(bias_data)

    west_pierside, east_pierside, flat_data = _get_pierside_info(files)
    master_flat = _make_master_flat(flat_data, master_bias)
    _image_reduce(west_pierside, east_pierside, master_bias, master_flat)
