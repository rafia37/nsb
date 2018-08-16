#!/usr/bin/env python3
import sys, os
import argparse
from astropy.io import fits
from astropy.io.fits import getheader, getval, getdata
import numpy as np
import shutil
from tqdm import tqdm

def _get_pierside_info(files, fil):
    """
    Get the pierside info to rotate the frame.
    """
    # create lists
    west_pierside = []
    east_pierside = []

    print('\n### Getting pierside info and matching filter {}'.format(fil))
    for filename in tqdm(files):
        pierside = getval(filename, 'PIERSIDE')
        image_fil = getval(filename, 'FILTER')

        if pierside =='WEST' and image_fil == fil:
            west_pierside.append(filename)

        elif pierside == 'EAST' and image_fil == fil:
            east_pierside.append(filename)

    with open('west_pierside-{}.txt'.format(fil), 'w') as file:
        for pierside in west_pierside:
            line = '{}\n'.format(pierside)
            file.write(line)

    with open('east_pierside-{}.txt'.format(fil), 'w') as file:
        for pierside in east_pierside:
            line = '{}\n'.format(pierside)
            file.write(line)

    print('Done...')
    return west_pierside, east_pierside

def _image_reduce(west_pierside, east_pierside, master_bias, master_flat):
    """
    Goes through the east and west pierside images and divides them by the
    master flat. If a west pierside image, it must be rotated by 180 degrees.
    """
    print('\n### Subtracting bias, dividing flat for west images')
    for filename in tqdm(west_pierside):
        with fits.open(filename, mode='update') as file:
            data = file[0].data
            science_data = data - master_bias

            #science_data = np.rot90(science_data, 0)
            science_header = file[0].header

            science_data = np.divide(science_data, master_flat)
            science_header['NSB_BIAS'] = 'True'
            science_header['NSB_FLAT'] = 'True'

            file[0].data = science_data

            file.flush()
    print('Done...')

    print('\n### Subtracting bias, dividing flat for east images')
    for filename in tqdm(east_pierside):
        with fits.open(filename, mode='update') as file:
            data = file[0].data

            #science_data = np.rot90(data, 0)
            science_header = file[0].header

            science_data = data - master_bias
            science_data = np.divide(science_data, master_flat)

            science_header['NSB_BIAS'] = 'True'
            science_header['NSB_FLAT'] = 'True'

            file[0].data = science_data

            file.flush()
    print('Done...\n')

if __name__ == '__main__':
    """
    Main loop that goes through, seperates files by pierside, and reduces the
    science images.
    """
    print('\n##### nsb_reduce.py')
    parser = argparse.ArgumentParser()
    parser.add_argument('files', nargs='+', help='files to reduce')
    parser.add_argument('-b', '--bias', help='the master bias')
    parser.add_argument('-f', '--flat', help='the master flat')

    args = parser.parse_args()
    files = args.files
    bias = args.bias
    flat = args.flat

    master_bias = getdata(bias)
    master_flat = getdata(flat)
    fil = getval(flat, 'FILTER')

    west_pierside, east_pierside = _get_pierside_info(files, fil)
    _image_reduce(west_pierside, east_pierside, master_bias, master_flat)
