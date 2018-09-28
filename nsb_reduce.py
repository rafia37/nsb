#!/usr/bin/env python3
import sys, os
import argparse
from astropy.io import fits
from astropy.io.fits import getheader, getval, getdata
import numpy as np
import shutil
from tqdm import tqdm

def _image_reduce(files, master_bias, master_flat, do_bias=True):
    """
    Goes through the east and west pierside images and divides them by the
    master flat. If a west pierside image, it must be rotated by 180 degrees.
    """
    print('\n#---------- Reducing images')
    print('Subtracting bias, dividing flat for images')
    for file in tqdm(files):
        with fits.open(file, mode='update') as f:
            data = f[0].data

            science_header = f[0].header

            if do_bias:
                data = np.subtract(data, master_bias)
                science_header['NSB_BIAS'] = 'True'
            else:
                data = np.subtract(data, 1100.)


            #science_data = np.rot90(science_data, 2)

            science_data = np.divide(data, master_flat)
            #science_data = np.clip(science_data, 0, 65535)

            science_header['NSB_FLAT'] = 'True'

            f[0].data = science_data

            f.flush()
    print('\nDone...')

if __name__ == '__main__':
    """
    Main loop that goes through, seperates files by pierside, and reduces the
    science images.
    """
    print('\n########## nsb_reduce.py')
    print('Reduces science frames by subtracting bias and dividing flat')
    parser = argparse.ArgumentParser()
    parser.add_argument('files', nargs='+', help='files to reduce')
    parser.add_argument('-b', '--bias', help='the master bias')
    parser.add_argument('-f', '--flat', help='the master flat')

    args = parser.parse_args()
    files = args.files
    bias = args.bias
    flat = args.flat


    if bias:
        master_bias = getdata(bias)
    else:
        do_bias = False
    master_flat = getdata(flat)

    _image_reduce(files, master_bias, master_flat)

    print('\nFINISHED\n')
