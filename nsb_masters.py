#!/usr/bin/env python3
import sys, os
import argparse
from astropy.io import fits
from astropy.io.fits import getheader, getval, getdata
import numpy as np
import shutil
from tqdm import tqdm
def _get_biases(files):
    """

    """
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

def _get_flats(files):
    """

    """
    flat_data_v = []
    flat_files_v = []
    flat_data_w = []
    flat_files_w = []
    flat_data_z = []
    flat_files_z = []
    for filename in files:
        type = getval(filename, 'IMAGETYP')

        if type == 'Flat Field':
            fil = getval(filename, 'FILTER')

            if fil == 'v':
                flat_data_v.append(getdata(filename))
                flat_files_v.append(filename)
            elif fil == 'w':
                flat_data_w.append(getdata(filename))
                flat_files_w.append(filename)
            elif fil == 'z':
                flat_data_z.append(getdata(filename))
                flat_files_z.append(filename)
            else: pass
        else: pass

    if len(flat_files_v) != 0:
        with open('flat_v.txt', 'w') as file:
            for flat in flat_files_v:
                line = '{}\n'.format(flat)
                file.write(line)
    if len(flat_files_w) != 0:
        with open('flat_w.txt', 'w') as file:
            for flat in flat_files_w:
                line = '{}\n'.format(flat)
                file.write(line)
    if len(flat_files_z) != 0:
        with open('flat_z.txt', 'w') as file:
            for flat in flat_files_z:
                line = '{}\n'.format(flat)
                file.write(line)

    return flat_data_v, flat_data_w, flat_data_z


def _make_master_bias(bias_data):
    """

    """
    print('\n### Making master bias')

    bias_amount = len(bias_data)
    print('Using {} bias images for master bias'.format(bias_amount))

    master_bias = np.median(bias_data, axis=0)
    biashead = fits.Header()
    biashead['IMAGETYP'] = 'Master Bias'
    biasfits = fits.PrimaryHDU(master_bias, header=biashead)

    print('Saving master bias ----> masterbias.fits')
    biasfits.writeto('./masterbias.fits', overwrite=True)

    return master_bias

def _make_master_flat(flat_data, master_bias, fil):
    """
    Makes the master flat. It normalizes the flat data then takes the median
    of the normalized flat to remove stars. Writes a file named
    "masterflat.fits".
    """
    if len(flat_data) == 0:
        return None

    print('\n### Creating master flat {}'.format(fil))

    flat_amount = len(flat_data)
    #print('Using {} images to make master flat'.format(flat_amount))

    norm_flat = []
    for flat in tqdm(flat_data):
        flat_bias_sub = flat - master_bias

        norm_flat.append(np.divide(flat_bias_sub, np.median(flat_bias_sub)))

    print('Saving master flat ----> masterflat-{}.fits'.format(fil))
    master_flat = np.median(norm_flat, axis=0)
    flathead = fits.Header()
    flathead['IMAGETYP'] = 'Master Flat'
    flathead['FILTER'] = fil.lower()
    newflat = fits.PrimaryHDU(master_flat, header=flathead)

    newflat.writeto('./masterflat-{}.fits'.format(fil), overwrite=True)

    print('Done...')
    #return master_flat

if __name__ == '__main__':
    """
    Main loop that goes through the fits files and creates master bias and master flats for each filter.
    """
    print('\n##### nsb_masters.py')
    print('Creates master bias and master flats for filters V, W, and z')
    parser = argparse.ArgumentParser()
    parser.add_argument('files', nargs='+', help='files to reduce')

    args = parser.parse_args()
    files = args.files

    bias_data = _get_biases(files)
    master_bias = _make_master_bias(bias_data)
    v, w, z = _get_flats(files)

    master_flat = _make_master_flat(v, master_bias, 'V')
    master_flat = _make_master_flat(w, master_bias, 'W')
    master_flat = _make_master_flat(z, master_bias, 'z')

    print('\nFINISHED\n')
