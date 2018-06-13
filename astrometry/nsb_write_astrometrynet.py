#!/usr/env/python3
from astropy.io import fits
import os, sys
DEST = './corrected_astrometry/'


with open('astrometry_results.txt', 'r') as file:
    for line in file:
        new_line = line.strip('\n').split(',')
        filename = new_line[0]
        ra = new_line[2]
        ra = ra.replace('h', '')
        ra = ra.replace('m', '')
        ra = ra.replace('s', '')

        dec = new_line[1]
        dec = dec.replace('d', '')
        dec = dec.replace('m', '')
        dec = dec.replace('s', '')

        with fits.open(filename, mode='update') as fits_file:
            print('Updating {}'.format(filename))
            header = fits_file[0].header
            header['RA'] = ra
            header['DEC'] = dec

            fits_file.flush()

        os.rename(filename, DEST+filename)
