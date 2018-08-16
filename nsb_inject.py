#!/usr/bin/env python3
from astropy.io import ascii
from astropy.io import fits
from tqdm import tqdm
import numpy as np

# constants that can be changed
BOLT_FILE = 'nsb_bolt.txt'
SQM_TEL_FILE = 'nsb_sqm_tel.txt'
SQM_ZENITH_FILE = 'nsb_sqm_zenith.txt'
NSB_FILE = 'nsb_data.txt'

def inject_headers_sqm(file, format):
    """
    Inject headers into the fits files with sqm measurements, telescope and
    zenith.

    Parameters
    file   | the sqm file (either nsb_sqm_tel.txt, nsb_sqm_zenith.txt)
    format | to know what header names to use
    """
    print('\n### Injecting sqm {} info into fits files'.format(format))
    data = ascii.read(file)

    if format == 'telescope':
        prefix = 'ST'
        do_az_elv = True
    elif format == 'zenith':
        prefix = 'SZ'
        do_az_elv = False
    for row in tqdm(data):
        with fits.open(row['filename'], mode='update') as fits_file:
            fits_header = fits_file[0].header

            fits_header['{}-OBS'.format(prefix)] = (
                row['sqm_ut'], 'UTC SQM expsoure'
            )
            fits_header['{}-TDEL'.format(prefix)] = (
                row['t_delta'], 'time diff (seconds) fits UTC to SQM UTC'
            )
            if do_az_elv:
                fits_header['{}-AZ'.format(prefix)] = (row['az'], 'degrees')
                fits_header['{}-ELV'.format(prefix)] = (row['elv'], 'degrees')

            fits_header['{}-COUNT'.format(prefix)] = (row['counts'], 'adu')
            fits_header['{}-NSB'.format(prefix)] = (
                row['sqm_nsb'], 'nsb mag/arcsec^2'
            )

            fits_file.flush()

    print('Done...')

def inject_headers_boltwood(file):
    """
    Inject headers into the fits files with boltwood measurements.

    Parameters
    file | the boltwood file (nsb_bolt.txt)
    """
    print('\n### Injecting boltwood info into fits files')
    data = ascii.read(file)

    for row in tqdm(data):
        with fits.open(row['filename'], mode='update') as fits_file:
            fits_header = fits_file[0].header

            fits_header['BT-OBS'] = (
                row['bolt_ut'], 'UTC of boltwood measure'
            )
            fits_header['BT-TDEL'] = (
                row['t_del'], 'time diff (seconds) fits UTC to bolt UTC'
            )
            fits_header['BT-SKYT'] = (row['sky_t'], 'sky temp (f)')
            fits_header['BT-AMBT'] = (row['amb_t'], 'ambient temp (f)')
            fits_header['BT-WIND'] = (row['wind'], 'wind (miles)')
            fits_header['BT-HUM'] = (row['hum'], 'humidity')
            fits_header['BT-DEW'] = (row['dew'], 'dew point')
            fits_header['BT-CLOUD'] = (
                row['cloud'], 'cloud flag 0=unk, 1=clr, 2=cloudy, 3=veryCloudy'
            )

            fits_file.flush()

    print('Done...')

def inject_headers_nsb(file):
    """
    Inject headers into the fits files with night sky brightness measurements.

    Parameters
    file | the nsb data file (nsb_data.txt)
    """
    print('\n### Injecting nsb measurements info fits files')
    data = ascii.read(file)

    for row in tqdm(data):
        with fits.open(row['filename'], mode='update') as fits_file:
            fits_header = fits_file[0].header

            fits_header['UA-MEAN'] = (
                row['mean_bkg'], 'mean background (counts/pixel)'
            )


            fits_header['UA-MED'] = (
                row['med_bkg'], 'median background (counts/pixel)'
            )
            fits_header['UA-STD'] = (row['std'], 'std of background')

            # check to see if value is nan, if so replace with 0
            if np.isnan(row['zp']) or np.isnan(row['nsb']):
                zp_value = (0, 'failed to find zeropoint')
                nsb_value = (0, 'failed to find nsb')
            else:
                zp_value = (row['zp'], 'zeropoint (mag)')
                nsb_value = (row['nsb'], 'nsb measurement (mag/arcsec^2)')

            fits_header['UA-ZP'] = zp_value
            fits_header['UA-NSB'] = nsb_value

            fits_file.flush()

    print('Done...')

if __name__ == '__main__':
    """
    Start of the program.
    """
    print('\n####### nsb_inject.py')

    inject_headers_sqm(SQM_TEL_FILE, 'telescope')

    inject_headers_sqm(SQM_ZENITH_FILE, 'zenith')

    inject_headers_boltwood(BOLT_FILE)

    inject_headers_nsb(NSB_FILE)

    print('\nFINISHED\n')
