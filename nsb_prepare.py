#!/usr/bin/env python3
from astropy.io import fits
import argparse
import ephem
from tqdm import tqdm
import math
from astropy import wcs
from astropy import units as u
from astropy.coordinates import SkyCoord
import os, sys

import warnings

def inject_headers(files):
    """
    Injects headers needed for the photometry pipeline. Adds the headers that
    LL would like for their processing. Updates the existing files.

    Parameters
    files : list of files from argparse

    Returns
    None
    """
    print('\n####### nsb_prepare.py')
    print('\n### Adding headers for PHOTOMETRY PIPELINE and NSB LL')
    for filename in tqdm(files):

        with fits.open(filename, mode='update') as image:
            image_header = image[0].header
            image_header['GAIN'] = ('1.4', 'ADU/electron')
            image_header['RDNOISE'] = ('23', 'electrons')
            image_header['TELESCOP'] = 'Pomenis'
            image_header['INSTRUME'] = 'Apogee F9000'
            image_header['LAT-OBS'] = ('32,442525', 'degrees')
            image_header['LON-OBS'] = ('-110.789161', 'degrees')
            image_header['ALT-OBS'] = ('2791', 'meters')
            image_header['BUNIT'] = ('ADU')
            image_header['BSCALE'] = 1
            image_header['PHOT-CO'] = 'na'
            image_header['RADESYSa'] = 'FK5'
            try:
                image_header.remove('RADECSYS')
            except KeyError:
                pass
            image_header['COORDSYS'] = ('FK5', 'coordinate system for ra,dec')
            datetime = image_header['DATE-OBS']
            telescope_az = math.radians(float(image_header['AZIMUTH']))
            telescope_alt = math.radians(float(image_header['ALTITUDE']))
            telescope_point = (telescope_az, telescope_alt)
            moon_distance = calculate_moon_distance(datetime, telescope_point)
            image_header['MOONDIST'] = (moon_distance, 'degrees')
            image_header['IFOV'] = ('10136', 'arcseconds')

            # get center ra, dec from image
            # use this for blocks in future
            try:
                w = wcs.WCS(image[0].header)
                lat, lon = w.wcs_pix2world(1528, 1528, 1)
                coord = SkyCoord(lat*u.degree, lon*u.degree, frame='fk5')
                coord_string = coord.to_string('hmsdms', sep=':')
                ra, dec = coord_string.split(' ')

                image_header['RA'] = ra
                image_header['DEC'] = dec

            except KeyError:
                print('ERROR {}'.format(filename))
                pass

            #ra_string = '{}:{}:{}'.format(ra)
            # comment out if you do not want to crop
            #science_data = image[0].data
            #image[0].data = science_data[500:2556,500:2556]

            image.flush()

def rename(files):
    """
    Renames the files to .fits.
    """

    for file in files:
        name, ext = os.path.splitext(file)

        if ext == '.fits':
            pass
        else:

            new_name = name + '.fits'
            print('{} ----. {}'.format(file, new_name))
            os.rename(file, new_name)

    return None


def calculate_moon_distance(datetime, telescope_point):
    """
    Uses pyephem to calculate the lunar distance from the telescope pointing.

    Parameters
    datetime        : date and time of the observation
    telescope_point : tuple of the azimuth and altitude of the telescope

    Returns
    degree_distance : distance from the moon in degrees, rounded to 2 points
    """
    # create pomenis observer for Lemmon
    # NOTE time is in UT of the observation
    pomenis_lemmon = ephem.Observer()
    pomenis_lemmon.lon, pomenis_lemmon.lat = '-110.789161', '32.442525'
    pomenis_lemmon.elevation = 2791
    datetime = '{} {}'.format(datetime[:10], datetime[12:])
    pomenis_lemmon.date = datetime

    # create moon object and calculate the az and alt for a time at pomenis
    moon = ephem.Moon()
    moon.compute(pomenis_lemmon)

    moon_az, moon_alt = moon.az, moon.alt
    moon_point = (moon_az, moon_alt)

    # calculate the seperation from the telescope pointing to moon
    distance = ephem.separation(telescope_point, moon_point)
    degree_distance = round(math.degrees(distance), 2)

    return degree_distance


if __name__ == '__main__':
    """
    Parses user input for list of files.
    """
    parser = argparse.ArgumentParser()
    parser.add_argument(dest='files', nargs='+')
    parser.add_argument('-r', '--rename', dest='rename', action='store_true')
    args = parser.parse_args()

    if args.rename:
        rename(args.files)
        print('\nDone...\n')

    else:
        inject_headers(args.files)
        print('\nDone...\n')
