#!/usr/bin/env python3
from astropy.io import fits
import argparse
import ephem
import math
from astropy import wcs
from astropy import units as u
from astropy.coordinates import SkyCoord
from astropy.coordinates import Angle
import os, sys
import numpy as np
from tqdm import tqdm

def _inject_headers(files, crop=False):
    """
    Injects headers needed for the photometry pipeline. Adds the headers that
    LL would like for their processing. Updates the existing files.

    Parameters
    files : list of files from argparse

    Returns
    None
    """
    print('\n#---------- Adding headers for PHOTOMETRY PIPELINE and NSB LL')
    print('Calculating moon distance and injecting headers')
    if crop:
        print('Cropping images')
    for filename in tqdm(files):

        with fits.open(filename, mode='update') as image:
            image_header = image[0].header
            image_header['GAIN'] = ('1.4', 'ADU/electron')
            image_header['RDNOISE'] = ('23', 'electrons')
            image_header['TELESCOP'] = 'Pomenis'
            image_header['INSTRUME'] = 'Apogee F9000'
            image_header['LAT-OBS'] = ('31.9583295', 'degrees')
            image_header['LON-OBS'] = ('-111.596664', 'degrees')
            image_header['ALT-OBS'] = ('2096', 'meters')
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
            site_latitude = float(image_header['LAT-OBS'])
            site_longitude = float(image_header['LONG-OBS'])
            site_altitude = float(image_header['ALT-OBS'])
            site = (site_latitude, site_longitude, site_altitude)
            telescope_point = (telescope_az, telescope_alt)
            moon_distance = _calculate_moon_distance(
                datetime,
                telescope_point,
                site
            )
            pixels = float(image_header['NAXIS1'])
            platescale = 4.93

            image_header['MOONDIST'] = (moon_distance, 'degrees')

            # get center ra, dec from image
            # use this for blocks in future
            try:
                w = wcs.WCS(image[0].header)
                lat, lon = w.wcs_pix2world(1528, 1528, 1)
                coord = SkyCoord(lat*u.degree, lon*u.degree, frame='fk5')
                coord_string = coord.to_string('hmsdms', sep=':')
                ra, dec = coord_string.split(' ')

                new_ra = Angle(ra, unit=u.hour)
                new_dec = Angle(dec, unit=u.degree)
                old_ra = Angle(image_header['RA'], unit=u.hour)
                old_dec = Angle(image_header['DEC'], unit=u.degree)

                """
                pointing_error = np.sqrt((((new_ra.degree - old_ra.degree)**2) +
                    ((new_dec.degree - old_dec.degree)**2))*3600.)

                print('Pointing error: {} arcseconds'.format(pointing_error))
                point1 = new_ra.degree * np.rad2deg(np.cos(new_dec.degree))
                point2 = old_ra.degree * np.rad2deg(np.cos(old_dec.degree))
                print((point1 - point2)*3600.)

                image_header['P-ERROR'] = (
                    pointing_error,
                    'Pointing Error arcsec'
                )
                """
                image_header['RA'] = ra
                image_header['DEC'] = dec


            except KeyError:
                print('Error, no WCS information, REMOVE {}'.format(filename))
                pass

            if crop:
                science_data = image[0].data
                image[0].data = science_data[500:2556,500:2556]
                pixels = len(image[0].data[0])

            image_header['IFOV'] = int(pixels * platescale)

            image.flush()

    print('\nDone...')

def _rename(files):
    """
    Renames the files to .fits.
    """

    for file in files:
        name, ext = os.path.splitext(file)

        if ext == '.fits':
            pass
        else:

            new_name = name + '.fits'
            print('{} ----> {}'.format(file, new_name))
            os.rename(file, new_name)
    print('\nDone...')

    return None

def _calculate_moon_distance(datetime, telescope_point, site):
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
    pomenis_lemmon.lon, pomenis_lemmon.lat = site[1], site[0]
    pomenis_lemmon.elevation = site[2]
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

def _parse_arguments():
    """
    Parse command line arguments.
    """
    parser = argparse.ArgumentParser()
    parser.add_argument(dest='files', nargs='+')
    parser.add_argument('-r', '--rename', dest='rename', action='store_true')
    parser.add_argument('-c', '--crop', dest='crop', action='store_true')
    args = parser.parse_args()

    return args.files, args.rename, args.crop

def main():
    print('\n########## nsb_prepare.py')
    print('Calculates moon distance from pointing and injects headers')

    files, rename, crop = _parse_arguments()

    if rename:
        _rename(files)

    else:
        _inject_headers(files, crop)
        print('\nFINISHED\n')

    sys.exit(0)

if __name__ == '__main__':
    main()

