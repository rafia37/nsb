#!/usr/env/python
from astropy.io import fits
import argparse
import ephem
import math

def inject_headers(files):
    print('### Adding headers for PHOTOMETRY PIPELINE')
    for filename in files:
        print('Adding headers to {}'.format(filename))

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
            image_header['COORDSYS'] = ('ICRS', 'coordinate system for ra,dec')
            datetime = image_header['DATE-OBS']
            telescope_az = math.radians(float(image_header['AZIMUTH']))
            telescope_alt = math.radians(float(image_header['ALTITUDE']))
            telescope_point = (telescope_az, telescope_alt)
            moon_distance = calculate_moon_distance(datetime, telescope_point)
            image_header['MOONDIST'] = (moon_distance, 'degrees')
            image_header['IFOV'] = ('15066', 'arcseconds')



            image.flush()

def calculate_moon_distance(datetime, telescope_point):
    moon = ephem.Moon()

    pomenis_lemmon = ephem.Observer()
    pomenis_lemmon.lon, pomenis_lemmon.lat = '-110.789161', '32.442525'
    pomenis_lemmon.elevation = 2791
    datetime = '{} {}'.format(datetime[:10], datetime[12:])
    # pass the ut time for observation, MST is -7
    pomenis_lemmon.date = datetime


    moon = ephem.Moon()

    moon.compute(pomenis_lemmon)

    moon_az, moon_alt = moon.az, moon.alt

    moon_point = (moon_az, moon_alt)

    distance = ephem.separation(telescope_point, moon_point)
    degree_distance = round(math.degrees(distance), 2)
    print('moon point = {}'.format(moon_point))
    print('telescope point = {}'.format(telescope_point))
    print('moon distance = {}'.format(degree_distance))
    return degree_distance


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument(dest='files', nargs='+')
    args = parser.parse_args()

    inject_headers(args.files)
