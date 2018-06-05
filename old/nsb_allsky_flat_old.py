#!/usr/env/python
import sys, os
from glob import glob
import argparse
from astropy.io import fits
from astropy.io.fits import getheader, getval, getdata
import numpy as np

# show the telescope side
def _get_pierside_info(files):
    """
    Get the pierside info. This will be moved to the setup section when bias
    and darks are introduced.
    """
    west_pierside = []
    east_pierside = []
    flat_data = []

    for filename in files:
        pierside = getval(filename, 'PIERSIDE')
        altitude = getval(filename, 'ALTITUDE')

        if pierside == 'WEST' and altitude >= 75:
            data = np.rot90(getdata(filename), 2)
            flat_data.append(data)
            west_pierside.append(filename)

        elif pierside =='WEST':
            west_pierside.append(filename)

        elif pierside == 'EAST' and altitude >= 75:
            data = getdata(filename)
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

    return west_pierside, east_pierside, flat_data

def _make_master_flat(flat_data):
    """
    Makes the master flat. It normalizes the flat data then takes the median
    of the normalized flat to remove stars. Writes a file named
    "masterflat.fits".
    """
    print('# Making master flat')

    flat_amount = len(flat_data)
    print('Using {} images to make master flat'.format(flat_amount))

    if flat_amount <= 4:
        print('Not enough data to make master flat!')
        sys.exit('Exiting...')

    norm_flat = []
    for flat in flat_data:
        norm_flat.append(np.divide(flat, np.median(flat)))

    master_flat = np.median(norm_flat, axis=0)
    flathead = fits.Header()
    flathead['IMAGETYP'] = 'Master Flat'
    newflat = fits.PrimaryHDU(master_flat, header=flathead)
    newflat.writeto('masterflat.fits', overwrite=True)

    return master_flat

def _image_reduce(west_pierside, east_pierside, master_flat):
    """
    Goes through the east and west pierside images and divides them by the
    master flat. If a west pierside image, it must be rotated by 180 degrees.
    """
    print('# Dividing WEST images')
    for filename in west_pierside:
        print('Updating {}'.format(filename))
        with fits.open(filename, mode='update') as file:
            science_data = np.rot90(file[0].data, 2)
            science_header = file[0].header
            science_header['NSB_FLAT'] = 'True'

            science_data = np.divide(science_data, master_flat)
            file[0].data = science_data

            file.flush()

    print('# Dividing EAST images')
    for filename in east_pierside:
        print('Updating {}'.format(filename))
        with fits.open(filename, mode='update') as file:
            science_data = file[0].data
            science_header = file[0].header
            science_header['NSB_FLAT'] = 'True'

            science_data = np.divide(science_data, master_flat)
            file[0].data = science_data

            file.flush()

def main():
    """
    Main loop that goes through, seperates files by pierside, and creates a
    master flat with images taken above 75 degrees altitude. This is because
    uniformity is perfect at the zenith and degrades to about two percent per
    degree at a zenith angle near 70 degrees.
    """
    parser = argparse.ArgumentParser()
    parser.add_argument(dest='files', nargs='+')
    files = parser.parse_args().files

    west_pierside, east_pierside, flat_data = _get_pierside_info(files)
    master_flat = _make_master_flat(flat_data)

    _image_reduce(west_pierside, east_pierside, master_flat)

if __name__ == '__main__':
    main()
