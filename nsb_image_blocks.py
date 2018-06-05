#!/usr/bin/env python3
import sys, os
from astropy.io.fits import getval
from astropy.io.fits import getdata, getheader
from astropy.io import fits
import argparse
import logging
from skimage.util.shape import view_as_blocks
import numpy as np
import matplotlib.pyplot as plt
from astropy import wcs

def create_blocks(files, side=191):
    """
    Read in the science images and create blocks of them.
    764 pixels x 764 pixels = 16 blocks per image
    382 pixels x 382 pixels = 64 blocks per image
    191 pixels x 191 pixels = 256 blocks per image

    If binned 2x2 then the amount of blocks returned are cut in half.
    """

    if not os.path.exists('./temp'):
        print('# Creating temp folder to store cubed fit files')
        os.makedirs('./temp')

    print('# Creating cubed fits files')
    for filename in files:
        # check if fits file
        science_data = getdata(filename)
        header = getheader(filename)
        #new_header = add_headers(header_data)
        print('{} ----> view as blocks in temp folder'.format(filename))
        images = view_as_blocks(science_data, block_shape=(side, side))

        for index, i in enumerate(images):
            new_science = fits.PrimaryHDU(
                i,
                header=header
            )
            new_science.writeto(
                './temp/temp{}-{}'.format(index+1, filename), overwrite=True)

def unpack_blocks():
    """
    Creates a final block directory, then goes over the cube fits files and
    unpacks them. Insert headers needed for photometry pipeline.
    """
    # make the small directory
    print('# Creating blocks folder to store unpacked blocks')
    if not os.path.exists('./blocks'):
        os.makedirs('./blocks')


    # open the fits files in multi
    print('# Unpacking blocks')
    index = 0
    for files in os.listdir('./temp'):
        print('{} ----> unpacking blocks'.format(files))

        with fits.open('./temp/{}'.format(files)) as filename:
            header = filename[0].header
            header['GAIN'] = 2
            header['RDNOISE'] = 15
            header['TELESCOP'] = 'Pomenis'
            header['INSTRUME'] = 'Apogee F9000'

            # go through the multi files and unpack
            for i in range(len(filename[0].data)):
                print('Unpacked block {}'.format(i+1))
                science_data = filename[0].data[i]
                new_image = fits.PrimaryHDU(
                    science_data,
                    header=header
                )

                # naming, need to work on this
                temp_name = files.split('-')[1:]
                temp_name = '-'.join(temp_name)
                temp_name = temp_name.split('.')[0]
                temp_name = '{}-block{}.fits'.format(temp_name, index+1)
                print(temp_name)

                unpacked_name = './blocks/{}'.format(
                    temp_name,
                    index+1
                )
                print('Saving as {}.fits'.format(unpacked_name))
                new_image.writeto(unpacked_name,overwrite=True)
                index += 1

def _cleanup():
    """
    Removes the cubed fits files to get rid of unneeded data.
    """
    os.system('rm -rf ./temp')

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument(dest='files', nargs='+')
    files = parser.parse_args().files
    create_blocks(files)
    unpack_blocks()
    _cleanup()


if __name__ == '__main__':

    main()
