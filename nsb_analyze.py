#!/usr/bin/env python3

import argparse
from astropy.io.fits import getval
import numpy as np
from astropy.time import TimeDelta
from astropy.time import Time
from astropy.stats import sigma_clipped_stats
from astropy.coordinates import Angle
import astropy.units as u
from tqdm import tqdm

def calculate_nsb(files):
    """
    calculate the nsb measurement using the equation below. Need to comment this
    code more haha.
    """
    for file in tqdm(files):
        zeropoint = []
        zeropoint_error = []
        with open(file, 'r') as f:
            pp_zeropoint = np.nan
            for line in f:
                split_line = line.strip('\n').split('\t')
                if split_line[0] == '#filename':

                    filename = split_line[1]
                    alt = getval(split_line[1], 'ALTITUDE')
                    az = getval(split_line[1], 'AZIMUTH')
                    date_obs = getval(split_line[1], 'DATE-OBS')
                    t_delta = TimeDelta(25200, format='sec')
                    local_date_obs = Time(date_obs) - t_delta
                    image_size = float(getval(split_line[1], 'NAXIS1')) ** 2
                elif split_line[0] == '#platescale_x':
                    platescale = float(split_line[1])
                elif split_line[0] == '#platescale_y':
                    pass
                elif split_line[0] == '#middle_jd':
                    pass
                elif split_line[0] == '#exp_time':
                    exposure = float(split_line[1])
                elif split_line[0] == '#mask_mean':
                    bkg_mean = float(split_line[1])
                elif split_line[0] == '#mask_med':
                    bkg = float(split_line[1])
                elif split_line[0] == '#std':
                    std = float(split_line[2])
                elif split_line[0] == '#mean':
                    nomask_mean = float(split_line[2])
                elif split_line[0] == '#med':
                    nomask_median = float(split_line[2])
                elif split_line[0] == '#amount':
                    pixels_used = float(split_line[2])
                    percent_used = pixels_used / image_size
                elif split_line[0] == '#center_ra':
                    ra = split_line[1]
                    ra = Angle(ra, u.degree)
                    ra = ra.to_string(u.hour, sep=':', precision=2, pad=True)
                elif split_line[0] == '#center_dec':
                    dec = split_line[1]
                    dec = Angle(dec, u.degree)
                    dec = dec.to_string(u.degree, sep=':', precision=1, alwayssign=True, pad=True)
                elif split_line[0] == '#ifov':
                    ifov = float(split_line[1]) * platescale
                elif split_line[0] == '#filter':
                    filter = split_line[1]


                elif split_line[0][0] == '#' and split_line[0] != '#pp zeropoint':
                    pass
                elif split_line[0] == 'id':
                    pass

                elif split_line[0] == '#pp zeropoint':
                    pp_zeropoint = float(split_line[1])
                else:
                    zeropoint.append(float(split_line[1]))
                    zeropoint_error.append(float(split_line[2]))


        med_zero = pp_zeropoint

        #print(pp_zero)
        #med_zero = np.median(zeropoint)
        #sigma_clipped_stats(zeropoint, sigma=1.0, iters=10)
        #print(med_zero)
        #med_zero_error = np.median(zeropoint_error)

        nsb = med_zero - (2.5*np.log10( ((bkg / exposure ) / (platescale ** 2)) ))


        phot_co = ''
        with open('nsb_data.csv', 'a') as f:
            line = '{},{},{},{:.2f},{:.2f},{},{},{},{},{},{:.2f},{:.2f},{:.2f},{:.2f},{:.2f},na\n'.format(
                filename,
                date_obs,
                str(local_date_obs)[:19],
                az,
                alt,
                ra,
                dec,
                exposure,
                filter,
                ifov,
                bkg_mean,
                bkg,
                std,
                med_zero,
                nsb
            )
            f.write(line)

if __name__ == '__main__':
    """

    """
    parser = argparse.ArgumentParser()
    parser.add_argument(dest='files', nargs='+')
    files = parser.parse_args().files
    with open('nsb_data.csv', 'w') as f:
        header = '#filename,ut,local,az,elv,ra,dec,exp,filt,ifov,mean_bkg,med_bkg,std,zp,nsb,phot_co\n'
        f.write(header)

    print('\n##### nsb_analyze\n')
    print('### calculating the nsb measurement')
    print('msky = Z - 2.5log10(med_sky / exposure / platescale^2)\n')
    calculate_nsb(files)
    print('\nresults ----> nsb_data.csv')
    print('\nDone...\n')
