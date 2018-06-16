import argparse
from astropy.io.fits import getval
import numpy as np
from astropy.time import TimeDelta
from astropy.time import Time
from astropy.stats import sigma_clipped_stats
from astropy.coordinates import Angle
import astropy.units as u

def calculate_nsb(files):
    """

    """
    for file in files:
        zeropoint = []
        zeropoint_error = []
        with open(file, 'r') as f:
            for line in f:
                split_line = line.strip('\n').split('\t')

                if split_line[0] == '#filename':

                    filename = split_line[1]
                    alt = getval(split_line[1], 'ALTITUDE')
                    az = getval(split_line[1], 'AZIMUTH')
                    date_obs = getval(split_line[1], 'DATE-OBS')
                    t_delta = TimeDelta(25200, format='sec')
                    local_date_obs = Time(date_obs) - t_delta
                    print(filename)
                elif split_line[0] == '#platescale_x':
                    platescale = float(split_line[1])
                elif split_line[0] == '#exp_time':
                    exposure = float(split_line[1])
                elif split_line[0] == '#mask_mean':
                    bkg_mean = float(split_line[1])
                elif split_line[0] == '#mask_med':
                    bkg = float(split_line[1])
                elif split_line[0] == '#mask_std':
                    std = float(split_line[1])
                elif split_line[0] == '#center_ra':
                    ra = split_line[1]
                    ra = Angle(ra, u.degree)
                    ra = ra.to_string(u.hour, sep=':', precision=2, pad=True)
                elif split_line[0] == '#center_dec':
                    dec = split_line[1]
                    dec = Angle(dec, u.degree)
                    dec = dec.to_string(u.degree, sep=':', precision=1, alwayssign=True, pad=True)
                elif split_line[0][0] == '#':
                    pass
                elif split_line[0] == 'id':
                    pass
                else:
                    zeropoint.append(float(split_line[1]))
                    zeropoint_error.append(float(split_line[2]))


        med_zero = np.median(zeropoint)
        sigma_clipped_stats(zeropoint, sigma=1.0, iters=10)
        print(med_zero)
        med_zero_error = np.median(zeropoint_error)

        #print(exposure, bkg, platescale)
        nsb = med_zero - (2.5*np.log( (bkg) / platescale ** 2))
        filter = 'g'
        ifov = '15066'
        phot_co = ''
        with open('nsb_data.txt', 'a') as f:
            line = '{}\t{}\t{}\t{:.2f}\t{:.2f}\t{}\t{}\t{}\t{}\t{}\t{:.2f}\t\t{:.2f}\t{:.2f}\t{:.2f}\t{:.2f}\tna\n'.format(
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
    with open('nsb_data.txt', 'w') as f:
        header = '#filename\t\t\t\t\tut\t\t\tlocal\t\t\taz\telv\tra\t\tdec\t\texp\tfilt\tifov\tmean_bkg\tmed_bkg\tstd\tzp\tnsb\tphot_co\n'
        f.write(header)
    calculate_nsb(files)
