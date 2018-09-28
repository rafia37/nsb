#!/usr/bin/env python3
import argparse
from astropy.time import Time
from astropy.io import fits
from astropy.time import TimeDelta
from tqdm import tqdm
import mmap
import warnings
warnings.filterwarnings('ignore')
# read the file

def get_num_lines(file_path):
    """

    """
    file_path = open(file_path, 'r+')
    buff = mmap.mmap(file_path.fileno(),0)
    lines = 0
    while buff.readline():
        lines+=1
    return lines

def parse_sqm(file):
    """

    """
    data = [[] for _ in range(3)]
    with open(file, 'r') as f:
        for line in f:
            if line[0] != '#':
                line = line.strip('\n').split(';')

                # 0 = ut date
                # 1 = counts
                # 2 = mag
                data[0].append(line[0])
                data[1].append(line[3])
                data[2].append(line[5])
        return data

def parse_bolt(file):
    """

    """
    data = [[] for _ in range(11)]
    with open(file, 'r') as f:
        for line in f:
            line = line.strip('\n').split()
            data[0].append('{}T{}'.format(line[0], line[1])) # date and time
            data[1].append(line[8]) # sky temp (f)
            data[2].append(line[9]) # ambient temp (f)
            data[3].append(line[4]) # wind speed (miles)
            data[4].append(line[10]) # humidity
            data[5].append(line[11]) # dew
            data[6].append(line[5]) # rain flag
            data[7].append(line[7]) # wet flag
            data[8].append(line[2]) # cloud flag
            data[9].append(line[3]) # wind flag
            data[10].append(line[6]) # rain flag
        return data

def calculate_tdelta_sqm(nsb_file, sqm_tel_data, sqm_zenith_data, filter_used):
    """

    """
    zenith_add = []
    tel_add = []
    print('\n### Finding zenith and telescopeSQM data for pointings')
    with open(nsb_file, 'r') as file:
        for line in tqdm(file, total=get_num_lines(nsb_file)):
            if line[0] != '#':
                line_split = line.strip('\n').split(',')
                min_delta = 10 # create just some min

                for index, x in enumerate(sqm_tel_data[0]):
                    delta = Time(line_split[1]) - Time(x)
                    if abs(float(delta.sec)) < min_delta:
                        min_delta = abs(float(delta.sec))
                        min_timestamp = x
                        min_counts = sqm_tel_data[1][index]
                        min_mag = sqm_tel_data[2][index]
                        az = line_split[3]
                        elv = line_split[4]
                        tel_col = '{},{},{},{},{},{},{:.1f}\n'.format(
                            line_split[0],
                            str(min_timestamp)[:19],
                            az,
                            elv,
                            min_counts,
                            min_mag,
                            min_delta
                        )
                    else:
                        pass
                tel_add.append(tel_col)

                min_delta = 10 # create just some min
                for index, x in enumerate(sqm_zenith_data[0]):
                    delta = Time(line_split[1]) - Time(x)
                    if abs(float(delta.sec)) < min_delta:
                        min_delta = abs(float(delta.sec))
                        min_timestamp = x
                        min_counts = sqm_zenith_data[1][index]
                        min_mag = sqm_zenith_data[2][index]
                        zenith_col = '{},{},{},{},{:.1f}\n'.format(
                            line_split[0],
                            str(min_timestamp)[:19],
                            min_counts,
                            min_mag,
                            min_delta
                        )

                    else:
                        pass
                zenith_add.append(zenith_col)
    print('Writing ----> nsb_sqm_tel.csv')
    with open('{}_sqm_tel.csv'.format(filter_used), 'w') as f:
        header = ('filename, sqm_ut,az,elv,counts,nsb,t_delta\n')
        f.write(header)
        for line in tel_add:
            f.write(line)
    print('Writing ----> nsb_sqm_zenith.csv')
    with open('{}_sqm_zenith.csv'.format(filter_used), 'w') as f:
        header = ('filename,sqm_ut,counts,nsb,t_delta\n')
        f.write(header)
        for line in zenith_add:
            f.write(line)
    print('    Done...')

def calculate_tdelta_bolt(nsb_file,bolt_data, filter_used):
    print('\n### Finding boltwood data for pointings')
    bolt_add = []
    with open(nsb_file, 'r') as file:
        for line in tqdm(file, total=get_num_lines(nsb_file)):

            if line[0] != '#':
                line_split = line.strip('\n').split(',')
                min_delta = 60 # create just some min

                for index, x in enumerate(bolt_data[0]):
                    t = Time(x) + TimeDelta(25200, format='sec')
                    #print(t) # convert UT
                    delta = Time(line_split[1]) - t
                    if abs(float(delta.sec)) < min_delta:
                        min_delta = abs(float(delta.sec))
                        min_timestamp = t
                        #min_counts = sqm_tel_data[1][index]
                        #min_mag = sqm_tel_data[2][index]
                        bolt_col = '{},{},{},{},{},{},{},{},{:.1f}\n'.format(
                        line_split[0],
                        str(min_timestamp)[:19],
                        bolt_data[1][index],
                        bolt_data[2][index],
                        bolt_data[3][index],
                        bolt_data[4][index],
                        bolt_data[5][index],
                        bolt_data[8][index],
                        min_delta
                        )

                    else:
                        pass
                """print('file time = {}'.format(line_split[1]))
                print('bolt time  = {}'.format(min_timestamp))

                print('t delta   = {:.2f}'.format(min_delta))
                print()"""
                bolt_add.append(bolt_col)
    print('    Writing ----> nsb_bolt.csv')

    with open('{}_bolt.csv'.format(filter_used), 'w') as f:
        header = ('filename,bolt_ut,sky_t,amb_t,wind,hum,dew,cloud,t_del\n')
        f.write(header)
        for line in bolt_add:
            f.write(line)
    print('    Done...')
    return bolt_add

if __name__ == '__main__':
    print('\n####### nsb_parser.py')
    parser = argparse.ArgumentParser()

    parser.add_argument('nsb')
    parser.add_argument('bolt')
    parser.add_argument('sqm_telescope')
    parser.add_argument('sqm_zenith')

    args = parser.parse_args()

    nsb_file = args.nsb
    filter_used = nsb_file.split('_')[0]
    bolt_data = parse_bolt(args.bolt)
    sqm_tel_data = parse_sqm(args.sqm_telescope)
    sqm_zenith_data = parse_sqm(args.sqm_zenith)

    calculate_tdelta_sqm(nsb_file, sqm_tel_data, sqm_zenith_data, filter_used)

    calculate_tdelta_bolt(nsb_file, bolt_data, filter_used)

    print('\n####### FINISHED\n')
