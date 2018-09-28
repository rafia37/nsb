#!/usr/bin/env python3
import subprocess
import argparse
import os, sys

def pp_run(bash_command):
    """

    """
    process = subprocess.Popen(bash_command.split(), stdout=subprocess.PIPE)
    output, error = process.communicate()

    print(output)

    return None

def _parse_arguments():
    parser = argparse.ArgumentParser()
    parser.add_argument(dest='files', nargs='+', help='files to run PP on')
    parser.add_argument(
        '-a',
        '--aperture_radius'
    )
    parser.add_argument(
        '-s',
        '--snr',
        type=float
    )

    args = parser.parse_args()

    return args.files, args.aperture_radius, args.snr

def main():
    """

    """
    print('\n########## nsb_process.py')
    print('Runs images through PHOTOMETRYPIPELINE with optimal settings')
    files, aperture_radius, snr = _parse_arguments()

    if not aperture_radius:
        aperture_radius = 1.2

    if not snr:
        snr = 10

    for file in files:
        commands = [
        'pp_prepare.py {}'.format(file),
        'pp_register.py {}'.format(file),
        'pp_photometry.py {} -aprad {} -snr {}'.format(
            file,
            aperture_radius,
            snr),
        'pp_calibrate.py {}'.format(file)
        ]
        for command in commands:
            print('\n#----------{}'.format(command.split()[0]))
            pp_run(command)
            print('\nDone...')

        print('\nFINISHED\n')

    sys.exit(0)

if __name__ == '__main__':
    main()
