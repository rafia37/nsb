#!/usr/bin/env python3

import argparse
import sys

def parse_boltwood(boltwood_file, log_name):
	"""
	Opens the boltwood file and parses the information.

	Parameters
	----------
	boltwood_file : the txt log file
	log_name  : the log file to write

	Returns
	-------
	None
	"""
	with open(boltwood_file, 'r') as file:
		for line in file:
			split_line = line.split()

			if len(split_line) <= 2:
				pass

			elif split_line[2] == 'M' and split_line[3] == '~D':

				date = split_line[0]
				time = split_line[1]
				cloud_condition = split_line[5]
				wind_condition = split_line[6]
				rain_condition = split_line[7]
				sky_temperature = split_line[10]
				ambient_temperature = split_line[11]
				wind_speed = split_line[12]
				wet_flag = split_line[13]
				rain_flag = split_line[14]
				humidity = split_line[15]
				dew = split_line[16]

				with open(boltwood_log, 'a+') as log:
					line = [
						date, # 0
						time, # 1
						cloud_condition, # 2
						wind_condition, # 3
						wind_speed, # 4
						rain_condition, # 5
						rain_flag, # 6
						wet_flag, # 7
						sky_temperature, # 8
						ambient_temperature, # 9
						humidity, # 10
						dew, # 11
						'\n'
					]
					line = '\t'.join(line)
					log.write(line)

	return None

def parse_arguments():
	"""
	Creates an argument parser for command line interface

	Parameters
	----------
	None

	Returns
	-------
	boltwood_file : the boltwood file read in from the command line
	"""
	parser = argparse.ArgumentParser()

	parser.add_argument('boltwood')
	args = parser.parse_args()

	boltwood_file = args.boltwood

	return boltwood_file

def main():
	"""
	Reads in the command line arguments and parses information from the
	boltwood file. Then writes it to a log in a cleaner and easier to
	use format.
	"""

	boltwood_log = 'boltwood.log'
	boltwood_file = parse_arguments()

	parse_boltwood(boltwood_file, boltwood_log)

	sys.exit(0)

if __name__ == '__main__':
	main()



