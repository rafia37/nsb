#Just a small shell wrapper to iterate over the fits files. 
#Ex. nsb_process.sh kp*.fits

for file in "$@"
	
	do
	    pp_prepare.py "$file"
	    pp_register.py "$file"
	    pp_photometry.py "$file" -aprad 1.2 -snr 10
	    pp_calibrate.py "$file"
	done
