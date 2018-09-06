for f in *.fits
do
    pp_prepare.py $f
    pp_register.py $f
    pp_photometry.py $f -aprad 1.2 -snr 10
    pp_calibrate.py $f
done
