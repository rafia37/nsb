for f in all*-w.fits
do
    pp_prepare.py $f
    pp_register.py $f
    pp_photometry.py $f -aprad 1.5
    pp_calibrate.py $f
done
