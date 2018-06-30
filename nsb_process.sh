for f in moon*-v.fits
do
    pp_prepare $f
    pp_register $f
    pp_photometry $f -aprad 1
    pp_calibrate $f
done
