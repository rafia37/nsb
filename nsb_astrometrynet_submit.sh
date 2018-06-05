"""
Shell script that submits each file to astrometry.net to find the center
RA and DEC.
"""
printf 'Sending images to astrometry.net'
for f in *.fits
do
    python2 nsb_astrometry.py --upload $f --wait
done
