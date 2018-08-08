"""
Shell script that submits each file to astrometry.net to find the center
RA and DEC.
"""
printf 'Sending images to astrometry.net'
for f in full*.fits
do
    python2 nsb_astrometrynet_client.py --upload $f --wait --apikey yiplytzetqmxeimx
done
