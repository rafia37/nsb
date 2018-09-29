"""
Personal Photometry Pipeline Configuation File
2016-11-01, michael.mommert@nau.edu
"""

# Photometry Pipeline
# Copyright (C) 2016  Michael Mommert, michael.mommert@nau.edu

# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.

# You should have received a copy of the GNU General Public License
# along with this program.  If not, see
# <http://www.gnu.org/licenses/>.

##### telescope/instrument configurations

# MYTELESCOPE setup parameters
mytelescope_param = {
    'telescope_instrument': 'Telescope/Instrument',  # telescope/instrument name
    'telescope_keyword': 'mytelescope',  # telescope/instrument keyword
    'observatory_code': '695',  # MPC observatory code
    'secpix': (0.1, 0.1),  # pixel size (arcsec) before binning

    # image orientation preferences
    'flipx': True,
    'flipy': False,
    'rotate': 0,

    # instrument-specific FITS header keywords
    'binning': ('CCDBIN1', 'CCDBIN2'),  # binning in x/y
    'extent': ('NAXIS1', 'NAXIS2'),  # N_pixels in x/y
    'ra': 'RA',  # telescope pointing, RA
    'dec': 'DEC',  # telescope pointin, Dec
    'radec_separator': ':',  # RA/Dec hms separator, use 'XXX'
    # if already in degrees
    'date_keyword': 'DATE-OBS',  # obs date/time
    # keyword; use
    # 'date|time' if
    # separate
    'obsmidtime_jd': 'MJD-OBS',  # obs midtime jd keyword
    # (usually provided by
    # pp_prepare
    'object': 'OBJECT',  # object name keyword
    'filter': 'FILTER',  # filter keyword
    'filter_translations': {'g stuff': 'g'},
    # filtername translation dictionary
    'exptime': 'EXPTIME',  # exposure time keyword (s)
    'airmass': 'AIRMASS',  # airmass keyword

    # source extractor settings
    'source_minarea': 8,  # default sextractor source minimum N_pixels
    'source_snr': 3,  # default sextractor source snr for registration
    'aprad_default': 5,  # default aperture radius in px
    'aprad_range': [2, 10],  # [minimum, maximum] aperture radius (px)
    'sex-config-file': rootpath + '/setup/mytelescope.sex',
    'mask_file': {},
    #                        mask files as a function of x,y binning

    # scamp settings
    'scamp-config-file': rootpath + '/setup/mytelescope.scamp',
    'reg_max_mag'          : 19,
    'reg_search_radius'    : 0.5, # deg
    'source_tolerance': 'high',

    # default catalog settings
    'astrometry_catalogs': ['GAIA'],
    'photometry_catalogs': ['SDSS-R9', 'APASS9', '2MASS']
}


pomenis_param = {
    'telescope_instrument': 'Pomenis/Apogee F9000',  # telescope/instrument name
    'telescope_keyword': 'Pomenis',  # telescope/instrument keyword
    'observatory_code': 'G96',  # MPC observatory code
    'secpix': (4.93, 4.93),  # pixel size (arcsec) before binning

    # image orientation preferences
    'flipx': False,
    'flipy': False,
    'rotate': 180,

    # instrument-specific FITS header keywords
    'binning': ('XBINNING', 'YBINNING'),  # binning in x/y
    'extent': ('NAXIS1', 'NAXIS2'),  # N_pixels in x/y
    'ra': 'RA',  # telescope pointing, RA
    'dec': 'DEC',  # telescope pointin, Dec
    'radec_separator': ':',  # RA/Dec hms separator, use 'XXX'
    # if already in degrees
    'distort': {'PV1_0': -0.0001027723533574,
            'PV1_1': 1.000022031272,
            'PV1_2': -9.204128949252e-05,
            'PV1_4': 9.054368603059e-06,
            'PV1_5': -4.00926084892e-05,
            'PV1_6': 2.544675964472e-05,
            'PV2_0': -7.418342081741e-05,
            'PV2_1': 1.000076590456,
            'PV2_2': 9.939637133313e-05,
            'PV2_4': 5.863234662142e-05,
            'PV2_5': 2.429569111788e-05,
            'PV2_6': 1.683944229526e-05},

    'date_keyword': 'DATE-OBS',  # obs date/time
            # obs date/time
    # keyword; use
    # 'date|time' if
    # separate
    'obsmidtime_jd': 'MJD-OBS',  # obs midtime jd keyword
    # (usually provided by
    # pp_prepare
    'object': 'OBJECT',  # object name keyword
    'filter': 'FILTER',  # filter keyword
    'filter_translations': {'g': 'g',
                            'r': 'r',
                            'v': 'V',
                            'w': 'W',
                            'z': 'z'
                            },
    # filtername translation dictionary
    'exptime': 'EXPOSURE',  # exposure time keyword (s)
    'airmass': 'AIRMASS',  # airmass keyword

    # source extractor settings
    'source_minarea': 3,  # default sextractor source minimum N_pixels
    'source_snr': 50,  # default sextractor source snr for registration
    'aprad_default': 1,  # default aperture radius in px
    'aprad_range': [1, 4],  # [minimum, maximum] aperture radius (px)
    'sex-config-file': rootpath + '/setup/pomenis.sex',
    'mask_file': {},
    #                        mask files as a function of x,y binning

    # scamp settings
    'scamp-config-file': rootpath + '/setup/pomenis.scamp',
    'reg_max_mag'          : 14,
    'reg_search_radius'    : .05, # deg
    'source_tolerance': 'high',

    # default catalog settings
    'astrometry_catalogs': ['TGAS'],
    'photometry_catalogs': ['APASS9']
}

implemented_telescopes.append('Pomenis')

### translate INSTRUME (or others, see _pp_conf.py) header keyword into
#   PP telescope keyword
# example: INSTRUME keyword in header is 'mytel'
instrument_identifiers['Apogee F9000'] = 'Pomenis'

### translate telescope keyword into parameter set defined here
telescope_parameters['Pomenis'] = pomenis_param
##### add telescope configurations to 'official' telescopes.py

implemented_telescopes.append('MYTELESCOPE')

### translate INSTRUME (or others, see _pp_conf.py) header keyword into
#   PP telescope keyword
# example: INSTRUME keyword in header is 'mytel'
instrument_identifiers['instrume_identifier'] = 'MYTELESCOPE'

### translate telescope keyword into parameter set defined here
telescope_parameters['MYTELESCOPE'] = mytelescope_param
