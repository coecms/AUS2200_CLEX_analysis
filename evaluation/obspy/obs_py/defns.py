from pathlib import Path

from util_helpers import logging_x as logging
from nwp_helpers import ops_params
logger = logging.getLogger('obs_py')

PKG_PATH = Path(__file__).parent.absolute()

code_version = "4.0.15"
institution  = 'BoM'

from .dt_fmt import dt_fmt

# since dtf_grid defined here (used in multiple routines)
# define dtf_obs for consistency
dtf_obs    = dt_fmt(date_fmt='%Y%m%d', time_fmt='%H%M%S')
dtf_grid   = dt_fmt(date_fmt='%Y%m%d', time_fmt='%H%M')
dtf_create = dt_fmt(date_fmt='%A %d-%B-%Y', time_fmt=' %H:%M:%S')

OPS = ops_params()

# field_dep is a list of lists as there are multiple options for
#    some variables
# this depends on which conversions are coded, and how
#    so source code dependent, and hence a constant
field_dep = { \
     'dew_point' : [['temperature','relative_humidity'], \
                    ['pressure','specific_humidity']], \
     'relative_humidity' : [['pressure','temperature','dew_point'], \
                            ['pressure','temperature','specific_humidity']], \
     'specific_humidity' : [['pressure','dew_point'], \
                            ['pressure','temperature','relative_humidity']], \
     'zonal_wind' : [['wind_speed','wind_direction']], \
     'meridional_wind' : [['wind_speed','wind_direction']], \
     'wind_speed' : [['zonal_wind', 'meridional_wind']], \
     'wind_direction' : [['zonal_wind', 'meridional_wind']] }

# ht_corr_dep also depends on which conversions are coded, and how
#    so source code dependent, and hence a constant
ht_corr_dep = { \
     'screen_pressure' : ['surface_height', \
                          'obs_surface_height', \
                          'obs_screen_temperature', \
                          'obs_screen_dew_point', 'screen_pressure'] , \
     '10m_pressure' : ['screen_temperature', \
                       'screen_dew_point','screen_pressure'] , \
     'screen_zonal_wind' : ['10m_zonal_wind','z_lev','z_lev_u','roughness_length'], \
     'screen_meridional_wind' : ['10m_meridional_wind','z_lev','z_lev_u','roughness_length'], \
     'screen_windgust' : ['10m_windgust','z_lev','z_lev_u','roughness_length'], \
     'screen_wind_speed' : ['10m_wind_speed','z_lev','z_lev_u','roughness_length'] }




# field_dep is a list of lists as there are multiple options for
#    some variables
# this depends on which conversions are coded, and how
#    so source code dependent, and hence a constant
