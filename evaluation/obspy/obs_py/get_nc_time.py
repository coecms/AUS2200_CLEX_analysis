from datetime import datetime as dt
from datetime import timedelta as del_t
import numpy as np

from .defns import logger

def get_nc_time(nc_time, ref_fmt='%Y-%m-%d %H:%M:%S', source=None):
    """
    Extracts time from netcdf variables
      - returns a numpy array of datetime objects
    Can pull out the date from
       1. source='units' the units string
             (this will be the reference time, which can be the basedate for a forecast file), or 
       2. add in the values from the CDMS variable nc_time
            source = anything else will just read from the variable

     nc_time: CDMS variable, e.g. cdmsfile[time_var] where time_var is the
                 netcdf name for the time variable
     ref_fmt: format for reference times
     source: string
                'units' : return reference date
                anything else: return values from nc_time with reference time
                                  included
    """

    ref_str = nc_time.units.split(' ')
    if ref_str[1] != 'since':
        logger.error(f"Unknown time description: {nc_time.units}")
        logger.error(f"required format: <time_periods> since <date>")
        raise ValueError

    ref_units = ref_str[0]

    dt_str = ' '.join(ref_str[2:])
    ref_date = dt.strptime(dt_str,ref_fmt)

    if source == 'units': return np.array([ref_date])

    # safest to use float 64 as may convert long periods to seconds
    # Just need to deal with oe of the annoying differences between CDMS
    #   axes and variables
    if hasattr(nc_time,'asdatetime'):
        dates = np.array( nc_time.asdatetime() )
        date_inc = None
    elif hasattr(nc_time,'axis'):
        date_inc = np.array([nc_time.getValue().astype(np.float64)])
    else:
        date_inc = np.array([nc_time.astype(np.float64)])

    if date_inc is not None:
       date_inc = date_inc.flatten()

       logger.debug(f"date_inc: {date_inc}")
       if ref_units == 'months':
          date_inc = [del_t(months=float(tt)) for tt in date_inc]
       elif ref_units == 'weeks':
          date_inc = [del_t(weeks=float(tt)) for tt in date_inc]
       elif ref_units == 'days':
          date_inc = [del_t(days=float(tt)) for tt in date_inc]
       elif ref_units == 'hours':
          date_inc = [del_t(hours=float(tt)) for tt in date_inc]
       elif ref_units == 'minutes':
          date_inc = [del_t(minutes=float(tt)) for tt in date_inc]
       elif ref_units == 'seconds':
          date_inc = [del_t(seconds=float(tt)) for tt in date_inc]
       elif ref_units == 'milliseconds':
          date_inc = [del_t(milliseconds=float(tt)) for tt in date_inc]
       else:
          logger.error(f"Unknown time period: {ref_units}")
          raise ValueError

       dates = np.array([ref_date + tt for tt in date_inc])
    return dates
