import unittest

import numpy as np
from zoneinfo import ZoneInfo
from datetime import datetime as dt

zone_utc = ZoneInfo('UTC')

def local_datetime_utc_date (local_dt, obs_tz, odt_fmt):
    # convert from local time to UTC 
    # local_dt : list of datetimes for each obs (datetime objects)
    # obs_tz : list of time zones for each obs (ZoneInfo objects)

    #  replace(tzinfo=) puts original time as local
    #  astimezone(..) converts to another timezone
    tt = [dd.replace(tzinfo=stz).astimezone(zone_utc) \
              for dd,stz in zip(local_dt,obs_tz) ]

    obdate = np.array([ int(dt.strftime(dd,odt_fmt)) for dd in tt] )
    return obdate

class test_local_datetime_utc_date(unittest.TestCase):
   def test_rtn(self):
      local_dt = [dt.strptime('202101040645','%Y%m%d%H%M'), 
                  dt.strptime('202106050615','%Y%m%d%H%M') ]
      obs_tz = [ZoneInfo('Australia/South')]*2
      odt_fmt = dt_fmt()
      od = local_datetime_utc_date(local_dt, obs_tz, odt_fmt.date_fmt() )
      self.assertEqual( od, ['20210103','20210604'] )
 
       
if __name__ == '__main__':
   unittest.main()
