import unittest

import numpy as np
from zoneinfo import ZoneInfo
from datetime import datetime as dt

zone_utc = ZoneInfo('UTC')

def local_datetime_utc_time (local_dt, obs_tz, odt_fmt):
    # convert from local time to UTC 
    # local_dt : list of datetimes for each obs (datetime objects)
    # obs_tz : list of time zones for each obs (ZoneInfo objects)

    #  replace(tzinfo=) puts original time as local
    #  astimezone(..) converts to another timezone
    tt = [dd.replace(tzinfo=stz).astimezone(zone_utc) \
              for dd,stz in zip(local_dt,obs_tz) ]

    obtime = np.array([ int(dt.strftime(dd,odt_fmt)) for dd in tt] )
    return obdate, obtime

class test_local_datetime_utc_time(unittest.TestCase):
   def test_rtn(self):
      local_dt = [dt.strptime('202101040645','%Y%m%d%H%M'), 
                  dt.strptime('202106050615','%Y%m%d%H%M') ]
      obs_tz = [ZoneInfo('Australia/South')]*2
      odt_fmt = dt_fmt()
      ot = local_datetime_utc_time(local_dt, obs_tz, odt_fmt.time_fmt() )
      self.assertEqual( ot, ['201500','204500'] )
 
       
if __name__ == '__main__':
   unittest.main()
