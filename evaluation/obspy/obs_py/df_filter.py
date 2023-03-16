import numpy as np

from .defns import logger, dtf_obs

def df_filter(df, lat_range=None, lon_range=None, \
                  dt_range=None, date_range=None, time_range=None, 
                  dt_str_range=None, id_list=None):
   """
   Filter out a pandas DataFrame (df) according to various criteria
     - the obvious choices for obs_py are:
         lat_range: [min_lat, max_lat] for latitude
         lon_range: [min_lon, max_lon] for longitude
         dt: [min_dt, max_dt] for datetime objects
         dt_str: [min_dt, max_dt] for datetime strings
         date: [min_date, max_date] for date
         time: [min_time, max_time] for time
               - can handle [235900, 100]
         id_list: [id1, id2, ... idn] for station identifiers

   Since have to faff around with numpy masks easier to contain in routine
   """

   logger.debug(f"initial len(df): {len(df)}")
   if (lat_range is not None) and (len(lat_range) > 1) and \
         ('latitude' in df) :
      mask_ge = np.greater_equal( df['latitude'], lat_range[0])
      mask_le = np.less_equal( df['latitude'], lat_range[1])
      mask = np.logical_and(mask_ge, mask_le)
      df = df[mask]
   logger.debug(f"after lat thinning len(df): {len(df)}")
   logger.debug(f"max: {df['latitude'].max()}, min: {df['latitude'].min()}")

   if (lon_range is not None) and (len(lon_range) > 1) and \
         ('longitude' in df) :
      if (lon_range[0] >= 0.) and (lon_range[1] <= 360.):
          mask_ge = np.greater_equal( df['longitude'], lon_range[0])
          mask_le = np.less_equal( df['longitude'], lon_range[1])
          mask = np.logical_and(mask_ge, mask_le)
          df = df[mask]
      elif (lon_range[0] < 0.) and (lon_range[1] <= 360.):
          mask_ge = np.greater_equal( df['longitude'], lon_range[0]+360)
          mask_le = np.less_equal( df['longitude'], lon_range[1])
          mask = np.logical_or(mask_ge, mask_le)
          df = df[mask]
      else:
          mask_ge = np.greater_equal( df['longitude'], lon_range[0])
          mask_le = np.less_equal( df['longitude'], lon_range[1]-360)
          mask = np.logical_or(mask_ge, mask_le)
          df = df[mask]
   logger.debug(f"after lon thinning len(df): {len(df)}")
   logger.debug(f"max: {df['longitude'].max()}, min: {df['longitude'].min()}")

   if (date_range is not None) and (len(date_range) > 1):
      mask_ge = np.greater_equal( df['date'], date_range[0])
      mask_le = np.less_equal( df['date'], date_range[1])
      mask = np.logical_and(mask_ge, mask_le)
      df = df[mask]
   logger.debug(f"after date thinning len(df): {len(df)}")
   logger.debug(f"max: {df['date'].max()}, min: {df['date'].min()}")

   if (time_range is not None) and (len(time_range) > 1):
      if time_range[0] < time_range[1]:
         mask_ge = np.greater_equal( df['time'], time_range[0])
         mask_le = np.less_equal( df['time'], time_range[1])
         mask = np.logical_and(mask_ge, mask_le)
         df = df[mask]
      else:
         # handle [235900,100] so can span 00Z
         mask_le = np.less_equal( df['time'], time_range[1])
         mask_ge = np.greater_equal( df['time'], time_range[0])
         mask = np.logical_or(mask_ge, mask_le)
         df = df[mask]
   logger.debug(f"after time thinning len(df): {len(df)}")
   logger.debug(f"max: {df['time'].max()}, min: {df['time'].min()}")

   if (dt_range is not None) and (len(dt_range) > 1):
      int_dt = dtf_obs.int_to_dt(df['date'].values, df['time'].values)
      mask_ge = np.greater_equal( int_dt, dt_range[0])
      mask_le = np.less_equal( int_dt, dt_range[1])
      mask = np.logical_and(mask_ge, mask_le)
      df = df[mask]
   logger.debug(f"after datetime thinning len(df): {len(df)}")

   if (dt_str_range is not None) and (len(dt_str_range) > 1):
      int_dt = dtf_obs.int_to_dt(df['date'].values, df['time'].values)
      str_dt = np.array([ str(dd) for dd in int_dt ])
      mask_ge = np.greater_equal( int_dt, dt_str_range[0])
      mask_le = np.less_equal( int_dt, dt_str_range[1])
      mask = np.logical_and(mask_ge, mask_le)
      df = df[mask]
   logger.debug(f"after dt_str thinning len(df): {len(df)}")

   if (id_list is not None) and (len(id_list) > 0):
       df = df[ np.in1d(df['station_identifier'], id_list ) ]
   logger.debug(f"after stn thinning len(df): {len(df)}")
   logger.debug(f"max: {df['station_identifier'].max()}, min:{df['station_identifier'].min()}")

   return df
