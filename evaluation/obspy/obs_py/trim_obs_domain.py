from time import time as timer

from .defns import logger

def trim_obs_domain(df, grid, lat_name='latitude', lon_name='longitude'):
   """
    Remove observations outside of a grid domain to avoid extrapolation
       df: Pandas DataFrame of observations
       grid: datastructure that provides grid latitudes and longitudes

    Currently drops rows, but an alternative option would be to return a mask
   """
   nobs_in = len(df)
   logger.debug(f"grid lat max: {grid.getLatitude().getValue().max()}")
   logger.debug(f"grid lat min: {grid.getLatitude().getValue().max()}")
   logger.debug(f"grid lon max: {grid.getLongitude().getValue().max()}")
   logger.debug(f"grid lon min: {grid.getLongitude().getValue().max()}")

   stime = timer()
   df = df[ df[lat_name] <= grid.getLatitude().getValue().max() ]
   df = df[ df[lat_name] >= grid.getLatitude().getValue().min() ]
   logger.debug(f"latitude check: number of obs reduced from {nobs_in} to {len(df)}")

   nobs_in = len(df)
   df = df[ df[lon_name] <= grid.getLongitude().getValue().max() ]
   df = df[ df[lon_name] >= grid.getLongitude().getValue().min() ]
   logger.debug(f"longitude check: number of obs reduced from {nobs_in} to {len(df)}")
   logger.info(f"number of obs reduced from {nobs_in} to {len(df)}")
   logger.info(f"time: {timer() - stime}")

   return df
