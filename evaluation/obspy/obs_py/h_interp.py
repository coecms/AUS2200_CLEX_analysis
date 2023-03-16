from numpy import array as np_array
from numpy import where as np_where

from .defns import logger
from .h_interp_getpoints import h_interp_getpoints

def h_interp( field, ind_lat, ind_lon, obs_ind, lat_wt, lon_wt, \
              ls_mask=None, \
              grid00=None, grid01=None, grid10=None, grid11=None, \
              lsm00=None, lsm01=None, lsm10=None, lsm11=None, \
              allow_all_sea=False ) :
       """
       Horizontal interpolation from grid to obs point
         field = grid field(lat, lon, time, level)
         ind_lat = index of latitude before obs
         ind_lon = index of longitude before obs
         obs_ind = list of indices in ind_lat, ind_lon
         lat_wt  = latitude weights
         lon_wt  = longitude weights
         grid00/lsm00 etc are the surrounding points (field & land_sea_mask)
             none or all must be None
         allow_all_sea : option to allow points to surrounded by sea to be used
       """
       eps       = 1.e-6

       if grid00 is None:
          grid00, grid01, grid10, grid11,  \
          lsm00, lsm01, lsm10, lsm11 =  \
            h_interp_getpoints( field, ind_lat, ind_lon, obs_ind, \
                                ls_mask=ls_mask )
     
       wt00 = (1. - lon_wt)*(1. - lat_wt)
       wt01 = (     lon_wt)*(1. - lat_wt)
       wt10 = (1. - lon_wt)*(     lat_wt)
       wt11 = (     lon_wt)*(     lat_wt)
       # adjust weights & renormalize
       if not(lsm00 is None):
             if allow_all_sea:
                all_sea = np_where( (lsm00+lsm01+lsm10+lsm11) < eps, True, False )
             else:
                all_sea = np_array( [False]*len(wt00) )
     
             #Multiply weights by land sea mask if at least one point is land
             #   If no land points then use original weight (if allow_all_sea)
             #   Else scale weights by land sea mask
             wtl00 = np_where( all_sea, wt00, wt00*lsm00 )
             wtl01 = np_where( all_sea, wt01, wt01*lsm01 )
             wtl10 = np_where( all_sea, wt10, wt10*lsm10 )
             wtl11 = np_where( all_sea, wt11, wt11*lsm11 )
             sum_wtl = wtl00 + wtl01 + wtl10 + wtl11
             
             # catch the case where obs is between two sea points, and so 
             #  the only land points nearby have zero weight
             # e.g. one point is land, but it has zero weight
             # If have thin finger of land - just take the original weights
             zero_land = np_where( sum_wtl < eps, True, False )
             wt00 = np_where( zero_land, wt00, wtl00 )
             wt01 = np_where( zero_land, wt01, wtl01 )
             wt10 = np_where( zero_land, wt10, wtl10 )
             wt11 = np_where( zero_land, wt11, wtl11 )
             sum_wt = wt00 + wt01 + wt10 + wt11

             logger.debug('sum_wt: {0}\n'.format(sum_wt))
             logger.debug('sum_wt min,arg: {0}\n'.format(sum_wt.min(),sum_wt.argmin() ))
             logger.debug('sum_wt max,arg: {0}\n'.format(sum_wt.max(),sum_wt.argmax() ))
             logger.debug('wt00: {0}\n'.format(wt00))
             logger.debug('wt01: {0}\n'.format(wt01))
             logger.debug('wt10: {0}\n'.format(wt10))
             logger.debug('wt11: {0}\n'.format(wt11))
        
             if (sum_wt < 1.e-6).any():
               logger.error(' All weights == 0 for a point,'+ \
                    ' allow_all_sea: {0}\n'.format(allow_all_sea))
               logger.error(' Point index: {0}  &shape: {1}\n'.\
                       format(sum_wt.argmin(),sum_wt.shape))
               raise ValueError

             wt00 = wt00 / sum_wt
             wt01 = wt01 / sum_wt
             wt10 = wt10 / sum_wt
             wt11 = wt11 / sum_wt

       # interpolate to obs location
       # works as long as regular longitude grid.

       npts = list(range( len( wt00 )))

       # grid00 etc can have time dimension so need 
       #    need to drop wt00 etc back to scalar and multiply
       #    each grid00[i] etc. which has dimension time

       ob_gval = [wt00[i]*grid00[i] + wt01[i]*grid01[i] + \
                  wt10[i]*grid10[i] + wt11[i]*grid11[i] for i in npts ]
       ob_gval = np_array( ob_gval )
       ob_lsm  = np_array( [ lsm00, lsm01, lsm10, lsm11 ] )

       logger.debug(f'ob_gval.shape {ob_gval.shape}')
       logger.debug(f'grid00.shape {grid00.shape},{grid11.shape}')
       logger.debug(f'ob_gval {ob_gval}')
       logger.debug(f'grid00  {grid00}')
       logger.debug(f'grid01  {grid01}')
       logger.debug(f'grid10  {grid10}')
       logger.debug(f'grid11  {grid11}')
       logger.debug(f'ind_lat {ind_lat}')
       logger.debug(f'ind_lon {ind_lon}')

       return ob_gval, ob_lsm

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
