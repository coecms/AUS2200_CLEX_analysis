import numpy as np
from .defns import logger
eps     = 1.e-8


def uv_to_ffddd( u, v ):
     """
     Convert wind components to wind speed and direction.
            u,v: wind components in the natural coord system (u to E, v to N)
            Output: wind speed and wind direction (in degrees).
     """
     # if v != 0: #   ddd = atan2(u,v)
     # else
     #   if u < 0  : ddd = 270 (will sobtract 180 so becomes 90)
     #   elif u==0 : ddd = 0
     #   else      : ddd = 90 (will subtract 180 and add 360 ... so becomes 270)

     ddd  = np.where( abs(v) > eps, np.degrees(np.ma.arctan2(u,v)), 270. )
     mask = np.logical_and( abs(v) < eps, abs(u) < eps )
     ddd  = np.where( mask, 0., ddd )
     mask = np.logical_and( abs(v) < eps,  u > eps )
     ddd  = np.where( mask, 90., ddd )

     ddd = np.where( ddd < 0., ddd+360., ddd )

     ddd = ddd - 180.
     ddd = np.where( ddd < 0., ddd+360., ddd )

     # NaN is not handled correctly
     ddd = np.where( np.logical_or(np.isnan(u), np.isnan(v)), np.NaN, ddd )

     ff = np.ma.sqrt(u*u+v*v)
      
     return ff,ddd

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

def uv_to_ff( u, v ):
    ff,ddd =  uv_to_ffddd( u, v )
    return ff

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

def uv_to_ddd( u, v ):
    ff,ddd =  uv_to_ffddd( u, v )
    return ddd

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
