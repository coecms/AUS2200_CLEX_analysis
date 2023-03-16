import numpy as np

from .defns import OPS, logger
from .dz_def_ok import dz_def_ok


def adjust_w10_by_dz( u10, v10, dz=None, zs=None, znew=None ):
       """
       Adjust 10m winds, u10, v10 [ ms-1 ] to a new height, as in OPS
       https://code.metoffice.gov.uk/doc/ops/ops-2021.03.0/doc/OSDP3.html
       from Howard & Clark (2003, 2007)
 
       Note that as this code takes model to obs (opposite of OPS) the definition
         of dz is reversed - hence the minus signs
       """
 
       dz,dz_check = dz_def_ok( dz, zs, znew )
       s_miss      = -9999999.
       s_miss_test = -9999998.
       logger.debug(f"adj_tol: {OPS.get_w10_adj_tol()}")
       logger.debug(f"dz range: {dz.min()} {dz.max()}")
       scale_fac = np.where( -dz <= OPS.get_w10_adj_tol()[0], OPS.get_w10_adj_min(), s_miss )
       scale_fac = np.where( -dz >= OPS.get_w10_adj_tol()[1], OPS.get_w10_adj_max(), scale_fac )
       # else linear ramp up scaling factor
       scale_fac = np.where( scale_fac < s_miss_test,  \
                                  OPS.get_w10_adj_min()-OPS.get_w10_adj_fac()*dz, scale_fac )
       logger.debug(f"count valid scale_fac: {(scale_fac > s_miss_test).sum()}")
       u10_new   = u10/scale_fac

       if v10 is not None:
           v10_new   = v10/scale_fac
           return u10_new, v10_new
       else:
           return u10_new

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
