from util_helpers import is_numpy_array

from .defns import logger
from .Set_Latind import Set_Latind
from .Set_Lonind import Set_Lonind
from .check_grid import check_grid
from .h_interp import h_interp
from .t_interp import t_interp


def interp( ob_df, grid, gfile_info, valid_dt, \
               ls_mask=None, extrap=True, allow_all_sea=False ):
       """
       Interpolate grid to obs point
       Adjust for land/sea differences if model land/sea mask available
       ob_df: pandas dataframe of all obs data
       grid:  grid to be interpolated 
       """

       check_grid( grid )
       ls_mask_input = is_numpy_array( ls_mask )

       obs_ind = list(range(len(ob_df)))
       # set up lat & lon  indices 
       ind_lat, lat_wt = Set_Latind( grid, ob_df['latitude'].values )
       ind_lon, lon_wt = Set_Lonind( grid, ob_df['longitude'].values )

       if not extrap :
          eps = 1.e-2
          mask_lt = np.logical_or( lat_wt< -eps, lat_wt > 1+eps )
          mask_ln = np.logical_or( lon_wt< -eps, lon_wt > 1+eps )
          mask_ll = np.logical_or( mask_lt, mask_ln )
          lat_wt  = np.ma.array( lat_wt, mask=mask_ll )
          lon_wt  = np.ma.array( lon_wt, mask=mask_ll )

       logger.debug(f'ind_lat: {ind_lat}')
       logger.debug(f'lat_wt: {lat_wt}')
       logger.debug(f'ind_lon: {ind_lon}')
       logger.debug(f'lon_wt: {lon_wt}')

       ob_val, ob_lsm  = h_interp( grid, ind_lat, ind_lon, \
                                   obs_ind, lat_wt, lon_wt, \
                                   ls_mask=ls_mask,
                                   allow_all_sea=allow_all_sea)

       if not extrap :
          ndim = np.array (ob_val[0].shape ).prod()
          ob_mask = np.array( [mask_ll]*ndim ).reshape( ob_val.shape )
          ob_val  = np.ma.array( ob_val, mask=ob_mask )

       logger.debug(f'grid.shape: {grid.shape}')

       if 't' in grid.getOrder():
           ob_val   = t_interp( ob_val, gfile_info, valid_dt, ob_df )
#      if 'z' in grid.getOrder():
#          ob_val   = v_interp( ob_val, gfile_info, vert_lev, ob_df )
       
#      remove any residual degenerate dimensions
       ob_val = ob_val.squeeze()

       if ls_mask_input:
           return ob_val, ob_lsm
       else :
           return ob_val


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
