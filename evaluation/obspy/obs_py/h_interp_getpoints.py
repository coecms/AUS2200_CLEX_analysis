from numpy import array as np_array

from .defns import logger


def h_interp_getpoints( field, ind_lat, ind_lon, obs_ind,  \
                        ls_mask=None ) :
       """
         Returns grid values for linear interpolation
            and land sea mask values
       """

       obs_vec = np_array(obs_ind)
       logger.debug(f'ind_lat max,min: {ind_lat.max()} {ind_lat.min()}')
       logger.debug(f'ind_lon max,min: {ind_lon.max()} {ind_lon.min()}')
       logger.debug(f'obs_ind max,min: {obs_vec.max()} {obs_vec.min()}')
       logger.debug(f'field.shape {field.shape}')
       if ls_mask is None:
           logger.debug(f'ls_mask.shape None')
       else:
           logger.debug(f'ls_mask.shape {ls_mask.shape}')

       # Don't know why, but CDMS2 does not allow indices of numpy.int64
       i_lt  = [int(ind_lat[k]) for k in obs_ind]
       j_ln  = [int(ind_lon[k]) for k in obs_ind]

       # No need to redo the zip command 8 times. Just do it once
       zip_ij = list(zip(i_lt, j_ln))

       grid00 = np_array( [ field[i  ,j  ] for i,j in zip_ij ] )
       grid01 = np_array( [ field[i  ,j-1] for i,j in zip_ij ] )
       grid10 = np_array( [ field[i-1,j  ] for i,j in zip_ij ] )
       grid11 = np_array( [ field[i-1,j-1] for i,j in zip_ij ] )

       # later on there is an assumed time dimension, even if degenerate
       if len(grid00.shape) == 1 :
          grid00 = grid00.reshape( len(grid00), 1)
          grid01 = grid01.reshape( len(grid01), 1)
          grid10 = grid10.reshape( len(grid10), 1)
          grid11 = grid11.reshape( len(grid11), 1)

       if ls_mask is None:
          lsm00 = None
          lsm01 = None
          lsm10 = None
          lsm11 = None

       else:
          lsm00 = np_array( [ ls_mask[i  ,j  ] for i,j in zip_ij ] )
          lsm01 = np_array( [ ls_mask[i  ,j-1] for i,j in zip_ij ] )
          lsm10 = np_array( [ ls_mask[i-1,j  ] for i,j in zip_ij ] )
          lsm11 = np_array( [ ls_mask[i-1,j-1] for i,j in zip_ij ] )

          logger.debug(f'lsm00.shape {lsm00.shape}')
          logger.debug(f'lsm01.shape {lsm01.shape}')
          logger.debug(f'lsm10.shape {lsm10.shape}')
          logger.debug(f'lsm11.shape {lsm11.shape}')
    
       return grid00, grid01, grid10, grid11, lsm00, lsm01, lsm10, lsm11
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
