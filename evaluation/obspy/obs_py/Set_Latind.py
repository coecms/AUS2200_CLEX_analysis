from util_helpers import vmax, vmin

from .defns import logger

def Set_Latind( grid, ob_val, debug=None ) :

       """
       Returns indices and weights for horizontal interpolation
          from grid to array of obs latitude (ob_val)
       Note that grid latitudes may be 90 to -90 so need to 
          flip (multiply by -1)
          Therefore need to do the same to the obs latitudes.
       """

       grid_val = grid.getLatitude().getValue()
       if grid_val[0] < grid_val[-1]:
          ind_val  = grid_val.searchsorted(ob_val[:])
       else:
          ind_val  = (-grid_val).searchsorted(-ob_val[:])
       ind_val  = vmax( ind_val, 1 )
       ind_val  = vmin( ind_val, len(grid_val)-1 )

       dval     = ( grid_val[ ind_val ] - ob_val[:] )  \
                   / ( grid_val[ ind_val ] - grid_val[ ind_val-1 ] )

       logger.debug(f'glat range: {grid_val[0]} {grid_val[-1]} {len(grid_val)}\n')
       logger.debug(f'grid_val: {grid_val[ind_val]}\n')
       logger.debug(f'grid_val-1: {grid_val[ind_val-1]}\n')
       logger.debug(f'ob_val: {ob_val}\n')
       logger.debug(f'ind_val: {ind_val}\n')
       logger.debug(f'dval: {dval}\n')

       return ind_val, dval

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
