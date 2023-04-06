from util_helpers import vmax, vmin
from .defns import logger

def Set_Lonind( grid, ob_val,  debug=None ) :

       """
       Returns indices and weights for horizontal interpolation
          from grid to obs longitude(ob_val)
       """

       grid_val = grid.getLongitude().getValue()
       ind_val  = grid_val.searchsorted(ob_val[:])
       ind_val  = vmax( ind_val, 1 )
       ind_val  = vmin( ind_val, len(grid_val)-1 ) 

       dval     = ( grid_val[ ind_val ] - ob_val[:] )  \
                   / ( grid_val[ ind_val ] - grid_val[ ind_val-1 ] )

       logger.debug(f'glon range: {grid_val[0]} {grid_val[-1]} {len(grid_val)}\n')
       logger.debug(f'grid_val: {grid_val[ind_val]}\n')
       logger.debug(f'grid_val-1: {grid_val[ind_val-1]}\n')
       logger.debug(f'ob_val: {ob_val}\n')
       logger.debug(f'ind_val: {ind_val}\n')
       logger.debug(f'dval: {dval}\n')

       return ind_val, dval
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
