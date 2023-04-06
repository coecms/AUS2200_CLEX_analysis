from .defns import logger

def order_grid(grid, ndim=4):
       """
       Put grid into appropriate order (lat,lon,time,level)
       """
      
       logger.warning("order_grid is to be deprecated, use transform_grid")

       if grid.getLevel() is None:
         order_str = 'yxt'
       else:
         order_str = 'yxtz'
       new_grid = grid.reorder(order_str)

       if (ndim == 2) and (len(grid.shape) > 2):
            new_grid = new_grid[:,:,0]
       elif (new_grid.shape[-1] == 1) and (len(new_grid.shape) == 4): 
            # e.g. file with just 1st model level 
            new_grid = new_grid[:,:,:,0]

       return new_grid
         
         
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
