import numpy as np
from .defns import logger

# from https://stackoverflow.com/questions/67509913/add-an-attribute-to-a-numpy-array-in-runtime

class mod_cdms_grid(np.ndarray):
    def __new__(cls, input_array, getLatitude=None, getLongitude=None, 
                    _order_=None):
       obj = np.asarray(input_array).view(cls)
       obj.getLatitude = getLatitude
       obj.getLongitude = getLongitude
       obj._order_ = _order_
       return obj
    def __array_finalize__(self, obj):
       if obj is None: return
       self.getLatitude = getattr(obj,'getLatitude',None)
       self.getLongitude = getattr(obj,'getLongitude',None)
       self._order_ = getattr(obj,'_order_',None)
    def getOrder(self):
       return self._order_

def transform_grid(grid, var_info, grid_dt, ndim=4, valid_dt=None):
       """
       apply various transforms to the grid
           - reorder grid into appropriate order (lat,lon,time,level)
           - apply any scalings etc.

       var_info is the dictionary from the relevant instance of archive.py, 
              e.g. archive.file_info()[variable]

       ndim: max number of dimensions of output grid (can truncate to 2d)
       valid_dt: list of valid datetimes of interest (could be obs or grid)
           - will trim down if warranted
            returns extra date time at start and end if available
       """

       logger.debug(f"var_info: {var_info}")
       axes = [ax.id for ax in grid.getAxisList()]
       order_str = 'yxtz'
       if 'level' not in axes: order_str = order_str.replace('z',"")
       if 'time' not in axes: order_str = order_str.replace('t',"")
       if ndim == 2: order_str = order_str[:2]
       logger.debug(f"order_str: {order_str}")

       new_grid = mod_cdms_grid(grid.reorder(order_str))
       new_grid.getLatitude = grid.getLatitude
       new_grid.getLongitude = grid.getLongitude
       new_grid._order_ = order_str

       if (ndim == 2) and (len(grid.shape) > 2):
            # this should handle ERA-5 data where topo etc. could be 4d
            ldim0 = new_grid.shape[0]
            ldim1 = new_grid.shape[1]
            if len(grid.shape) == 3:
              new_grid = new_grid[:,:,0]
            elif len(grid.shape) == 4:
              new_grid = new_grid[:,:,0,0]
       elif (new_grid.shape[-1] == 1) and (len(new_grid.shape) == 4): 
            # e.g. file with just 1st model level 
            new_grid = new_grid[:,:,:,0]

       # if ndim == 2 then time dimension has just been removed 
       #    ..and grid_dt will not be accessed by interp
       if ("t" in order_str) and (ndim > 2):
          # grid_dt >=  test should give first index > min valid_dt
          # start seems to be out by one for COtL interpolation test
          #    - not sure why, and will take quite some time to work out 
          #    - adding one extra time field should not have significant impact
          i_start = max(0, (grid_dt >= min(valid_dt)).argmax()-1)
          if min(valid_dt) < grid_dt[0]:
             logger.warning(f"valid datetime out of range {min(valid_dt)}")
             logger.warning(f"grid_dt min, max: {grid_dt[0]}, {grid_dt[-1]}")
             # grid_dt >= min_valid_dt will be True everywhere, argmax = 0
             #    so works by accident. Just to make safe!
             i_start = 0

          i_end   = min(len(grid_dt), (grid_dt >= max(valid_dt)).argmax()+1)
          if max(valid_dt) > grid_dt[-1]:
             logger.warning(f"valid datetime out of range {max(valid_dt)}")
             logger.warning(f"grid_dt min, max: {grid_dt[0]}, {grid_dt[-1]}")
             # grid_dt <= max_valid_dt will be True everywhere, argmax = 0
             i_end = len(grid_dt)

          new_grid = new_grid[:,:,i_start: i_end]
          newgrid_dt = grid_dt[i_start: i_end]
          logger.info(f" i_start: {i_start}, {grid_dt[i_start]}")
          # last element used in a slice is i_end-1
          #  if this is the end of list then it is wrong. Tough - it's Python.
          logger.info(f" i_end: {i_end}, {grid_dt[i_end-1]}")
          logger.info(f" valid_dt: {min(valid_dt)}, {max(valid_dt)}")
       else:
          newgrid_dt = grid_dt

       if 'rescale_mul' in var_info: 
            new_grid = new_grid*var_info['rescale_mul']
       if 'rescale_div' in var_info: 
            new_grid = new_grid/var_info['rescale_div']
       if 'rescale_add' in var_info: 
            new_grid = new_grid+var_info['rescale_add']
       if 'rescale_sub' in var_info: 
            new_grid = new_grid-var_info['rescale_sub']

       return new_grid, newgrid_dt
         
         
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
