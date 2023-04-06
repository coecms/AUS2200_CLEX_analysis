from numpy import array as np_array
from numpy import where as np_where
from numpy import ones as np_ones
from numpy import nan as np_nan
from numpy import arange as np_arange
from datetime import datetime as dt

from .defns import logger
from .dbg_parm import dbg_parm

def t_interp( ob_grid, gfile_info, grid_dt, o_df,
              time_wt=None ) :
       """
       Time interpolation from grid to obs point
         (This allows for different formats for obs & grid date times)
         ob_grid   = grid field(lat, lon, time, level) interpolated to obs(lat,lon)
         grid_dt   = grid validity datetime object/string
         o_df      = obs dataframe
         time_wt   = gridpoint weights
 
         Assumes grid_dt and ob_dt are similar types of variables
       """
       nprt = dbg_parm.get_dbg_parm('np_dbg_len')

       lt = len(gfile_info.var_dict()['time']['attributes']['units'])
       ob_dt = [dt.strptime('{}{}'.format(od,str(ot).zfill(lt)), 
                 gfile_info.var_dict()['date']['attributes']['units']+
                 gfile_info.var_dict()['time']['attributes']['units'] ) \
                 for od,ot in zip(o_df['date'].values,o_df['time'].values) ]
       
       ind_time = grid_dt.searchsorted(ob_dt)
       logger.debug(f'ind_time: {ind_time}')

       return_wt = (time_wt == [])

       logger.debug(f'grid_dt: {grid_dt}')
       logger.debug(f'ob_dt: {ob_dt}')

       if len(grid_dt) == 1 :
          if return_wt:
             return ob_grid, [1.]*len(o_df)
          else: 
             return ob_grid

       logger.debug(f'ob_dt[:10]: {ob_dt[:nprt]}')

#      Mask for observations beyond time window defined by grids
       ob_dt_arr = np_array(ob_dt)
       fc_mask   = np_where(ob_dt_arr > grid_dt[-1],False,True)
       
       if (time_wt is None) or return_wt :
          gdt  = list(grid_dt[ind_time[fc_mask]])
          gdt1 = list(grid_dt[ind_time[fc_mask]-1])
          
          time_wt = [  (g-o).total_seconds() \
                      /(g-g1).total_seconds() for g,g1,o in zip(gdt,gdt1,ob_dt_arr[fc_mask]) ]
   
          time_wt  = np_array(time_wt)

       logger.debug(f'time_wt: {time_wt.max()}, {time_wt.min()}')
       gti_max = np_array(ind_time[fc_mask]).max()
       gti_min = np_array(ind_time[fc_mask]).min()
       logger.debug(f'ind_time: {gti_max}, {gti_min}')
       logger.debug(f'grid_dt max: {grid_dt[gti_max]}, {grid_dt[gti_max-1]}')
       logger.debug(f'grid_dt min: {grid_dt[gti_min]}, {grid_dt[gti_min-1]}')
       logger.debug(f'valid fclens: {fc_mask.sum()}')
       logger.debug(f'time_wt (shape): {time_wt.shape}')
       logger.debug(f'time_wt: {time_wt}')
       logger.debug(f'obs date/time: {ob_dt}')
       logger.debug(f'grid date/time: {grid_dt[ind_time[fc_mask]]}')
       logger.debug(f'grid date/time -1: {grid_dt[ind_time[fc_mask]-1]}')

       iobs = np_arange(len(ind_time))
       grid0  = [ ob_grid[i,ind_time[i]  ] for i in iobs[fc_mask] ]
       grid0  = np_array( grid0 )
       grid1  = [ ob_grid[i,ind_time[i]-1] for i in iobs[fc_mask] ]
       grid1  = np_array( grid1 )

       logger.debug(f'grid shape: {grid0.shape}')

       # This is to catch 3d fields. It works, but not sure of its generality
       #   for 4d
       if ( grid0.shape != time_wt.shape ):
            time_wt = time_wt.repeat(grid0.shape[1]).reshape(grid0.shape)

       ob_val = (1-time_wt)*grid0 + time_wt*grid1

       logger.debug(f'grid0: {grid0}')
       logger.debug(f'grid1: {grid1}')
       logger.debug(f'ob_val: {ob_val}')

#      revert to full shape
       ob_val_full = np_ones(fc_mask.shape,float)*np_nan
       ob_val_full[fc_mask] = ob_val
       time_wt_full = np_ones(fc_mask.shape,float)*np_nan
       time_wt_full[fc_mask] = time_wt

       if return_wt:
          return ob_val_full, time_wt
       else: 
          return ob_val_full


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
