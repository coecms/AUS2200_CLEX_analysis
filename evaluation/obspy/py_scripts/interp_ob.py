import os
import argparse
import cdms2 as cdms
import numpy as np
from datetime import datetime as dt
from datetime import timedelta as delt
from pathlib import Path

import obs_py as opy
from util_helpers import yaml_cfg, miss_val

def update_val(var,f_info,indx,upd,ob_vals):
    # return vector of updates scattered back to full vector of obs
    if var in ob_vals:
       try:
         upd_vec = ob_vals[var].values
       except:
         upd_vec = ob_vals[var]
    else:
       if f_info.var_dict()[var]['type'][:6] == 'np.int':
            upd_vec = np.array( [miss_val().imiss()]*len(ob_vals) )
       else:
            upd_vec = np.array( [np.nan]*len(ob_vals) )
    np.put(upd_vec, indx, upd)
    return upd_vec
   
def trim_cfg_lists(cfg):
    # delete unwanted items from standard configurations
    for cfg_k,cfg_i in cfg.items():
      if ('del_'+cfg_k in cfg) and isinstance(cfg_i,list):
         for i in range(len(cfg_i),0,-1):
            if cfg_i[i-1] in cfg['del_'+cfg_k]:
                del( cfg_i[i-1] )
    return cfg
     

def check_grid_files(arc_info, base_dt, max_iter, da_window, var_list, file_dt_check, sfc_fields):
   # Check grid files are present, and find possible substitutes
   # will search around for suitable data from other datetimes for topography and land-sea mask
   #   only out to a certain distance (max_iter jumps across da_window)
   # da_window = width of da-window = time between successive analyses (minutes)

   file_d = {}
   for bdt in base_dt:
      file_d[bdt] = {}
      field_checked = []

      if sfc_fields:
         for field in ['surface_height', 'lsm']:
             vdt = bdt
             try:
                if 'no_date_check' in arc_info.file_info()[field]:
                   opy.logger.debug(f"field: {field},    bdt: {bdt}")
                   gfile, g_fname, g_var = arc_info.grid_file(field, bdt )
                else:
                   if file_dt_check[:3] == 'val':
                       gfile, g_fname, g_var = arc_info.grid_file(field, bdt, valid_dt=[bdt])
                   else:
                       gfile, g_fname, g_var = arc_info.grid_file(field, bdt, base_dt=[bdt])
             except ValueError:
                gfile, g_fname, g_var, vdt = find_close_file(field, bdt, arc_info,
                                                            da_window, 0, max_iter, file_dt_check)
   
             if not os.access(g_fname, os.R_OK) :
                 gfile, g_fname, g_var, vdt = find_close_file(field, vdt, arc_info,
                                                              da_window, 0, max_iter, file_dt_check)
             file_d[bdt][field] = {'gfile':gfile, 'g_fname':g_fname,
                                   'g_var':g_var, 'valid_dt':vdt}
             field_checked.append(field)

      for field in var_list:
         if field in field_checked: continue
         try:
            if file_dt_check[:3] == 'val':
               gfile, g_fname, g_var = arc_info.grid_file(field, bdt, valid_dt=[bdt])
            else:
               gfile, g_fname, g_var = arc_info.grid_file(field, bdt, base_dt=[bdt])
         except ValueError:
            gfile = g_fname = g_var = None
         file_d[bdt][field] = {'gfile':gfile, 'g_fname':g_fname, 'g_var':g_var}
         if (g_fname is not None) and os.access(g_fname, os.R_OK) :
            file_d[bdt][field]['valid_dt'] = bdt
         else:
            file_d[bdt][field]['valid_dt'] = None
    
   return file_d


def find_close_file(field, bdt, arc_info, da_window, iter_val, max_iter, file_dt_check):
      # search around looking for possible alternative base datetimes that contain useful fino
      #   for topography and land sea mask only

      # check bdt+da_window, bdt-da_window, bdt+2*da_window, bdt-2*da_window, ...
      # if iter_val < 1 : check basetime ahead in time, otherwise jump back in time relative to 
      #    original
      if iter_val <= 1:
         iter_val = -iter_val
         bdt = bdt + delt(seconds=(2*iter_val+1)*da_window*60)  # bdt updated from previous passes
         iter_val += 1
      else:
         bdt = bdt - delt(seconds=(2*iter_val+1)*da_window*60)  # bdt updated from previous passes
         iter_val  = -iter_val

      if iter_val <= max_iter:
         opy.logger.info(f"Searching for file for {field} valid at {bdt}")
         if file_dt_check[:3] == 'val':
            gfile, g_fname, g_var = arc_info.grid_file(field, bdt,valid_dt=[bdt])
         else:
            gfile, g_fname, g_var = arc_info.grid_file(field, bdt, base_dt=[bdt])
         vdt = bdt
         if not os.access(g_fname, os.R_OK) :
             gfile, g_fname, g_var, vdt = find_close_file(field, bdt, arc_info,
                                                          da_window, iter_val,
                                                          max_iter, file_dt_check)
         opy.logger.info(f"Substitute file {g_fname}")
      else:
         gfile = g_fname = g_var = vdt = None

      return gfile, g_fname, g_var, vdt


if __name__ == '__main__':

#  Initialization

   cfg = yaml_cfg('interp_ob', Path(__file__).parent.absolute() )
   
   parser = argparse.ArgumentParser()
   parser.add_argument('-i','--input', 
                       help="input netCDF file of observations")
   parser.add_argument('-o','--output', 
                       help="output netCDF file of observations and model values")
   parser.add_argument('-f','--forecast', 
                       help="forecast step (hours)")
   parser.add_argument('-v','--verbose', 
                       help="verbosity (logging) levels", default=20)
   cmd_args = parser.parse_args()
   opy.logger.setLevel(cmd_args.verbose)
   obs_f  = cmd_args.input
   out_f  = cmd_args.output
   base_dt_fmt = '%Y%m%d%H%M%S'
   
   for k,i in cfg['dbg_parms'].items() :
       opy.dbg_parm.set_dbg_parm(k,i)
   nprint = opy.dbg_parm.get_dbg_parm('np_dbg_len') 

   # set up timedelta for for forecast from analysis
   if cfg['fc_step_unit'].lower()[0] == 'h':
      fc_step = delt(seconds=float(cmd_args.forecast)*60*60)
   elif cfg['fc_step_unit'].lower()[0] == 'm':
      fc_step = delt(seconds=float(cmd_args.forecast)*60)
   elif cfg['fc_step_unit'].lower()[0] == 's':
      fc_step = delt(seconds=int(cmd_args.forecast))
   
   if cfg['gfile_meta']['dt_fmt'] is not None:
       dtf_grid = opy.dt_fmt(date_fmt=cfg['gfile_meta']['dt_fmt']['date_fmt'],
                             time_fmt=cfg['gfile_meta']['dt_fmt']['time_fmt'])
   else:
       dtf_grid = opy.dtf_grid
   
   arc_info = opy.archive(dir_type=cfg['gfile_meta']['fmt_type'], 
                          root=cfg['gfile_meta']['arch_root'], 
                          dtfmt=dtf_grid, 
                          da_win=cfg['da_window'],
                          fcst_freq=cfg['fc_frequency'], 
                          fc_type=cfg['gfile_meta']['fc_type'], 
                          lev_type=cfg['gfile_meta']['lev_type'], 
                          info_items=cfg['gfile_meta']['var_info'])

   # open obs file
   df, f_info, nw_info = opy.nc_df(obs_f)
   
   # select according to lat/lon domain
   df = opy.df_filter(df, lat_range=cfg['lat_range'], 
                          lon_range=cfg['lon_range'], 
                          time_range=cfg['time_range'], 
                          date_range=cfg['date_range'], 
                          id_list=cfg['stn_id'])

   # delete items from lists
   cfg = trim_cfg_lists(cfg)
   opy.logger.debug(f'configuration: {cfg}')
   
   # generate list of datetimes of obs
   ovd = f_info.var_dict()
   odt_fmt = ovd['date']['attributes']['units']+\
             ovd['time']['attributes']['units']
   len_t = len(ovd['time']['attributes']['units'])
   
   obs_dt = [dt.strptime(str(dfd)+str(dft).zfill(len_t),odt_fmt) for dfd,dft in zip(df['date'][:],df['time'][:])]
   
   # generate list of model basedates and valid_datetimes
   base_dt  = [arc_info.da_time_from_dt(odt-fc_step, use_basedt=True) for odt in obs_dt ]

   # can now select obs with same base_dt
   df['base_dt'] = [odt.strftime(base_dt_fmt) for odt in base_dt]
   
   # unique list of base_dt, valid_dt
   base_dt = list(set(base_dt))
   base_dt.sort()
   opy.logger.info(f'basedate-times: {base_dt}')

   # debugging option - truncate number of base datetimes
   #  but [:-1] gives an empty list if only a single item in list
   if len(base_dt) > 1:
      base_dt = base_dt[:cfg['max_base_dt'] ]
  
   # get list of known variables
   var_list = []
   for k in df.keys():
       if k[:4] == 'obs_' and \
          k[4:] in arc_info.file_info() :
               var_list.append(k[4:])

   if cfg['auxiliary_field'] is not None :
     for f in cfg['auxiliary_field']:
       if f not in var_list: var_list.append(f)

   # generate file_names etc
   # required for checking which data is available
   file_d = check_grid_files(arc_info, base_dt,
                             cfg['max_check_file_iter'],
                             cfg['da_window'], var_list, cfg['file_dt_check'], 
                             cfg['sfc_fields'])
     
# loop over base dates generating interpolated model values

   obsv = {'indx': np.arange(len(df)), 
           'dz':np.array( [miss_val().xmiss()]*len(df) ) }
   
   # lists of variables processed at each stage should be only set once
   #    otherwise on 2nd pass everything is already defined
   skip_var = None

# Grab surface height (required for everything) and adjust latitudes and
#   make sure all obs are in the domain

# check that obs & grid longitudes have the same sense [180,180] or [0,360]
# should also extract missing value and generate appropriate class for
#    each grid
      
   bdt    = base_dt[0]
   if cfg['sfc_fields']:
      field  = 'surface_height' 
      ndim=2
   else:
      field  = 'surface_height' 
      ndim=3

   gfile  = file_d[bdt][field]['gfile']
   g_var  = file_d[bdt][field]['g_var']
   g_val, g_val_dt = opy.transform_grid(gfile(g_var),
                                 arc_info.file_info()[field],
                                 arc_info.valid_dt(gfile, field),
                                 ndim=ndim)
   df, adj_obs_type = opy.match_longitudes(df, g_val)
   df = opy.trim_obs_domain(df, g_val)

   for bdt in base_dt:
   
      opy.logger.debug(f"bdt: {bdt}")

      # Extract infor for height correction if osb_surface_height (=stattion elevation)
      # is in data frame
      if cfg['sfc_fields']:
          if file_d[bdt]['surface_height']['valid_dt'] is None :
              opy.logger.WARN(f'Could not find surface height field for {bdt}')
              continue
          if file_d[bdt]['lsm']['valid_dt'] is None :
              opy.logger.WARN(f'Could not find land-sea-ice mask for {bdt}')
              continue
   
          gfile = file_d[bdt]['surface_height']['gfile']
          g_fname = file_d[bdt]['surface_height']['g_fname']
          g_var = file_d[bdt]['surface_height']['g_var']
    
          f_lsm = file_d[bdt]['lsm']['gfile']
          lsm_fname = file_d[bdt]['lsm']['g_fname']
          lsm_var = file_d[bdt]['lsm']['g_var']
    
          obs_mask = np.equal( df['base_dt'].values, bdt.strftime(base_dt_fmt) )
          if len(obs_mask) == 0: continue
          obs_indx = obsv['indx'][obs_mask]
          opy.logger.debug(f"obs_indx: {obs_indx[:nprint]}")
          opy.logger.debug(f"obs_dt: {obs_dt[:nprint]}")
          opy.logger.debug(f"obs_mask: {obs_mask[:nprint]}")

          # can't mask a list of datetimes - convert to numpy array
          obs_val_dt = list(set(np.array(obs_dt)[obs_mask]))
    
          field = 'surface_height' 
          bkg_f = 'bkg_'+field
          g_val, g_val_dt = opy.transform_grid(gfile(g_var), 
                                     arc_info.file_info()[field],
                                     arc_info.valid_dt(file_d[bdt][field]['gfile'], field),
                                     ndim=2)
          g_lsm, g_lsm_vdt = opy.transform_grid(f_lsm(lsm_var), 
                                     arc_info.file_info()['lsm'],
                                     arc_info.valid_dt(file_d[bdt]['lsm']['gfile'], 'lsm'),
                                     ndim=2)
    

          opy.logger.info(f'Interpolating file: {lsm_fname} {lsm_var}')
          # Land/Sea Mask shouldn't have time dimension so grid_valid_dt is a dummy
          g_ob, lsm = opy.interp(df[obs_mask], g_val, f_info, g_val_dt, 
                                 ls_mask=g_lsm,
                                 allow_all_sea=cfg['allow_all_sea'])
       
          df[bkg_f] = update_val(bkg_f,f_info,obs_indx,g_ob,df)
          obsv['dz'] = update_val('dz',f_info,obs_indx,
                                  df['obs_surface_height'][obs_mask]-g_ob,obsv)
    
          gfile.close()
          f_lsm.close()

      #end if surface_fields:

      if skip_var is None:
         skip_var = list(df.keys())
   
      opy.logger.info(f"var_list: {var_list}")
      opy.logger.info(f"skip_var: {skip_var}")

      for field in var_list:
         # catch any fields used by other fields such as surface_height
         if file_d[bdt][field]['valid_dt'] is None :
             opy.logger.warning(f'Could not find {field} grid for {bdt}')
             continue
         bkg_f = 'bkg_'+field

         opy.logger.debug(f"{bkg_f} in var_dict: {bkg_f in f_info.var_dict()}")

         if bkg_f in skip_var: continue
         if bkg_f not in f_info.var_dict(): continue

         g_var = file_d[bdt][field]['g_var']
         opy.logger.info(f"Interpolating file: {file_d[bdt][field]['g_fname']} {g_var}")

         g_val, g_val_dt = opy.transform_grid(file_d[bdt][field]['gfile'](g_var),
                                    arc_info.file_info()[field],
                                    arc_info.valid_dt(file_d[bdt][field]['gfile'], field),
                                    valid_dt=obs_val_dt)
         # archive should return some info on which variables
         #  are on same grid as land sea mask
         g_ob, lsm = opy.interp(df[obs_mask], g_val, f_info, g_val_dt, 
                                ls_mask=g_lsm,
                                allow_all_sea=cfg['allow_all_sea'])
         opy.logger.debug(f"interpolated values: {g_ob[:nprint]}")
         df[bkg_f] = update_val(bkg_f,f_info,obs_indx,g_ob,df)
         file_d[bdt][field]['gfile'].close()
   

   df = opy.derive_fields_df(df, 'obs_', f_info)
   df = opy.height_correct_df(df, 'obs_')
   df = opy.derive_fields_df(df, 'bkg_', f_info)
   df = opy.height_correct_df(df, 'bkg_')
   df = opy.derive_fields_df(df, 'adj_', f_info)
   df = opy.height_correct_df(df, 'adj_')

   if adj_obs_type != 0:
        df = opy.restore_longitudes(df, adj_obs_type)
    
   opy.df_nc(df, out_f, nw_info, f_info, list_var=list(df.keys()) )
