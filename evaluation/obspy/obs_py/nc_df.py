import pandas as pd
import os
import sys
import json   # use loads to go from string to dictionary
from netCDF4 import Dataset

from .obfile_info import obs_var_info
from .ob_network import ob_network

def nc_df(fname):
   """
   Read netcdf file and set up pandas dataframe
      ... and dictionary of attributes
   """
   # get list of variables
   # build dictionary of type, attributes & type 
 
   ncd = Dataset(fname,'r')

   # following would be easier BUT
   # user may re-define individual elements outside yaml files
   if hasattr(ncd,'obfile_base_cfg') :
      obfile_yaml = ncd.obfile_base_cfg
      if not os.access(obfile_yaml, os.R_OK): obfile_yaml = None
   else:
      obfile_yaml = os.path.join(os.getcwd(),'obfile_info.yml')

   if hasattr(ncd,'obfile_user_cfg'):
      user_ob_yaml = ncd.obfile_user_cfg
      if not os.access(user_ob_yaml, os.R_OK): user_ob_yaml = None
   else:
      user_ob_yaml = None
    
   # Assume only one index for obs
   # read data from file, in case of modifications

   nc_info = {}
   ix = list(ncd.dimensions.keys())[0]
   i_nc = ncd.variables[ix]
   df_index = pd.Index(i_nc[:],name=ix)
   df_val = {}
  
   for k,i_nc in ncd.variables.items()  :
      if hasattr(i_nc,'units') and (i_nc.units[:6] == 'string'):
         nc_info[k] = { 'type': '<U'+i_nc.units.split('-')[1] }
      else:
         nc_info[k] = { 'type': 'np.'+str(i_nc.dtype) }

      if hasattr(i_nc,'units') :
          nc_info[k]['attributes'] =  { 'units': i_nc.units } 
      else :
          nc_info[k]['attributes'] =  { 'units': 'Unknown' } 

      if hasattr(i_nc,'long_name'): 
         nc_info[k]['attributes']['long_name'] = i_nc.long_name
      if k == ix : 
         nc_info[k] = {'nc_axis': None}  # just needs to have the attribute
      df_val[k] = i_nc[:]

   nc_hdr = {}
   for k in ncd.ncattrs():
       nc_hdr[k] = ncd.getncattr(k)

   fob_info = obs_var_info(var_dict=nc_info, 
                           nc_axis=ix,
                           file_attr=nc_hdr,
                           yaml_base_cfg=obfile_yaml, 
                           yaml_user_cfg=user_ob_yaml)

   if hasattr(ncd,'network_user_cfg'):
       network_yaml = ncd.network_base_cfg
       if not os.access(network_yaml, os.R_OK): network_yaml = None
   else:
       network_yaml = os.path.join(os.getcwd(),'ob_network.yml')

   user_nw_yaml = None
   if hasattr(ncd,'network_user_cfg'):
      user_nw_yaml = ncd.network_user_cfg
      if not os.access(user_nw_yaml, os.R_OK): user_nw_yaml = None

   network_name = 'odb_surface'
   if hasattr(ncd,'network_name'):
        network_name = ncd.getncattr('network_name')
        if network_name == 'None': network_name = None

   odb_obs_type = None
   if hasattr(ncd,'odb_obs_type'):
       odb_obs_type = ncd.getncattr('odb_obs_type')
       if odb_obs_type != 'None':
          odb_obs_type = json.loads(odb_obs_type)

   varno_translate = None
   if hasattr(ncd,'varno_translate'):
      varno_translate = ncd.getncattr('varno_translate')
      if varno_translate != 'None':
         varno_translate = json.loads(varno_translate)

   odb_columns = None
   if hasattr(ncd,'odb_columns'):
      odb_columns = ncd.getncattr('odb_columns')
      if odb_columns != 'None':
         odb_columns = json.loads(odb_columns)

   odb_regionfilter = None
   if hasattr(ncd,'odb_regionfilter'):
       odb_regionfilter = ncd.getncattr('odb_regionfilter')

   nw_info = ob_network(name=network_name,
                        yaml_user_cfg=user_nw_yaml,
                        varno_translate=varno_translate,
                        odb_obs_type=odb_obs_type,
                        odb_columns=odb_columns,
                        odb_regionfilter=odb_regionfilter)

   df = pd.DataFrame(data=df_val,index=df_index)

   ncd.close()

   return df, fob_info, nw_info

