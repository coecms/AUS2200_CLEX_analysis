import sys
import os 
import json   # use json.dumps to convert strucutres to string
import numpy as np
from netCDF4 import Dataset
from datetime import datetime as dt

from util_helpers import miss_val

from .defns import PKG_PATH, dtf_create, code_version, institution, logger
from .ob_convert import ob_convert


def df_nc(df,out_file,netw,fob_info,obs_type='surface',
          file_attr={},list_var=None,
          xmiss=miss_val()):
    """
    Write pandas dataframe to netcdf file
    df: pandas dataframe
    out_file: netcdf filename
    netw: instance of ob_network class
    fob_info: instance of obs_var_info class
    obs_type: category of observations (string)
    file_attr: dictionary of attributes for file
    """

    file_name=out_file
    ncdataset=Dataset(file_name,'w')
   
    file_attr['CreationDate'] = dt.strftime(dt.utcnow(),dtf_create.dt_fmt())
    file_attr['OBSPY_location'] = str(PKG_PATH)
    file_attr['obfile_base_cfg'] = str(fob_info.yaml_base_cfg)
    file_attr['obfile_user_cfg'] = str(fob_info.yaml_user_cfg)
    file_attr['network_base_cfg'] = str(netw.yaml_base_cfg)
    file_attr['network_user_cfg'] = str(netw.yaml_user_cfg)
    file_attr['network_name'] = str(netw.name)
    file_attr['varno_translate'] = json.dumps(netw.varno_translate())
    file_attr['odb_obs_type'] = json.dumps(netw.odb_obs_type())
    file_attr['odb_columns'] = json.dumps(netw.odb_columns())
    file_attr['odb_regionfilter'] = netw.odb_regionfilter()
    file_attr['ObsStartDate'] = str(df['date'].min())
    file_attr['ObsEndDate'] = str(df['date'].max())
    file_attr['ObsStartTime'] = str(df['time'].min())
    file_attr['ObsEndTime'] = str(df['time'].max())
    file_attr['OBSPY_version'] = code_version
    file_attr['institution'] = institution
    if 'FileHistory' in file_attr: 
        file_attr['FileHistory'] = file_atttr['FileHistory']+'\n'+' '.join(sys.argv)
    else:
        file_attr['FileHistory'] = ' '.join(sys.argv)
    file_attr['CurrentDirectory'] = os.getcwd()

    for k,i in fob_info.var_dict().items():
       if 'nc_axis' in i: 
          nc_axis = k
          break
    ncdataset.createDimension(nc_axis,size=len(df))

    # add a copy of the index as a variable so automatically
    #  will create netcdf variable corresponding to the dimension

    if nc_axis not in list(df.keys()):
      df[nc_axis] = df.index.get_level_values(level=nc_axis).values

    logger.debug(f'df.keys = {df.keys()}')

    if list_var is None:
       if 'inv_dict' in netw.info:
          list_var = list(netw.info['inv_dict'].keys())
          list_var.append(fob_info.nc_axis())
       elif 'obs_dict' in netw.info:
          list_var = list(netw.info['obs_dict'].keys())

    # ensure index is in list_var
    list_var.append(nc_axis)
    # make list_var a list of unique elements
    list_var  = list(set(list_var))   

    logger.debug(f'list_var: {list_var}')
    for k,i in fob_info.var_dict().items():
         
        if k not in list_var: continue
        logger.debug(f'k: {k}')
        if (k[:4]=='obs_') and (i['obs_type'] is None): continue

        tipe=i['type']
        if tipe is None : continue
        if k not in df: continue

        if (obs_type in i['obs_type']) or (i['obs_type'][0] == 'all') \
            or (k[:4] in ['bkg_','adj_']) :
           logger.debug(f'k,tipe = {k},{tipe}')
           logger.debug(f'df,shape =  {df[k].shape}')

           if len(df[k]) == 0 : continue

           if tipe[:6]=='np.int':
              fill_value=xmiss.imiss()
              ncdataset.createVariable(k,
                                       ob_convert(tipe,1), nc_axis,
                                       fill_value=fill_value)
           elif tipe[:8]=='np.float':
              fill_value=xmiss.xmiss()
              ncdataset.createVariable(k,
                                       ob_convert(tipe,1), nc_axis,
                                       fill_value=fill_value)
           elif tipe[:2]=='<U':
              fill_value=xmiss.cmiss()
              ncdataset.createVariable(k, tipe, nc_axis,
                                       fill_value=fill_value)
           else: 
              continue

           var=ncdataset[k]
           #var.setncattr('_FillValue',fill_value)

           try:
             var[:]=df[k].values.filled(fill_value=fill_value)
             logger.debug(df[k].values.filled()[0:5])
           except:
#            This gives floats correctly masked in netcdf file
             if tipe[:8] == 'np.float':
                var[:] = np.where(np.isnan(df[k].values), fill_value, df[k].values)
             else:
                var[:]=df[k].values

           logger.debug(f'var_dict k: {k},  i: {i}')
           var_attr = {}
           for ka,ia in i['attributes'].items():
               var_attr[ka] = ia
               logger.debug(f"attribute {ka}: {ia}")
           var_attr['grid_type'] = "None"
           var_attr['fill_value'] = fill_value
           var_attr['missing_value'] = fill_value
           var.setncatts(var_attr)
           print('written variable ',k,' type: ',tipe)
    
    ncdataset.setncatts(file_attr)
    ncdataset.close()
