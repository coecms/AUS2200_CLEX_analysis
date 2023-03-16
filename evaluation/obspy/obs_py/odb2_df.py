import numpy as np
import pandas as pd
import obsio as odb
import sys

from .make_odb_typefilter import make_odb_typefilter
from .ob_convert import ob_convert
from .defns import logger

def odb2_df(odb_file,netwk,fob_info,odb_filter=None,obs_type='surface'):

    """
     Generate pandas DataFrame from ODB2
     odb_file : name of ODB2 file
     netwk : instance of ob_network class
     fob_info : instance of obs_var_info (observation file info)
     odb_filter : output from make_odb_typefilter
     obs_type : type of observations to use
    """

    if 'obs_dict' not in netwk.info:
       logger.error('obs_dict undefined in ob_network configuration')
       raise KeyError

    if odb_filter is None:
        odb_obs_types = tuple(netwk.odb_obs_type()[obs_type])
        odb_filter = odb.set_odb_filter(ops_obstype=odb_obs_types)

    logger.debug(f"odb_filter : {odb_filter}")
    obsdata = odb.fromodb(odb_file, netwk.odb_columns(), where=(odb_filter))
    logger.debug(f"obsdata keys: {odb.list_odb(odb_file,table=False)}")
    
    varno=obsdata['varno']
    repmask = (varno == netwk.info['report_varno'])
    nrep = np.count_nonzero(repmask)
    
    dfd_index = {}
    dfd_val = {}

    # extract information for dataframe axes
    # kv & iv will be obs_py parameter names from ob_network.yml (latitude, station_identifier etc.)
    #    - from ob_network.info['var_dict']
    #    - the description of these parameters is in obs_var_info.var_dict (obfile_info.py/yml)
    for kv in netwk.info['df_index']:
        iv = netwk.info['obs_dict'][kv]
        logger.debug(f"build dataframe axis, kv: {kv}, iv: {iv}")

        if 'varno' in iv:
           dfd_index[kv] = obsdata[iv['odbcol']][varno==iv['varno']]
           logger.debug(f"varno k,i = {kv},{iv['varno']}")

        elif 'odbcol' in iv:
           # This is some code from Susan Rennie &/or Andy Smith. Not 100% sure how it actually works
           dfd_index[kv] = obsdata.pop(iv['odbcol'])[repmask]
           logger.debug(f"odbcol k,i = {kv},{iv['odbcol']}")

        dfd_index[kv] = ob_convert( fob_info.var_dict()[kv]['type'], 
                                    dfd_index[kv] )

        logger.debug(f"dfd_index.keys = {dfd_index.keys()}")
        dfd_val[kv] = dfd_index[kv]

    # set up MultiIndex indices
    dfd_val[fob_info.nc_axis()]=np.arange(len(dfd_index[kv]))
    logger.debug(f"dfd_index keys: {dfd_index.keys()}")
    logger.debug(f"kv: {kv}, len obs axis: {len(dfd_index[kv])}")
   

    indx_arrays = [ va for kv,va in dfd_index.items() ]
    indx_names = list( dfd_index.keys() )
    df_index = pd.MultiIndex.from_arrays(indx_arrays,names=indx_names)

    # extract rest of info
    skip_obstype = []
    for kv,iv in netwk.info['obs_dict'].items():
        if kv in dfd_val:
            continue
        
        obs_val = []
        if 'varno' in iv:
            obs_val = obsdata[iv['odbcol']][varno==iv['varno']]
            logger.debug(f"varno k,i = {kv},{iv['varno']}")
        elif 'odbcol' in iv:
            obs_val = obsdata.pop(iv['odbcol'])[repmask]
            logger.debug(f"odbcol k,i = {kv},{iv['odbcol']}")

        # Need to use temporary storage as field may not be in ODB, and returns data length 0
        if len(obs_val) > 1:
            dfd_val[kv] = obs_val
            logger.debug(f"building dfd_val.keys: {dfd_val.keys()}, len: {len(dfd_val[kv])}")
        else:
            skip_obstype.append(kv)

    logger.debug(f"End building dfd_val directly from obsdata")

    #Needed because old OPS doesn't have these both values included for some obs (u,v etc.)
    mask_copy = {}

    for kv,iv in netwk.info['obs_dict'].items():
        if ('mask_copy' in iv) and (iv['mask_copy'] not in skip_obstype) :
           mask_copy[kv] = dfd_val[iv['mask_copy']].mask.copy()

    for kv,iv in netwk.info['obs_dict'].items():
        if ('mask_copy' in iv) and (iv['mask_copy'] not in skip_obstype) :
           dfd_val[kv][mask_copy[kv]] = \
                    obsdata['obsvalue'][varno==iv['varno']][mask_copy[kv]]
           logger.debug(f"mask_copy kv: {kv}, iv: {iv}, len: {len(dfd_val[kv])}")

    # variable conversions
    for kv,iv in netwk.info['obs_dict'].items():
      if 'transform' in iv:
         var1 = iv['transform']['var1'] 
         var2 = iv['transform']['var2'] 
         if (var1 in dfd_val) and (var2 in dfd_val):
            dfd_val[kv] = ob_convert( iv['transform']['func'],  \
                                      dfd_val[var1], dfd_val[var2] )
            logger.debug(f"transform kv: {kv}, len: {len(dfd_val[kv])}")

    # type conversions (int32, float32, float64 etc.)
    for kv in dfd_val:
       logger.debug(f"type convert kv: {kv}")
       dfd_val[kv] = ob_convert( fob_info.var_dict()[kv]['type'], 
                                 dfd_val[kv] )
       logger.debug(f"len: {len(dfd_val[kv])}")
        
    return pd.DataFrame(data=dfd_val, index=df_index)
