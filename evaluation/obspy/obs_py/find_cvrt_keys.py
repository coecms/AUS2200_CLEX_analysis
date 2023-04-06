import copy
import numpy as np
import thermo
import os

from .defns import logger
from .ffddd_to_uv import ffddd_to_u, ffddd_to_v
from .uv_to_ffddd import uv_to_ddd, uv_to_ff

def find_cvrt_keys(df_vars, cvrt_info, var_dep,  \
                   npass=0, pfx_input=None, pfx_derive=None):
   """
   Returns a dictionary of information (cvrt_info) on variables that can be derived from 
     a list of existing vars (df_vars)
   Uses recursion as first pass may allow one variable to be derived that
     can then be used to derive another
   So order of names in upd_df_vars is important
     upd_df_vars is the list of dataframe keys once derivaed variables are calculated
     cvrt_info holds info on which conversion to use
         - keys of cvrt_info are list of variables to convert

   var_dep: dictionary of derived variables and their dependencies
               one of defns.field_dep or defns.ht_corr_dep
               If use var_dep=ht_corr_dep, then should also use
                   pfx_input='bkg_' and pfx_derive='adj_'

   npass is belt and braces for recursion limit

   pfx_input: prefix for input variables in dataframe to be considered (None --> all)
   pfx_output: prefix for output variables (None --> same as pfx_input)
   """

   all_found = False
   upd_df_vars = copy.deepcopy(df_vars)
   add_cvrt_info = False   # Assume no more variables found that can be converted
   if (pfx_input is not None) and (pfx_derive is None): pfx_derive = pfx_input
   npass += 1

   # search through variables in dataframe
   #    1. working out which variables are already there
   #    2. and which variables can be derived

   for df_v in df_vars:
      # check if the variable is in the dict of field dependencies
      #     and extract all the stuff before the dictionary name (pfx)
      #
      #  pfx could be 'obs_10m_', 'bkg_screen_' etc.
      #
      #  pfx_input = None ==> do all
      #
      pfx = None
      for k in var_dep.keys():
          if (df_v[-len(k):] == k)  and \
             ((pfx_input is None) or (pfx[:len(pfx_input)] == pfx_input)):
              pfx = df_v[:-len(k)]
              break

      # variable name is not one that can be used to derive another
      #   move on to the next
      if pfx is None: continue
      if pfx_derive is not None:
          pfx_d = pfx_derive + pfx[len(pfx_derive):]
      else:
          pfx_d = pfx
      
      # generate list of variables with same prefix
      #   related_var (and therefore upd_df_vars) refers to input variables, so use pfx
      related_var = []
      for k in var_dep.keys():
         if pfx+k in upd_df_vars: related_var.append(k)

      logger.debug(f"prefix: {pfx}, related_var: {related_var}")

      # update dictionary of variables that can be derived
      for k,i in var_dep.items():
     
          if k not in related_var:
             # have found variable (pfx+k) that may need to be derived
             # check if it can be
             #   use_vars refers to input variables, so use pfx
             for il_dep,l_dep in enumerate(i):  # l_dep is a list of dependencies
                dep_ok = True
                use_vars = []
                for v in l_dep:
                   if pfx+v not in upd_df_vars: dep_ok = False
                   use_vars.append(pfx+v)
                if dep_ok:
                    # can generate this variable (pfx+k)
                    logger.debug(f"adding variable to be derived: {pfx_d + k}")
                    cvrt_info[pfx_d + k] = {}
                    cvrt_info[pfx_d + k]['l_indx'] = il_dep
                    cvrt_info[pfx_d + k]['var'] = k
                    cvrt_info[pfx_d + k]['input'] = use_vars
                    upd_df_vars.append(pfx_d + k)
                    add_cvrt_info = True
                    break
                 # else cannot derive variable, as some variables are missing
          # else related variable already exists
                 
   if npass > 10:
      logger.error("Too many recursions")
      raise ValueError
   elif add_cvrt_info:
      # go round again seeing if anything new can be derived
      upd_df_vars, cvrt_info, npass =  \
           find_cvrt_keys(upd_df_vars, cvrt_info, var_dep, 
                          npass, pfx_input, pfx_derive)

   # else: return

   return upd_df_vars, copy.deepcopy(cvrt_info), npass

if __name__ == '__main__':
  from .defns import logger, field_dep, ht_corr_dep
  from .nc_df import nc_df
  from .derive_fields_df import derive_fields_df
  indir = os.path.expandvars('${kgo_topdir}/odb')
  #opy.logger.setLevel(10)

  #pfx = 'obs_'
  #df,fi,nw = nc_df(os.path.join(indir,'surface.nc'))
  pfx = 'bkg_'
  df,fi,nw = nc_df(os.path.join(indir,'interp_surface_g3.nc'))

  dfk = []
  dfk0 = df.keys()
  for v in df.keys():
     if v[:4] == pfx:
        dfk.append(v)
  ci = {}

  list_var, cv_i, ip = find_cvrt_keys(dfk,ci,field_dep)

  df2 = derive_fields_df(df, pfx, fi)

  dfk2 = list(df2.keys())
  der_v = []
  for v in dfk2:
     if v not in dfk0: der_v.append(v)
  print(f" new fields: {der_v}")
  
