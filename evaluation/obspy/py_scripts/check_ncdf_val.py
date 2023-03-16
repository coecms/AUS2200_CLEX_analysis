import numpy as np
import obs_py as opy
from glob import glob
import os
import argparse

if __name__ == '__main__':
   parser = argparse.ArgumentParser()
   parser.add_argument('-i','--input',help='Input file')
   parser.add_argument('-r','--reference',help='Reference file for comparison')
   parser.add_argument('-t','--tol',help='Tolerance for comparison', 
                        default='0.0001')
   parser.add_argument('-m','--min_ref',help='Minimum value for reference', 
                        default='0.00001')
   parser.add_argument('-v','--verbose',help='verbosity/debug level', \
                       default=20)
   cfg = parser.parse_args()

   opy.logger.setLevel(cfg.verbose)
   tol = float(cfg.tol)
   min_ref = float(cfg.min_ref)

   f_ex  = glob(os.path.expandvars(cfg.input))
   if len(f_ex) == 1:
      f_ex = f_ex[0]
      print(f"test file: {f_ex}") 
   else:
      opy.logger.error(f" error finding test file: {os.path.expandvars(cfg.input)}")
      raise ValueError

   f_ref = os.path.expandvars(cfg.reference)
   if f_ref[-3:] != '.nc.': f_ref = os.path.join(f_ref,f_ex)
   f_ref = glob(f_ref)
   if len(f_ref) == 1:
      f_ref = f_ref[0]
      print(f"reference file: {f_ref}") 
   else:
      opy.logger.error(f" error finding reference file: {os.path.expandvars(cfg.reference)}")
      raise ValueError

   df_ex,fex_info,nw_ex = opy.nc_df(f_ex)
   df_ref,fref_info,nw_ref = opy.nc_df(f_ref)

   var_ex = df_ex.keys()
   var_ref = df_ref.keys()

   if len(df_ex) == len(df_ref):
      print(f"PASS: Both files have same length {len(df_ex)}")
   else:
      print(f"FAIL:  files have different lengths")
      print(f"     Experiment DataFrame: {len(df_ex)}")
      print(f"     Reference DataFrame: {len(df_ref)}")

   list_var = []
   for kx in var_ex:
      if kx not in var_ref:
         print(f"FAIL: variable {kx} in test file, missing from reference file")
      else:
         print(f"PASS: variable {kx} in test file also in reference")
         list_var.append(kx)

   for kr in var_ref:
      if kr not in var_ex:
          print(f"FAIL: variable {kr} in reference file, missing from test file")
         
   for k in list_var:
       try:
          nan_ex  = np.isnan(df_ex[k]).sum()
          # only want to catch problems with comparing with NaN, no others
       except:
          ndiff = np.not_equal(df_ref[k], df_ex[k]).sum()
          if ndiff == 0:
             print(f'PASS: zero differences for column {k}')
          else:
             print(f'FAIL: Number of differences for {k} is {ndiff}')
          continue

       nan_ref = np.isnan(df_ref[k]).sum()
       if nan_ex == nan_ref:
           print(f"PASS: column {k} has same number of NaN {nan_ex}")
       else:
           print(f"FAIL: column {k} has mismatched number of NaN")
           print(f"      test count NaN: {nan_ex}, reference count NaN: {nan_ref}")
       av = abs(df_ref[k])
       norm = np.where(av > min_ref, av, min_ref)
       rel_diff = abs(df_ex[k] - df_ref[k])/norm
       rel_diff_max = rel_diff.max()
       if rel_diff_max < tol:
           print(f"PASS: column {k} relative differences within {tol}")
       else:
           imax = rel_diff.argmax()
           print(f'FAIL: Maximum relative difference for {k} is {rel_diff_max} at {imax}')
           for k in list_var:
              try:
                 abs_diff = abs(df_ex[k] - df_ref[k]) 
                 if abs_diff > tol : 
                      print(f"{k} {df_ex[k].values[imax]} {df_ref[k].values[imax]} {abs_diff}")
              except:
                 print(f"{k} {df_ex[k].values[imax]} {df_ref[k].values[imax]}")
            
