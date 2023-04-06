from thermo import C2K

def df_std_units(df, nw_info, fvar_dict):
   """
   Convert data to standard units
   df: pandas dataframe
   nw_info : instance of ob_network().info
             info on reported variables
   fvar_dict : instance of obfile_info().var_dict()
             info on standard variables

   This could be done as read in data from obs file
     - but wil need repeated code in multiple places
   """

   for v in df.keys():
      if v in nw_info['inv_dict']:
         ob_units = nw_info['inv_dict'][v]['units'] 
         f_units = fvar_dict[v]['attributes']['units']
         if ob_units != f_units:
            units = [ ob_units, f_units ]
            if units == ['C','K']:
              df[v] = df[v] + C2K
            elif units == ['hPa','Pa']:
              df[v] = df[v]*100
            else:
              print(' unknown pair of units: ',units)
              raise IndexError

   return df
