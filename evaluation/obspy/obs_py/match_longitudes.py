def match_longitudes(df, grid, lat_name='latitude', lon_name='longitude'):
   """
    Remove observations outside of a grid domain to avoid extrapolation
       df: Pandas DataFrame of observations
       grid: datastructure that provides grid latitudes and longitudes

    Currently drops rows, but an alternative option would be to return a mask
   """

   if (df[lon_name].max() > 180.)  \
        and (grid.getLongitude().getValue().min() < 0.):
       df[lon_name] = np.where(df[lon_name] > 180., \
                               df[lon_name]-360, \
                               df[lon_name])
       trim_type=-180

   elif (df[lon_name].min() < 0.)  \
        and (grid.getLongitude().getValue().max() > 180.):
       df[lon_name] = np.where(df[lon_name] < 0., \
                               df[lon_name]+360, \
                               df[lon_name])
       trim_type=360
   else:
       trim_type = 0

   return df, trim_type
    
