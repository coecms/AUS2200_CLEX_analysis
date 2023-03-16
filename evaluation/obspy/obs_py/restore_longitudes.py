def restore_longitudes(df, adj_type, lat_name='latitude', lon_name='longitude'):
   """
    Remove observations outside of a grid domain to avoid extrapolation
       df: Pandas DataFrame of observations
       adj_type: type fo adjustment made originally
                 -180 --> was remapped from [   0,360] to [-180,180]
                  360 --> was remapped from [-180,180] to [   0,360]
       grid: datastructure that provides grid latitudes and longitudes

    Currently drops rows, but an alternative option would be to return a mask
   """

   if adj_type == -180:
       df[lon_name] = np.where(df[lon_name] < 0., \
                               df[lon_name]+360, \
                               df[lon_name])
   elif adj_type == 360:
       df[lon_name] = np.where(df[lon_name] > 180., \
                               df[lon_name]-360, \
                               df[lon_name])
   return df
    
