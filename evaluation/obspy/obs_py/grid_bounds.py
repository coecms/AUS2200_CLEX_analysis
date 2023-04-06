
def grid_bounds( field, buff_width):
   """
   returns latitude and longtiude bounds that exclude an
     outer buffer zone of a CDMS grid

   required to guard against wind grids being greater than land-sea grids
       in any dimension

   field: netcdf Land-sea-ice mask
   buff_frac: fractional value giving width of buffer relative to domain size
   """

   lat = field.getLatitude()
   lon = field.getLongitude()

   if lat[0] < lat[-1]:
      lat_min = lat[buff_width]
      lat_max = lat[-buff_width-1]
   else:
      lat_max = lat[buff_width]
      lat_min = lat[-buff_width-1]

   if lon[0] < lon[-1]:
      lon_min = lon[buff_width]
      lon_max = lon[-buff_width-1]
   else:
      lon_max = lon[buff_width]
      lon_min = lon[-buff_width-1]

   return [lat_min,lat_max], [lon_min, lon_max]
