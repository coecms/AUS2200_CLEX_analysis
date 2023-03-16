from .defns import logger

def check_grid(grid):
       """
       Check that a grid dimensions are OK
       Should be lat=dim[0], lon=dim[1], time=dim[2]
       """

       idim_lat = 0
       idim_lon = 1

       if grid.shape[idim_lat] != grid.getLatitude().shape[0] :
          logger.error('Latitudes must be dimension {0}\n'.format(idim_lat))
          logger.error('grid shape: {} \n'.format(grid.shape))
          logger.error('latitude shape: {}\n'.format(grid.getLatitude().shape))
          raise ValueError

       if grid.shape[idim_lon] != grid.getLongitude().shape[0] :
          logger.error('Longitudes must be dimension {0}\n'.format(idim_lon))
          logger.error('grid shape: {}\n'.format(grid.shape))
          logger.error('longitude shape: {}\n'.format(grid.getLongitude().shape))
          raise ValueError
       
       return 

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
