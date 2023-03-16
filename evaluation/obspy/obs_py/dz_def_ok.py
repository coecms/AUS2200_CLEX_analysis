from .defns import logger

def dz_def_ok( dz, zs, znew):
       """
       Check that enough info provided to calculate height difference
       https://code.metoffice.gov.uk/doc/ops/ops-2021.03.0/doc/OSDP3.html
       Note: For dz has a "positive sense" i.e. dz > 0 implies go upwards
           This is different from some of the from the OPS where adjust screen obs to model levels.
       """
 
       if not( dz is None ): 
          flag = True
       elif not(zs is None) and not(znew is None):
          dz = znew-zs 
          flag = True
       else:
          flag = False

       if not flag:
            logger.error('could not calc dz. dz, zs, znew: {} {} {}\n'.format(dz, zs, znew))
            raise ValueError

       return dz,flag

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
