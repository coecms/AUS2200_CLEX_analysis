from thermo import lowLvlLapse as lowLvlLapse

from .dz_def_ok import dz_def_ok

def adjust_T_by_dz( T, dz=None, zs=None, znew=None ):
       """
       Adjust temperature(Celsius or Kelvin) for changes in height(m).
       Really on valid for small dz
       Can provide zs = height of input temperature and znew = new height
       Same as for OPS:
       https://code.metoffice.gov.uk/doc/ops/ops-2021.03.0/doc/OSDP3.html
       """
 
       dz,dz_check = dz_def_ok( dz, zs, znew )
       Tnew = T - dz*lowLvlLapse
       return Tnew

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
