from numpy import where as np_where

from .defns import OPS
from .dz_def_ok import dz_def_ok

def adjust_RH_by_dz( RH, dz=None, zs=None, znew=None ):
       """
       Adjust relative humidity (%) for changes in height(m).
       Really on valid for small dz
       Can provide zs = height of input temperature and znew = new height
       Opposite sign to OPS.
       https://code.metoffice.gov.uk/doc/ops/ops-2021.03.0/doc/OSDP3.html

       """

       dz,dz_check = dz_def_ok( dz, zs, znew )
       RHnew = RH + dz*OPS.get_lowLvlLapse_RH()

     # Ensure 0 <= RHnew  <= 100.
       RHnew = np_where( RHnew < 0.,     0., RHnew )
       RHnew = np_where( RHnew > 100., 100., RHnew )
       return RHnew
       
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
