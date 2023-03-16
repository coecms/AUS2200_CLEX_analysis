from thermo import mix_rto_from_dewpt, mix_rto_from_specificHumidity
from thermo import virtualtemp
from thermo import gravAccel, Rd

from .defns import logger
from .dbg_parm import dbg_parm

def adjust_p_HY(T_obs, Ps_bkg, obs_dz, Td_obs=None, Q_obs=None):
   """
   Adjust surface pressure TO model elevation using
      hydrostatic approx (quick and dirty)
   Choice of obs -> model is to be consistent with OPS and adjust_p_ND
   Note:
      adjust_p_ND is a better option, but requires T850 and Q850
      This routine is quick, dirty and simpler
   
   obs    : Obs data structure for lat, lon & time
   T_obs  : Temperature of obs (K)
   Td_obs : Dewpoint of obs (need q to calculate virtual temp) (K)
   Q_obs  : Specific Humidty of obs (need q to calculate virtual temp) (g/kg)
             (require Td_obs or Q_obs to calculate virtual temp )
   Ps_bkg : Model Surface Pressure (will return in same units, 
                 UNLESS use Q_obs NOT specified
                 Conversion from Td_obs to mixing ratio requires hPa
   obs_dz : station height - model orography (m)
   """
   nprt = dbg_parm.get_dbg_parm('np_dbg_len')

   if Q_obs is None:
      mr = mix_rto_from_dewpt( Ps_bkg, Td_obs)
   else:
      mr = mix_rto_from_specificHumidity(Q_obs)

   tv = virtualtemp( T_obs, mr )
   pob_zs = Ps_bkg -gravAccel*Ps_bkg*obs_dz / Rd / tv

   logger.debug(f'Ps_bkg[:ndp]: {Ps_bkg[:nprt]}')
   logger.debug(f'tv[:ndp]: {tv[:nprt]}')
   logger.debug(f'obs_dz[:ndp]: {obs_dz[:nprt]}')
   logger.debug(f'pob_zs[:ndp]: {pob_zs[:nprt]}')

   return pob_zs

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

