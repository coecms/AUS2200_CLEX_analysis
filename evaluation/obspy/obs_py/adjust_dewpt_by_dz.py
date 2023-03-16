from thermo import mix_rto_from_dewpt, relativeHumidity, mix_rto_from_RH
from thermo import water_vap_pres, dewpt

from .adjust_RH_by_dz import adjust_RH_by_dz
from .adjust_T_by_dz import adjust_T_by_dz

def adjust_dewpt_by_dz( Td,P,T,Pnew, Tnew=None, dz=None, zs=None, znew=None, logger=None ):
       """
       Adjust dew-point (Td) for changes in height
          - convert to relative humidity(RH),
            adjust RH for height as in OPS
            convert back to dewpoint
       Td      = dewpoint [ K ]    at original height
       P,Pnew  = pressure [ Pa ]  at original,new height
       T,Tnew  = temperature [ K ] at original,new height
       dz      = change in height [ m ]. Alternatively provide heights
       zs,znew = original,new heights [ m ]
       """
     
       mix_ratio = mix_rto_from_dewpt(P, Td)
       rel_hum   = relativeHumidity( P, T, mix_ratio )
       rel_hum2  = adjust_RH_by_dz( rel_hum, dz, zs, znew )
       if (Tnew is None):
         Tnew = adjust_T_by_dz( T, dz, zs, znew )
       mix_ratio2 = mix_rto_from_RH( Pnew, Tnew, rel_hum2 )
       vap_press  = water_vap_pres( Pnew, mix_ratio2 )
       Td_new     = dewpt( vap_press )

       if logger is not None:
          logger.debug(f"mix_ratio: {mix_ratio}")
          logger.debug(f"mix_ratio2: {mix_ratio2}")
          logger.debug(f"rel_hum: {rel_hum}")
          logger.debug(f"rel_hum2: {rel_hum2}")
          logger.debug(f"Tnew: {Tnew}")
          logger.debug(f"vap_press: {vap_press}")
          logger.debug(f"Td_new: {Td_new}")
       return Td_new


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
