def make_odb_typefilter(odb_obs_type, region=None):
     """
     given an obs type (list of OPS types) construct a filter for use with ODB
       must use | for "or" and & for "and"
     Need to encase each sub expression in () due to precedence of operators
       i.e. '|' precedes '==' normally
     """
     odb_typefilter=['(ops_obstype=={:})'.format(x) for x in odb_obs_type]

     if (region is not None) and (region.strip()[0] == '('):
          odb_typefilter='( ('+' | '.join(odb_typefilter)+' ) & ( '+region+') )'

     return odb_typefilter



