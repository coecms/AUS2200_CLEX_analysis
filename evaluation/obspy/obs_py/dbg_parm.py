from .defns import logger

class dbg_parm_class(object):
   def __init__(self):
      self.list_dir = [ 'get_dbg_parm', 'set_dbg_parm' ]
      return

   def __dir__(self):
       return self.list_dir

   def set_dbg_parm(self, var, val):
      if not hasattr(self,var):
         self.list_dir.append(var)
      setattr(self,var,val)
      return
      
   def get_dbg_parm(self, var):
      if hasattr(self,var):
         return getattr(self,var)
      else:
         logger.error(f"Unknown parameter name: {var}")
         raise ValueError

dbg_parm = dbg_parm_class()
