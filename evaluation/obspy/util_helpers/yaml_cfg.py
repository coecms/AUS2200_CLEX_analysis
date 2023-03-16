from pathlib import Path
import os
import yaml

def yaml_cfg(rtn_name,rtn_path):
    """
    General routine to read config from yaml files

    rtn_name: name of routine, will search for file rtn_name.yml
 
    rtn_path: full path to calling script/package. 
          default settings assumed to be in this directory
          pathlib.Path(__file__).parent.absolute() will suffice
    """

    cfg_def = os.path.join(rtn_path,rtn_name+'.yml')

    print('Default configuration for ',rtn_name,'.py from: ',cfg_def)

    with open(cfg_def,'r') as stream:
       cfg = yaml.safe_load(stream)
           
    cfg_user = rtn_name+'.yml'
    if os.access(cfg_user,os.R_OK):
      with open(cfg_user,'r') as stream:
         user_cfg = yaml.safe_load(stream)

#   following could probably be done better with recursion
#      but this already taken much longer than planned

         for k,i in user_cfg.items():
           if k in cfg:
             if cfg[k] is None:
               cfg[k] = i
             elif isinstance(i,list):
               for j in i:
                  if not j in cfg[k]: cfg[k].append(j)
             elif isinstance(i,dict):
                for j in i:
                    if isinstance(user_cfg[k][j],dict):
                      cfg[k][j] = user_cfg[k][j].copy()
                    else:
                      cfg[k][j] = user_cfg[k][j]
             else:
                cfg[k] = i
           else:
              cfg[k] = i
       
    return cfg
