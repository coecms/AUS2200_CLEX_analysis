# Sets up meta data for observations

import os
import yaml
from pathlib import Path

from .defns import PKG_PATH, logger
from .defns import dtf_obs

class obs_var_info( object ) :
    """
       Class to set up meta data for observations
    """
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # class constructor
    #    initialization of variables

    def __init__(self, var_dict=None, 
                       file_attr=None,
                       nc_axis=None, 
                       yaml_base_cfg=None,
                       yaml_user_cfg=None):

        if yaml_base_cfg is None:
            yaml_base_cfg = os.path.join(PKG_PATH,'obfile_info.yml')
        self.yaml_base_cfg = yaml_base_cfg
        print('obfile meta-data read from ',yaml_base_cfg)
        with open(yaml_base_cfg,"r") as stream:
            obfile_yaml = yaml.safe_load(stream)

        if yaml_user_cfg is not None:
            yaml_user_cfg = Path(yaml_user_cfg).absolute()
            print('user defined obfile meta-data read from ',yaml_base_cfg)
            with open(yaml_base_cfg,"r") as stream:
                obfile_yaml_user = yaml.safe_load(stream)
            for k,i in obfile_yaml_user.items():
                obfile_yaml[k] = i.copy()
        self.yaml_user_cfg = yaml_user_cfg
        
        if var_dict is not None:
           # copy updated items - but keep existing items from yaml file
           #   that are not in netcdf (e.g. obs_type)
           for k,i in var_dict.items():
              if isinstance(i,dict):
                 for kk,ii in i.items():
                    obfile_yaml['var_dict'][k][kk] = ii
              else:
                 obfile_yaml['var_dict'][k] = i

        if file_attr is not None:
           for k,i in file_attr.items():
                 obfile_yaml['file_attr'][k] = i
              
        self.set_var_dict(obfile_yaml['var_dict'])
        self.set_file_attr(obfile_yaml['file_attr'])
        self.set_nc_axis(val=nc_axis)
        
        return
  
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    def __dir__(self):

        list_dir = [ '__init__', '__dir__', \
                     '__var_dict', '__file_attr', \
                     '__nc_axis', \
                     'file_attr', 'var_dict', 'nc_axis', \
                     'add_var_dict_item', \
                     'add_file_attr_key', \
                     'set_file_attr', \
                     'set_nc_axis', \
                     'set_var_dict', 'set_var_dict_key', \
                   ]

        return list_dir

    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    def var_dict(self):

        return self.__var_dict

    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    def file_attr(self):

        return self.__file_attr

    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    def nc_axis(self):

        return self.__nc_axis

    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    def set_var_dict(self, val):

        self.__var_dict = val
        self.set_nc_axis( fdict=self.var_dict() )

        return

    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    def set_file_attr(self, val):

        self.__file_attr = val

        return

    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    def set_nc_axis(self, val=None, fdict=None):

        if val is not None:
           self.__nc_axis = val
        if fdict is not None:
           for k,i in fdict.items():
              if 'nc_axis' in i:
                 self.__nc_axis = k
        return

    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    def add_var_dict_item(self, new_key, item, key_list=None):
        """
        Extend var_dict, the dictionary of meta-data about 
        variables in obs files
          add_var_dict_item(new_key, item) will add a new item 
                to every key in var_dict
          add_var_dict_item(new_key, item, key_list) will add a new item 
                to specific keys in var_dict

          NOTE: can be used to create entirely new keys in var_dict, by
                having key_list being a list of only new keys
        """

        var_dict = self.var_dict()
        if key_list is None:
            key_list = list(var_dict.keys())

        for k in key_list:
            var_dict[k][new_key] = item

        self.set_var_dict(var_dict)
            
        return


    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    def add_file_attr_key(self, key, item):

       self.__file_attr[key] = item
       return

    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    def set_var_dict_key(self, key, itemd):
        """
        Set one key of var_dict, the dictionary of meta-data about 
        variables in obs files
          key : the key in var_dict that it is to be modified
          itemd : a dictionary, that will be added to/modify var_dict element by element
        """

        for k,i in itemd:
           self.__var_dict[key][k] = i

        return
