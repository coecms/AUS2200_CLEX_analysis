"""

**** This does not work for DataFrames *****
    Dunno why

Generalize JSON dump & load to handle numpy arrays as Pandas tseries

Example usage:
   expected = np.arange(100, dtype=np.float)
   dumped = json.dumps(expected, cls=GenEncoder)
   result = json.loads(dumped, object_hook=json_gen_obj_hook)

   # None of the following assertions will be broken.
   assert result.dtype == expected.dtype, "Wrong Type"
   assert result.shape == expected.shape, "Wrong Shape"
   assert np.allclose(expected, result), "Wrong Values"
"""
import base64
import json
import numpy as np
import pandas as pd
from sys import exit

use_b64 = False # For some reason have problems with writing base64 encoded numpy arrays
                # to JSON files in Python3

mod_name = 'jsonUtilities'

class GenEncoder(json.JSONEncoder):
    def default(self, obj):
        """
        if input object is a ndarray it will be converted into a dict holding dtype, shape and the data base64 encoded

        Extended to handle pandas.Series
        """

        # Order of tests is important here as pd.Series can be
        #  mistaken for np.ndarray
        if isinstance(obj, pd.DataFrame):

            print(mod_name+"GenEncoder does not handle DataFrames")
            exit(1)

            return dict(__pandas__=obj.to_json(),
                        pdtype='frame')

        elif isinstance(obj, pd.Series):

            return dict(__pandas__=obj.to_json(),
                        pdtype='series')


        elif isinstance(obj, np.ndarray):
            if use_b64:
               data_64 = base64.b64encode(obj.data)
            else:
               data_64 = obj.data.tolist()

            return dict(__ndarray__=data_64,
                        dtype=str(obj.dtype),
                        shape=obj.shape)

        # Let the base class default method raise the TypeError
        return json.JSONEncoder(self, obj)


def json_gen_obj_hook(dct):
    """
    Decodes a previously encoded numpy ndarray
    with proper shape and dtype
    :param dct: (dict) json encoded ndarray
    :return: (ndarray) if input was an encoded ndarray

    Extended to handle pandas.Series
    """
    if isinstance(dct, dict) and '__pandas__' in dct:

        return pd.read_json(dct['__pandas__'],typ=dct['pdtype'])

    elif isinstance(dct, dict) and '__ndarray__' in dct:

        if use_b64 : 
           np_data = base64.b64decode(dct['__ndarray__'])
           return np.frombuffer(np_data, dct['dtype']).reshape(dct['shape'])
        else:
           return np.array(dct['__ndarray__'],dtype=dct['dtype']).reshape(dct['shape'])

    return dct

