import unittest

def encode_id( info, loc_ind ):
    """
    Encode network station index
      - does not work for lists, vectors etc.
    """
    len_id = info['len']-len(info['add'])
    if 'delim' in info:
        loc_ind = loc_ind.split(info['delim'])[info['delim_indx']]
    id = f"{int(loc_ind)}".zfill(len_id)
    id = f"{info['add']}{id}"

    return id

class test_encode(unittest.TestCase):
  def test_rtn(self):
     nw_info = {'add':'cotl', 'scale':10} 
     self.assertEqual( encode_id(nw_info, 5), ['cotl0005'] )


if __name__ == '__main__':
   unittest.main()

