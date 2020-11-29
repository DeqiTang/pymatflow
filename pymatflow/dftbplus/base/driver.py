"""

"""

import os 

from pymatflow.dftbplus.base.hsd import hsd_block

def new_driver():
    out = hsd_block(name="Driver", val=None, block_type="method", level=0)
     
    out.status = True
    return out
