"""

"""

import os 

from pymatflow.dftbplus.base.hsd import hsd_block

def new_reks():
    out = hsd_block(name="REKS", val=None, block_type="method", level=0)
    
    out.status = False
    

    return out
