"""

"""

import os 

from pymatflow.dftbplus.base.hsd import hsd_block

def new_options():
    out = hsd_block(name="Options", val=None, block_type="method", level=0)
    
    out.status = True
    
    out.scalar["WriteDetailedOut"] = "Yes"

    return out
