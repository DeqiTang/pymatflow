"""

"""

import os 

from pymatflow.dftbplus.base.hsd import HsdBlock

def new_options():
    out = HsdBlock(name="Options", val=None, block_type="method", level=0)
    
    out.status = True
    
    out.scalar["WriteDetailedOut"] = "Yes"

    return out
