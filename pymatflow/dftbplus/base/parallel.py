"""

"""

import os 

from pymatflow.dftbplus.base.hsd import hsd_block

def new_parallel():
    out = hsd_block(name="Parallel", val=None, block_type="method", level=0)
    
    out.status = False
    
    out.scalar["Groups"] = 1

    return out
