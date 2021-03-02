"""

"""

import os 

from pymatflow.dftbplus.base.hsd import HsdBlock

def new_parallel():
    out = HsdBlock(name="Parallel", val=None, block_type="method", level=0)
    
    out.status = False
    
    out.scalar["Groups"] = 1

    return out
