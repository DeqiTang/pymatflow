"""

"""

import os 

from pymatflow.dftbplus.base.hsd import hsd_block

def new_analysis():
    out = hsd_block(name="Analysis", val=None, block_type="method", level=0)
    
    out.status = False
    
    out.scalar["CalculateForces"] = "No"

    return out
