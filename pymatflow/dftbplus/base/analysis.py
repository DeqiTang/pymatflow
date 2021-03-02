"""

"""

import os 

from pymatflow.dftbplus.base.hsd import HsdBlock

def new_analysis():
    out = HsdBlock(name="Analysis", val=None, block_type="method", level=0)
    
    out.status = False
    
    out.scalar["CalculateForces"] = "No"

    return out
