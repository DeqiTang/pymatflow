"""

"""

import os 

from pymatflow.dftbplus.base.hsd import HsdBlock

def new_electrondynamics():
    out = HsdBlock(name="ElectronDynamics", val=None, block_type="method", level=0)
    
    out.status = False
    
    out.scalar["WriteRestart"] = "Yes"

    return out
