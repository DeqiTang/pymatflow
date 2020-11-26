"""

"""

import os 

from pymatflow.dftbplus.base.hsd import hsd_block

def new_electrondynamics():
    out = hsd_block(name="ElectronDynamics", val=None, block_type="method", level=0)
    
    out.status = False
    
    out.scalar["WriteRestart"] = "Yes"

    return out
