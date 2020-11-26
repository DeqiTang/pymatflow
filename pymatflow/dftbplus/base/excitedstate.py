"""

"""

import os 

from pymatflow.dftbplus.base.hsd import hsd_block

def new_excitedstate():
    out = hsd_block(name="ExcitedState", val=None, block_type="method", level=0)
    
    out.status = False
    
    #out.scalar[""] = 

    return out
