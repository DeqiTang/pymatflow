"""

"""

import os 

from pymatflow.dftbplus.base.hsd import HsdBlock

def new_excitedstate():
    out = HsdBlock(name="ExcitedState", val=None, block_type="method", level=0)
    
    out.status = False
    
    #out.scalar[""] = 

    return out
