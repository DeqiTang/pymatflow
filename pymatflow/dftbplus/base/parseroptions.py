"""

"""

import os 

from pymatflow.dftbplus.base.hsd import HsdBlock

def new_parseroptions():
    out = HsdBlock(name="ParserOptions", val=None, block_type="method", level=0)
    
    out.status = False
    
    out.scalar["ParserVersion"] = 7

    return out
