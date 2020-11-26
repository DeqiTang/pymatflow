"""

"""

import os 

from pymatflow.dftbplus.base.hsd import hsd_block

def new_parseroptions():
    out = hsd_block(name="ParserOptions", val=None, block_type="method", level=0)
    
    out.status = False
    
    out.scalar["ParserVersion"] = 7

    return out
