"""

"""

import os 

from pymatflow.dftbplus.base.hsd import HsdBlock

def new_driver():
    out = HsdBlock(name="Driver", val=None, block_type="method", level=0)

    out.method["Thermostat"] = HsdBlock(name="Thermostat", val=None, block_type="method", level=1)

    out.status = True
    return out
