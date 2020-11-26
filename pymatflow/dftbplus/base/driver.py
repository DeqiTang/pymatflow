"""

"""

import os 

from pymatflow.dftbplus.base.hsd import hsd_block

def new_driver():
    out = hsd_block(name="Driver", val=None, block_type="method", level=0)
    out.method["SteepestDescent"] = hsd_block(name="SteepestDescent", val=None, block_type="method", level=1)
    out.method["ConjugateGradient"] = hsd_block(name="ConjugateGradient", val=None, block_type="method", level=1)
    out.method["gDIIS"] = hsd_block(name="gDIIS", block_type="method", level=1) 
    out.method["LBFGS"] = hsd_block(name="LBFGS", block_type="method", level=1)
    out.method["SecondDerivatives"] = hsd_block(name="SecondDerivatives", block_type="method", level=1)
    out.method["VelocityVerlet"] = hsd_block(name="VelocityVerlet", block_type="method", level=1)
    out.method["Socket"] = hsd_block(name="Socket", block_type="method", level=1)

    out.status = True
    out.method["LBFGS"].status = True
    return out
