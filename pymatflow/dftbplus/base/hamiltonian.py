"""

"""

import os 

from pymatflow.dftbplus.base.hsd import hsd_block

def new_hamiltonian():
    out = hsd_block(name="Hamiltonian", val="DFTB", block_type="method", level=0)
    out.method["Mixer"] = hsd_block(name="Mixer", val="Broyden", block_type="method", level=1)
    out.method["SpinPolarisation"] = hsd_block(name="SpinPolarisation", val=None, block_type="method", level=1)
    out.method["SpinOrbit"] = hsd_block(name="SpinOrbit", block_type="method", level=1) 

    out.method["SlaterKosterFiles"] = hsd_block(name="SlaterKosterFiles", val="Type2FileNames", block_type="method", level=1)
    out.method["SlaterKosterFiles"].scalar["Prefix"] = "./"
    out.method["SlaterKosterFiles"].scalar["Separator"] = "-"
    out.method["SlaterKosterFiles"].scalar["Suffix"] = ".skf"
    out.method["SlaterKosterFiles"].scalar["LowerCaseTypeName"] = "No"

    out.list_of_property["MaxAngularMomentum"] = hsd_block(name="MaxAngularMomentum", block_type="property", level=1)

    out.scalar["Scc"] = "Yes"
    out.scalar["SCCTolerance"] = 1.0e-5
    out.scalar["Charge"] = 0

    out.method["KPointsAndWeights"] = hsd_block(name="KPointsANdWeights", val="SupercellFolding", block_type="method", level=1)
    

    out.status = True
    out.method["Mixer"].status =True
    out.method["SlaterKosterFiles"].status = True
    out.list_of_property["MaxAngularMomentum"].status = True

    out.method["KPointsAndWeights"].status = True

    return out
