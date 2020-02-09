from pymatflow.qe.pwscf import pwscf
from pymatflow.qe.dfpt import dfpt_run 

class qe:
    def __init__(self):
        self.pwscf = pwscf()
        self.dfpt = dfpt_run()

