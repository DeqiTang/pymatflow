from pymatflow.qe.pwscf import PwScf
from pymatflow.qe.dfpt import DfptRun 

class Qe:
    def __init__(self):
        self.pwscf = PwScf()
        self.dfpt = DfptRun()

