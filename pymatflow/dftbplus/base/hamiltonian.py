"""

"""

import os 

class dftb_mixer:
    def __init__(self):
        pass
    def to_string(self):
        pass


class dftb:
    """
    """
    def __init__(self):
        self.params = {
                "SCC": "Yes",
                "SCCTolerance": 1.0e-5,
                "MaxSCCIterations": 100,
        }
        self.mixer = dftb_mixer()
        self.maxangularmomentum = dftb_maxangularmomentum()
        self.spinpolarisation = dftb_spinpolarisation()
        self.spinconstants = dftb_spin_constants()
        self.spinorbit = dftb_spinorbit()
        self.solver = dftb_solver()
        self.filling = dftb_filling()
        self.nonaufbau = dftb_nonaufbau()
        self.slaterkosterfiles = dftb_slaterkosterfiles()
        self.polynomialrepulsive = dftb_polynoialrepulsive()
        self.kpointsandweights = dftb_kpointsandweights()
        self.orbitalpotential = dftb_orbitalpotential()
        self.initialcharge = dftb_initialcharge()
        self.electricfield = dftb_electricfield()
        self.dispersion = dftb_dispersion()

    def to_string(self):
        out = ""



class hamiltonian:
    """
    """
    def __init__(self):
        self.method = "DFTB" # currently only dftb implemented in dftb+ 
        self.methods = {}
        self.dftb = dftb()

    def to_string(self):
        pass
        
