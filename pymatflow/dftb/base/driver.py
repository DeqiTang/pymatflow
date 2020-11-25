"""

"""

import os 

class Driver:
    """
    The driver is responsible for changing the geometry of the input structure during the calculation.
    """
    def __init__(self):
        self.method = "Static" # SteepestDescent, ConjugateGradient, gDIIS, LBFGS, FIRE, SecondDerivatives, VelocityVerlet, Socket
        self.methods = {}

    def to_string(self):
        pass
        
