from pymatflow.abinit.namelist import OpticNamelist


class Optic:
    def __init__(self):
        super().__init__()
        self.files = OpticNamelist("FILES")
        self.parameters = OpticNamelist("PARAMETERS")
        self.computations = OpticNamelist("COMPUTATIONS")
        
        self._initialize()

    def _initialize(self):
        self.files.set_params({
            "ddkfile_1" : 'optic-runo_DS4_1WF16',
            "ddkfile_2" : 'optic-runo_DS5_1WF17',
            "ddkfile_3" : 'optic-runo_DS6_1WF18',
            "wfkfile": 'optic-runo_DS3_WFK'
        })

        self.parameters.set_params({
            "broadening": 0.002,
            "domega": 0.0003,
            "maxomega": 0.3,
            "scissor": 0.000,
            "tolerance": 0.002
        })

        self.computations.set_params({
            "num_lin_comp": 1,
            "lin_comp": 11,
            "num_nonlin_comp": 2,
            "nonlin_comp": [123, 222],
            "num_nonlin2_comp": 0
        })

    def to_string(self):
        out = ""
        out += self.files.to_string() +"\n"
        out += self.parameters.to_string() +"\n"
        out += self.computations.to_string() +"\n"
        return out