

from pymatflow.vasp.base.start import start_incharge
from pymatflow.vasp.base.electrons import electrons_incharge
from pymatflow.vasp.base.ions import ions_incharge
from pymatflow.vasp.base.write import write_incharge
from pymatflow.vasp.base.dipolecorrection import dipolecorrection_incharge
from pymatflow.vasp.base.xc import xc_incharge
from pymatflow.vasp.base.lr import lr_incharge
from pymatflow.vasp.base.orbitalmagnet import orbitalmagnet_incharge


from pymatflow.vasp.group import VaspVariableGroup

from . import vasp_variable_to_string

neb_incharge =  ["IOPT", "SPRING", "IMAGES", "LCLIMB", "ICHAIN", "LNEBCELL",
    "MAXMOVE", "LGLOBAL", "LAUTOSCALE", "INVCURV", "LLINEOPT", "FDSTEP"]

"""
in misc now, waiting for further classification:
    ISPIN, MAGMOM, LNONCOLLINEAR, LSORBIT,
    ALGO,
"""

class VaspIncar(VaspVariableGroup):
    """
    """
    def __init__(self):
        """
        """
        super().__init__()

        self.runtype = None


    def set_runtype(self, runtype="static"):
        """
        self.runtype: static | opt | md | dfpt | phonon | neb | custom
        """
        self.runtype = runtype
        self.basic_setting()

    def to_incar(self, fout):
        fout.write(self.to_string())
        
    def to_string(self):
        incar_out = ""

        start = "# =============================\n# Start parameter\n# =============================\n"
        electrons = "# =============================\n# Electrons related parameters\n# ======================================\n"
        ions = "# =============================\n# Ions related parameters\n# =============================\n"
        write = "# =============================\n# Write flags\n# =============================\n"
        dipolecorrection = "# =============================\n# Dipole correction related parameters\n# =============================\n"
        xc = "# =============================\n# Exchange correlation related parameters\n# =============================\n"
        lr = "# =============================\n# Linear response parameters\n# =============================\n"
        orbitalmagnet = "# =============================\n# Orbitalmagnet related parameters\n# =============================\n"
        misc ="# =============================\n# Miscellaneous parameters\n# =============================\n"
        neb = "# =============================\n# Neb related parameters\n# ================================\n"


        for item in self.params:
            if self.params[item].as_val() == None:
                continue
            if item in start_incharge:
                start = start + vasp_variable_to_string(self.params[item]) + "\n\n"
            elif item in electrons_incharge:
                electrons = electrons + vasp_variable_to_string(self.params[item]) + "\n\n"
            elif item in ions_incharge:
                ions = ions + vasp_variable_to_string(self.params[item]) + "\n\n"
            elif item in write_incharge:
                write = write + vasp_variable_to_string(self.params[item]) + "\n\n"
            elif item in dipolecorrection_incharge:
                dipolecorrection = dipolecorrection + vasp_variable_to_string(self.params[item]) + "\n\n"
            elif item in xc_incharge:
                xc = xc + vasp_variable_to_string(self.params[item]) + "\n\n"
            elif item in lr_incharge:
                lr = lr + vasp_variable_to_string(self.params[item]) + "\n\n"
            elif item in orbitalmagnet_incharge:
                orbitalmagnet = orbitalmagnet + vasp_variable_to_string(self.params[item]) + "\n\n"
            elif item in neb_incharge:
                neb = neb + vasp_variable_to_string(self.params[item]) + "\n\n"
            else:
                misc = misc + vasp_variable_to_string(self.params[item]) + "\n\n"
                
        if self.runtype in ["static", "opt", "md", "dfpt", "neb", "phonon", "custom"]:
            incar_out += start
            incar_out += electrons
            incar_out += ions

        incar_out += xc
        incar_out += write
        incar_out += misc

        if self.runtype in ["dfpt", "custom"]:
            incar_out += lr
        if self.runtype in ["neb", "custom"]:
            incar_out += neb

        print(incar_out)
        return incar_out


    def basic_setting(self):
        """
        self.runtype: static | opt | md | dfpt | phonon | neb | custom
        Note:
            basic_setting() will do the basic setting for incar
            according to the value of self.runtype
        """
        self.set_params({
            "ENCUT": 300,
            "EDIFF": 1.0E-6,
            "ISMEAR": 0,
            "NWRITE": 2,
            })
        if self.runtype == "static":
            self.set_params({
                "IBRION": -1,
                })
        elif self.runtype == "opt":
            self.set_params({
                "EDIFFG": 1.0E-3,
                "IBRION": 1, # or 2
                "ISIF": 0, # default is non variable cell
                "NSW": 100,
                })
        elif self.runtype == "md":
            self.set_params({
                "EDIFFG": None,
                "IBRION": 0,
                "POTIM": 0.5, # in md, POTIM means tiem step in unit of fs while in relaxation it means step width
                "NSW": 1000,
                })
        elif self.runtype == "phonon":
            self.set_params({
                "IBRION": 7, # 计算声子: 5, 6, 7, 8
                "NFREE": 4, # 设置为2或4, 推荐不要设置为1
                "POTIM": 0.015, # 默认值为: 0.015
                "NSW": 1, # 如何设置NSW?
                })
        elif self.runtype == "dfpt":
            self.set_params({
                "LEPSILON": "T",
                "LRPA": "T",
                })
        elif self.runtype == "neb":
            self.set_params({
            "ISIF": 2,
            "EDIFFG": -0.01, # when you set IOPT > 0, must set EDIFFG < 0
            "IBRION": 3, # use VTST optimizer
            "POTIM": 0, # use VTST optimizer
            "LCLIMB": "T",
            "IOPT": 1,
            "SPRING": -5,
            })
        elif self.runtype == "neb":
            pass

    def set_properties_calculation(self, option=[]):
        """
        option:
        """
        pass

    def set_md(self, ensemble=0, thermostat=0):
        """
        set md running

        ensemble:
            0: NVE
            1: NVT
            2: NPT
            3: NPH
        thermostat:
            0: Anderson
            1: Nose-Hoover
            2: Langevin
            3: Multiple Anderson
        """
        if ensemble == 0:
            self.set_param("MDALGO", 0)
            self.set_param("SMASS", -3)
        elif ensemble == 1:
            self.set_param("ISIF", 2)
            if thermostat == 0:
                self.set_param("MDALGO", 1)
            elif thermostat == 1:
                self.set_param("MDALGO", 2)
            elif thermostat == 2:
                self.set_param("MDALGO", 3)
            elif thermostat == 3:
                self.set_param("MDALGO", 13)
        elif ensemble == 2:
            self.set_param("MDALGO", 3)
            self.set_param("ISIF", 3)
        elif ensemble == 3:
            self.set_param("MDALGO", 3)
            self.set_param("ISIF", 3)
            self.set_param("LANGEVIN_GAMMA_L", 0.0)
