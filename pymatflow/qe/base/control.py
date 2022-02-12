"""
in control of &control /
"""
import sys

from pymatflow.qe.group import QeVariableGroup
from . import qe_variable_to_string

"""
usage:
"""

class QeControl(QeVariableGroup):
    """

    """
    def __init__(self):
        super().__init__()
        #self.incharge 
        self.set_params({
            "calculation": None,
            "title": None,
            "verbosity": None,
            "restart_mode": None,
            "wf_collect": None,
            "nstep": None,
            "iprint": None,
            "tstress": None,
            "tprnfor": None,
            "dt": None,
            "outdir": None,
            "wfcdir": None,
            "prefix": None,
            "lkpoint_dir": None,
            "max_seconds": None,
            "etot_conv_thr": None,
            "forc_conv_thr": None,
            "disk_io": None,
            "pseudo_dir": None,
            "tefield": None,
            "dipfield": None,
            "lelfield": None,
            "nberrycyc": None,
            "lorbm": None,
            "lberry": None,
            "gdir": None,
            "nppstr": None,
            "lfcpopt": None,
            "gate": None,
        })

    def to_in(self, fout):
        """
        :param fout: a file stream for writing
        """
        fout.write(self.to_string())
    
    def to_string(self):
        out = ""
        out += "&control\n"
        for item in self.params:
            if self.params[item].as_val() == None:
                continue
            out += qe_variable_to_string(self.params[item])
            out += "\n"
        out += "/\n"
        out += "\n"
        return out
    
    def calculation(self, calc="scf"):
        self.set_param("calculation", calc)

    def basic_setting(self, calc="scf"):
        """
            do a basic setting for all kinds of calculation
        """
        self.set_param("prefix", 'pwscf')
        self.calculation(calc)
        if calc == "scf":
            self.set_param("outdir", "./tmp")
            self.set_param("pseudo_dir", "./")
            self.set_param("wf_collect", ".true.")
        elif calc == "nscf":
            self.set_param("outdir", "./tmp")
            self.set_param("pseudo_dir", "./tmp")
            self.set_param("wf_collect", ".true.")
        elif calc == "bands":
            self.set_param("outdir", "./tmp")
            self.set_param("pseudo_dir", "./")
            self.set_param("wf_collect", ".true.")
        elif calc == "relax":
            self.set_param("outdir", "./tmp")
            self.set_param("pseudo_dir", "./")
            self.set_param("wf_collect", ".true.")
        elif calc == "md":
            self.set_param("outdir", "./tmp")
            self.set_param("pseudo_dir", "./")
            self.set_param("wf_collect", ".true.")
        elif calc == "vc-relax":
            self.set_param("outdir", "./tmp")
            self.set_param("pseudo_dir", "./")
            self.set_param("wf_collect", ".true.")
        elif calc == "vc-md":
            self.set_param("outdir", "./tmp")
            self.set_param("pseudo_dir", "./")
            self.set_param("wf_collect", ".true.")
