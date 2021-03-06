"""
in control of &control /
"""
import sys

"""
usage:
"""

class qe_control:
    """

    """
    def __init__(self):
        self.params = {
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
                }

    def to_in(self, fout):
        """
        :param fout: a file stream for writing
        """
        fout.write("&control\n")
        for item in self.params:
            if self.params[item] is not None:
                if type(self.params[item]) == str:
                    if self.params[item] == ".true." or self.params[item] == ".false.":
                        fout.write("%s = %s\n" % (item, str(self.params[item])))
                    else:
                        fout.write("%s = '%s'\n" % (item, str(self.params[item])))
                else:
                    fout.write("%s = %s\n" % (item, str(self.params[item])))
        fout.write("/\n")
        fout.write("\n")
    
    def calculation(self, calc="scf"):
        self.params["calculation"] = calc

    def basic_setting(self, calc="scf"):
        """
            do a basic setting for all kinds of calculation
        """
        self.params["prefix"] = 'pwscf'
        self.calculation(calc)
        if calc == "scf":
            self.params["outdir"] = "./tmp"
            self.params["pseudo_dir"] = "./"
            self.params["wf_collect"] = True
        elif calc == "nscf":
            self.params["outdir"] = "./tmp"
            self.params["pseudo_dir"] = "./tmp"
            self.params["wf_collect"] = ".true."
        elif calc == "bands":
            self.params["outdir"] = "./tmp"
            self.params["pseudo_dir"] = "./"
            self.params["wf_collect"] = ".true."
        elif calc == "relax":
            self.params["outdir"] = "./tmp"
            self.params["pseudo_dir"] = "./"
            self.params["wf_collect"] = ".true."
        elif calc == "md":
            self.params["outdir"] = "./tmp"
            self.params["pseudo_dir"] = "./"
            self.params["wf_collect"] = ".true."
        elif calc == "vc-relax":
            self.params["outdir"] = "./tmp"
            self.params["pseudo_dir"] = "./"
            self.params["wf_collect"] = ".true."
        elif calc == "vc-md":
            self.params["outdir"] = "./tmp"
            self.params["pseudo_dir"] = "./"
            self.params["wf_collect"] = ".true."

    def set_params(self, params):
        """
        :param params: a dict storing the parameters and values
        """
        for item in params:
            self.params[item] = params[item]
