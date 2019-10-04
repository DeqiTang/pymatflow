#!/usr/bin/env python
# _*_ coding: utf-8 _*_

import os
import shutil
import pymatgen as mg
import matplotlib.pyplot as plt

from emuhelper.qe.base.control import qe_control
from emuhelper.qe.base.system import qe_system
from emuhelper.qe.base.electrons import qe_electrons
from emuhelper.qe.base.ions import qe_ions
from emuhelper.qe.base.arts import qe_arts


class static_run:
    """
    """
    def __init__(self, xyz_f):
        self.control = qe_control()
        self.system = qe_system()
        self.electrons = qe_electrons()
        #self.ions = qe_ions()
        self.arts = qe_arts(xyz_f)

        self.control.params["calculation"] = 'scf'
        self.control.params["outdir"] = "./tmp"
        self.control.params["pseudo_dir"] = "./"
        self.control.params["wf_collect"] = ".true."
        self.system.params["ibrav"] = 0
        self.system.params["nat"] = self.arts.xyz.natom
        self.system.params["ntyp"] = self.arts.xyz.nspecies
        self.system.params["ecutwfc"] = 100
        self.system.params["input_DFT"] = 'PBE'
        self.system.params["occupations"] = 'smearing'
        self.system.params["smearing"] = "gaussian"
        self.system.params["degauss"] = 0.0001

    def gen_input(self, directory="tmp-static-qe", inpname="static.in"):
        """
        directory: a place for all the generated files
        """
        if os.path.exists(directory):
            shutil.rmtree(directory)
        os.mkdir(directory)
        
        os.system("cp *.UPF %s/" % directory)

        with open(os.path.join(directory, inpname), 'w') as fout:
            self.control.to_in(fout)
            self.system.to_in(fout)
            self.electrons.to_in(fout)
            self.arts.to_in(fout)
    
    def run(self, directory="tmp-static-qe", inpname="static.in", output="static.out"):
        os.chdir(directory)
        os.system("pw.x < %s | tee %s" % (inpname, output))
        os.chdir("../")
    
    def converge_ecutwfc(self, emin, emax, step, directory="tmp-ecutwfc-qe"):
        if os.path.exists(directory):
            shutil.rmtree(directory)
        os.mkdir(directory)
        os.system("cp *.UPF %s/" % directory)
    
        os.chdir(directory)
        n_test = int((emax - emin) / step)
        for i in range(n_test + 1):
            ecut_wfc = int(emin + i * step)
            ecut_rho = ecut_wfc * 4 # using default value for ecut_rho: 4 * ecutwfc
            inp_name = "ecutwfc-%d.in" % ecut_wfc
            self.control.params['outdir'] = './tmp-' + str(ecut_wfc)
            self.system.params['ecutwfc'] = ecut_wfc
            with open(inp_name, 'w') as fout:
                self.control.to_in(fout)
                self.system.to_in(fout)
                self.electrons.to_in(fout)
                self.arts.to_in(fout)
        # run the simulation
        for i in range(n_test + 1):
            ecut_wfc = int(emin + i * step)
            inp_name = "ecutwfc-%d.in" % ecut_wfc
            out_f_name = "ecutwfc-%d.out" % ecut_wfc
            os.system("pw.x < %s | tee %s" % (inp_name, out_f_name))

        # analyse the result
        for i in range(n_test + 1):
            ecut_wfc = int(emin + i * step)
            out_f_name = "ecutwfc-%d.out" % ecut_wfc
            os.system("cat %s | grep '!    total energy' >> energy-ecutwfc.data" % out_f_name)

        ecut_wfc_all = [emin + i * step for i in range(n_test + 1)]
        energy_all = []
        with open("energy-ecutwfc.data", 'r') as fin:
            for line in fin:
                energy_all.append(float(line.split()[4]))

        plt.plot(ecut_wfc_all, energy_all)
        plt.show()
        os.chdir("../")

        
    def converge_ecutrho(self, emin, emax, step, ecutwfc, directory="tmp-ecutrho-qe"):
        if os.path.exists(directory):
            shutil.rmtree(directory)
        os.mkdir(directory)
        os.system("cp *.UPF %s/" % directory)
    
        os.chdir(directory)
        n_test = int((emax - emin) / step)
        for i in range(n_test + 1):
            ecut_rho = int(emin + i * step)
            inp_name = "ecutrho-%d.in" % ecut_rho
            self.control.params['outdir'] = './tmp-' + str(ecut_rho)
            self.system.params['ecutwfc'] = ecutwfc
            self.system.params["ecutrho"] = ecut_rho
            with open(inp_name, 'w') as fout:
                self.control.to_in(fout)
                self.system.to_in(fout)
                self.electrons.to_in(fout)
                self.arts.to_in(fout)
        # run the simulation
        for i in range(n_test + 1):
            ecut_rho = int(emin + i * step)
            inp_name = "ecutrho-%d.in" % ecut_rho
            out_f_name = "ecutrho-%d.out" % ecut_rho
            os.system("pw.x < %s | tee %s" % (inp_name, out_f_name))
        # analyse the result
        for i in range(n_test + 1):
            ecut_rho = int(emin + i * step)
            out_f_name = "ecutrho-%d.out" % ecut_rho
            os.system("cat %s | grep '!    total energy' >> energy-ecutrho.data" % out_f_name)

        ecut_rho_all = [emin + i * step for i in range(n_test + 1)]
        energy_all = []
        with open("energy-ecutrho.data", 'r') as fin:
            for line in fin:
                energy_all.append(float(line.split()[4]))

        plt.plot(ecut_rho_all, energy_all)
        plt.show()
        os.chdir("../")
    #
