#!/usr/bin/env python
# _*_ coding: utf-8 _*_

import os
import re
import shutil
import pymatgen as mg
import matplotlib.pyplot as plt

from emuhelper.qe.base.control import qe_control
from emuhelper.qe.base.system import qe_system
from emuhelper.qe.base.electrons import qe_electrons
from emuhelper.qe.base.arts import qe_arts


class neb_run:
    """
    Reference:
        http://www.quantum-espresso.org/Doc/INPUT_NEB.html
    """
    def __init__(self, image1, image2, image3):
        self.control = qe_control()
        self.system = qe_system()
        self.electrons = qe_electrons()
        self.arts1 = qe_arts(image1)
        self.arts2 = qe_arts(image2)
        self.arts3 = qe_arts(image3)
        
        self.control.basic_setting("scf") 
        self.system.basic_setting(self.arts1)
        self.electrons.basic_setting()
        
    def neb(self, directory="tmp-qe-neb", inpname="neb.in", output="neb.out", 
            mpi="", runopt="gen", control={}, system={}, electrons={}, kpoints_mp=[1, 1, 1, 0, 0, 0]):
        """
        directory: a place for all the generated files
        """
        if runopt == "gen" or runopt == "genrun":
            if os.path.exists(directory):
                shutil.rmtree(directory)
            os.mkdir(directory)
            os.system("cp *.UPF %s/" % directory)
            
            self.control.set_params(control)
            self.system.set_params(system)
            self.electrons.set_params(electrons)
            self.arts1.set_kpoints(kpoints_mp) # use arts1 to set kpoints and cells and species
            # must set "wf_collect = False", or it will come across with erros in davcio
            # Error in routine davcio (10): 
            # error while reading from file ./tmp/pwscf_2/pwscf.wfc1
            self.control.params["wf_collect"] = False

            with open(os.path.join(directory, inpname), 'w') as fout:
                fout.write("BEGIN\n")
                fout.write("BEGIN_PATH_INPUT\n")
                fout.write("&PATH\n")
                fout.write("string_method = 'neb'\n")
                fout.write("nstep_path = 100\n")
                fout.write("opt_scheme = 'broyden'\n")
                fout.write("num_of_images = 5\n")
                fout.write("k_max = 0.3\n")
                fout.write("k_min = 0.2\n")
                fout.write("CI_scheme = 'auto'\n")
                fout.write("path_thr = %f\n" % 0.05e0)
                fout.write("/\n")
                fout.write("END_PATH_INPUT\n")
                fout.write("BEGIN_ENGINE_INPUT\n")
                self.control.to_in(fout)
                self.system.to_in(fout)
                self.electrons.to_in(fout)
                self.arts_to_neb(fout)
                fout.write("END_ENGINE_INPUT\n")
                fout.write("END\n")

        if runopt == "run" or runopt == "genrun":
            os.chdir(directory)
            os.system("%s neb.x -inp %s | tee %s" % (mpi, inpname, output))
            os.chdir("../")

    def arts_to_neb(self, fout):
        # fout: a file stream for writing
        fout.write("ATOMIC_SPECIES\n")
        for element in self.arts1.xyz.specie_labels:
            tmp = os.listdir("./")
            pseudo_file = ""
            for f in tmp:
                match_string = "%s\." % element
                match = re.match(match_string, f)
                if match is not None and match.string.split(".")[-1] == 'UPF':
                    pseudo_file = match.string
                    break
            fout.write("%s %f %s\n" % (element, mg.Element(element).atomic_mass, pseudo_file))
        fout.write("\n")
        cell = self.arts1.xyz.cell
        fout.write("CELL_PARAMETERS angstrom\n")
        fout.write("%f %f %f\n" % (cell[0], cell[1], cell[2]))
        fout.write("%f %f %f\n" % (cell[3], cell[4], cell[5]))
        fout.write("%f %f %f\n" % (cell[6], cell[7], cell[8]))
        fout.write("\n")
        # writing KPOINTS to the fout
        self.arts1.write_kpoints(fout)
        fout.write("\n")
        # =========================
        fout.write("BEGIN_POSITIONS\n")
        fout.write("FIRST_IMAGE\n")
        fout.write("ATOMIC_POSITIONS angstrom\n")
        for atom in self.arts1.xyz.atoms:
            fout.write("%s\t%f\t%f\t%f\n" % (atom.name, atom.x, atom.y, atom.z))
        fout.write("\n")
        fout.write("INTERMEDIATE_IMAGE\n")
        fout.write("ATOMIC_POSITIONS angstrom\n")
        for atom in self.arts2.xyz.atoms:
            fout.write("%s\t%f\t%f\t%f\n" % (atom.name, atom.x, atom.y, atom.z))
        fout.write("\n")
        fout.write("LAST_IMAGE\n")
        fout.write("ATOMIC_POSITIONS angstrom\n")
        for atom in self.arts3.xyz.atoms:
            fout.write("%s\t%f\t%f\t%f\n" % (atom.name, atom.x, atom.y, atom.z))
        fout.write("\n")
        fout.write("END_POSITIONS\n")
