#!/usr/bin/env python
# _*_ coding: utf-8 _*_

import os
import re
import shutil
import pymatgen as mg
import matplotlib.pyplot as plt

from pymatflow.qe.pwscf import pwscf
from pymatflow.qe.base.arts import qe_arts

class neb_run(pwscf):
    """
    Reference:
        http://www.quantum-espresso.org/Doc/INPUT_NEB.html
    Note:
        check the official manual for neb for some knowledge 
        about using neb.x (which is helpfule!)

        A gross estimate of the required number of iterations is 
        (number of images) * (number of atoms) * 3. 
        Atoms that do not move should not be counted

        the neb calculation is usually difficult to converge,
        and experience is needed.
    Q&A:
        can we fix some atoms when doing neb calculation?

    """
    def __init__(self):
        """
        """
        super().__init__()
        self.images = []
        self.path = {} # Namelist: &PATH
        self.set_basic_path()
        
    def get_images(self, images):
        """
        images:
            ["first.xyz", "intermediate-1.xyz", "intermediate-2.xyz", ...,"last.xyz"]
        self.images containe all the images while self.arts only contains the first image.
        self.arts is provided by base class pwscf, and is used to set kpoints.
        """
        self.images = []
        for image in images:
            arts = qe_arts()
            arts.xyz.get_xyz(image)
            self.images.append(arts)
        self.system.basic_setting(self.images[0])
        for image in self.images:
            image.basic_setting(ifstatic=True)
        # self.arts is actually the same as self.images[0]
        # it is used to set kpoints convinently
        self.arts.xyz.get_xyz(self.images[0].xyz.file)

    def set_params(self, control={}, system={}, electrons={}):

        super().set_params(control=control, system=system, electrons=electrons)

        # must set "wf_collect = False", or it will come across with erros in davcio
        # Error in routine davcio (10): 
        # error while reading from file ./tmp/pwscf_2/pwscf.wfc1
        self.control.params["wf_collect"] = False


    def set_path(self, path={}):
        for item in path:
            self.path[item] = path[item]

    def neb(self, directory="tmp-qe-neb", inpname="neb.in", output="neb.out", mpi="", runopt="gen", restart_mode="from_scratch"):
        """
        directory: a place for all the generated files
        """
        if runopt == "gen" or runopt == "genrun":
            if restart_mode == "from_scratch":
                self.path["restart_mode"] = restart_mode
                if os.path.exists(directory):
                    shutil.rmtree(directory)
                os.mkdir(directory)
                os.system("cp *.UPF %s/" % directory)
                for art in self.images:
                    os.system("cp %s %s/" % (art.xyz.file, directory))
            elif restart_mode == "restart":
                self.path["restart_mode"] = restart_mode
                # first check whether there is a previous neb running
                if not os.path.exists(directory):
                    print("===================================================\n")
                    print("                 Warning !!!\n")
                    print("===================================================\n")
                    print("restart neb calculation:\n")
                    print("  directory of previous neb calculattion not found!\n")
                    sys.exit(1)
                # this assumes the previous neb run and the current neb run are using the same inpname
                os.chdir(directory)
                os.system("mv %s %s.old" % (inpname, inpname)) 
                os.chdir("../")

            with open(os.path.join(directory, inpname), 'w') as fout:
                fout.write("BEGIN\n")
                fout.write("BEGIN_PATH_INPUT\n")
                fout.write("&PATH\n")
                for item in self.path:
                    if self.path[item] is not None:
                        if type(self.path[item]) == str:
                            fout.write("%s = '%s'\n" % (item, self.path[item]))
                        else:
                            fout.write("%s = %s\n" % (item, str(self.path[item])))

                fout.write("/\n")
                fout.write("END_PATH_INPUT\n")
                fout.write("BEGIN_ENGINE_INPUT\n")
                self.control.to_in(fout)
                self.system.to_in(fout)
                self.electrons.to_in(fout)
                self.images_to_neb(fout)
                fout.write("END_ENGINE_INPUT\n")
                fout.write("END\n")
            # gen yhbatch script
            self.gen_yh(directory=directory, inpname=inpname, output=output)

        if runopt == "run" or runopt == "genrun":
            os.chdir(directory)
            os.system("%s neb.x -inp %s | tee %s" % (mpi, inpname, output))
            os.chdir("../")

    def images_to_neb(self, fout):
        # fout: a file stream for writing
        fout.write("ATOMIC_SPECIES\n")
        for element in self.arts.xyz.specie_labels:
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
        cell = self.arts.xyz.cell
        fout.write("CELL_PARAMETERS angstrom\n")
        fout.write("%.9f %.9f %.9f\n" % (cell[0][0], cell[0][1], cell[0][2]))
        fout.write("%.9f %.9f %.9f\n" % (cell[1][0], cell[1][1], cell[1][2]))
        fout.write("%.9f %.9f %.9f\n" % (cell[2][0], cell[2][1], cell[2][2]))
        fout.write("\n")
        # writing KPOINTS to the fout
        self.arts.write_kpoints(fout)
        fout.write("\n")
        # =========================
        fout.write("BEGIN_POSITIONS\n")
        fout.write("\n")
        fout.write("FIRST_IMAGE\n")
        fout.write("ATOMIC_POSITIONS angstrom\n")
        for atom in self.images[0].xyz.atoms:
            fout.write("%s\t%.9f\t%.9f\t%.9f\n" % (atom.name, atom.x, atom.y, atom.z))
        fout.write("\n")
        for i in range(1, len(self.images) - 1):
            fout.write("INTERMEDIATE_IMAGE\n")
            fout.write("ATOMIC_POSITIONS angstrom\n")
            for atom in self.images[i].xyz.atoms:
                fout.write("%s\t%.9f\t%.9f\t%.9f\n" % (atom.name, atom.x, atom.y, atom.z))
            fout.write("\n")
        fout.write("LAST_IMAGE\n")
        fout.write("ATOMIC_POSITIONS angstrom\n")
        for atom in self.images[-1].xyz.atoms:
            fout.write("%s\t%.9f\t%.9f\t%.9f\n" % (atom.name, atom.x, atom.y, atom.z))
        fout.write("\n")
        fout.write("END_POSITIONS\n")


    def gen_yh(self, directory, inpname, output):
        """
        generating yhbatch job script for calculation
        """
        with open(os.path.join(directory, inpname.split(".in")[0]+".sub"), 'w') as fout:
            fout.write("#!/bin/bash\n")
            fout.write("yhrun -N 1 -n 24 neb.x -inp %s > %s\n" % (inpname, output))

    def set_basic_path(self):
        self.path["string_method"] = 'neb'
        self.path["nstep_path"] = 100
        self.path["opt_scheme"] = 'broyden'
        self.path["num_of_images"] = 5
        self.path["k_max"] = 0.3
        self.path["k_min"] = 0.2
        self.path["CI_scheme"] = 'auto'
        self.path["path_thr"] = 0.05e0
        self.path["ds"] = 1.0e0


