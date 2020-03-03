"""
Phonopy calculation
"""
import sys
import re
import os
import shutil
import numpy as np


from pymatflow.remote.server import server_handle
from pymatflow.qe.pwscf import pwscf

"""
Note:
    At present, only crystal type ATOMIC_POSITIONS is supported
    by phonopy:
    reference:
    https://atztogo.github.io/phonopy/qe.html
"""




class phonopy_run(pwscf):
    """
    Note:
        kpoints as well as energy cutoff both will have a significant
        influence on the precision of phonopy running.
        like when you have a small cell and use a small kpoints might
        result in a really useless phonon band. this is just as what
        we do in scf or opt running, generally a small structure requires
        large kpoint sets, while the larger one requires smaller kpoint
        sets.
    """
    def __init__(self):
        super().__init__()

        self.control.basic_setting("scf")

        self.supercell_n = [1, 1, 1]

        #must print print out forces and stress after scf
        # so that phonopy can parse the scf output file and
        # construct the FORCE CONSTANT MATIRX
        self.control.params["tprnfor"] = True
        self.control.params["tstress"] = True

    def phonopy(self, directory="tmp-qe-phonopy", pos_inpname="pos.in", head_inpname="head.in", mpi="", runopt="gen", auto=0):
        """
        :param directory: a place for all the generated files
        """
        if runopt == "gen" or runopt == "genrun":
            if os.path.exists(directory):
                shutil.rmtree(directory)
            os.mkdir(directory)

            #os.system("cp *.UPF %s/" % directory)
            #os.system("cp %s %s/" % (self.arts.xyz.file, directory))

            # do not copy too many files at the same time or it will be slow
            # so we do not copy all UPF files in the directory but just copy
            # those used in the calculation.
            shutil.copyfile(self.arts.xyz.file, os.path.join(directory, os.path.basename(self.arts.xyz.file)))
            all_upfs = [s for s in os.listdir() if s.split(".")[-1] == "UPF"]
            for element in self.arts.xyz.specie_labels:
                for upf in all_upfs:
                    if upf.split(".")[0] == element:
                        shutil.copyfile(upf, os.path.join(directory, upf))
                        break
            self.arts.pseudo.dir = os.path.abspath(directory)
            self.control.pseudo_dir = os.path.abspath(directory)
            #


            with open(os.path.join(directory, head_inpname), 'w') as fout:
                self.control.to_in(fout)
                self.system.to_in(fout)
                self.electrons.to_in(fout)
                self.arts.write_kpoints(fout)

            # set up the Phonopy calculation
            os.chdir(directory)
            os.system("cat %s > %s" % (head_inpname, pos_inpname))
            with open(pos_inpname, 'a') as fout:
                self.arts.to_in(fout, coordtype="crystal")
            os.system("phonopy --qe -d --dim='%d %d %d' -c %s" % (self.supercell_n[0], self.supercell_n[1], self.supercell_n[2], pos_inpname))
            os.system("ls | grep 'supercell-' > pos.data")
            disp_dirs = []
            with open("pos.data", 'r') as fin:
                for line in fin:
                    disp_dirs.append(line.split(".")[0].split("-")[1])
            # IMPORTANT:
            # we need to overwrite the head_inpname file, reasons are bellow
            # here we must pay attention to system.params["nat"]
            # as system.params["nat"] is the value of the original structure
            # but after phonopy build the supercell the number of atoms changes
            # accordingly. so we must change the value of system.params["nat"]
            # according to self.supercell_n
            self.system.params["nat"] = self.system.params["nat"] * self.supercell_n[0] * self.supercell_n[1] * self.supercell_n[2]
            with open(head_inpname, 'w') as fout:
                self.control.to_in(fout)
                self.system.to_in(fout)
                self.electrons.to_in(fout)
                self.arts.write_kpoints(fout)
            #
            for disp in disp_dirs:
                os.system("cat %s supercell-%s.in > supercell-%s-full.in" % (head_inpname, disp, disp))
                os.system("rm supercell-%s.in" % disp)
            os.chdir("../")
            # end build the phonopy

            # gen yhbatch script
            with open(os.path.join(directory, "phonopy-job.sub"), 'w') as fout:
                fout.write("#!/bin/bash\n")
                fout.write("\n")
                for disp in disp_dirs:
                    fout.write("yhrun -N 1 -n 24 pw.x < supercell-%s-full.in > supercell-%s.out\n" % (disp, disp))
            # gen pbs script
            with open(os.path.join(directory, "phonopy-job.pbs"), 'w') as fout:
                fout.write("#!/bin/bash\n")
                fout.write("#PBS -N %s\n" % self.run_params["jobname"])
                fout.write("#PBS -l nodes=%d:ppn=%d\n" % (self.run_params["nodes"], self.run_params["ppn"]))
                fout.write("cd $PBS_O_WORKDIR\n")
                fout.write("NP=`cat $PBS_NODEFILE | wc -l`\n")
                fout.write("\n")
                for disp in disp_dirs:
                    fout.write("mpirun -np $NP -machinefile $PBS_NODEFILE pw.x < supercell-%s-full.in > supercell-%s.out\n" % (disp, disp))

        if runopt == "run" or runopt == "genrun":
            os.chdir(directory)
            # run the dft
            disp_dirs = []
            with open("pos.data", 'r') as fin:
                for line in fin:
                    disp_dirs.append(line.split(".")[0].split("-")[1])
            for disp in disp_dirs:
                os.system("%s pw.x < supercell-%s-full.in | tee supercell-%s.out" % (self.run_params["mpi"], disp, disp))
            os.chdir("../")
        server_handle(auto=auto, directory=directory, jobfilebase="phonopy-job", server=self.run_params["server"])
    #
    #

    def phonopy_qha(self, directory="tmp-qe-phonopy-qha", pos_inpname="pos.in", head_inpname="head.in", mpi="", runopt="gen", auto=0):
        """
        :param directory: a place for all the generated files
        Refer to https://phonopy.github.io/phonopy/qha.html
        """
        pass
