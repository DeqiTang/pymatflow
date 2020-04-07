"""
Phonopy calculation
"""
import sys
import os
import shutil
import seekpath
import numpy as np

from pymatflow.remote.server import server_handle

from pymatflow.cp2k.base.xyz import cp2k_xyz

from pymatflow.cp2k.cp2k import cp2k
#from pymatflow.cp2k.base.glob import cp2k_glob
#from pymatflow.cp2k.base.force_eval import cp2k_force_eval

"""
Usage:
    python phonon_cp2k.py xxx.xyz
    xxx.xyz is the input structure file

    make sure the xyz structure file and pseudopotential file
    for all the elements of the system is in the directory.

Dependencies:
    pip install --user phonopy
    pip install --user cp2k_tools
    pip install --user cp2k_input_tools

Note:
    phonopy read the xxx.inp and it can only read the system structure
    by COORD specified in SUBSYS. So I can not use TOPOLOGY.
    PLUS: only scaled coordinates are currently supported!

    phonopy and seekpath both use spglib to decide the space group.

References:
    https://www.cp2k.org/exercises:2018_uzh_cmest:phonon_calculation
"""


class phonopy_run(cp2k):
    """
    """
    def __init__(self):
        super().__init__()
        #self.glob = cp2k_glob()
        #self.force_eval = cp2k_force_eval()

        self.glob.basic_setting(run_type="ENERGY_FORCE")
        self.force_eval.basic_setting()

        self.supercell_n = [1, 1, 1]


    def phonopy(self, directory="tmp-cp2k-phonopy", runopt="gen", auto=0):
        self.force_eval.check_spin()
        if runopt == "gen" or runopt == "genrun":
            if os.path.exists(directory):
                shutil.rmtree(directory)
            os.mkdir(directory)

            os.chdir(directory)
            shutil.copyfile(self.force_eval.subsys.xyz.file, "%s" % os.path.basename(self.force_eval.subsys.xyz.file))

            inp_name = "phonon.inp"
            with open(inp_name, 'w') as fout:
                self.glob.to_input(fout)
                #fout.write("\n")
                fout.write("&FORCE_EVAL\n")
                fout.write("\tMETHOD Quickstep\n")
            # subsys
            self.to_subsys_phonopy(inp_name)
            # end subsys
            with open(inp_name, 'a') as fout:
                # dft
                self.force_eval.dft.to_input(fout)
                # end dft
                fout.write("&END FORCE_EVAL\n")


            # build the phonopy running files
            os.system("phonopy --cp2k -c %s -d --dim='%d %d %d'" % (inp_name, self.supercell_n[0], self.supercell_n[1], self.supercell_n[2]))
            # now phonon-supercell-00x.inp is generated which will be used to construct input for cp2k
            # in the past the above command will generate the supercell-00x.inp which can not be the
            # input file of cp2k directly, and need us to construct the final input for cp2k for every
            # displacement. but after an update of phonopy to version 2.4.2. the enerated phonon-supercell-00x.inp
            # can be run by cp2k directly! so we comment those old code for constructing the final input file.
            os.system("ls | grep 'phonon-supercell-' > geo.data")
            disps = []
            with open("geo.data", 'r') as fin:
                for line in fin:
                    disps.append(line.split(".")[0].split("-")[2])

            #for disp in disps:
            #    in_name = "phonon-supercell-%s.inp" % disp
            #    if os.path.exists(in_name) is not True:
            #        break
            #    tmp_file = "phonon-supercell-%s.tmp.txt" % disp
            #    shutil.copyfile(in_name, tmp_file)
            #    # important: different disp calculation should have different PROJECT name
            #    self.glob.params["PROJECT"] = "abinitio" + "-supercell-" + disp
            #    with open(in_name, 'w') as fout:
            #        self.glob.to_input(fout)
            #        fout.write("\n")
            #        fout.write("&FORCE_EVAL\n")
            #        fout.write("\tMETHOD Quickstep\n")
            #        fout.write("\t&SUBSYS\n")
            #    self.print_kinds(in_name)
            #    os.system("cat %s | sed '1d;2d;3d;4d;5d;6d;7d' | sed '$d' | sed '$d' | sed '$d' | sed '$d' | sed '$d' >> %s" % (tmp_file, in_name))
            #    with open(in_name, 'a') as fout:
            #        fout.write("\t&END SUBSYS\n")
            #        # dft
            #        self.force_eval.dft.to_input(fout)
            #        # end dft
            #        fout.write("\t&PRINT\n")
            #        fout.write("\t\t&FORCES\n")
            #        fout.write("\t\t\tFILENAME forces\n")
            #        fout.write("\t\t&END FORCES\n")
            #        fout.write("\t&END PRINT\n")
            #        fout.write("&END FORCE_EVAL\n")
            os.chdir("../")


            #
            # generate yhbatch file
            with open(os.path.join(directory, "phonopy-job.slurm"), 'w') as fout:
                fout.write("#!/bin/bash\n")
                fout.write("#SBATCH -p %s\n" % self.run_params["partition"])
                fout.write("#SBATCH -N %d\n" % self.run_params["nodes"])
                fout.write("#SBATCH -n %d\n" % self.run_params["ntask"])
                fout.write("#SBATCH -J %s\n" % self.run_params["jobname"])
                fout.write("#SBATCH -o %s\n" % self.run_params["stdout"])
                fout.write("#SBATCH -e %s\n" % self.run_params["stderr"])
                for disp in disps:
                    fout.write("yhrun -N 1 -n 24 $PMF_CP2K -in phonon-supercell-%s.inp > phonon-supercell-%s.inp.out\n" % (disp, disp))

            # generate pbs file
            with open(os.path.join(directory, "phonopy-job.pbs"), 'w') as fout:
                fout.write("#!/bin/bash\n\n")
                fout.write("#PBS -N %s\n" % self.run_params["jobname"])
                fout.write("#PBS -l nodes=%d:ppn=%d\n" % (self.run_params["nodes"], self.run_params["ppn"]))
                if "queue" in self.run_params and self.run_params["queue"] != None:
                    fout.write("#PBS -q %s\n" %self.run_params["queue"])                
                fout.write("\n")
                fout.write("cd $PBS_O_WORKDIR\n")
                fout.write("NP=`cat $PBS_NODEFILE | wc -l`\n")
                for disp in disps:
                    fout.write("mpirun -np $NP -machinefile $PBS_NODEFILE $PMF_CP2K -in phonon-supercell-%s.inp > phonon-supercell-%s.inp.out\n" % (disp, disp))


        if runopt == "run" or runopt == "genrun":
            os.chdir(directory)
            disps = []
            with open("geo.data", 'r') as fin:
                for line in fin:
                    disps.append(line.split(".")[0].split("-")[2])
            for disp in disps:
                in_name = "supercell-%s.inp" % disp
                os.system("%s $PMF_CP2K -in phonon-supercell-%s.inp | tee phonon-supercell-%s.inp.out" % (self.run_params["mpi"], disp, disp))
            os.chdir("../")
        server_handle(auto=auto, directory=directory, jobfilebase="phonopy-job", server=self.run_params["server"])


    def to_subsys_phonopy(self, fname):
        cell = self.force_eval.subsys.xyz.cell
        with open(fname, 'a') as fout:
            fout.write("\t&SUBSYS\n")
            for element in self.force_eval.subsys.xyz.specie_labels:
                fout.write("\t\t&KIND %s\n" % element)
                fout.write("\t\t\tBASIS_SET DZVP-MOLOPT-SR-GTH\n")
                fout.write("\t\t\tPOTENTIAL GTH-PBE\n")
                fout.write("\t\t&END KIND\n")
            fout.write("\t\t&CELL\n")
            #fout.write("\t\t\tABC %f %f %f\n" % (cell[0], cell[4], cell[8]))
            # unit cell here can only be specified via ABC when doing phonopy calculation
            fout.write("\t\t\tABC %f %f %f\n" % (cell[0][0], cell[1][1], cell[2][2]))
            fout.write("\t\t&END CELL\n")
            #fout.write("\t\t&TOPOLOGY\n")
            #fout.write("\t\t\tCOORD_FILE_FORMAT xyz\n")
            #fout.write("\t\t\tCOORD_FILE_NAME %s\n" % sys.argv[1])
            #fout.write("\t\t&END TOPOLOGY\n")
            fout.write("\t\t&COORD\n")
            fout.write("\t\t\tSCALED .TRUE.\n")
            for atom in self.force_eval.subsys.xyz.atoms:
                #fout.write("\t\t\t%s\t%f\t%f\t%f\n" % (atom.name, atom.x/cell[0], atom.y/cell[4], atom.z/cell[8]))
                fout.write("\t\t\t%s\t%f\t%f\t%f\n" % (atom.name, atom.x/cell[0][0], atom.y/cell[1][1], atom.z/cell[2][2]))
            fout.write("\t\t&END COORD\n")
            fout.write("\t&END SUBSYS\n")
            fout.write("\n")

    def print_kinds(self, fname):
        with open(fname, 'a') as fout:
            for element in self.force_eval.subsys.xyz.specie_labels:
                fout.write("\t\t&KIND %s\n" % element)
                fout.write("\t\t\tBASIS_SET DZVP-MOLOPT-SR-GTH\n")
                fout.write("\t\t\tPOTENTIAL GTH-PBE\n")
                fout.write("\t\t&END KIND\n")
