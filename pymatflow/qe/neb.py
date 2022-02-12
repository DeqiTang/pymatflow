"""
control Nudged Elastic Band(NEB)calculation
"""
import os
import re
import sys
import shutil
import pymatflow.base as base

from pymatflow.remote.server import server_handle
from pymatflow.qe.pwscf import PwScf
from pymatflow.qe.base.arts import QeArts

class NebRun(PwScf):
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
        surely we can through ATOMIC_POSITIONS
    """
    def __init__(self):
        """
        """
        super().__init__()
        self.images = []
        self.path = {} # Namelist: &PATH
        self.set_basic_path()

        self.control.basic_setting("scf")

    def get_images(self, images):
        """
        :param images:
            ["first.xyz", "intermediate-1.xyz", "intermediate-2.xyz", ...,"last.xyz"]
            self.images containe all the images while self.arts only contains the first image.
            self.arts is provided by base class pwscf, and is used to set kpoints.
        """
        self.images = []
        for image in images:
            arts = QeArts()
            arts.xyz.get_xyz(image)
            self.images.append(arts)
        self.system.basic_setting(self.images[0])
        for image in self.images:
            #image.basic_setting(ifstatic=True)
            image.basic_setting(ifstatic=False)
        # self.arts is actually the same as self.images[0]
        # it is used to set kpoints convinently
        self.arts.xyz.get_xyz(self.images[0].xyz.file)
        self.arts.basic_setting(ifstatic=True)

    def set_params(self, control={}, system={}, electrons={}):

        super().set_params(control=control, system=system, electrons=electrons)

        # must set "wf_collect = False", or it will come across with erros in davcio
        # Error in routine davcio (10):
        # error while reading from file ./tmp/pwscf_2/pwscf.wfc1
        self.control.set_param("wf_collect", "False")


    def set_path(self, path={}):
        for item in path:
            self.path[item] = path[item]

    def neb(self, directory="tmp-qe-neb", inpname="neb.in", output="neb.out", runopt="gen", restart_mode="from_scratch", auto=0):
        """
        :param directory: a place for all the generated files
        """
        if runopt == "gen" or runopt == "genrun":
            if restart_mode == "from_scratch":
                self.path["restart_mode"] = restart_mode
                if os.path.exists(directory):
                    shutil.rmtree(directory)
                os.mkdir(directory)

                #os.system("cp *.UPF %s/" % directory)
                #os.system("cp %s %s/" % (self.arts.xyz.file, directory))

                # do not copy too many files at the same time or it will be slow
                # so we do not copy all UPF files in the directory but just copy
                # those used in the calculation.
                for art in self.images:
                    shutil.copyfile(art.xyz.file, os.path.join(directory, os.path.basename(art.xyz.file)))
                #all_upfs = [s for s in os.listdir() if s.split(".")[-1] == "UPF"]
                all_file = os.listdir()
                for element in self.arts.xyz.specie_labels:
                    for item in all_file:
                        #if re.match("(%s)(.*)(upf)" % element, item, re.IGNORECASE):
                        if item.split(".")[0].lower() == element.lower() or item.split("_")[0].lower() == element.lower():
                            shutil.copyfile(item, os.path.join(directory, item))
                            break
                self.arts.pseudo.dir = os.path.abspath(directory)
                self.control.set_params({"pseudo_dir": os.path.abspath(directory)})
                #
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
                            if self.path[item] == ".true." or self.path[item] == ".false.":
                                fout.write("%s = %s\n" % (item, str(self.path[item])))
                            else:
                                fout.write("%s = '%s'\n" % (item, str(self.path[item])))
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
            self.gen_llhpc(directory=directory, inpname=inpname, output=output, cmd="$PMF_NEBX")
            # gen pbs script
            self.gen_pbs(directory=directory, inpname=inpname, output=output, cmd="$PMF_NEBX", jobname=self.run_params["jobname"], nodes=self.run_params["nodes"], ppn=self.run_params["ppn"], queue=self.run_params["queue"])
            # gen cdcloud script
            self.gen_cdcloud(directory=directory, inpname=inpname, output=output, cmd="$PMF_NEBX")

        if runopt == "run" or runopt == "genrun":
            os.chdir(directory)
            os.system("%s $PMF_NEBX -inp %s | tee %s" % (self.run_params["mpi"], inpname, output))
            os.chdir("../")
        server_handle(auto=auto, directory=directory, jobfilebase="neb", server=self.run_params["server"])

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
            fout.write("%s %f %s\n" % (element, base.element[element].mass, pseudo_file))
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
        if self.images[0].ifstatic == True:
            for atom in self.images[0].xyz.atoms:
                fout.write("%s\t%.9f\t%.9f\t%.9f\n" % (atom.name, atom.x, atom.y, atom.z))
        elif self.images[0].ifstatic == False:
            for atom in self.images[0].xyz.atoms:
                fout.write("%s\t%.9f\t%.9f\t%.9f" % (atom.name, atom.x, atom.y, atom.z))
                for fix in atom.fix:
                    if fix == True:
                        fout.write("\t0")
                    elif fix == False:
                        fout.write("\t1")
                fout.write("\n")
        else:
            print("===============================================\n")
            print("warning: qe.neb.neb_run.images_to_neb():\n")
            print("self.images[i].ifstatic could only be True or False\n")
            sys.exit(1)
        fout.write("\n")

        for i in range(1, len(self.images) - 1):
            fout.write("INTERMEDIATE_IMAGE\n")
            fout.write("ATOMIC_POSITIONS angstrom\n")
            if self.images[i].ifstatic == True:
                for atom in self.images[i].xyz.atoms:
                    fout.write("%s\t%.9f\t%.9f\t%.9f\n" % (atom.name, atom.x, atom.y, atom.z))
            elif self.images[i].ifstatiwc == False:
                for atom in self.images[i].xyz.atoms:
                    fout.write("%s\t%.9f\t%.9f\t%.9f" % (atom.name, atom.x, atom.y, atom.z))
                    for fix in atom.fix:
                        if fix == True:
                            fout.write("\t0")
                        elif fix == False:
                            fout.write("\t1")
                    fout.write("\n")
            else:
                print("===============================================\n")
                print("warning: qe.neb.neb_run.images_to_neb():\n")
                print("self.images[i].ifstatic could only be True or False\n")
                sys.exit(1)
            fout.write("\n")
                        
        fout.write("LAST_IMAGE\n")
        fout.write("ATOMIC_POSITIONS angstrom\n")
        if self.images[-1].ifstatic == True:
            for atom in self.images[-1].xyz.atoms:
                fout.write("%s\t%.9f\t%.9f\t%.9f\n" % (atom.name, atom.x, atom.y, atom.z))
        elif self.images[-1].ifstatic == False:
            for atom in self.images[-1].xyz.atoms:
                fout.write("%s\t%.9f\t%.9f\t%.9f" % (atom.name, atom.x, atom.y, atom.z))
                for fix in atom.fix:
                    if fix == True:
                        fout.write("\t0")
                    elif fix == False:
                        fout.write("\t1")
                fout.write("\n")
        else:
            print("===============================================\n")
            print("warning: qe.neb.neb_run.images_to_neb():\n")
            print("self.images[i].ifstatic could only be True or False\n")
            sys.exit(1)
        fout.write("\n")
        fout.write("END_POSITIONS\n")


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

    def gen_yh(self, inpname, output, directory, cmd="neb.x"):
        """
        generating yhbatch job script for calculation
        """
        with open(os.path.join(directory, inpname.split(".in")[0]+".sub"), 'w') as fout:
            fout.write("#!/bin/bash\n")
            fout.write("yhrun -N 1 -n 24 %s -inp %s > %s\n" % (cmd, inpname, output))

    def gen_pbs(self, inpname, output, directory, cmd="neb.x", jobname="neb.x", nodes=1, ppn=32, queue=None):
        """
        generating pbs job script for calculation
        """
        with open(os.path.join(directory, inpname.split(".in")[0]+".pbs"), 'w') as fout:
            fout.write("#!/bin/bash\n")
            fout.write("#PBS -N %s\n" % jobname)
            fout.write("#PBS -l nodes=%d:ppn=%d\n" % (nodes, ppn))
            if queue != None:
                fout.write("#PBS -q %s\n" % queue)
            fout.write("\n")
            fout.write("cd $PBS_O_WORKDIR\n")
            fout.write("NP=`cat $PBS_NODEFILE | wc -l`\n")
            fout.write("mpirun -np $NP -machinefile $PBS_NODEFILE %s -inp %s > %s\n" % (cmd, inpname, output))
