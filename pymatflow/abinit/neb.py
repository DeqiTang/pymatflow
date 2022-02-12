"""
Nudged Elastic Band calculation
"""
import os
import shutil


from pymatflow.remote.server import server_handle
from pymatflow.abinit.abinit import Abinit
from pymatflow.abinit.base.system import AbinitSystem

class NebRun(Abinit):
    """
    """
    def __init__(self):
        """
        """
        super().__init__()
        self.images = []
        self.dataset[0].electrons.basic_setting()

        self.dataset[0].guard.set_queen(queen="neb")

        self.params = {
                "imgmov": 5,
                "nimage": 7,
                "ntimimage": 30,
                "mep_solver": None,
                "tolimg": 5.0e-05,
                "imgwfstor": None,
                "mep_mxstep": None,
                "neb_algo": None,
                "neb_string": None,
                "npimage": None,
                "string_algo": None,
                "cineb_start": 7,
                "dynimage": "0 5*1 0",
                "fxcartfactor": 1.0,
                }

    def get_images(self, images):
        """
        :param images:
            ["first.xyz", "last.xyz"]
        """
        self.images = []
        for image in images:
            system = abinit_system()
            system.xyz.get_xyz(image)
            self.images.append(system)

    def neb(self, directory="tmp-abinit-neb", runopt="gen", auto=0):

        self.files.name = "neb.files"
        self.files.main_in = "neb.in"
        self.files.main_out = "neb.out"
        self.files.wavefunc_in = "neb-i"
        self.files.wavefunc_out = "neb-o"
        self.files.tmp = "tmp"

        if runopt == "gen" or runopt == "genrun":
            if os.path.exists(directory):
                shutil.rmtree(directory)
            os.mkdir(directory)
            os.system("cp *.psp8 %s/" % directory)
            os.system("cp *.GGA_PBE-JTH.xml %s/" % directory)
            for image in self.images:
                os.system("cp %s %s/" % (image.xyz.file, directory))

            self.dataset[0].electrons.set_scf_nscf("scf")
            #

            # generate llhpc submit script
            script="neb.slurm"
            with open(os.path.join(directory, script),  'w') as fout:
                fout.write("#!/bin/bash\n")
                fout.write("#SBATCH -p %s\n" % self.run_params["partition"])
                fout.write("#SBATCH -N %d\n" % self.run_params["nodes"])
                fout.write("#SBATCH -n %d\n" % self.run_params["ntask"])
                fout.write("#SBATCH -J %s\n" % self.run_params["jobname"])
                fout.write("#SBATCH -o %s\n" % self.run_params["stdout"])
                fout.write("#SBATCH -e %s\n" % self.run_params["stderr"])
                fout.write("cat > %s<<EOF\n" % self.files.main_in)
                #self.dataset[0].electrons.to_dataset[0](fout)
                fout.write(self.dataset[0].electrons.to_string())
                self.images_to_input(fout)
                fout.write("\n")
                fout.write("# ===============================\n")
                fout.write("# neb related setting\n")
                fout.write("# ===============================\n")
                for item in self.params:
                    if self.params[item] is not None:
                        fout.write("%s %s\n" % (item, str(self.params[item])))
                fout.write("EOF\n")
                fout.write("cat > %s<<EOF\n" % self.files.name)
                #self.files.to_files(fout, system=self.dataset[0].system)
                #self.files.to_files(fout, system=self.images[0])
                fout.write(self.files.to_string(system=self.images[0]))
                fout.write("EOF\n")
                fout.write("yhrun %s < %s\n" % ("$PMF_ABINIT", self.files.name))



            # generate pbs submit script
            script="neb.pbs"
            with open(os.path.join(directory, script),  'w') as fout:
                fout.write("#!/bin/bash\n")
                fout.write("#PBS -N %s\n" % self.run_params["jobname"])
                fout.write("#PBS -l nodes=%d:ppn=%d\n" % (self.run_params["nodes"], self.run_params["ppn"]))
                if "queue" in self.run_params and self.run_params["queue"] != None:
                    fout.write("#PBS -q %s\n" % self.run_params["queue"])
                fout.write("\n")
                fout.write("cd $PBS_O_WORKDIR\n")
                fout.write("NP=`cat $PBS_NODEFILE | wc -l`\n")
                fout.write("cat > %s<<EOF\n" % self.files.main_in)
                #self.dataset[0].electrons.to_dataset[0](fout)
                fout.write(self.dataset[0].electrons.to_string())
                self.images_to_input(fout)
                fout.write("\n")
                fout.write("# ===============================\n")
                fout.write("# neb related setting\n")
                fout.write("# ===============================\n")
                for item in self.params:
                    if self.params[item] is not None:
                        fout.write("%s %s\n" % (item, str(self.params[item])))
                fout.write("EOF\n")
                fout.write("cat > %s<<EOF\n" % self.files.name)
                #self.files.to_files(fout, system=self.dataset[0].system)
                #self.files.to_files(fout, system=self.images[0])
                fout.write(self.files.to_string(system=self.images[0]))
                fout.write("EOF\n")
                fout.write("mpirun -np $NP -machinefile $PBS_NODEFILE %s < %s\n" % ("$PMF_ABINIT", self.files.name))

            # generate cdcloud submit script
            script="neb.slurm_cd"
            with open(os.path.join(directory, script),  'w') as fout:
                fout.write("#!/bin/bash\n")
                fout.write("#SBATCH -p %s\n" % self.run_params["partition"])
                fout.write("#SBATCH -N %d\n" % self.run_params["nodes"])
                fout.write("#SBATCH -n %d\n" % self.run_params["ntask"])
                fout.write("#SBATCH -J %s\n" % self.run_params["jobname"])
                fout.write("#SBATCH -o %s\n" % self.run_params["stdout"])
                fout.write("#SBATCH -e %s\n" % self.run_params["stderr"])
                fout.write("#\n")
                fout.write("export I_MPI_PMI_LIBRARY=/opt/gridview/slurm/lib/libpmi.so\n")
                fout.write("export FORT_BUFFERED=1\n")
                fout.write("cat > %s<<EOF\n" % self.files.main_in)
                #self.dataset[0].electrons.to_dataset[0](fout)
                fout.write(self.dataset[0].electrons.to_string())
                self.images_to_input(fout)
                fout.write("\n")
                fout.write("# ===============================\n")
                fout.write("# neb related setting\n")
                fout.write("# ===============================\n")
                for item in self.params:
                    if self.params[item] is not None:
                        fout.write("%s %s\n" % (item, str(self.params[item])))
                fout.write("EOF\n")
                fout.write("cat > %s<<EOF\n" % self.files.name)
                #self.files.to_files(fout, system=self.dataset[0].system)
                #self.files.to_files(fout, system=self.images[0])
                fout.write(self.files.to_string(system=self.images[0]))
                fout.write("EOF\n")
                fout.write("srun --mpi=pmix_v3 %s < %s\n" % ("$PMF_ABINIT", self.files.name))


        if runopt == "run" or runopt == "genrun":
            os.chdir(directory)
            os.system("$PMF_ABINIT < %s" % inpname.split(".")[0]+".files")
            os.chdir("../")
        server_handle(auto=auto, directory=directory, jobfilebase="neb", server=self.run_params["server"])

    def images_to_input(self, fout):
        #self.images[0].to_dataset[0](fout)
        fout.write(self.images[0].to_string())
        fout.write("\n")
        fout.write("xangst_lastimg\n")
        for atom in self.images[-1].xyz.atoms:
            fout.write("%.9f %.9f %.9f\n" % (atom.x, atom.y, atom.z))
        fout.write("\n")
