"""
Static calculation
"""
import numpy as np
import sys
import os
import shutil
import matplotlib.pyplot as plt

from pymatflow.remote.server import server_handle
from pymatflow.siesta.siesta import siesta
#from pymatflow.siesta.base.system import siesta_system
#from pymatflow.siesta.base.electrons import siesta_electrons
#from pymatflow.siesta.base.properties import siesta_properties

class static_run(siesta):
    """
    """
    def __init__(self):
        super().__init__()
        #self.system = siesta_system()
        #self.electrons = siesta_electrons()
        #self.properties = siesta_properties()

        self.electrons.basic_setting()


    def scf(self, directory="tmp-siesta-static", inpname="static-scf.fdf", output="static-scf.out", runopt="gen", auto=0, properties=[]):
        if runopt == "gen" or runopt == "genrun":
            if os.path.exists(directory):
                shutil.rmtree(directory)
            os.mkdir(directory)

            for element in self.system.xyz.specie_labels:
                shutil.copyfile("%s.psf" % element, os.path.join(directory, "%s.psf" % element))

            # use self.properties.options to contorl the calculation of properties
            self.properties.options = properties
            with open(os.path.join(directory, inpname), 'w') as fout:
                self.system.to_fdf(fout)
                self.electrons.to_fdf(fout)
                self.properties.to_fdf(fout)

            # gen yhbatch script
            self.gen_yh(directory=directory, inpname=inpname, output=output, cmd="siesta")
            # gen pbs script
            self.gen_pbs(directory=directory, inpname=inpname, output=output, cmd="siesta", jobname=self.run_params["jobname"], nodes=self.run_params["nodes"], ppn=self.run_params["ppn"])

        if runopt == "run" or runopt == "genrun":
            os.chdir(directory)
            os.system("%s siesta < %s | tee %s" % (self.run_params["mpi"], inpname, output))
            os.chdir("../")
        server_handle(auto=auto, directory=directory, jobfilebase="static-scf", server=self.params["server"])

    def scf_restart(self, directory="tmp-siesta-static", inpname="static-scf-restart.fdf", output="static-scf-restart.out",
        runopt="gen", auto=0, properties=[]):

        # first check whether there is a previous scf running
        if not os.path.exists(directory):
            print("===================================================\n")
            print("                 Warning !!!\n")
            print("===================================================\n")
            print("scf(restart) calculation:\n")
            print("  directory of previous scf calculattion not found!\n")
            sys.exit(1)
        if runopt == "gen" or runopt == "genrun":

            self.electrons.dm["UseSaveDM"] = "true"
            # use self.properties.option to contorl the calculation of properties
            self.properties.options = properties
            with open(os.path.join(directory, inpname), 'w') as fout:
                self.system.to_fdf(fout)
                self.electrons.to_fdf(fout)
                self.properties.to_fdf(fout)

            # gen yhbatch script
            self.gen_yh(directory=directory, inpname=inpname, output=output, cmd="siesta")
            # gen pbs script
            self.gen_pbs(directory=directory, inpname=inpname, output=output, cmd="siesta", jobname=self.run_params["jobname"], nodes=self.run_params["nodes"], ppn=self.run_params["ppn"])

        if runopt == "run" or runopt == "genrun":
            os.chdir(directory)
            os.system("%s siesta < %s | tee %s" % (self.run_params["mpi"], inpname, output))
            os.chdir("../")
        server_handle(auto=auto, directory=directory, jobfilebase="static-scf-restart", server=self.params["server"])


    def converge_cutoff(self, emin, emax, step, directory="tmp-siesta-cutoff", runopt="gen", auto=0):
        if runopt == "gen" or runopt == "genrun":
            if os.path.exists(directory):
                shutil.rmtree(directory)
            os.mkdir(directory)

            for element in self.system.xyz.specie_labels:
                shutil.copyfile("%s.psf" % element, os.path.join(directory, "%s.psf" % element))

            self.electrons.dm["UseSaveDM"] = "false"

            n_test = int((emax - emin) / step)
            for i in range(n_test + 1):
                meshcutoff = int(emin + i * step)
                self.electrons.params["MeshCutoff"] = meshcutoff
                self.system.label = "siesta-" + str(meshcutoff)
                with open(os.path.join(directory, "cutoff-%d.fdf" % meshcutoff), 'w') as fout:
                    self.system.to_fdf(fout)
                    self.electrons.to_fdf(fout)
                    #self.properties.to_fdf(fout)

            # gen yhbatch running script
            with open("converge-cutoff.sub", 'w') as fout:
                fout.write("#!/bin/bash\n")
                for i in range(n_test + 1):
                    meshcutoff = int(emin + i * step)
                    inp_name = "cutoff-%d.fdf" % meshcutoff
                    out_f_name = "cutoff-%d.out" % meshcutoff
                    fout.write("yhrun -N 1 -n 24 siesta < %s > %s\n" % (inp_name, out_f_name))

            # gen pbs running script
            with open("converge-cutoff.pbs", 'w') as fout:
                fout.write("#!/bin/bash\n")
                fout.write("#PBS -N %s\n" % self.run_params["jobname"])
                fout.write("#PBS -l nodes=%d:ppn=%d\n" % (self.run_params["nodes"], self.run_params["ppn"]))
                fout.write("\n")
                fout.write("cd $PBS_O_WORKDIR\n")
                fout.write("NP=`cat $PBS_NODEFILE | wc -l`\n")
                for i in range(n_test + 1):
                    meshcutoff = int(emin + i * step)
                    inp_name = "cutoff-%d.fdf" % meshcutoff
                    out_f_name = "cutoff-%d.out" % meshcutoff
                    fout.write("mpirun -np $NP -machinefile $PBS_NODEFILE siesta < %s > %s\n" % (inp_name, out_f_name))

        if runopt == "run" or runopt == "genrun":
            # run
            os.chdir(directory)
            for i in range(n_test + 1):
                meshcutoff = int(emin + i * step)
                os.system("%s siesta < cutoff-%d.fdf | tee cutoff-%d.out" % (self.run_params["mpi"], meshcutoff, meshcutoff))
            os.chdir("../")
        server_handle(auto=auto, directory=directory, jobfilebase="converge-cutoff", server=self.params["server"])

    #
