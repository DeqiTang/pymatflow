#!/usr/bin/env python
# _*_ coding: utf-8 _*_

import os
import shutil
import matplotlib.pyplot as plt

from pymatflow.abinit.abinit import abinit


class static_run(abinit):
    """
    GOAL: support for both single dataset and multi-dataset mode in abinit,
          currently, only for single dataset mode
    """
    def __init__(self):
        super().__init__()
        #self.system = abinit_system()
        #self.electrons = abinit_electrons()
        #self.properties = abinit_properties()
        self.input.guard.set_queen(queen="static", electrons=self.input.electrons, system=self.input.system)

        self.input.electrons.basic_setting()


    def scf(self, directory="tmp-abinit-static", inpname="static-scf.in", mpi="", runopt="gen",
        jobname="abinit-scf", nodes=1, ppn=32):
        self.files.name = "static-scf.files"
        self.files.main_in = "static-scf.in"
        self.files.main_out = "static-scf.out"
        self.files.wavefunc_in = "static-scf-i"
        self.files.wavefunc_out = "static-scf-o"
        self.files.tmp = "tmp"
        if runopt == "gen" or runopt == "genrun":
            if os.path.exists(directory):
                shutil.rmtree(directory)
            os.mkdir(directory)
            os.system("cp *.psp8 %s/" % directory)
            os.system("cp *.GGA_PBE-JTH.xml %s/" % directory)
            os.system("cp %s %s/" % (self.input.system.xyz.file, directory))

            self.input.electrons.set_scf_nscf("scf")
            #
            self.input.guard.check_all()



            # generate pbs job submit script
            self.gen_pbs(directory=directory, script="static-scf.pbs", cmd="abinit", jobname=jobname, nodes=nodes, ppn=ppn)

            # generate local bash job run script
            self.gen_bash(directory=directory, script="static-scf.sh", cmd="abinit", mpi=mpi)

        if runopt == "run" or runopt == "genrun":
            os.chdir(directory)
            #os.system("abinit < %s" % inpname.split(".")[0]+".files")
            os.system("bash %s" % "static-scf.sh")
            os.chdir("../")

    def nscf(self, directory="tmp-abinit-static", mpi="", runopt="gen",
        jobname="staic-nscf", nodes=1, ppn=32):
        # first check whether there is a previous scf running

        self.files.name = "static-nscf.files"
        self.files.main_in = "static-nscf.in"
        self.files.main_out = "static-nscf.out"
        self.files.wavefunc_in = "static-scf-o"
        self.files.wavefunc_out = "static-nscf-o"
        self.files.tmp = "tmp"
        if not os.path.exists(directory):
            print("===================================================\n")
            print("                 Warning !!!\n")
            print("===================================================\n")
            print("non-scf calculation:\n")
            print("  directory of previous scf calculattion not found!\n")
            sys.exit(1)
        if runopt == "gen" or runopt == "genrun":

            self.input.electrons.set_scf_nscf("nscf")
            self.input.electrons.params["irdwfk"] = 1
            self.input.electrons.params["irdden"] = 1
            # dos
            self.input.electrons.params["nband"] = 10
            self.input.electrons.params["dosdeltae"] = 0.00005
            self.input.electrons.params["occopt"] = 7
            self.input.electrons.params["tsmear"] = 0.0001
            # end dos

            #
            self.input.guard.check_all()



            # generate pbs job submit script
            self.gen_pbs(directory=directory, script="static-nscf.pbs", cmd="abinit", jobname=jobname, nodes=nodes, ppn=ppn)

            # generate local bash job run script
            self.gen_bash(directory=directory, script="static-nscf.sh", cmd="abinit", mpi=mpi)

        if runopt == "run" or runopt == "genrun":
            os.chdir(directory)
            #os.system("abinit < %s" % inpname.split(".")[0]+".files")
            os.system("bash %s" % "static-nscf.sh")
            os.chdir("../")


    def bands(self, directory="tmp-abinit-static", inpname="static-band.in", mpi="", runopt="gen",
        jobname="abinit-band", nodes=1, ppn=32):
        """
        we can use abiopen.py static-band-output_GSR.nc --expose -sns=talk to view the band structure.
        """
        self.files.name = "static-bands.files"
        self.files.main_in = "static-bands.in"
        self.files.main_out = "static-bands.out"
        self.files.wavefunc_in = "static-nscf-o"
        self.files.wavefunc_out = "static-bands-o"
        self.files.tmp = "tmp"
        # first check whether there is a previous scf running
        if not os.path.exists(directory):
            print("===================================================\n")
            print("                 Warning !!!\n")
            print("===================================================\n")
            print("band structure calculation:\n")
            print("  directory of previous scf calculattion not found!\n")
            sys.exit(1)
        if runopt == "gen" or runopt == "genrun":

            self.input.electrons.params["iscf"] = -2
            self.input.electrons.params["nband"] = 8
            self.input.electrons.params["tolwfr"] = 1.0e-12 # when kptopt < 0 namely band structure calculatin, we can only use tolwfr
            self.input.electrons.params["tolvrs"] = None
            self.input.electrons.params["toldfe"] = None
            #self.input.electrons.params["irdden"] = 1 # actually irdden will be 1 by default if iscf < 0

            #with open(os.path.join(directory, inpname), 'w') as fout:
                #self.electrons.to_in(fout)
                #self.system.to_in(fout)

            self.input.guard.check_all()

            # generate pbs job submit script
            self.gen_pbs(directory=directory, script="static-bands.pbs", cmd="abinit", jobname=jobname, nodes=nodes, ppn=ppn)

            # generate local bash job run script
            self.gen_bash(directory=directory, script="static-bands.sh", cmd="abinit", mpi=mpi)

        if runopt == "run" or runopt == "genrun":
            os.chdir(directory)
            #os.system("abinit < %s" % inpname.split(".")[0]+".files")
            os.system("bash %s" % "static-bands.sh")
            os.chdir("../")

    def converge_ecut(self, emin, emax, step, directory="tmp-abinit-ecut", mpi="", runopt="gen",
        jobname="converge-ecut", nodes=1, ppn=32):

        if runopt == "gen" or runopt == "genrun":
            if os.path.exists(directory):
                shutil.rmtree(directory)
            os.mkdir(directory)
            os.system("cp *.psp8 %s/" % directory)
            os.system("cp %s %s/" % (self.input.system.xyz.file, directory))

            os.chdir(directory)
            n_test = int((emax - emin) / step)
            for i in range(n_test + 1):
                cutoff = int(emin + i * step)
                inp_name = "ecut-%d.in" % cutoff
                files_name = "ecut-%d.files" % cutoff
                with open(files_name, 'w') as fout:
                    fout.write(inp_name)
                    fout.write("\n")
                    fout.write("ecut-%d.out\n" % cutoff)
                    fout.write("ecut-%d-input\n" % cutoff)
                    fout.write("ecut-%d-output\n" % cutoff)
                    fout.write("temp\n")
                    for element in self.input.system.xyz.specie_labels:
                        fout.write("%s\n" % (element + ".psp8"))
                #
                self.input.electrons.params["ecut"] = cutoff
                with open(inp_name, 'w') as fout:
                    self.input.electrons.to_input(fout)
                    self.input.system.to_input(fout)
            os.chdir("../")

            # generate pbs script files
            with open(os.path.join(directory, "converge-ecut.pbs"), 'w') as fout:
                fout.write("#!/bin/bash\n")
                fout.write("#PBS -N %s\n" % jobname)
                fout.write("#PBS -l nodes=%d:ppn=%d\n" % (nodes, ppn))
                fout.write("\n")
                fout.write("cd $PBS_O_WORKDIR\n")
                fout.write("NP=`cat $PBS_NODEFILE | wc -l`\n")
                for i in range(n_test + 1):
                    cutoff = int(emin + i * step)
                    #inp_name = "ecut-%d.in" % cutoff
                    files_name = "ecut-%d.files" % cutoff
                    fout.write("mpirun -np $NP -machinefile $PBS_NODEFILE %s < %s\n" % ("abinit", files_name))
            # generate the result analsysis scripts
            os.system("mkdir -p %s/post-processing" % directory)
            with open(os.path.join(directory, "post-processing/analysis-ecut.sh"), 'w') as fout:
                fout.write("#!/bin/bash\n")
                fout.write("\n")
                for i in range(n_test + 1):
                    cutoff = int(emin + i * step)
                    out_f_name = "ecut-%d.out" % cutoff
                    #fout.write("cat ../%s | grep 'Etotal=' >> energy-ecut.data" % out_f_name)
                    fout.write("energy=`cat ../%s | grep \'Etotal=\' | cut -d \"=\" -f 2 `\n" % out_f_name)
                    fout.write("cat >> energy-ecut.data<<EOF\n")
                    fout.write("%d ${energy}\n" % cutoff)
                    fout.write("EOF\n")
                fout.write("\n")
                fout.write("cat >> ecut-energy.gp<<EOF\n")
                fout.write("set term gif\n")
                fout.write("set output 'energy-ecut.gif'\n")
                fout.write("set title 'Ecut Converge Test'\n")
                fout.write("set xlabel 'Ecut()'\n")
                fout.write("set ylabel 'Total Energy()'\n")
                fout.write("plot 'energy-ecut.data' w l\n")
                fout.write("EOF\n")
                fout.write("\n")
                fout.write("gnuplot ecut-energy.gp")

            #
        if runopt == "run" or runopt == "genrun":
            # run the simulation
            os.chdir(directory)
            for i in range(n_test + 1):
                cutoff = int(emin + i * step)
                files_name = "ecut-%d.files" % cutoff
                os.system("abinit < %s" % (files_name))
            os.chdir("../")
