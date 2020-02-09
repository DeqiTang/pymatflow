
import os
import sys
import shutil
import matplotlib.pyplot as plt

from pymatflow.qe.base.control import qe_control
from pymatflow.qe.base.system import qe_system
from pymatflow.qe.base.electrons import qe_electrons
from pymatflow.qe.base.ions import qe_ions
from pymatflow.qe.base.cell import qe_cell
from pymatflow.qe.base.arts import qe_arts


class pwscf:
    """
    About:
        this class provide a base representation of pwscf
    """
    def __init__(self):
        self.control = qe_control()
        self.system = qe_system()
        self.electrons = qe_electrons()
        self.ions = qe_ions()
        self.cell = qe_cell()
        self.arts = qe_arts()

        #self.control.basic_setting("scf")
        self.electrons.basic_setting()
        self.set_kpoints() # default kpoint setting

    def get_xyz(self, xyzfile):
        """
        xyz_f:
            a modified xyz formatted file(the second line specifies the cell of the
            system).
        """
        self.arts.xyz.get_xyz(xyzfile)
        self.system.basic_setting(self.arts)
        self.arts.basic_setting(ifstatic=True)


    def set_params(self, control={}, system={}, electrons={}, ions={}, cell={}):
        # check if user try to set occupations and smearing and degauss
        # through system. if so, use self.set_occupations() which uses
        # self.system.set_occupations() to set them, as self.system.set_params()
        # is suppressed from setting occupations related parameters
        self.set_occupations(system)
        self.control.set_params(control)
        self.system.set_params(system)
        self.electrons.set_params(electrons)
        self.ions.set_params(ions)
        self.cell.set_params(cell)

    def set_kpoints(self, kpoints_option="automatic", kpoints_mp=[2, 2, 2, 0, 0, 0], crystal_b=None):
        # about the format of crystal_b
        # see corresponding comment for function self.arts.set_kpoints()
        self.arts.set_kpoints(option=kpoints_option, kpoints_mp=kpoints_mp, crystal_b=crystal_b)

    def run(self, directory="tmp-pwscf", inpname="pwscf.in", output="pwscf.out", mpi="", runopt="gen",
            jobname="pwscf", nodes=1, ppn=32):
        """
        directory: the place where all the magic happening

        parameters:
            directory: the overall static calculation directory

        runopt: determine whether the calculation is executed.
                there are three values: 'gen', 'genrun', 'run'
                'gen': only generate the input files
                'genrun': generate input files and run
                'run': run from the previously generated input files
                P.S.  run is run by local command directly rather than
                      commit the job through job manager
        Note:
            two mode of generating the input files: (1) a brand new calculation
            remove the directory if it already exists and create a brand new
            directory for calculation. (2) a new calculation but if there exists
            the directory already, will not remove it but just generate the input
            files inside.
        """
        if runopt == 'gen' or runopt == 'genrun':
            if os.path.exists(directory):
                shutil.rmtree(directory)
            os.mkdir(directory)

            #os.system("cp *.UPF %s/" % directory)
            #os.system("cp %s %s/" % (self.arts.xyz.file, directory))

            # do not copy too many files at the same time or it will be slow
            # so we do not copy all UPF files in the directory but just copy
            # those used in the calculation.
            shutil.copyfile(self.arts.xyz.file, os.path.join(directory, self.arts.xyz.file))
            all_upfs = [s for s in os.listdir() if s.split(".")[-1] == "UPF"]
            for element in self.arts.xyz.specie_labels:
                for upf in all_upfs:
                    if upf.split(".")[0] == element:
                        shutil.copyfile(upf, os.path.join(directory, upf))
                        break
            #

            with open(os.path.join(directory, inpname), 'w') as fout:
                self.control.to_in(fout)
                self.system.to_in(fout)
                self.electrons.to_in(fout)
                self.arts.to_in(fout)

            # gen yhbatch script
            self.gen_yh(directory=directory, inpname=inpname, output=output, cmd="pw.x")
            # gen pbs scripts
            self.gen_pbs(directory=directory, inpname=inpname, output=output, cmd="pw.x", jobname=jobname, nodes=nodes, ppn=ppn)

        if runopt == 'genrun' or runopt == 'run':
            os.chdir(directory)
            os.system("%s pw.x < %s | tee %s" % (mpi, inpname, output))
            os.chdir("../")

    def set_atomic_forces(self, pressure=None, pressuredir=None):
        self.arts.set_atomic_forces(pressure=pressure, direction=pressuredir)

    def set_occupations(self, system):
        """
            # check if user try to set occupations and smearing and degauss
            # through system. if so, use self.system.set_occupations() to
            # set them, as self.system.set_params() is suppressed from setting
            # occupations related parameters
            # if occupations == None, use default smearing occupation. and
            # if occupations == "tetrahedra" the value set for smearing and degauss is ignored.
            # if occupations == "smearing", the value of smearing and degauss
            # should be legacy, not None or other illegal values.
        """
        if "occupations" in system:
            if system["occupations"] == None: # user default setting of set_occupations()
                self.system.set_occupations()
            elif system["occupations"] == "tetrahedra":
                self.system.set_occupations(occupations="tetrahedra")
            elif system["occupations"] == "smearing":
                if "smearing" in system and "degauss" in system:
                    self.system.set_occupations(occupations="smearing", smearing=system["smearing"], degauss=system["degauss"])
                elif "smearing" in system:
                    self.system.set_occupations(occupations="smearing", smearing=system["smearing"])
                elif "degauss" in system:
                    self.system.set_occupations(occupations="smearing", degauss=system["degauss"])
                else:
                    self.system.set_occupations(occupations="smearing")
            elif system["occupations"] == "tetrahedra_lin":
                self.system.set_occupations(occupations="tetrahedra_lin")
            elif system["occupations"] == "tetrahedra_opt":
                self.system.set_occupations(occupations="tetrahedra_opt")
            elif system["occupations"] == "fixed":
                self.system.set_occupations(occupations="fixed")
            elif system["occupations"] == "from_input":
                self.system.set_occupations(occupations="from_input")


    def set_spin(self):
        pass

    def gen_yh(self, inpname, output, directory, cmd="pw.x"):
        """
        generating yhbatch job script for calculation
        """
        with open(os.path.join(directory, inpname.split(".in")[0]+".sub"), 'w') as fout:
            fout.write("#!/bin/bash\n")
            fout.write("yhrun -N 1 -n 24 %s < %s > %s\n" % (cmd, inpname, output))

    def gen_pbs(self, inpname, output, directory, cmd="pw.x", jobname="pwscf", nodes=1, ppn=32):
        """
        generating pbs job script for calculation
        """
        with open(os.path.join(directory, inpname.split(".in")[0]+".pbs"), 'w') as fout:
            fout.write("#!/bin/bash\n")
            fout.write("#PBS -N %s\n" % jobname)
            fout.write("#PBS -l nodes=%d:ppn=%d\n" % (nodes, ppn))
            fout.write("\n")
            fout.write("cd $PBS_O_WORKDIR\n")
            fout.write("NP=`cat $PBS_NODEFILE | wc -l`\n")
            fout.write("mpirun -np $NP -machinefile $PBS_NODEFILE %s < %s > %s\n" % (cmd, inpname, output))
