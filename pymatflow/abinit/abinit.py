"""
Overall representation of Abinit
"""

import os
import shutil

from pymatflow.abinit.base.input import abinit_input
from pymatflow.abinit.base.files import abinit_files



class abinit:
    """
    Note: support for both single dataset mode and multiple data set mode now
        if you want to use only the single dataset mode use ndtset = 0
        otherwise use ndtset to the number of dataset yout want to use
    """
    def __init__(self, ndtset=0):
        """
        :parma ndtset: -> self.ndtset. it is uesed to decide how much dataset
            will be output to the input string. like even when the length of
            self.dataset if 7, it will only output the number of datasets
            defined by self.ndtset. default behavior is output all the datasets
            Note: we should not directly set value of ndtset, but should do it
            via self.set_ndtset()
        """
        self.ndtset = ndtset

        self.dataset = []
        for i in range(self.ndtset+1):
            self.dataset.append(abinit_input())
            self.dataset[-1].n = i

        self.dataset[0].electrons.basic_setting()
        self.files = abinit_files()
        self._initialize()

    def _initialize(self):
        """ initialize the object, do some default setting
        """
        self.run_params  = {}
        self.set_run()

    def get_xyz(self, xyzfile):
        # only get structure for the default dataset 0
        self.dataset[0].system.xyz.get_xyz(xyzfile)

    def set_params(self, params={}):
        # only set the parameters for default dataset 0
        self.dataset[0].set_params(params)

    def set_kpoints(self, kpoints={}):
        # only set the kpoints for default dataset 0
        self.dataset[0].electrons.kpoints.set_params(kpoints)

    def set_properties(self, properties=[]):
        # only set the parameters for default dataset 0
        self.dataset[0].properties.get_option(option=properties)
    #

    def dft_plus_u(self):
        # only set the parameters for default dataset 0
        self.dataset[0].electrons.dft_plus_u()

    def set_run(self, mpi="", server="pbs", jobname="abinit", nodes=1, ppn=32):
        self.run_params["mpi"] = mpi
        self.run_params["server"] = server
        self.run_params["jobname"] = jobname
        self.run_params["nodes"] = nodes
        self.run_params["ppn"] = ppn

    def to_string(self):
        """
        :return the input as a string
        """

        inp_str = "ndtset %d\n" % self.ndtset
        inp_str += "# -----------------Overall Default Dataset: 0 ------------------\n"
        inp_str += self.dataset[0].to_string()

        for i in range(1, self.ndtset+1):
            # do not return the system related setting for dataset other than 0
            self.dataset[i].system.status = False
            head = "# -----------------------Dataset: %d ---------------------------\n" % (i)
            inp_str += head + self.dataset[i].to_string()
        return inp_str


    def set_ndtset(self, ndtset=0):
        """
        :param ndtset: setting how many dadaset will be used, default is 0 which
            means ndtset equals to 0, and multidataset mode is not used.
            if there are not enough dataset if self.dataset, it will automatically
            add them.
            if there are more dataset in self.dataset than ndtset, they will
            not be removed, but self.ndtset will be set to ndtset, and that will
            affect how many dataset will be output to the input string.
        """
        self.ndtset = ndtset
        while len(self.dataset) - 1 < ndtset:
            self.dataset.append(abinit_input())
            self.dataset[-1].n = len(self.dataset) - 1

    def __getitem__(self, item):
        """
        Note: if item is within self.dataset, it can be returned.
            if item is equal to the length of self.dataset, it is over
            the self.dataset list by 1 element, this way, another abinit_input()
            is append to self.dataset.
            but if item is larger than length of self.dataset, it will
            print out warning, and exit the program.
        """
        if item < len(selfl.dataset):
            self.ndtset = len(self.dataset) - 1
            return self.dataset[item]
        elif item == len(self.dataset):
            self.dataset.append(abinit_input())
            self.dataset[-1] = len(self.dataset) - 1
            self.ndtset = len(self.dataset) - 1
            return self.dataset[item]
        else:
            print("=============================================================\n")
            print("                      Warning !!!\n")
            print("-------------------------------------------------------------\n")
            print("pymatflow.abinit.abinit.__getitem__:\n")
            print("trying to get the unexisted element from self.dataset\n")
            sys.exit()

    def gen_yh(self, inpname, output, directory, cmd="abinit"):
        """
        generating yhbatch job script for calculation
        """
        with open(os.path.join(directory, inpname.split(".in")[0]+".sub"), 'w') as fout:
            fout.write("#!/bin/bash\n")
            fout.write("yhrun -N 1 -n 24 %s < %s > %s\n" % (cmd, inpname, output))

    def gen_pbs(self, directory, script="abinit.pbs", cmd="abinit", jobname="abinit", nodes=1, ppn=32):
        """
        generating pbs job script for calculation
        """
        with open(os.path.join(directory, script),  'w') as fout:
            fout.write("#!/bin/bash\n")
            fout.write("#PBS -N %s\n" % jobname)
            fout.write("#PBS -l nodes=%d:ppn=%d\n" % (nodes, ppn))
            fout.write("\n")
            fout.write("cd $PBS_O_WORKDIR\n")
            fout.write("NP=`cat $PBS_NODEFILE | wc -l`\n")
            fout.write("cat > %s<<EOF\n" % self.files.main_in)
            #self.input.to_input(fout)
            fout.write(self.to_string())
            fout.write("EOF\n")
            fout.write("cat > %s<<EOF\n" % self.files.name)
            #self.files.to_files(fout, system=self.input.system)
            fout.write(self.files.to_string(system=self.dataset[0].system))
            fout.write("EOF\n")
            fout.write("mpirun -np $NP -machinefile $PBS_NODEFILE %s < %s\n" % (cmd, self.files.name))

    def gen_bash(self, directory, script="abinit.sh", cmd="abinit", mpi=""):
        """
        generating pbs job script for calculation
        """
        with open(os.path.join(directory, script),  'w') as fout:
            fout.write("#!/bin/bash\n")
            fout.write("cat > %s<<EOF\n" % self.files.main_in)
            #self.input.to_input(fout)
            fout.write(self.to_string())
            fout.write("EOF\n")
            fout.write("cat > %s<<EOF\n" % self.files.name)
            #self.files.to_files(fout, system=self.input.system)
            fout.write(self.files.to_string(system=self.dataset[0].system))
            fout.write("EOF\n")
            fout.write("%s %s < %s\n" % (mpi, cmd, self.files.name))
