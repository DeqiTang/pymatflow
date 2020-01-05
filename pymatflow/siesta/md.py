#!/usr/bin/env python
# _*_ coding: utf-8 _*_

import numpy as np
import sys
import os
import shutil
import pymatgen as mg


from pymatflow.siesta.siesta import siesta
#from pymatflow.siesta.base.system import siesta_system
#from pymatflow.siesta.base.electrons import siesta_electrons
#from pymatflow.siesta.base.ions import siesta_ions


class md_run(siesta):
    """
    """
    def __init__(self):
        super().__init__()
        #self.system = siesta_system()
        #self.electrons = siesta_electrons()
        #self.ions = siesta_ions()

        self.electrons.basic_setting()
        self.ions.basic_setting(option="md")

    def md(self, directory="tmp-siesta-md", inpname="molecular-dynamics.fdf", output="molecular-dynamics.out",
            mpi="", runopt="gen"):
        """
        """
        if runopt == "gen" or runopt == "genrun":
            if os.path.exists(directory):
                shutil.rmtree(directory)
            os.mkdir(directory)
        
            for element in self.system.xyz.specie_labels:
                shutil.copyfile("%s.psf" % element, os.path.join(directory, "%s.psf" % element))

            with open(os.path.join(directory, inpname), 'w') as fout:
                self.system.to_fdf(fout)
                self.electrons.to_fdf(fout)
                self.ions.to_fdf(fout)

            # gen yhbatch script
            self.gen_yh(directory=directory, inpname=inpname, output=output, cmd="siesta")

        if runopt == "run" or runopt == "genrun":
            # run the simulation
            os.chdir(directory)
            os.system("siesta < %s | tee %s" % (inpname, output))
            os.chdir("../")

    #
