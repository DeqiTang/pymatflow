#!/usr/bin/env python
# _*_ coding: utf-8 _*_

import os
import sys
import numpy as np
import shutil
import pymatgen as mg
import matplotlib.pyplot as plt

from pymatflow.siesta.base.system import siesta_system
from pymatflow.siesta.base.electrons import siesta_electrons
from pymatflow.siesta.base.properties import siesta_properties

from pymatflow.siesta.base.transiesta import siesta_transiesta
from pymatflow.siesta.base.tbtrans import siesta_tbtrans

class ts_run:
    """
    """
    def __init__(self, electrodes, scattering):
        self.system_electrodes = []
        self.properties = []
        for xyz_f in electrodes:
            self.system_electrodes.append(siesta_system(xyz_f))

        for system in self.system_electrodes:
            self.properties.append(siesta_properties(system.xyz))
        
        self.system_scattering = siesta_system(scattering)

        self.electrons = siesta_electrons()

        self.transiesta = siesta_transiesta()

        self.tbtrans = siesta_tbtrans()
        
        self.electrons.basic_setting()
        self.electrons.params["SolutionMethod"] = "transiesta"
                

    def ts(self, directory="tmp-siesta-ts", inpname="transiesta.fdf", output="transiesta.out",
            mpi="", runopt="gen", electrons={}, properties=[], kpoints_mp=[1, 1, 1]):
        if runopt == "gen" or runopt == "genrun":
            if os.path.exists(directory):
                shutil.rmtree(directory)
            os.mkdir(directory)
        
            for element in self.system_electrodes[0].xyz.specie_labels:
                shutil.copyfile("%s.psf" % element, os.path.join(directory, "%s.psf" % element))
       
            self.electrons.kpoints_mp = kpoints_mp
            self.electrons.set_params(electrons)
            # use self.properties.options to contorl the calculation of properties
            #self.properties.options = properties

            # calculate the electrondes first
            os.mkdir(os.path.join(directory, "electrodes"))
            for i in range(len(self.system_electrodes)):
                os.mkdir(os.path.join(directory, "electrodes", "electrode-%d" % i))
                with open(os.path.join(directory, "electrodes", "electrode-%d" % i, "electrode.fdf"), 'w') as fout:
                    self.system_electrodes[i].to_fdf(fout)
                    self.electrons.to_fdf(fout)
                    #self.properties.to_fdf(fout)
                for element in self.system_electrodes[0].xyz.specie_labels:
                    shutil.copyfile("%s.psf" % element, os.path.join(directory, "electrodes", "electrode-%d" % i, "%s.psf" % element))

            #self.transiesta.to_fdf(fout)
            #self.tbtrans.to_fdf(fout)


            # gen yhbatch script
            #self.gen_yh(directory=directory, inpname=inpname, output=output, cmd="siesta")

        if runopt == "run" or runopt == "genrun":
            os.chdir(directory)
            # run the electrode calculation
            for i in range(len(self.system_electrodes)):
                os.chdir(os.path.join("electrodes", "electrode-%d" % i))
                os.system("%s siesta --electrode < %s | tee %s" % (mpi, "electrode.fdf", "electrode.fdf.out"))
                os.chdir("../../")
            #
            # run the scattering zone
            #os.system("%s transiesta < %s | tee %s" % (mpi, inpname, output))
            #os.system('%s tbtrans < %s | tee %s' % (mpi, inpname, output))
            os.chdir("../")


    def gen_yh(self, inpname, output, directory="tmp-siesta-static", cmd="siesta"):
        """
        generating yhbatch job script for calculation
        """
        with open(os.path.join(directory, inpname+".sub"), 'w') as fout:
            fout.write("#!/bin/bash\n")
            fout.write("yhrun -N 1 -n 24 %s < %s > %s\n" % (cmd, inpname, output))

