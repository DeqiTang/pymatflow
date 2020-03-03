#!/usr/bin/env python
# _*_ coding: utf-8 _*_

import os
import sys
import numpy as np
import shutil

from pymatflow.siesta.base.system import siesta_system
from pymatflow.siesta.base.electrons import siesta_electrons
from pymatflow.siesta.base.properties import siesta_properties

from pymatflow.siesta.base.transiesta import siesta_transiesta
from pymatflow.siesta.base.tbtrans import siesta_tbtrans


"""
Note:
    Gamma-only calculation on electrode is not feasible in transiesta.
"""

class ts_run:
    """
    """
    def __init__(self):
        self.system_electrodes = []
        self.system_device = None


        self.electrons = siesta_electrons()

        self.transiesta = siesta_transiesta()

        self.tbtrans = siesta_tbtrans()

        self.electrons.basic_setting()
        #self.electrons.params["SolutionMethod"] = "transiesta"

    def get_electrodes_device(self, electrodes, device):
        self.system_electrodes = []
        #self.properties = []
        for xyz_f in electrodes:
            electrode = siesta_system()
            electrode.xyz.get_xyz(xyz_f)
            self.system_electrodes.append(electrode)

        self.system_device = siesta_system()
        self.system_device.xyz.get_xyz(device)


    def ts(self, directory="tmp-siesta-ts", inpname="transiesta.fdf", output="transiesta.out",
            mpi="", runopt="gen", electrons={}, properties=[], kpoints_mp=[1, 1, 1], bias=[0, 1, 0.2]):
        if runopt == "gen" or runopt == "genrun":
            if os.path.exists(directory):
                shutil.rmtree(directory)
            os.mkdir(directory)

            for electrode in self.system_electrodes:
                for element in electrode.xyz.specie_labels:
                    shutil.copyfile("%s.psf" % element, os.path.join(directory, "%s.psf" % element))
                shutil.copyfile(electrode.xyz.file, os.path.join(directory, electrode.xyz.file))
            for element in self.system_device.xyz.specie_labels:
                shutil.copyfile("%s.psf" % element, os.path.join(directory, "%s.psf" % element))
            shutil.copyfile(self.system_device.xyz.file, os.path.join(directory, self.system_device.xyz.file))


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
                    self.electrons.set_params({
                        "SolutionMethod": None,
                        #"TS.WriteHS": "true",
                        #"TS.HS.Save": "true",
                        #"TS.DE.Save": "true",
                        })
                    self.electrons.to_fdf(fout)
                    #self.properties.to_fdf(fout)
                for element in self.system_electrodes[i].xyz.specie_labels:
                    shutil.copyfile("%s.psf" % element, os.path.join(directory, "electrodes", "electrode-%d" % i, "%s.psf" % element))
                    shutil.copyfile(self.system_electrodes[i].xyz.file, os.path.join(directory, "electrodes", "electrode-%d/%s" % (i, self.system_electrodes[i].xyz.file)))

            # now set the transiesta calculation of the device(including scattering zone)
            os.mkdir(os.path.join(directory, "device"))
            for v in np.arange(bias[0], bias[1], bias[2]):
                os.mkdir(os.path.join(directory, "device", "bias-%.6f" % v))
                with open(os.path.join(directory, "device", "bias-%.6f" % v, "device.fdf"), 'w') as fout:
                    self.transiesta.ts["Voltage"] = v
                    self.transiesta.set_params(ts={
                        "TS.WriteHS": "true",
                        "TS.HS.Save": "true",
                        "TS.DE.Save": "true",
                    })
                    self.transiesta.to_fdf(fout)
                    self.tbtrans.to_fdf(fout)
                    self.system_device.to_fdf(fout)
                    self.electrons.set_params({
                        "SolutionMethod": "transiesta",
                    })
                    self.electrons.to_fdf(fout)
                for element in self.system_device.xyz.specie_labels:
                    shutil.copyfile("%s.psf" % element, os.path.join(directory, "device/bias-%.6f/%s.psf" % (v, element)))
                shutil.copyfile(self.system_device.xyz.file, os.path.join(directory, "device/bias-%.6f/" % v, self.system_device.xyz.file))

            # gen yhbatch script
            with open(os.path.join(directory, "job.sub"), 'w') as fout:
                fout.write("#!/bin/bash\n")
                fout.write("\n")
                fout.write("# ----------------------------------------------------------------------\n")
                fout.write("# run the electrode calculation task\n")
                fout.write("# ----------------------------------------------------------------------\n")
                fout.write("\n")
                for i in range(len(self.system_electrodes)):
                    fout.write("# running on electrode: %d\n" % i)
                    fout.write("cd electrodes/electrode-%d\n" % i)
                    fout.write("%s siesta --electrode < %s > %s\n" % (mpi, "electrode.fdf", "electrode.fdf.out"))
                    fout.write("cd ../../\n\n")
                fout.write("\n")
                fout.write("# ----------------------------------------------------------------------\n")
                fout.write("# run the scattering calculation tasks(different bias)\n")
                fout.write("# ----------------------------------------------------------------------\n")
                fout.write("cd device\n\n\n")
                bias_all = np.arange(bias[0], bias[1], bias[2])
                fout.write("# *************************************************first bias calculation\n")
                fout.write("cd bias-%.6f\n" % bias_all[0])
                fout.write("%s transiesta < %s > %s\n" % (mpi, "device.fdf", "device.fdf.transiesta.out"))
                fout.write("%s tbtrans < %s > %s\n" % (mpi, "device.fdf", "device.fdf.tbtrans.out"))
                fout.write("cat device.fdf.tbtrans.out | grep \" Voltage, Current(A) =\" | awk \'{print $4\"   \"$5}\' >> ../IV.dat\n")
                fout.write("cd ../\n\n")
                fout.write("# *************************************************other bias calculation\n")
                fout.write("# copying previous bias output SystemLabel.TSDE\n")
                fout.write("# to current bias calculation this can drastically\n")
                fout.write("# improve convergence, recommended by transiesta manual\n")
                fout.write("\n")
                for i in range(1, len(bias_all)):
                    fout.write("cd bias-%.6f\n" % bias_all[i])
                    fout.write("cp ../bias-%.6f/siesta.TSDE ./\n" % bias_all[i-1])
                    fout.write("%s transiesta < %s > %s\n" % (mpi, "device.fdf", "device.fdf.transiesta.out"))
                    fout.write("%s tbtrans < %s > %s\n" % (mpi, "device.fdf", "device.fdf.tbtrans.out"))
                    fout.write("cat device.fdf.tbtrans.out | grep \" Voltage, Current(A) =\" | awk \'{print $4\"   \"$5}\' >> ../IV.dat\n")
                    fout.write("cd ../\n\n")

        if runopt == "run" or runopt == "genrun":
            os.chdir(directory)
            # run the electrode calculation
            for i in range(len(self.system_electrodes)):
                os.chdir(os.path.join("electrodes", "electrode-%d" % i))
                os.system("%s siesta --electrode < %s | tee %s" % (mpi, "electrode.fdf", "electrode.fdf.out"))
                os.chdir("../../")

            # run the scattering zone
            os.chdir("device")
            for v in np.arange(bias[0], bias[1], bias[2]):
                os.chdir("bias-%.6f" % v)
                # this can drastically improve convergence, recommended by transiesta manual
                os.system("cp ../bias-%.6f/siesta.TSDE ./" % (v-bias[2]))
                os.system("%s transiesta < %s | tee %s" % (mpi, "device.fdf", "device.fdf.transiesta.out"))
                os.system('%s tbtrans < %s | tee %s' % (mpi, "device.fdf", "device.fdf.tbtrans.out"))
                os.system("cat device.fdf.tbtrans.out | grep \" Voltage, Current(A) =\" | awk \'{print $4\"   \"$5}\' >> ../IV.dat")
                os.chdir("../")
            os.chdir("../")
