#!/usr/bin/env python
# _*_ coding: utf-8 _*_

import os
import sys
import shutil

from pymatflow.remote.server import server_handle
from pymatflow.vasp.vasp import vasp
from pymatflow.vasp.base.poscar import vasp_poscar

"""
usage:
"""

class neb_run(vasp):
    """
    Note:
        in vasp running with vtst, the number of cpu cores used must be equal to
        the number of images to be calculated(initial and final will not be calculated). or there might be images that will not be calculated.
    """
    def __init__(self):
        """
        """
        super().__init__()

        self.incar.set_runtype(runtype="neb")
        self.poscars = []

        self.nimage = 5 # default value


    def get_images(self, images):
        """
        images:
            ["first.xyz", "intermediate-1.xyz", "intermediate-2.xyz", ..., "last.xyz"]
        """
        self.poscars = []
        for image in images:
            poscar = vasp_poscar()
            poscar.xyz.get_xyz(image)
            self.poscars.append(poscar)


    def neb(self, directory="tmp-vasp-neb-vtst", runopt="gen", properties={}, auto=0):
        """
        directory: a place for all the generated files
        """
        if runopt == "gen" or runopt == "genrun":
            if os.path.exists(directory):
                shutil.rmtree(directory)
            os.mkdir(directory)
            shutil.copyfile("POTCAR", os.path.join(directory, "POTCAR"))
            for image in self.poscars:
                os.system("cp %s %s/" % (image.xyz.file, directory))

            #self.incar.properties.set_params(properties)
            # =======================================
            # Constructing the input file for VASP
            # =======================================

            # using nebmake.pl to build images
            os.chdir(directory)
            os.mkdir("./is")
            os.mkdir("./fs")
            with open("./is/POSCAR", 'w') as fout:
                self.poscars[0].to_poscar(fout)
            with open("./fs/POSCAR", 'w') as fout:
                self.poscars[-1].to_poscar(fout)

            os.system("nebmake.pl is/POSCAR fs/POSCAR %d" % self.nimage)
            os.chdir("../")

            # 对每个镜像进行几何优化时可以用vasp的优化器也可以用vtst的
            # 设置IOPT为0就使用vasp的优化器, 需要通过IBRION和POTIM来指定
            # 否者就可以通过设置IBRION = 3, POTIM = 0来关闭vasp的优化器
            # 然后设置IOPT= 1或者2来使用VTST的优化器
            # 当设置IOPT > 0，的时候也要注意需要明确设置EDIFFG < 0

            # gen llhpc script
            self.gen_llhpc(directory=directory, cmd="$PMF_VASP_STD_NEB", scriptname="neb.slurm")
            # gen pbs script
            self.gen_pbs(directory=directory, cmd="$PMF_VASP_STD_NEB", scriptname="neb.pbs", jobname=self.run_params["jobname"], nodes=self.run_params["nodes"], ppn=self.run_params["ppn"])
            # gen local bash script
            self.gen_bash(directory=directory, cmd="$PMF_VASP_STD_NEB", scriptname="neb.sh")


        if runopt == "run" or runopt == "genrun":
            # run the neb calculation
            # each image on one core
            os.chdir(directory)
            #os.system("%s vasp" % mpi)
            os.system("bash neb.sh")
            os.chdir("../")
        server_handle(auto=auto, directory=directory, jobfilebase="neb", server=self.run_params["server"])

    #
