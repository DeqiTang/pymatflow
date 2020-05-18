"""
Overall manger for surface PES calculation

Prospects:
we might use a special scanning matter(moving atoms), like those used in AFM, and get the 
image of the surface potential energy, which may help build the dataset for traning of the
surface recoginization and classification for microscopy picture and help use research on
the surface with image processing technique.
"""
import numpy as np
import sys
import os
import shutil
from pymatflow.remote.server import server_handle

import pymatflow.cp2k as cp2k
import pymatflow.qe as qe
import pymatflow.vasp as vasp

"""
"""

class cp2k_run(cp2k.opt_run):
    """
    Note:
        calculate the surface potential energy surface via CP2K.
    """
    def __init__(self):
        """
        """
        super().__init__()
        self.set_geo_opt()
        self.pes_params = {}
        self.set_pes() # set default value

    def set_pes(self, move_atom=[-1], xrange=[0, 1.5, 0.1], yrange=[0, 1.5, 0.5], zshift=0.0, fix_z=1):
        """
        :parma move_atom: the atoms that will move in the calculation, list start from 0.
        :param xrange: x range for moving the specified moving atoms.
        :param: yrange: y range for moving the specified moving atoms
        :param: zshift: z shift for the moving atoms, will shift the z of specified moving atoms by value of zshift
        :param: fix_z: 0 -> do not fix any z of the atoms, 1 -> only fix z of the buttom atoms, 2: fix z of both the buttom and the moving atoms. 
                note x y are all fixed           
        """
        self.pes_params["move_atom"] = move_atom
        self.pes_params["xrange"] = xrange
        self.pes_params["yrange"] = yrange
        self.pes_params["zshift"] = zshift
        self.pes_params["fix_z"] = fix_z

    def run(self, directory="tmp-cp2k-pes-opt", runopt="gen", auto=0):
        if runopt == "gen" or runopt == "genrun":
            if os.path.exists(directory):
                shutil.rmtree(directory)
            os.mkdir(directory)
            shutil.copyfile(self.force_eval.subsys.xyz.file, os.path.join(directory, os.path.basename(self.force_eval.subsys.xyz.file)))

            # use &COORD to input structure rather than &TOPOLOGY
            self.force_eval.subsys.coord.status = True
            self.force_eval.subsys.topology.status = False

            xrange = self.pes_params["xrange"]
            yrange = self.pes_params["yrange"]
            zshift = self.pes_params["zshift"]
            fix_z = self.pes_params["fix_z"]

            os.chdir(directory)
            # generate the input files and the initial trajectory
            os.system("mkdir -p post-processing")
            for deltay in np.arange(yrange[0], yrange[1], yrange[2]):
                for deltax in np.arange(xrange[0], xrange[1], xrange[2]):
                    os.mkdir("_%.3f_%.3f_" % (deltax if np.abs(deltax) >= 0.001 else 0.0, deltay if np.abs(deltay) >= 0.001 else 0.0))

                    for i in self.pes_params["move_atom"]:
                        self.force_eval.subsys.xyz.atoms[i].x += deltax
                        self.force_eval.subsys.xyz.atoms[i].y += deltay
                        # shift z of the specified atoms by self.pes_params["zshfit"]
                        #----------------------------------------------
                        self.force_eval.subsys.xyz.atoms[i].z += zshift

                    # fix xy of no move atoms and move atoms
                    # first fix all atoms
                    for i in range(len(self.poscar.xyz.atoms)):
                        self.poscar.xyz.atoms[i].fix = [True, True, True]
                    # unfix z of moving atoms or z of no moving atoms
                    if fix_z == 0:
                        for i in range(len(self.poscar.xyz.atoms)):
                            self.poscar.xyz.atoms[i].fix[2] = False
                    elif fix_z == 1:
                        for i in self.pes_params["move_atom"]:
                            self.poscar.xyz.atoms[i].fix[2] = False
                    elif fix_z == 2:
                        # nothing need to do
                        pass

                    with open("_%.3f_%.3f_/geo-opt.inp" % (deltax if np.abs(deltax) >= 0.001 else 0.0, deltay if np.abs(deltay) >= 0.001 else 0.0), 'w') as fout:
                        self.glob.to_input(fout)
                        self.force_eval.to_input(fout)
                        self.motion.to_input(fout)

                    with open("post-processing/trajectory-initial.xyz", 'a') as fout:
                        # generate the xyz trajectory file -> (unrelaxed original traj)
                        fout.write("%d\n" % self.force_eval.subsys.xyz.natom)
                        fout.write("deltax: %.3f | deltay: %.3f\n" % (deltax if np.abs(deltax) >= 0.001 else 0.0, deltay if np.abs(deltay) >= 0.001 else 0.0))
                        for atom in self.force_eval.subsys.xyz.atoms:
                            fout.write("%s %.9f %.9f %.9f\n" % (atom.name, atom.x, atom.y, atom.z))

                    for i in self.pes_params["move_atom"]:
                        # now we move the x y z back to the original value
                        self.force_eval.subsys.xyz.atoms[i].x -= deltax
                        self.force_eval.subsys.xyz.atoms[i].y -= deltay
                        self.force_eval.subsys.xyz.atoms[i].z -= zshift


            # write pbs job control script
            with open("pes-relax.pbs", 'w') as fout:
                fout.write("#!/bin/bash\n")
                fout.write("#PBS -N %s\n" % self.run_params["jobname"])
                fout.write("#PBS -l nodes=%d:ppn=%d\n" % (self.run_params["nodes"], self.run_params["ppn"]))
                fout.write("\n")
                fout.write("cd $PBS_O_WORKDIR\n")
                fout.write("NP=`cat $PBS_NODEFILE | wc -l`\n")

                fout.write("for deltay in `seq %.3f %.3f %.3f`\n" % (yrange[0], yrange[2], yrange[1]))
                fout.write("do\n")
                fout.write("for deltax in `seq %.3f %.3f %.3f`\n" % (xrange[0], xrange[2], xrange[1]))
                fout.write("do\n")
                fout.write("  # run the calculation\n")
                fout.write("  cd _${deltax}_${deltay}_\n")
                fout.write("  mpirun -np $NP -machinefile $PBS_NODEFILE %s -in %s > %s\n" % ("$PMF_CP2K", "geo-opt.inp", "geo-opt.out"))
                fout.write("  cd ../\n")
                fout.write("done\n")
                fout.write("done\n")
                # write bash script to generate the xyz trajectory file -> (relaxed)
                #with open("get_traj_relaxed.sh", 'w') as fout:
                fout.write("\n\n")
                fout.write("# code to extract final structure for each combination of deltax and deltay\n")
                #fout.write("#!/bin/bash\n")
                #fout.write("\n")
                #fout.write("\n")
                fout.write("output_trajfile=./post-processing/trajectory-relaxed.xyz\n")
                fout.write("natom=%d\n" % self.force_eval.subsys.xyz.natom)
                fout.write("for deltay in `seq %.3f %.3f %.3f`\n" % (yrange[0], yrange[2], yrange[1]))
                fout.write("do\n")
                fout.write("for deltax in `seq %.3f %.3f %.3f`\n" % (xrange[0], xrange[2], xrange[1]))
                fout.write("do\n")
                fout.write("  echo ${natom} >> ${output_trajfile}\n")
                fout.write("  cat >> ${output_trajfile}<<EOF\n")
                fout.write("deltax: ${deltax} | deltay: ${deltay}\n")
                fout.write("EOF\n")
                fout.write("  cat _${deltax}_${deltay}_/ab-initio-pos-1.xyz | tail -n -${natom} >> ${output_trajfile}\n")
                fout.write("done\n")
                fout.write("done\n")
                # write result analysis file
                #with open("get_pes.sh", 'w') as fout:
                fout.write("\n\n")
                #fout.write("#!/bin/bash\n")
                fout.write("cat > post-processing/pes.data<<EOF\n")
                fout.write("# format: x y energy(Ry)\n")
                fout.write("EOF\n")
                fout.write("\n")
                fout.write("for deltay in `seq %.3f %.3f %.3f`\n" % (yrange[0], yrange[2], yrange[1]))
                fout.write("do\n")
                fout.write("for deltax in `seq %.3f %.3f %.3f`\n" % (xrange[0], xrange[2], xrange[1]))
                fout.write("do\n")
                fout.write("  energy=`cat _${deltax}_${deltay}_/geo-opt.out | grep 'ENERGY| Total FORCE_EVAL ( QS ) energy (a.u.):' | tail -1`\n")
                fout.write("  cat >> post-processing/pes.data<<EOF\n")
                fout.write("${deltax} ${deltay} ${energy:32:-2}\n")
                fout.write("EOF\n")
                fout.write("done\n")
                fout.write("done\n")
                fout.write("\n")
                fout.write("cat > post-processing/plot.gnuplot<<EOF\n")
                fout.write("set term png\n")
                fout.write("set output 'pes.png'\n")
                fout.write("set xlabel 'x'\n")
                fout.write("set ylabel 'y'\n")
                fout.write("splot 'pes.data'\n")
                fout.write("EOF\n")
                fout.write("cd post-processing; gnuplot plot.gnuplot; cd ../\n")

            # write local bash run script
            with open("pes-relax.sh", 'w') as fout:
                fout.write("#!/bin/bash\n")
                fout.write("\n")
                fout.write("for deltay in `seq %.3f %.3f %.3f`\n" % (yrange[0], yrange[2], yrange[1]))
                fout.write("do\n")
                fout.write("for deltax in `seq %.3f %.3f %.3f`\n" % (xrange[0], xrange[2], xrange[1]))
                fout.write("do\n")
                fout.write("  # run the calculation\n")
                fout.write("  cd _${deltax}_${deltay}_\n")
                fout.write("  %s %s -in %s | tee %s\n" % (self.run_params["mpi"], "$PMF_CP2K", "geo-opt.inp", "geo-opt.out"))
                fout.write("  cd ../\n")
                fout.write("done\n")
                fout.write("done\n")
                # write bash script to generate the xyz trajectory file -> (relaxed)
                #with open("get_traj_relaxed.sh", 'w') as fout:
                fout.write("\n\n")
                fout.write("# code to extract final structure for each combination of deltax and deltay\n")
                #fout.write("#!/bin/bash\n")
                #fout.write("\n")
                #fout.write("\n")
                fout.write("output_trajfile=./post-processing/trajectory-relaxed.xyz\n")
                fout.write("natom=%d\n" % self.force_eval.subsys.xyz.natom)
                fout.write("for deltay in `seq %.3f %.3f %.3f`\n" % (yrange[0], yrange[2], yrange[1]))
                fout.write("do\n")
                fout.write("for deltax in `seq %.3f %.3f %.3f`\n" % (xrange[0], xrange[2], xrange[1]))
                fout.write("do\n")
                fout.write("  echo ${natom} >> ${output_trajfile}\n")
                fout.write("  cat >> ${output_trajfile}<<EOF\n")
                fout.write("deltax: ${deltax} | deltay: ${deltay}\n")
                fout.write("EOF\n")
                fout.write("  cat _${deltax}_${deltay}_/ab-initio-pos-1.xyz | tail -n -${natom} >> ${output_trajfile}\n")
                fout.write("done\n")
                fout.write("done\n")
                # write result analysis file
                #with open("get_pes.sh", 'w') as fout:
                fout.write("\n\n")
                #fout.write("#!/bin/bash\n")
                fout.write("cat > post-processing/pes.data<<EOF\n")
                fout.write("# format: x y energy(Ry)\n")
                fout.write("EOF\n")
                fout.write("\n")
                fout.write("for deltay in `seq %.3f %.3f %.3f`\n" % (yrange[0], yrange[2], yrange[1]))
                fout.write("do\n")
                fout.write("for deltax in `seq %.3f %.3f %.3f`\n" % (xrange[0], xrange[2], xrange[1]))
                fout.write("do\n")
                fout.write("  energy=`cat _${deltax}_${deltay}_/geo-opt.out | grep 'ENERGY| Total FORCE_EVAL ( QS ) energy (a.u.):' | tail -1`\n")
                fout.write("  cat >> post-processing/pes.data<<EOF\n")
                fout.write("${deltax} ${deltay} ${energy:32:-2}\n")
                fout.write("EOF\n")
                fout.write("done\n")
                fout.write("done\n")
                fout.write("\n")
                fout.write("cat > post-processing/plot.gnuplot<<EOF\n")
                fout.write("set term png\n")
                fout.write("set output 'pes.png'\n")
                fout.write("splot 'pes.data'\n")
                fout.write("set xlabel 'x'\n")
                fout.write("set ylabel 'y'\n")                
                fout.write("EOF\n")
                fout.write("cd post-processing; gnuplot plot.gnuplot; cd ../\n")

            os.chdir("../")
        if runopt == "genrun" or runopt == "run":
            os.chdir(directory)
            os.system("bash pes-relax.sh")
            oschdir("../")
        # server handle
        server_handle(auto=auto, directory=directory, jobfilebase="pes-relax", server=self.run_params["server"])


class qe_run(qe.opt_run):
    """
    Note:
        calculate the surface potential energy surface via Quantum ESPRESSO.
    """
    def __init__(self):
        """
        """
        super().__init__()

        self.pes_params = {}
        self.set_pes() # set default value

    def set_pes(self, move_atom=[-1], xrange=[0, 1.5, 0.1], yrange=[0, 1.5, 0.5], zshift=0.0, fix_z=1):
        """
        :parma move_atom: the atoms that will move in the calculation, list start from 0.
        :param xrange: x range for moving the specified moving atoms.
        :param: yrange: y range for moving the specified moving atoms
        :param: zshift: z shift for the moving atoms, will shift the z of specified moving atoms by value of zshift
        :param: fix_z: 0 -> do not fix any z of the atoms, 1 -> only fix z of the buttom atoms, 2: fix z of both the buttom and the moving atoms. 
                note x y are all fixed   
        """
        self.pes_params["move_atom"] = move_atom
        self.pes_params["xrange"] = xrange
        self.pes_params["yrange"] = yrange
        self.pes_params["zshift"] = zshift
        self.pes_params["fix_z"] = fix_z

    def run(self, directory="tmp-qe-pes-opt", runopt="gen", auto=0):
        if runopt == "gen" or runopt == "genrun":
            if os.path.exists(directory):
                shutil.rmtree(directory)
            os.mkdir(directory)

            shutil.copyfile(self.arts.xyz.file, os.path.join(directory, os.path.basename(self.arts.xyz.file)))
            all_upfs = [s for s in os.listdir() if s.split(".")[-1] == "UPF"]
            for element in self.arts.xyz.specie_labels:
                for upf in all_upfs:
                    if upf.split(".")[0] == element:
                        shutil.copyfile(upf, os.path.join(directory, upf))
                        break
            self.arts.pseudo.dir = os.path.abspath(directory)
            self.control.set_params({"pseudo_dir": os.path.abspath(directory)})
            #

            xrange = self.pes_params["xrange"]
            yrange = self.pes_params["yrange"]
            zshift = self.pes_params["zshift"]
            fix_z = self.pes_params["fix_z"]
            os.chdir(directory)
            # generate the input files and the initial trajectory
            os.system("mkdir -p post-processing")
            # first iterate y and iterate x which is good for post processing to get the imgage
            for deltay in np.arange(yrange[0], yrange[1], yrange[2]):
                for deltax in np.arange(xrange[0], xrange[1], xrange[2]): 
                    # to avoid float -0.000 be translated to string -0.000 we use 0.0 when value ==0 whether it is 0.0 or -0.0
                    os.mkdir("_%.3f_%.3f_" % (deltax if np.abs(deltax) >= 0.001 else 0.0, deltay if np.abs(deltay) >= 0.001 else 0.0))

                    for i in self.pes_params["move_atom"]:
                        self.arts.xyz.atoms[i].x += deltax
                        self.arts.xyz.atoms[i].y += deltay
                        # shift z of the specified atoms by self.pes_params["zshfit"]
                        #----------------------------------------------
                        self.arts.xyz.atoms[i].z += zshift

                    # fix xy of no move atoms and move atoms
                    # first fix all atoms
                    for i in range(len(self.arts.xyz.atoms)):
                        self.arts.xyz.atoms[i].fix = [True, True, True]
                    # unfix z of moving atoms or z of no moving atoms
                    if fix_z == 0:
                        for i in range(len(self.arts.xyz.atoms)):
                            self.arts.xyz.atoms[i].fix[2] = False
                    elif fix_z == 1:
                        for i in self.pes_params["move_atom"]:
                            self.arts.xyz.atoms[i].fix[2] = False
                    elif fix_z == 2:
                        # nothing need to do
                        pass

                    with open("_%.3f_%.3f_/relax.in" % (deltax if np.abs(deltax) >= 0.001 else 0.0, deltay if np.abs(deltay) >= 0.001 else 0.0), 'w') as fout:
                        self.control.to_in(fout)
                        self.system.to_in(fout)
                        self.electrons.to_in(fout)
                        self.ions.to_in(fout)
                        self.arts.to_in(fout)

                    with open("post-processing/trajectory-initial.xyz", 'a') as fout:
                        # generate the xyz trajectory file -> (unrelaxed original traj)
                        fout.write("%d\n" % self.arts.xyz.natom)
                        fout.write("deltax: %.3f | deltay: %.3f\n" % (deltax if np.abs(deltax) >= 0.001 else 0.0, deltay if np.abs(deltay) >= 0.001 else 0.0))
                        for atom in self.arts.xyz.atoms:
                            fout.write("%s %.9f %.9f %.9f\n" % (atom.name, atom.x, atom.y, atom.z))

                    for i in self.pes_params["move_atom"]:
                        # now we move the x y z back to the original value
                        self.arts.xyz.atoms[i].x -= deltax
                        self.arts.xyz.atoms[i].y -= deltay
                        self.arts.xyz.atoms[i].z -= zshift


            # write pbs job control script
            with open("pes-relax.pbs", 'w') as fout:
                fout.write("#!/bin/bash\n")
                fout.write("#PBS -N %s\n" % self.run_params["jobname"])
                fout.write("#PBS -l nodes=%d:ppn=%d\n" % (self.run_params["nodes"], self.run_params["ppn"]))
                if "queue" in self.run_params and self.run_params["queue"] != None:
                    fout.write("#PBS -q %s\n" %self.run_params["queue"])                      
                fout.write("\n")
                fout.write("cd $PBS_O_WORKDIR\n")
                fout.write("NP=`cat $PBS_NODEFILE | wc -l`\n")

                # do not add -w to seq
                fout.write("for deltay in `seq %.3f %.3f %.3f`\n" % (yrange[0], yrange[2], yrange[1]))
                fout.write("do\n")
                fout.write("for deltax in `seq %.3f %.3f %.3f`\n" % (xrange[0], xrange[2], xrange[1]))
                fout.write("do\n")
                fout.write("  # run the calculation\n")
                fout.write("  cd _${deltax}_${deltay}_\n")
                fout.write("  mpirun -np $NP -machinefile $PBS_NODEFILE %s < %s > %s\n" % ("$PMF_PWX", "relax.in", "relax.out"))
                fout.write("  cd ../\n")
                fout.write("done\n")
                fout.write("done\n")
                # write bash script to generate the xyz trajectory file -> (relaxed)
                #with open("get_traj_relaxed.sh", 'w') as fout:
                fout.write("\n\n")
                fout.write("# code to extract final structure for each combination of deltax and deltay\n")
                #fout.write("#!/bin/bash\n")
                #fout.write("\n")
                #fout.write("\n")
                fout.write("output_trajfile=./post-processing/trajectory-relaxed.xyz\n")
                fout.write("natom=%d\n" % self.arts.xyz.natom)
                fout.write("for deltay in `seq %.3f %.3f %.3f`\n" % (yrange[0], yrange[2], yrange[1]))
                fout.write("do\n")
                fout.write("for deltax in `seq %.3f %.3f %.3f`\n" % (xrange[0], xrange[2], xrange[1]))
                fout.write("do\n")
                fout.write("  post-qe-relax.py -d _${deltax}_${deltay}_\n")
                fout.write("  echo ${natom} >> ${output_trajfile}\n")
                fout.write("  cat >> ${output_trajfile}<<EOF\n")
                fout.write("deltax: ${deltax} | deltay: ${deltay}\n")
                fout.write("EOF\n")
                fout.write("  cat _${deltax}_${deltay}_/post-processing/trajectory.xyz | tail -n -${natom} >> ${output_trajfile}\n")
                fout.write("done\n")
                fout.write("done\n")
                # write result analysis file
                #with open("get_pes.sh", 'w') as fout:
                fout.write("\n\n")
                #fout.write("#!/bin/bash\n")
                fout.write("cat > post-processing/pes.data<<EOF\n")
                fout.write("# format: x y energy(Ry)\n")
                fout.write("EOF\n")
                fout.write("\n")
                fout.write("for deltay in `seq %.3f %.3f %.3f`\n" % (yrange[0], yrange[2], yrange[1]))
                fout.write("do\n")
                fout.write("for deltax in `seq %.3f %.3f %.3f`\n" % (xrange[0], xrange[2], xrange[1]))
                fout.write("do\n")
                fout.write("  energy=`cat _${deltax}_${deltay}_/relax.out | grep '!    total energ' | tail -1`\n")
                fout.write("  cat >> post-processing/pes.data<<EOF\n")
                fout.write("${deltax} ${deltay} ${energy:32:-2}\n")
                fout.write("EOF\n")
                fout.write("done\n")
                fout.write("done\n")
                fout.write("\n")
                fout.write("cat > post-processing/plot.gnuplot<<EOF\n")
                fout.write("set term png\n")
                fout.write("set output 'pes.png'\n")
                fout.write("splot 'pes.data'\n")
                fout.write("set xlabel 'x'\n")
                fout.write("set ylabel 'y'\n")                
                fout.write("EOF\n")
                fout.write("cd post-processing; gnuplot plot.gnuplot; cd ../\n")

            # write local bash run script
            with open("pes-relax.sh", 'w') as fout:
                fout.write("#!/bin/bash\n")
                fout.write("#\n")
                fout.write("for deltay in `seq %.3f %.3f %.3f`\n" % (yrange[0], yrange[2], yrange[1]))
                fout.write("do\n")
                fout.write("for deltax in `seq %.3f %.3f %.3f`\n" % (xrange[0], xrange[2], xrange[1]))
                fout.write("do\n")
                fout.write("  # run the calculation\n")
                fout.write("  cd _${deltax}_${deltay}_\n")
                fout.write("  %s %s < %s > %s\n" % (self.run_params["mpi"], "$PMF_PWX", "relax.in", "relax.out"))
                fout.write("  cd ../\n")
                fout.write("done\n")
                fout.write("done\n")
                # write bash script to generate the xyz trajectory file -> (relaxed)
                #with open("get_traj_relaxed.sh", 'w') as fout:
                fout.write("\n\n")
                fout.write("# code to extract final structure for each combination of deltax and deltay\n")
                #fout.write("#!/bin/bash\n")
                #fout.write("\n")
                #fout.write("\n")
                fout.write("output_trajfile=./post-processing/trajectory-relaxed.xyz\n")
                fout.write("natom=%d\n" % self.arts.xyz.natom)
                fout.write("for deltay in `seq %.3f %.3f %.3f`\n" % (yrange[0], yrange[2], yrange[1]))
                fout.write("do\n")
                fout.write("for deltax in `seq %.3f %.3f %.3f`\n" % (xrange[0], xrange[2], xrange[1]))
                fout.write("do\n")
                fout.write("  post-qe-relax.py -d _${deltax}_${deltay}_\n")
                fout.write("  echo ${natom} >> ${output_trajfile}\n")
                fout.write("  cat >> ${output_trajfile}<<EOF\n")
                fout.write("deltax: ${deltax} | deltay: ${deltay}\n")
                fout.write("EOF\n")
                fout.write("  cat _${deltax}_${deltay}_/post-processing/trajectory.xyz | tail -n -${natom} >> ${output_trajfile}\n")
                fout.write("done\n")
                fout.write("done\n")
                # write result analysis file
                #with open("get_pes.sh", 'w') as fout:
                fout.write("\n\n")
                #fout.write("#!/bin/bash\n")
                fout.write("cat > post-processing/pes.data<<EOF\n")
                fout.write("# format: x y energy(Ry)\n")
                fout.write("EOF\n")
                fout.write("\n")
                fout.write("for deltay in `seq %.3f %.3f %.3f`\n" % (yrange[0], yrange[2], yrange[1]))
                fout.write("do\n")
                fout.write("for deltax in `seq %.3f %.3f %.3f`\n" % (xrange[0], xrange[2], xrange[1]))
                fout.write("do\n")
                fout.write("  energy=`cat _${deltax}_${deltay}_/relax.out | grep '!    total energ' | tail -1`\n")
                fout.write("  cat >> post-processing/pes.data<<EOF\n")
                fout.write("${deltax} ${deltay} ${energy:32:-2}\n")
                fout.write("EOF\n")
                fout.write("done\n")
                fout.write("done\n")
                fout.write("\n")
                fout.write("cat > post-processing/plot.gnuplot<<EOF\n")
                fout.write("set term png\n")
                fout.write("set output 'pes.png'\n")
                fout.write("splot 'pes.data'\n")
                fout.write("set xlabel 'x'\n")
                fout.write("set ylabel 'y'\n")                
                fout.write("EOF\n")
                fout.write("cd post-processing; gnuplot plot.gnuplot; cd ../\n")

            os.chdir("../")
        if runopt == "genrun" or runopt == "run":
            os.chdir(directory)
            os.system("bash pes-relax.sh")
            oschdir("../")
        # server handle
        server_handle(auto=auto, directory=directory, jobfilebase="pes-relax", server=self.run_params["server"])


class vasp_run(vasp.opt_run):
    """
    Note:
        calculate the surface potential energy surface via VASP.
    """
    def __init__(self):
        """
        """
        super().__init__()

        self.pes_params = {}
        self.set_pes() # set default value

        self.batch_x_y = None

    def set_pes(self, move_atom=[-1], xrange=[0, 1.5, 0.1], yrange=[0, 1.5, 0.5], zshift=0.0, fix_z=1):
        """
        :parma move_atom: the atoms that will move in the calculation, list start from 0.
        :param xrange: x range for moving the specified moving atoms.
        :param: yrange: y range for moving the specified moving atoms
        :param: zshift: z shift for the moving atoms, will shift the z of specified moving atoms by value of zshift
        :param: fix_z: 0 -> do not fix any z of the atoms, 1 -> only fix z of the buttom atoms, 2: fix z of both the buttom and the moving atoms. 
                note x y are all fixed            
        """
        self.pes_params["move_atom"] = move_atom
        self.pes_params["xrange"] = xrange
        self.pes_params["yrange"] = yrange
        self.pes_params["zshift"] = zshift
        self.pes_params["fix_z"] = fix_z

    def run(self, directory="tmp-vasp-pes-opt", runopt="gen", auto=0):
        """
        """
        xrange = self.pes_params["xrange"]
        yrange = self.pes_params["yrange"]
        zshift = self.pes_params["zshift"]
        fix_z = self.pes_params["fix_z"]
        
        nx = len(np.arange(xrange[0], xrange[1], xrange[2]))
        ny = len(np.arange(yrange[0], yrange[1], yrange[2]))
        if self.batch_x_y == None:
            # namely all in one batch
            self.batch_x_y = [nx, ny]
        else:
            pass
        
        if nx % self.batch_x_y[0] == 0:
            n_batch_x = int(nx / self.batch_x_y[0])
        else:
            n_batch_x = int(nx / self.batch_x_y[0]) + 1
        
        if ny % self.batch_x_y[1] == 0:
            n_batch_y = int(ny / self.batch_x_y[1])
        else:
            n_batch_y = int(ny / self.batch_x_y[1]) + 1         
        
        #

        if runopt == "gen" or runopt == "genrun":
            if os.path.exists(directory):
                shutil.rmtree(directory)
            os.mkdir(directory)

            shutil.copyfile("POTCAR", os.path.join(directory, "POTCAR"))
            os.system("cp %s %s/" % (self.poscar.xyz.file, directory))
            #
            #xrange = self.pes_params["xrange"]
            #yrange = self.pes_params["yrange"]
            #zshift = self.pes_params["zshift"]
            #fix_z = self.pes_params["fix_z"]
            os.chdir(directory)
            # generate the input files and the initial trajectory
            os.system("mkdir -p post-processing")
            # first iterate y and iterate x which is good for post processing to get the imgage
            for deltay in np.arange(yrange[0], yrange[1], yrange[2]):
                for deltax in np.arange(xrange[0], xrange[1], xrange[2]): 
                    # to avoid float -0.000 be translated to string -0.000 we use 0.0 when value ==0 whether it is 0.0 or -0.0
                    os.mkdir("_%.3f_%.3f_" % (deltax if np.abs(deltax) >= 0.001 else 0.0, deltay if np.abs(deltay) >= 0.001 else 0.0))

                    for i in self.pes_params["move_atom"]:
                        self.poscar.xyz.atoms[i].x += deltax
                        self.poscar.xyz.atoms[i].y += deltay
                        # shift z of the specified atoms by self.pes_params["zshfit"]
                        #----------------------------------------------
                        self.poscar.xyz.atoms[i].z += zshift

                    # fix xy of no move atoms and move atoms
                    # first fix all atoms
                    for i in range(len(self.poscar.xyz.atoms)):
                        self.poscar.xyz.atoms[i].fix = [True, True, True]
                    # unfix z of moving atoms or z of no moving atoms
                    if fix_z == 0:
                        for i in range(len(self.poscar.xyz.atoms)):
                            self.poscar.xyz.atoms[i].fix[2] = False
                    elif fix_z == 1:
                        for i in self.pes_params["move_atom"]:
                            self.poscar.xyz.atoms[i].fix[2] = False
                    elif fix_z == 2:
                        # nothing need to do
                        pass

                    with open("_%.3f_%.3f_/POSCAR" % (deltax if np.abs(deltax) >= 0.001 else 0.0, deltay if np.abs(deltay) >= 0.001 else 0.0), 'w') as fout:
                        self.poscar.to_poscar(fout)                    

                    with open("post-processing/trajectory-initial.xyz", 'a') as fout:
                        # generate the xyz trajectory file -> (unrelaxed original traj)
                        fout.write("%d\n" % self.poscar.xyz.natom)
                        fout.write("deltax: %.3f | deltay: %.3f\n" % (deltax if np.abs(deltax) >= 0.001 else 0.0, deltay if np.abs(deltay) >= 0.001 else 0.0))
                        for atom in self.poscar.xyz.atoms:
                            fout.write("%s %.9f %.9f %.9f\n" % (atom.name, atom.x, atom.y, atom.z))

                    for i in self.pes_params["move_atom"]:
                        # now we move the x y z back to the original value
                        self.poscar.xyz.atoms[i].x -= deltax
                        self.poscar.xyz.atoms[i].y -= deltay
                        self.poscar.xyz.atoms[i].z -= zshift
            
            with open("INCAR", 'w') as fout:
                self.incar.to_incar(fout)

            with open("KPOINTS", 'w') as fout:
                self.kpoints.to_kpoints(fout)            

            # write local bash run script
            with open("pes-relax.sh", 'w') as fout:
                fout.write("#!/bin/bash\n")
                fout.write("#\n")
                fout.write("for deltay in `seq %.3f %.3f %.3f`\n" % (yrange[0], yrange[2], yrange[1]))
                fout.write("do\n")
                fout.write("for deltax in `seq %.3f %.3f %.3f`\n" % (xrange[0], xrange[2], xrange[1]))
                fout.write("do\n")
                fout.write("  # run the calculation\n")
                fout.write("  cd _${deltax}_${deltay}_\n")
                fout.write("  cp ../INCAR .; cp ../POTCAR .; cp ../KPOINTS .;\n")
                fout.write("  mpirun -np $NP -machinefile $PBS_NODEFILE %s\n" % ("$PMF_VASP_STD"))
                fout.write("  cd ../\n")
                fout.write("done\n")
                fout.write("done\n")
                fout.write("cd post-processing; bash get_pes.sh; bash get_trajectory.sh; cd ../\n")

            # result collection bash script
            with open("post-processing/get_pes.sh", 'w') as fout:
                fout.write("#!/bin/bash\n")
                # write pes analysis file
                fout.write("\n\n")
                #fout.write("#!/bin/bash\n")
                fout.write("cat > ./pes.data<<EOF\n")
                fout.write("# format: x y energy(Ry)\n")
                fout.write("EOF\n")
                fout.write("\n")
                fout.write("for deltay in `seq %.3f %.3f %.3f`\n" % (yrange[0], yrange[2], yrange[1]))
                fout.write("do\n")
                fout.write("for deltax in `seq %.3f %.3f %.3f`\n" % (xrange[0], xrange[2], xrange[1]))
                fout.write("do\n")
                fout.write("  energy=`cat ../_${deltax}_${deltay}_/OUTCAR | grep 'entropy=' | tail -1 | cut -d \"=\" -f 2 | cut -d \"e\" -f 1`\n")
                fout.write("  cat >> ./pes.data<<EOF\n")
                fout.write("${deltax} ${deltay} ${energy}\n")
                fout.write("EOF\n")
                fout.write("done\n")
                fout.write("done\n")
                fout.write("\n")
                fout.write("cat > ./plot.gnuplot<<EOF\n")
                fout.write("set term png\n")
                fout.write("set output 'pes.png'\n")
                fout.write("splot 'pes.data'\n")
                fout.write("set xlabel 'x'\n")
                fout.write("set ylabel 'y'\n")                
                fout.write("EOF\n")
                fout.write("gnuplot plot.gnuplot\n")

            with open("post-processing/get_trajectory.sh", 'w') as fout:
                fout.write("#!/bin/bash\n")
                # write bash script to generate the xyz trajectory file -> (relaxed)
                #with open("get_traj_relaxed.sh", 'w') as fout:
                fout.write("\n\n")
                fout.write("# code to extract final structure for each combination of deltax and deltay\n")
                #fout.write("#!/bin/bash\n")
                #fout.write("\n")
                #fout.write("\n")
                fout.write("output_trajfile=./trajectory-relaxed.xyz\n")
                fout.write("natom=%d\n" % self.poscar.xyz.natom)
                fout.write("for deltay in `seq %.3f %.3f %.3f`\n" % (yrange[0], yrange[2], yrange[1]))
                fout.write("do\n")
                fout.write("for deltax in `seq %.3f %.3f %.3f`\n" % (xrange[0], xrange[2], xrange[1]))
                fout.write("do\n")
                fout.write("  sflow convert -i ../_${deltax}_${deltay}_/CONTCAR -o ../_${deltax}_${deltay}_/optimized.xyz\n")
                fout.write("  echo ${natom} >> ${output_trajfile}\n")
                fout.write("  cat >> ${output_trajfile}<<EOF\n")
                fout.write("deltax: ${deltax} | deltay: ${deltay}\n")
                fout.write("EOF\n")
                fout.write("  cat ../_${deltax}_${deltay}_/optimized.xyz | tail -n -${natom} >> ${output_trajfile}\n")
                fout.write("done\n")
                fout.write("done\n")


            # batch submitting script
            # dividing structures into groups, each group has one sub script   
            
            # generate job script for each batch
            for i_batch_y in range(n_batch_y):
                for i_batch_x in range(n_batch_x):
                    # write pbs job control script
                    with open("pes-relax-batch-%d-%d.pbs" % (i_batch_x, i_batch_y), 'w') as fout:
                        fout.write("#!/bin/bash\n")
                        fout.write("#PBS -N %s-%d-%d\n" % (self.run_params["jobname"], i_batch_x, i_batch_y))
                        fout.write("#PBS -l nodes=%d:ppn=%d\n" % (self.run_params["nodes"], self.run_params["ppn"]))
                        if "queue" in self.run_params and self.run_params["queue"] != None:
                            fout.write("#PBS -q %s\n" %self.run_params["queue"])                      
                        fout.write("\n")
                        fout.write("cd $PBS_O_WORKDIR\n")
                        fout.write("NP=`cat $PBS_NODEFILE | wc -l`\n")

                        y_start = yrange[0] + i_batch_y * self.batch_x_y[1] * yrange[2]
                        y_end = yrange[0] + (i_batch_y+1) * self.batch_x_y[1] * yrange[2] - yrange[2] / 2
                        # - yrange[2] / 2, so that the last value is ignored which is actually the begining of next batch
                        if y_end  > yrange[1]:
                            y_end = yrange[1]
                        x_start = xrange[0] + i_batch_x * self.batch_x_y[0] * xrange[2]
                        x_end = xrange[0] + (i_batch_x+1) * self.batch_x_y[0] * xrange[2] - xrange[2] / 2
                        # - xrange[2] / 2, so that the last value is ignored which is actually the begining of next batch
                        if x_end > xrange[1]:
                            x_end = xrange[1]

                        # do not add -w to seq
                        fout.write("for deltay in `seq %.3f %.3f %.3f`\n" % (y_start, yrange[2], y_end))
                        fout.write("do\n")
                        fout.write("for deltax in `seq %.3f %.3f %.3f`\n" % (x_start, xrange[2], x_end))
                        fout.write("do\n")
                        fout.write("  # run the calculation\n")
                        fout.write("  cd _${deltax}_${deltay}_\n")
                        fout.write("  cp ../INCAR .; cp ../POTCAR .; cp ../KPOINTS .;\n")
                        fout.write("  mpirun -np $NP -machinefile $PBS_NODEFILE %s\n" % ("$PMF_VASP_STD"))
                        fout.write("  cd ../\n")
                        fout.write("done\n")
                        fout.write("done\n")


            os.chdir("../")
        if runopt == "genrun" or runopt == "run":
            os.chdir(directory)
            os.system("bash pes-relax.sh")
            oschdir("../")
        # server handle
        for i_batch_y in range(n_batch_y):
            for i_batch_x in range(n_batch_x):
                #print("i_batch_x: %d, ibatch_y: %d\n" % (i_batch_x, i_batch_y))
                server_handle(auto=auto, directory=directory, jobfilebase="pes-relax-batch-%d-%d" % (i_batch_x, i_batch_y), server=self.run_params["server"])

