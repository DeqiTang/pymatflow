"""
Overall manger for surface PES calculation
"""
import numpy as np
import sys
import os
import shutil
from pymatflow.remote.server import server_handle

import pymatflow.cp2k as cp2k
import pymatflow.qe as qe

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

    def set_pes(self, move_atom=[-1], xrange=[0, 1.5, 0.1], yrange=[0, 1.5, 0.5], zshift=0.0):
        """
        :parma move_atom: the atoms that will move in the calculation, list start from 0.
        :param xrange: x range for moving the specified moving atoms.
        :param: yrange: y range for moving the specified moving atoms
        :param: zshift: z shift for the moving atoms, will shift the z of specified moving atoms by value of zshift
        """
        self.pes_params["move_atom"] = move_atom
        self.pes_params["xrange"] = xrange
        self.pes_params["yrange"] = yrange
        self.pes_params["zshift"] = zshift

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
            os.chdir(directory)
            # generate the input files and the initial trajectory
            os.system("mkdir -p post-processing")
            for deltax in np.arange(xrange[0], xrange[1], xrange[2]):
                for deltay in np.arange(yrange[0], yrange[1], yrange[2]):
                    os.mkdir("_%.3f_%.3f_" % (deltax if np.abs(deltax) >= 0.001 else 0.0, deltay if np.abs(deltay) >= 0.001 else 0.0))

                    for i in self.pes_params["move_atom"]:
                        self.force_eval.subsys.xyz.atoms[i].x += deltax
                        self.force_eval.subsys.xyz.atoms[i].y += deltay
                        # shift z of the specified atoms by self.pes_params["zshfit"]
                        #----------------------------------------------
                        self.force_eval.subsys.xyz.atoms[i].z += zshift

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

                fout.write("for deltax in `seq %.3f %.3f %.3f`\n" % (xrange[0], xrange[2], xrange[1]))
                fout.write("do\n")
                fout.write("for deltay in `seq %.3f %.3f %.3f`\n" % (yrange[0], yrange[2], yrange[1]))
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
                fout.write("for deltax in `seq %.3f %.3f %.3f`\n" % (xrange[0], xrange[2], xrange[1]))
                fout.write("do\n")
                fout.write("for deltay in `seq %.3f %.3f %.3f`\n" % (yrange[0], yrange[2], yrange[1]))
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
                fout.write("for deltax in `seq %.3f %.3f %.3f`\n" % (xrange[0], xrange[2], xrange[1]))
                fout.write("do\n")
                fout.write("for deltay in `seq %.3f %.3f %.3f`\n" % (yrange[0], yrange[2], yrange[1]))
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
                fout.write("EOF\n")
                fout.write("cd post-processing; gnuplot plot.gnuplot; cd ../\n")

            # write local bash run script
            with open("pes-relax.sh", 'w') as fout:
                fout.write("#!/bin/bash\n")
                fout.write("\n")
                fout.write("for deltax in `seq %.3f %.3f %.3f`\n" % (xrange[0], xrange[2], xrange[1]))
                fout.write("do\n")
                fout.write("for deltay in `seq %.3f %.3f %.3f`\n" % (yrange[0], yrange[2], yrange[1]))
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
                fout.write("for deltax in `seq %.3f %.3f %.3f`\n" % (xrange[0], xrange[2], xrange[1]))
                fout.write("do\n")
                fout.write("for deltay in `seq %.3f %.3f %.3f`\n" % (yrange[0], yrange[2], yrange[1]))
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
                fout.write("for deltax in `seq %.3f %.3f %.3f`\n" % (xrange[0], xrange[2], xrange[1]))
                fout.write("do\n")
                fout.write("for deltay in `seq %.3f %.3f %.3f`\n" % (yrange[0], yrange[2], yrange[1]))
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
        calculate the surface potential energy surface via CP2K.
    """
    def __init__(self):
        """
        """
        super().__init__()

        self.pes_params = {}
        self.set_pes() # set default value

    def set_pes(self, move_atom=[-1], xrange=[0, 1.5, 0.1], yrange=[0, 1.5, 0.5], zshift=0.0):
        """
        :parma move_atom: the atoms that will move in the calculation, list start from 0.
        :param xrange: x range for moving the specified moving atoms.
        :param: yrange: y range for moving the specified moving atoms
        :param: zshift: z shift for the moving atoms, will shift the z of specified moving atoms by value of zshift
        """
        self.pes_params["move_atom"] = move_atom
        self.pes_params["xrange"] = xrange
        self.pes_params["yrange"] = yrange
        self.pes_params["zshift"] = zshift

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
            os.chdir(directory)
            # generate the input files and the initial trajectory
            os.system("mkdir -p post-processing")
            for deltax in np.arange(xrange[0], xrange[1], xrange[2]):
                for deltay in np.arange(yrange[0], yrange[1], yrange[2]):
                    # to avoid float -0.000 be translated to string -0.000 we use 0.0 when value ==0 whether it is 0.0 or -0.0
                    os.mkdir("_%.3f_%.3f_" % (deltax if np.abs(deltax) >= 0.001 else 0.0, deltay if np.abs(deltay) >= 0.001 else 0.0))

                    for i in self.pes_params["move_atom"]:
                        self.arts.xyz.atoms[i].x += deltax
                        self.arts.xyz.atoms[i].y += deltay
                        # shift z of the specified atoms by self.pes_params["zshfit"]
                        #----------------------------------------------
                        self.arts.xyz.atoms[i].z += zshift

                    # fix no move atoms and xy of move atoms
                    # first fix all atoms
                    for i in range(len(self.arts.xyz.atoms)):
                        self.arts.xyz.atoms[i].fix = [True, True, True]
                    # unfix z of moving atoms
                    for i in self.pes_params["move_atom"]:
                        self.arts.xyz.atoms[i].fix[2] = False

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
                fout.write("for deltax in `seq %.3f %.3f %.3f`\n" % (xrange[0], xrange[2], xrange[1]))
                fout.write("do\n")
                fout.write("for deltay in `seq %.3f %.3f %.3f`\n" % (yrange[0], yrange[2], yrange[1]))
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
                fout.write("for deltax in `seq %.3f %.3f %.3f`\n" % (xrange[0], xrange[2], xrange[1]))
                fout.write("do\n")
                fout.write("for deltay in `seq %.3f %.3f %.3f`\n" % (yrange[0], yrange[2], yrange[1]))
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
                fout.write("for deltax in `seq %.3f %.3f %.3f`\n" % (xrange[0], xrange[2], xrange[1]))
                fout.write("do\n")
                fout.write("for deltay in `seq %.3f %.3f %.3f`\n" % (yrange[0], yrange[2], yrange[1]))
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
                fout.write("EOF\n")
                fout.write("cd post-processing; gnuplot plot.gnuplot; cd ../\n")

            # write local bash run script
            with open("pes-relax.sh", 'w') as fout:
                fout.write("#!/bin/bash\n")
                fout.write("#\n")
                fout.write("for deltax in `seq %.3f %.3f %.3f`\n" % (xrange[0], xrange[2], xrange[1]))
                fout.write("do\n")
                fout.write("for deltay in `seq %.3f %.3f %.3f`\n" % (yrange[0], yrange[2], yrange[1]))
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
                fout.write("for deltax in `seq %.3f %.3f %.3f`\n" % (xrange[0], xrange[2], xrange[1]))
                fout.write("do\n")
                fout.write("for deltay in `seq %.3f %.3f %.3f`\n" % (yrange[0], yrange[2], yrange[1]))
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
                fout.write("for deltax in `seq %.3f %.3f %.3f`\n" % (xrange[0], xrange[2], xrange[1]))
                fout.write("do\n")
                fout.write("for deltay in `seq %.3f %.3f %.3f`\n" % (yrange[0], yrange[2], yrange[1]))
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
                fout.write("EOF\n")
                fout.write("cd post-processing; gnuplot plot.gnuplot; cd ../\n")

            os.chdir("../")
        if runopt == "genrun" or runopt == "run":
            os.chdir(directory)
            os.system("bash pes-relax.sh")
            oschdir("../")
        # server handle
        server_handle(auto=auto, directory=directory, jobfilebase="pes-relax", server=self.run_params["server"])

