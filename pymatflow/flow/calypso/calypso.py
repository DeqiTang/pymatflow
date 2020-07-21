"""
providing interface to calypso software
"""
import os
import sys
import shutil
import numpy as np

from pymatflow.remote.server import server_handle

from pymatflow.vasp.opt import opt_run as vasp_opt_run

class calypso:
    """
    """
    def __init__(self):
        self._initialize()
        self.vasp = vasp_opt_run()

    def _initialize(self):
        self.params = {
            "SystemName": None, "NumberOfSpecies": None, "NameOfAtoms": None, "NumberOfFormula": None,
            "Volume": None, "Ialgo": None, "ICode": None, "NumberOfLocalOptim": None, "PsoRatio": None,
            "PopSize": None, "Kgrid": None, "Command": None, "MaxStep": None, "PickUp": None, "PickStep": None,
            "MaxTime": None, "LMC": None, "2D": None, "Areal": None, "DeltaZ": None, "MultiLayer": None,
            "LayerGap": None, "VacuumGap": None, "LAtom_Dis": None, "Cluster": None, "Vacancy": None,
            "MOL": None, "NumberOfTypeMolecule": None, "NumberOfMolecule": None, "DistOfMol": None,
            "SpeSpaceGroup": None, "FixCell": None, "FixAtom": None, "VSC": None, "MaxNumAtom": None,
            "LSurface": None, "CalSubstrate": None, "SurfaceThickness": None, "SubstrateThickness": None,
            "Substrate": None, "ECR": None, "CifFilePath": None, "MillerIndex": None, "SlabVacuumThick": None,
            "SlabDepth": None, "SlabNumLayers": None, "NumRelaxedLayers": None, "CapBondsWithH": None,
            "Hardness": None, "Band_edge": None, "TarBandGap": None, "Adsorption": None, "AdsorptionStyle": None,
            "NumberOfTypeAtom": None, "SuperCell": None, "RangeOfBondLenth": None, "AdsorptionSymmetry": None,
            "TypeOfSubstrate": None, "BothSide": None,
        }
        
        self.run_params = {}
        self.set_run()

    def run_vasp(self, directory="tmp-calypso-vasp", runopt="gen", auto=0):
        """
        directory: a place for all the generated files
        """
        
        self.params["Command"] = "sh submit.sh" 
        
        if runopt == "gen" or runopt == "genrun":
            if os.path.exists(directory):
                shutil.rmtree(directory)
            os.mkdir(directory)
            shutil.copyfile("POTCAR", os.path.join(directory, "POTCAR"))

            #os.system("cp %s %s/" % (self.vasp.poscar.xyz.file, directory))

            #os.system("sflow convert -i %s -o %s" % (self.vasp.poscar.xyz.file, os.path.join(directory, "input_structure.cif")))
            #self.set_params({"CifFilePath": "input_structure.cif"})
            
            # Parallel and NumberOfParallel are controlled internally here
            if self.run_params["nodes"] > 1:
                self.params["Parallel"] = "T"
                self.params["NumberOfParallel"] = self.run_params["nodes"]
            else:
                self.params["Parallel"] = "F"
                self.params["NumberOfParallel"] = None
            with open(os.path.join(directory, "submit.sh"), 'w') as fout:
                fout.write("#!/bin/bash\n")
                #fout.write("mpirun -machinefile snodefile -np %d $PMF_VASP_STD > log 2>/dev/null\n" % (self.run_params["ppn"]))
                fout.write("mpiexec -machinefile snodefile -n %d $PMF_VASP_STD > log 2>/dev/null\n" % (self.run_params["ppn"]))
        
            #for i in range(slef.gen_incar_n):
            #    self.vasp.set_params(self.multi_incar_params[i], runtype="opt")
            #    with open(os.path.join(directory, "INCAR_%d" % (i+1)), 'w') as fout:
            #        self.vasp.incar.to_incar(fout)
                
            #with open(os.path.join(directory, "POSCAR"), 'w') as fout:
            #    self.vasp.poscar.to_poscar(fout)
            
            #with open(os.path.join(directory, "input.dat"), 'w') as fout:
            #    self.to_input_dat(fout)
            
            with open(os.path.join(directory, "calypso.pbs"), 'w') as fout:
                fout.write("#!/bin/bash\n")
                fout.write("#PBS -N %s\n" % self.run_params["jobname"])
                fout.write("#PBS -l nodes=%d:ppn=%d\n" % (self.run_params["nodes"], self.run_params["ppn"]))
                if self.run_params["queue"] != None:
                    fout.write("#PBS -q %s\n" % self.run_params["queue"])            
                fout.write("\n")
                fout.write("cd $PBS_O_WORKDIR\n")
                for i in range(self.gen_incar_n):
                    fout.write("cat > 'INCAR_%d'<<EOF\n" % (i+1))
                    self.vasp.set_params(self.multi_incar_params[i], runtype="opt")
                    self.vasp.incar.to_incar(fout)
                    fout.write("EOF\n")
                fout.write("cat > KPOINTS<<EOF\n")
                self.vasp.kpoints.to_kpoints(fout)
                fout.write("EOF\n")
                fout.write("cat > input.dat<<EOF\n")
                self.to_input_dat(fout)
                fout.write("EOF\n")
                fout.write("NP=`cat $PBS_NODEFILE | wc -l`\n")
                fout.write("cat $PBS_NODEFILE > machinefile\n")
                fout.write("$PMF_CALYPSOX > pso.log\n")
            
                
        if runopt == "run" or runopt == "genrun":
            # run the neb calculation
            # each image on one core
            os.chdir(directory)
            #os.system("%s vasp" % mpi)
            os.system("bash calypso.sh")
            os.chdir("../")
        server_handle(auto=auto, directory=directory, jobfilebase="calypso", server=self.run_params["server"])

    #
    
    def set_run(self, mpi="", server="pbs", jobname="calypso", nodes=1, ppn=32, queue=None):
        """ used to set  the parameters controlling the running of the task
        :param mpi: you can specify the mpi command here, it only has effect on native running
        """
        self.run_params["server"] = server
        self.run_params["mpi"] = mpi
        self.run_params["jobname"] = jobname
        self.run_params["nodes"] = nodes
        self.run_params["ppn"] = ppn
        self.run_params["queue"] = queue

    def set_llhpc(self, partition="free", nodes=1, ntask=24, jobname="matflow_job", stdout="slurm.out", stderr="slurm.err"):
        self.run_params["partition"] = partition
        self.run_params["jobname"] = jobname
        self.run_params["nodes"] = nodes
        self.run_params["ntask"] = ntask
        self.run_params["stdout"] = stdout
        self.run_params["stderr"] = stderr
        
    def set_params(self, params={}):
        for item in params:
            self.params[item] = params[item]
    
    def to_input_dat(self, fout):
        for item in self.params:
            if self.params[item] == None:
                continue
            if type(self.params[item]) != list:
                fout.write("%s = %s\n" % (item, self.params[item]))
                fout.write("\n")
            elif type(self.params[item][0]) != list:
                # self.params[item] is vec
                fout.write("%s =" % item)
                for value in self.params[item]:
                    fout.write(" %s" % value)
                fout.write("\n")
                fout.write("\n")
            else:
                # self.params[item][0] is n x n matrix
                fout.write("@%s\n" % (item))
                for row in self.params[item]:
                    fout.write(" ")
                    for value in row:
                        fout.write(" %s" % (value))
                    fout.write("\n")
                fout.write("@End\n")
                
                
    def run_vasp_split(self, directory="tmp-calypso-vasp", runopt="gen", auto=0, split_batch=None):
        """
        directory: a place for all the generated files
        """
        
        self.params["Split"] = "T"
        
        pop_size = int(self.params["PopSize"]) if "PopSize" in self.params and self.params["PopSize"] != None else 30
        max_step = int(self.params["MaxStep"]) if "MaxStep" in self.params and self.params["MaxStep"] != None else 50
        if split_batch == None:
            split_batch = pop_size // 2
        
        self.params["Command"] = "sh submit.sh"
          
        
        if runopt == "gen" or runopt == "genrun":
            if os.path.exists(directory):
                shutil.rmtree(directory)
            os.mkdir(directory)
            shutil.copyfile("POTCAR", os.path.join(directory, "POTCAR"))

            #os.system("cp %s %s/" % (self.vasp.poscar.xyz.file, directory))

            #os.system("sflow convert -i %s -o %s" % (self.vasp.poscar.xyz.file, os.path.join(directory, "input_structure.cif")))
            #self.set_params({"CifFilePath": "input_structure.cif"})
            
            # Parallel and NumberOfParallel are controlled internally here
            if self.run_params["nodes"] > 1:
                self.params["Parallel"] = "T"
                self.params["NumberOfParallel"] = self.run_params["nodes"]
            else:
                self.params["Parallel"] = "F"
                self.params["NumberOfParallel"] = None
            with open(os.path.join(directory, "submit.sh"), 'w') as fout:
                fout.write("#!/bin/bash\n")
                #fout.write("mpirun -machinefile snodefile -np %d $PMF_VASP_STD > log 2>/dev/null\n" % (self.run_params["ppn"]))
                fout.write("mpiexec -machinefile snodefile -n %d $PMF_VASP_STD > log 2>/dev/null\n" % (self.run_params["ppn"]))
        
            
            n_batch = int(np.ceil(pop_size / split_batch))
            for i_batch in range(n_batch):
                with open(os.path.join(directory, "batch-%d.pbs" % i_batch), 'w') as fout:
                    fout.write("#!/bin/bash\n")
                    fout.write("#PBS -N %s\n" % self.run_params["jobname"])
                    fout.write("#PBS -l nodes=%d:ppn=%d\n" % (self.run_params["nodes"], self.run_params["ppn"]))
                    if self.run_params["queue"] != None:
                        fout.write("#PBS -q %s\n" % self.run_params["queue"])            
                    fout.write("\n")
                    fout.write("cd $PBS_O_WORKDIR\n")
                    fout.write("NP=`cat $PBS_NODEFILE | wc -l`\n")                 
                    if i_batch == 0:
                        for i in range(self.gen_incar_n):
                            fout.write("cat > 'INCAR_%d'<<EOF\n" % (i+1))
                            self.vasp.set_params(self.multi_incar_params[i], runtype="opt")
                            self.vasp.incar.to_incar(fout)
                            fout.write("EOF\n")
                        fout.write("cat > KPOINTS<<EOF\n")
                        self.vasp.kpoints.to_kpoints(fout)
                        fout.write("EOF\n")
                        fout.write("cat > input.dat<<EOF\n")
                        self.to_input_dat(fout)
                        fout.write("EOF\n")
                        
                        fout.write("mkdir -p running\n")
                        fout.write("$PMF_CALYPSOX 2&> pso.1.log\n")
                        fout.write("maxstep=%s\n" % int(self.params["MaxStep"]))
                        fout.write("popsize=%s\n" % int(self.params["PopSize"]))
                        fout.write("istep=2\n")
                        fout.write("while [ $istep -lt ${maxstep} ]\n")
                        fout.write("do\n")
                        #fout.write("  #start=`echo \"(${istep}-2) * ${popsize} + 1\" | bc`\n")
                        #fout.write("  #end=`echo \"(${istep}-1) * ${popsize}\" | bc`\n")
                        #fout.write("  start=1\n")
                        #fout.write("  end=${popsize}\n")
                        fout.write("  start=%s\n" % int(i_batch * pop_size / n_batch + 1))
                        if ((i_batch+1) * pop_size / n_batch) > pop_size:
                            fout.write("  end=%s\n" % pop_size)
                        else:
                            fout.write("  end=%s\n" % int((i_batch+1)*pop_size/n_batch))
                        fout.write("  for i in `seq -f \"%g\" $start 1 $end`\n")
                        fout.write("  do\n")
                        fout.write("    if [ ! -f OUTCAR_${i} ]\n")
                        fout.write("    then\n")
                        fout.write("      running_dir=running/step_`echo \"${istep}-1\" | bc`/${i}\n")
                        fout.write("      mkdir -p ${running_dir}\n")
                        fout.write("      cp POSCAR_${i} ${running_dir}/POSCAR\n")
                        fout.write("      cp POTCAR ${running_dir}/\n")
                        fout.write("      cp INCAR_* ${running_dir}/\n")
                        fout.write("      cp KPOINTS_${i} ${running_dir}/KPOINTS\n")
                        fout.write("      cd ${running_dir}\n")
                        fout.write("      for incar in INCAR_*\n")
                        fout.write("      do\n")
                        fout.write("        cp ${incar} INCAR\n")
                        fout.write("        if [[ -f CONTCAR && `cat CONTCAR | wc -l` -gt 0 ]]\n")
                        fout.write("        then\n")
                        fout.write("          cp CONTCAR POSCAR\n")
                        fout.write("        fi\n")
                        fout.write("        mpirun -np $NP -machinefile $PBS_NODEFILE $PMF_VASP_STD\n")
                        fout.write("      done\n")
                        fout.write("      cp OUTCAR ../../../OUTCAR_${i}\n")
                        fout.write("      contcar=`cat CONTCAR | wc -l`\n")
                        fout.write("      if [ $contcar == 0 ]\n")
                        fout.write("      then\n")
                        fout.write("        cp POSCAR ../../../CONTCAR_${i}\n")
                        fout.write("      else\n")
                        fout.write("        cp CONTCAR ../../../CONTCAR_${i}\n")
                        fout.write("      fi\n")
                        fout.write("      cd ../../../\n")
                        fout.write("    fi\n")
                        fout.write("  done\n")
                        fout.write("  istep_done=T\n")
                        fout.write("  for j in `seq  -f \"%%g\" 1 1 %s`\n" % pop_size)
                        fout.write("  do\n")
                        fout.write("    if [ ! -f running/step_`echo \"${istep}-1\" | bc`/${j}/OUTCAR ]\n")
                        fout.write("    then\n")
                        fout.write("      istep_done=F\n")
                        fout.write("      break\n")
                        fout.write("    fi\n")
                        fout.write("    # copy CONTCAR and OUTCAR (which may be already copied)\n")
                        fout.write("    cp running/step_`echo \"${istep}-1\" | bc`/${j}/OUTCAR OUTCAR_${j}\n")
                        fout.write("    cp running/step_`echo \"${istep}-1\" | bc`/${j}/CONTCAR CONTCAR_${j}\n")
                        fout.write("  done\n")
                        fout.write("  if [ ${istep_done} == \"T\" ]\n")
                        fout.write("  then\n")
                        fout.write("    $PMF_CALYPSOX 2&> pso.${istep}.log\n")
                        fout.write("  fi\n")
                        fout.write("  istep=`cat step`\n")                 
                        fout.write("done\n")
                    else:
                        fout.write("mkdir -p running\n")
                        #fout.write("istep=2\n")
                        fout.write("while [ ! -e step ]\n")
                        fout.write("do\n")
                        fout.write("  sleep 10s\n")
                        fout.write("done\n")
                        fout.write("\n")
                        fout.write("maxstep=%s\n" % int(self.params["MaxStep"]))
                        fout.write("istep=`cat step`\n")    
                        #fout.write("  #start=`echo \"(${istep}-2) * ${popsize} + 1\" | bc`\n")
                        #fout.write("  #end=`echo \"(${istep}-1) * ${popsize}\" | bc`\n")
                        #fout.write("  start=1\n")
                        #fout.write("  end=${popsize}\n")
                        fout.write("while [ $istep -lt ${maxstep} ]\n")
                        fout.write("do\n")                        
                        
                        fout.write("  start=%s\n" % (i_batch * pop_size / n_batch + 1))
                        if ((i_batch+1) * pop_size / n_batch) > pop_size:
                            fout.write("  end=%s\n" % pop_size)
                        else:
                            fout.write("  end=%s\n" % ((i_batch+1)*pop_size/n_batch))
                        fout.write("  for i in `seq -f \"%g\" $start 1 $end`\n")
                        fout.write("  do\n")
                        fout.write("    if [ ! -f OUTCAR_${i} ]\n")
                        fout.write("    then\n")
                        fout.write("      running_dir=running/step_`echo \"${istep}-1\" | bc`/${i}\n")
                        fout.write("      mkdir -p ${running_dir}\n")
                        fout.write("      cp POSCAR_${i} ${running_dir}/POSCAR\n")
                        fout.write("      cp POTCAR ${running_dir}/\n")
                        fout.write("      cp INCAR_* ${running_dir}/\n")
                        fout.write("      cp KPOINTS_${i} ${running_dir}/KPOINTS\n")
                        fout.write("      cd ${running_dir}\n")
                        fout.write("      for incar in INCAR_*\n")
                        fout.write("      do\n")
                        fout.write("        cp ${incar} INCAR\n")
                        fout.write("        if [[ -f CONTCAR && `cat CONTCAR | wc -l` -gt 0 ]]\n")
                        fout.write("        then\n")
                        fout.write("          cp CONTCAR POSCAR\n")
                        fout.write("        fi\n")                        
                        fout.write("        mpirun -np $NP -machinefile $PBS_NODEFILE $PMF_VASP_STD\n")
                        fout.write("      done\n")
                        fout.write("      cp OUTCAR ../../../OUTCAR_${i}\n")
                        fout.write("      contcar=`cat CONTCAR | wc -l`\n")
                        fout.write("      if [ $contcar == 0 ]\n")
                        fout.write("      then\n")
                        fout.write("        cp POSCAR ../../../CONTCAR_${i}\n")
                        fout.write("      else\n")
                        fout.write("        cp CONTCAR ../../../CONTCAR_${i}\n")
                        fout.write("      fi\n")
                        fout.write("      cd ../../../\n")
                        fout.write("    fi\n")
                        fout.write("    istep=`cat step`\n")
                        fout.write("    if [ ! ${istep} -lt %s ]; then exit; fi\n" % max_step)
                        fout.write("  done\n")
                        fout.write("done\n")
                        
            
        if runopt == "run" or runopt == "genrun":
            # run the neb calculation
            # each image on one core
            os.chdir(directory)
            #os.system("%s vasp" % mpi
            #os.system("bash calypso.sh")
            os.chdir("../")
        for i_batch in range(n_batch):
            server_handle(auto=auto, directory=directory, jobfilebase="batch-%d" % i_batch, server=self.run_params["server"])
                