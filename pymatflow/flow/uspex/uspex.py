"""
providing interface to USPEX software
"""
import os
import sys
import shutil

from pymatflow.remote.server import server_handle

from pymatflow.vasp.opt import opt_run as vasp_opt_run

class uspex:
    """
    """
    def __init__(self):
        self._initialize()
        self.vasp = vasp_opt_run()

    def _initialize(self):
        self.params = {
            "calculationMethod": None, "calculationType": None, "optType": None, "atomType": None,
            "numSpecies": None, "magRatio": None, "ldaU": None, "ExternalPressure": None, "valences": None,
            "goodBonds": None, "checkMolecules": None, "checkConnectivity": None, "fitLimit": None, "populationSize": None, "initialPopSize": None,
            "numGenerations": None, "stopCrit": None, "bestFrac": None, "keepBestHM": None, "reoptOld": None, "symmetries": None,
            "splitInto": None, "fracGene": None, "fracRand": None, "fracTopRand": None, "fracPerm": None,
            "fracAtomsMut": None, "fracRotMut": None, "fracLatMut": None, "fracSpinMut": None,
            "howManySwaps": None, "specificSwaps": None, "AutoFrac": None, "minVectorLength": None, "IonDistances": None,
            "constraint enhancement": None, "MolCenters": None, "Latticevalues": None, "abinitioCode": None,
            "KresolStart": None, "vacuumSize": None, "numParallelCalcs": None, "commandExecutable": None, "whichCluster": None,
            "remoteFolder": None, "PhaseDiagram": None, "RmaxFing": None, "deltaFing": None,
            "sigmaFing": None, "antiSeedsMax": None, "antiSeedsSigma": None, "doSpaceGroup": None, "SymTolerance": None,
            "repeatForStatistics": None, "stopFitness": None, "fixRndSeed": None, "collectForces": None,
            "mutationRate": None, "mutationDegree": None, "ordering active": None, "symmetrize": None, "valenceElectr": None,
            "percSliceShift": None, "maxDistHeredity": None, "manyParents": None, "minSlice": None, "maxSlice": None,
            "numberparents": None, "firstGeneMax": None, "minAt", None, "maxAt": None, "fracTrans": None, "howManyTrans": None,
            "specificTrans": None, "BoltzTraP_T_max": None, "BoltzTraP_T_delta": None, "BoltzTraP_T_efcut": None, "TE_T_interest": None,
            "TE_threshold": None, "TE_goal": None, "thicknessS": None, "thicknessB": None, "reconstruct": None, "GaussianWidth": None,
            "GaussianHeight": None, "FullRelax": None, "maxVectorLength": None, "minVectorLength": None, "PSO_softMut": None, "PSO_BestStruc": None,
            "PSO_BestEver": None, "vcnebType": None, "numImages": None, "numSteps": None, "optReadImages": None, "optimizerType": None, "optRelaxType": None, 
            "dt": None, "ConvThreshold": None, "VarPathLength": None, "K min": None, "K max": None, "Kconstant": None, "optFreezing": None, "optMethodCIDI": None, 
            "startCIDIStep": None, "pickupImages": None, "FormatType": None, "PrintStep": None, "numIterations": None, "speciesSymbol": None, "mass": None, 
            "amplitudeShoot": None, "magnitudeShoot": None, "shiftRatio": None, "orderParaType": None, "opCriteria": None, "cmdOrderParameter": None, 
            "cmdEnthalpyTemperature": None, "orderParameterFile": None, "enthalpyTemperatureFile": None, "trajectoryFile": None, "MDrestartFile": None,
        }
        
        self.run_params = {}
        self.set_run()

    def run_vasp(self, directory="tmp-uspex-vasp", runopt="gen", auto=0):
        """
        directory: a place for all the generated files
        """
        
        self.params["numSpecies"] = self.vasp.poscar.xyz.nspecies
        
        if runopt == "gen" or runopt == "genrun":
            if os.path.exists(directory):
                shutil.rmtree(directory)
            os.mkdir(directory)
            shutil.copyfile("POTCAR", os.path.join(directory, "POTCAR"))

            os.system("cp %s %s/" % (self.vasp.poscar.xyz.file, directory))

            #os.system("sflow convert -i %s -o %s" % (self.vasp.poscar.xyz.file, os.path.join(directory, "input_structure.cif")))
            self.set_params({"CifFilePath": "input_structure.cif"})
            
                
            with open(os.path.join(directory, "POSCAR"), 'w') as fout:
                self.vasp.poscar.to_poscar(fout)
            
            #with open(os.path.join(directory, "input.dat"), 'w') as fout:
            #    self.to_input_dat(fout)
            
            with open(os.path.join(directory, "uspex.pbs"), 'w') as fout:
                fout.write("#!/bin/bash\n")
                fout.write("#PBS -N %s\n" % self.run_params["jobname"])
                fout.write("#PBS -l nodes=%d:ppn=%d\n" % (self.run_params["nodes"], self.run_params["ppn"]))
                if self.run_params["queue"] != None:
                    fout.write("#PBS -q %s\n" % self.run_params["queue"])            
                fout.write("\n")
                fout.write("cd $PBS_O_WORKDIR\n")
                fout.write("cat > INCAR<<EOF\n")
                self.vasp.incar.to_incar(fout)
                fout.write("EOF\n")
                fout.write("cat > KPOINTS<<EOF\n")
                self.vasp.kpoints.to_kpoints(fout)
                fout.write("EOF\n")
                fout.write("cat > INPUT.txt<<EOF\n")
                self.to_input_txt(fout)
                fout.write("EOF\n")
                fout.write("NP=`cat $PBS_NODEFILE | wc -l`\n")
                fout.write("cat $PBS_NODEFILE > machinefile\n")
                
            
                
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
    
    def to_input_txt(self, fout):
        pass