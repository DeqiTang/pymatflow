"""
providing interface to calypso software
"""
from pymatflow.vasp.opt import opt_run as vasp_opt_run

class calypso:
    """
    """
    def __init__(self):
        self._initialize()
        self.vasp = vasp_opt_run

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

        if runopt == "gen" or runopt == "genrun":
            if os.path.exists(directory):
                shutil.rmtree(directory)
            os.mkdir(directory)
            shutil.copyfile("POTCAR", os.path.join(directory, "POTCAR"))
            
            with open(os.path.join(directory, "input.dat")) as fout:
                self.to_input_dat(fout)
            
            with open(os.path.join(directory, "calypso.pbs"), 'w') as fout:
                fout.write("#!/bin/bash\n")
                fout.write("#PBS -N %s\n" % self.run_params["jobname"])
                fout.write("#PBS -l nodes=%d:ppn=%d\n" % (self.run_params["nodes"], self.run_params["ppn"]))
                if self.run_params["queue"] != None:
                    fout.write("#PBS -q %s\n" % self.run_params["queue"])            
                fout.write("\n")
                fout.write("cd $PBS_O_WORKDIR\n")
                fout.write("cat > INCAR<<EOF\n")
                self.incar.to_incar(fout)
                fout.write("EOF\n")
                fout.write("cat > KPOINTS<<EOF\n")
                self.kpoints.to_kpoints(fout)
                fout.write("EOF\n")
                fout.write("NP=`cat $PBS_NODEFILE | wc -l`\n")
                fout.write("cat $PBS_NODEFILE > machinefile\n")
                fout.write("./$PMF_CALYPSOX > pso.log\n")
            
            with open(os.path.join(directory, "submit.sh"), 'w') as fout:
                fout.write("#!/bin/bash\n")
                fout.write("mpiexec -n %d $PMF_VASP_STD > log 2>/dev/null\n" % (self.run_params["ppn"]))
            self.params["Command"] = "sh submit.sh"
                
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
            if type(self.params[item]) != list:
                fout.write("%s = %s\n" % (item, self.params))
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