
def dftbPlusSubparser(subparsers):
    # --------------------------------------------------------------------------
    # DFTB+
    # --------------------------------------------------------------------------
    subparser = subparsers.add_parser("dftb+", help="using DFTB+ as calculator")

    #                      run params
    # -----------------------------------------------------------------
    gp = subparser.add_argument_group(title="overall running control", description="control the task running parameters")

    gp.add_argument("-d", "--directory", type=str, default="matflow-running",
        help="directory to generate all the files, do not specify the current directory")

    gp.add_argument("-r", "--runtype", type=int, default=0,
        choices=[0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10],
        help="choices of runtype. 0->static_run; 1->optimization; 2->cubic-cell; 3->hexagonal-cell; 4->tetragonal-cell; 5->neb; 6->vasp-phonon; 7->phonopy; 8->surf pes; 9->abc; 10->AIMD")

    # run option
    gp.add_argument("--runopt", type=str, default="gen",
        choices=["gen", "run", "genrun"],
        help="Generate or run or both at the same time.")

    gp.add_argument("--auto", type=int, default=3,
        choices=[0, 1, 2, 3],
        help="auto:0 nothing, 1: copying files to server, 2: copying and executing, 3: pymatflow run inserver with direct submit,  in order use auto=1, 2, you must make sure there is a working ~/.pymatflow/server_[pbs|llhpc].conf")

    gp.add_argument("--mpi", type=str, default="",
        help="MPI command, used in single node running, namely --auto 0 --runopt genrun")

    gp.add_argument("--server", type=str, default="pbs",
        choices=["pbs", "llhpc", "lsf_sz", "tianhe2", "lsf_sustc"],
        help="type of remote server, can be pbs or llhpc or lsf_sz or lsf_sustc")

    gp.add_argument("--jobname", type=str, default="matflow-running",
        help="jobname on the pbs server")

    gp.add_argument("--nodes", type=int, default=1,
        help="Nodes used in server")

    gp.add_argument("--ppn", type=int, default=32,
        help="ppn of the server")
    
    gp.add_argument("--queue", type=str, default=None,
        help="the queue to submit to job, default is not set")    
    
    # llhpc
    gp.add_argument("--partition", type=str, default="free",
        help="choose partition to submit job, now only apply for llhpc")

    gp.add_argument("--ntask", type=int, default=24,
        help="choose task number, now only apply for llhpc")

    gp.add_argument("--stdout", type=str, default="slurm.out",
        help="set standard out, now only apply for llhpc")

    gp.add_argument("--stderr", type=str, default="slurm.err",
        help="set standard err, now only apply for llhpc")

    # actually this can be put in the main subparser, but it will make the command not like git sub-cmmand
    # so we put them in every subsubparser
    gp = subparser.add_mutually_exclusive_group(required=True) # at leaset one of cif and xyz is provided
    # argparse will make sure only one of argument in structfile(xyz, cif) appear on command line
    gp.add_argument("--xyz", type=str, default=None,
        help="The xyz structure file with the second line specifying the cell parameter")

    gp.add_argument("--cif", type=str, default=None,
        help="The cif structure file")

    gp.add_argument("--xsd", type=str, default=None,
        help="The xsd structure file")

    gp.add_argument("--xsf", type=str, default=None,
        help="The xsf structure file")

    gp.add_argument("--images", type=str, nargs="+",
        help="the image stucture file(--images first.cif final.xsd), can only be cif, xsd, xsd, or xyz(second line is cell parameter) format")

    # Slater Koster file
    gp = subparser.add_argument_group(title="Slater Koster files")

    gp.add_argument("--slako-dir", type=str, default="./",
        help="specify the prefix path of the slako filess, default is ./.")
        
    # --------------------------------------------------------
    #                   DFTB+ PARAMETERS
    # --------------------------------------------------------

    # Hamiltonian

    gp = subparser.add_argument_group(title="Hamiltonian")

    gp.add_argument("--maxangularmomentum", type=str, default=None, required=True,
        help="Specifies the highest angular momentum for each atom type.")


    gp.add_argument("--kpoints-mp", type=int, nargs=6,
        default=[1, 1, 1, 0, 0, 0],
        help="set kpoints like -k 1 1 1 0 0 0")

    #gp.add_argument("--kpoints-mp-scf", type=int, nargs="+",
    #    default=[1, 1, 1, 0, 0, 0],
    #    help="set kpoints like -k 1 1 1 0 0 0")

    gp.add_argument("--kpoints-mp-nscf", type=int, nargs="+",
    default=None,
    help="set kpoints like -k 1 1 1 0 0 0")

    gp.add_argument("--kpath-manual", type=str, nargs="+", default=None,
        help="set kpoints for band structure calculation manually")

    gp.add_argument("--kpath-file", type=str, default="kpath.txt",
        help="set kpoints for band structure calculation manually from file")

    # Driver
    gp = subparser.add_argument_group(title="Driver")

    gp.add_argument("--steps", type=int, default=None,
        help="Number of MD steps to perform")

    gp.add_argument("--timestep", type=float, default=None,
        help="Time interval between two MD steps")

    #                     neb related PARAMETERS
    # --------------------------------------------------------------------------


    # PHONOPY parameters
    # ----------------------------------------
    gp = subparser.add_argument_group(title="phonopy parameters",
        description="setting of parameters needed by phonopy")

    gp.add_argument("--supercell-n", type=int, nargs="+",
        default=[1, 1, 1],
        help="supercell for phonopy, like [2, 2, 2]")


    # range_a range_b range_c
    # ----------------------------------------------
    gp = subparser.add_argument_group(title="cell optimization",
        description="cubic, hexagonal, tetragonal cell or general abc optimization parameters")

    gp.add_argument("--range-a", type=float, nargs=3,
        default=[-0.1, 0.1, 0.01],
        help="test range for a")

    gp.add_argument("--range-b", type=float, nargs=3,
        default=[-0.1, 0.1, 0.01], 
        help="test range for b")
            
    gp.add_argument("--range-c", type=float, nargs=3,
        default=[-0.1, 0.1, 0.01],
        help="test range for c")

    gp.add_argument("--batch-a", type=int,
        default=None,
        help="number of structure each batch a")
    
    gp.add_argument("--batch-b", type=int, 
        default=None,
        help="number of structure each batch b")
            
    gp.add_argument("--batch-c", type=int,
        default=None,
        help="number of structure each batch c")


    # surf pes
    gp = subparser.add_argument_group(title="surf pes",
        description="surf pes")

    gp.add_argument("--move-atom", type=int, nargs="+",
        default=[-1],
        help="specify the atoms to move, index starts from 0")

    gp.add_argument("--xrange", type=float, nargs="+",
        default=[1, 3, 0.5],
        help="x range for moving the specified moving atoms.")

    gp.add_argument("--yrange", type=float, nargs="+",
        default=[3, 5, 0.5],
        help="y range for moving the specified moving atoms.")

    gp.add_argument("--zshift", type=float,
        default=0.0,
        help="z shift for the moving atoms, will shift the z of specified moving atoms by value of zshift")

    gp.add_argument("--fix-z", type=int, default=1,
        choices=[0, 1, 2],
        help="0 -> do not fix any z of the atoms, 1: only fix z of the buttom atoms, 2: fix z of both the buttom and the moving atoms. default is 1")
    
    gp.add_argument("--fix-y", type=int, default=2,
        choices=[0, 1, 2],
        help="0 -> do not fix any z of the atoms, 1: only fix z of the buttom atoms, 2: fix z of both the buttom and the moving atoms. default is 2")
    
    gp.add_argument("--fix-x", type=int, default=2,
        choices=[0, 1, 2],
        help="0 -> do not fix any z of the atoms, 1: only fix z of the buttom atoms, 2: fix z of both the buttom and the moving atoms. default is 2")

    gp.add_argument("--batch-x-y", type=int, nargs=2, default=None,
        help="number of structures to calculate each batch x and y, default is all in one batch")

    # fix atoms
    gp = subparser.add_argument_group(title="fix atoms",
        description="specify atoms to fix in optimization, only used when --runtype=1")

    gp.add_argument("--fix", help="list of fixed atoms, index start from 1, have privilege over --fix-around-z", nargs='+', type=int)

    gp.add_argument("--fix-around-z", type=float, nargs=3, default=None,
        help="select atoms around specified z in Angstrom with tolerance, like this --fix-around-z 10 -0.5 0.5")            
            
    gp.add_argument("--color-fix", type=str, default="white",
        choices=["red", "green", "blue", "white"],
        help="select color to color the fixed atoms in xsd file, can be: red green blue and white")


def dftbPlusDriver(args):
    #params = {}
    #params["xxx"] = "xxx"

    mam = {}
    for item in args.maxangularmomentum.split("|"):
        mam[item.split(":")[0]] = item.split(":")[1]
    
    if args.runtype == 0:
        # static
        pass
                    
    elif args.runtype == 1:
        # optimization
        from pymatflow.dftbplus.opt import OptRun
        #  
        #            
        task = OptRun()
        task.get_xyz(xyzfile)

        kp = []
        kp.append(args.kpoints_mp[0])
        kp.append(0)
        kp.append(0)
        kp.append(0)
        kp.append(args.kpoints_mp[1])
        kp.append(0)
        kp.append(0)
        kp.append(0)
        kp.append(args.kpoints_mp[2])
        kp.append(args.kpoints_mp[3])
        kp.append(args.kpoints_mp[4])
        kp.append(args.kpoints_mp[5])
        task.hamiltonian.method["KPointsAndWeights"].list_of_scalar[""] = kp

        task.hamiltonian.method["SlaterKosterFiles"].scalar["Prefix"] = args.slako_dir
        for item in mam:
            task.hamiltonian.list_of_property["MaxAngularMomentum"].scalar[item] = mam[item]
        
        task.driver.val = "ConjugateGradient"  # "LBFGS"
        task.driver.scalar["ConvergentForcesOnly"] = "No"
        task.driver.scalar["AppendGeometries"] = "Yes"

        task.set_run(mpi=args.mpi, server=server, jobname=args.jobname, nodes=args.nodes, ppn=args.ppn, queue=args.queue) 
        task.set_llhpc(partition=args.partition, nodes=args.nodes, ntask=args.ntask, jobname=args.jobname, stdout=args.stdout, stderr=args.stderr)
        task.optimize(directory=args.directory, runopt=args.runopt, auto=args.auto)
    elif args.runtype == 2:
        # cubic
        pass
    elif args.runtype == 3:
        # hexagonal
        pass
    elif args.runtype == 4:
        # tetragonal
        pass
    elif args.runtype == 5:
        from pymatflow.dftbplus.md import MdRun
        #  
        #            
        task = MdRun()
        task.get_xyz(xyzfile)

        kp = []
        kp.append(args.kpoints_mp[0])
        kp.append(0)
        kp.append(0)
        kp.append(0)
        kp.append(args.kpoints_mp[1])
        kp.append(0)
        kp.append(0)
        kp.append(0)
        kp.append(args.kpoints_mp[2])
        kp.append(args.kpoints_mp[3])
        kp.append(args.kpoints_mp[4])
        kp.append(args.kpoints_mp[5])
        task.hamiltonian.method["KPointsAndWeights"].list_of_scalar[""] = kp

        task.hamiltonian.method["SlaterKosterFiles"].scalar["Prefix"] = args.slako_dir
        for item in mam:
            task.hamiltonian.list_of_property["MaxAngularMomentum"].scalar[item] = mam[item]
        
        task.driver.val = "VelocityVerlet"
        task.driver.scalar["MDRestartFrequency"] = 1
        task.driver.scalar["Steps"] = args.steps
        task.driver.scalar["TimeStep"] = args.timestep
        task.driver.scalar["ConvergentForcesOnly"] = "No"

        task.driver.method["Thermostat"].status = True
        task.driver.method["Thermostat"].val = "NoseHoover"
        task.driver.method["Thermostat"].scalar["Temperature"] = 400
        task.driver.method["Thermostat"].scalar["CouplingStrength"] = 3200

        task.set_run(mpi=args.mpi, server=server, jobname=args.jobname, nodes=args.nodes, ppn=args.ppn, queue=args.queue) 
        task.set_llhpc(partition=args.partition, nodes=args.nodes, ntask=args.ntask, jobname=args.jobname, stdout=args.stdout, stderr=args.stderr)
        task.md(directory=args.directory, runopt=args.runopt, auto=args.auto)
