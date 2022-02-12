import os

def abinitSubparser(subparsers):
    # --------------------------------------------------------------------------
    # Abinit
    # --------------------------------------------------------------------------
    subparser = subparsers.add_parser("abinit", help="using abinit as calculator")

    gp = subparser.add_argument_group(title="overall running control",
        description="control the overall running parameters")

    gp.add_argument("-r", "--runtype", type=int, default=0,
        choices=[0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10],
        help="choices of runtype. 0->static_run; 1->optimization; 2->cubic-opt; 3->hexagonal-opt; 4->tetragonal-opt; 5->dfpt-elastic-piezo-dielec; 6->dfpt-phonon; 7->phonopy; 8->abc; 9->Optic; 10->BSE")

    gp.add_argument("-d", "--directory", type=str, default="matflow-running",
        help="Directory to do the calculation")
    # run params
    # -----------------------------------------------------------------
    # run option
    gp.add_argument("--runopt", type=str, default="gen",
        choices=["gen", "run", "genrun"],
        help="Generate or run or both at the same time.")

    gp.add_argument("--auto", type=int, default=3,        
        choices=[0, 1, 2, 3],
        help="auto:0 nothing, 1: copying files to server, 2: copying and executing, 3: pymatflow run inserver with direct submit,  in order use auto=1, 2, you must make sure there is a working ~/.pymatflow/server_[pbs|llhpc].conf")

    gp.add_argument("--mpi", type=str, default="",
        help="MPI command: like 'mpirun -np 4'")

    gp.add_argument("--server", type=str, default="pbs",
        choices=["pbs", "llhpc", "tianhe2", "cdcloud"],
        help="type of remote server, can be pbs or llhpc or lsf_sz, or tianhe2, cdcloud")

    gp.add_argument("--jobname", type=str, default="matflow-job",
        help="jobname on the pbs server")

    gp.add_argument("--nodes", type=int, default=1,
        help="Nodes used in server")

    gp.add_argument("--ppn", type=int, default=32,
        help="ppn of the server")

    gp.add_argument("--queue", type=str, default=None,
        help="the queue to submit to job, default is not set")

    # llhpc
    gp.add_argument("--partition", type=str, default="free",
        help="choose partition to submit job")

    gp.add_argument("--ntask", type=int, default=24,
        help="choose task number")

    gp.add_argument("--stdout", type=str, default="slurm.out",
        help="set standard out")

    gp.add_argument("--stderr", type=str, default="slurm.err",
        help="set standard err")

    # structure file: either xyz or cif. they are exclusive
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

    # potential file
    gp = subparser.add_argument_group(title="pseudopotential")

    gp.add_argument("--pot", type=str, default="./",
        help="specify the path to the directory containing all the needed pseudopotential, default behavior is find them in the current directory automatically. if you pass 'auto' to it, matflow will get the pots automatically(need simple configuration, see manual)")

    # parallelism
    gp = subparser.add_argument_group(title="parallelism",
        description="set parameters for a parallel calculation with the ABINIT package.")
    
    gp.add_argument("--autoparal", type=int, default=None,
        choices=[0, 1, 2, 3, 4],
        help="AUTOmatisation of the PARALlelism")

    # electronic structure
    gp = subparser.add_argument_group(title="electronic structure:",
        description="setting of electronic structure related parameters")

    gp.add_argument("--chkprim", type=int, default=1,
        choices=[0, 1],
        help="check whether the input cell is primitive. if your cell is not primitive, set chkprim to 0. for more information, refer to https://docs.abinit.org/variables/gstate/#chkprim")

    gp.add_argument("--kpath-manual", type=str, nargs="+", default=None,
        help="manual input kpath for band structure calculation")

    gp.add_argument("--kpath-file", type=str,
        help="file to read the kpath for band structure calculation")

    gp.add_argument("--iscf", type=int, default=7,
        choices=[0, 1, 2, 3, 4, 5, 7, 12, 13, 14, 15, 17],
        help="set scf or nscf type. for more information, refer to https://docs.abinit.org/variables/basic/#iscf")

    gp.add_argument("--ecut", type=int, default=15,
        help="Kinetic energy cutoff for wave functions in unit of Hartree, default value: 15 Hartree. for more information, refer to https://docs.abinit.org/variables/basic/#ecut")

    gp.add_argument("--nstep", type=int, default=30,
        help="Number of (non-)self-consistent field STEPS. for more information, refer to https://docs.abinit.org/variables/basic/#nstep")

    gp.add_argument("--ixc", type=int, default=11,
        choices=[1, 2, 3 ,4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 26, 27, 28, 40, 41, 42],
        help="type of exchage-correlation functional. for more information, refer to https://docs.abinit.org/variables/basic/#ixc")

    gp.add_argument("--kptopt", type=int, default=1,
        choices=[1],
        help="Kpoints Generation scheme option: 0, 1, 2, 3, 4 or a negative value. for more information, refer to https://docs.abinit.org/variables/basic/#kptopt")

    gp.add_argument("--ngkpt", nargs=3, type=int,
        default=[1, 1, 1],
        help="number of grid points for kpoints generation. for more information, refer to https://docs.abinit.org/variables/basic/#ngkpt")

    #gp.add_argument("--kpoints-mp", type=int, nargs=6,
    #   default=[1, 1, 1, 0, 0, 0],
    #   help="monkhorst-pack type k mesh generation using ngkpt.")

    # electron occupation
    gp = subparser.add_argument_group(title="electronic structure: occupation",
        description="setting of occupation parameters in electronic structure calc")

    gp.add_argument("--occopt", type=int, default=3,
        choices=[0, 1, 2, 3, 4, 5, 6, 7],
        help="Controls how input parameters nband, occ, and wtk are handled. for more information, refer to https://docs.abinit.org/variables/basic/#occopt")

    gp.add_argument("--nband", type=int, nargs="+", default=None,
        help="Gives number of bands, occupied plus possibly unoccupied, for which wavefunctions are being computed along with eigenvalues. for more information, refer to https://docs.abinit.org/variables/basic/#nband")

    gp.add_argument("--occ", nargs="+", type=float, default=None,
        help="Gives occupation numbers for all bands in the problem. Needed if occopt == 0 or occopt == 2. Ignored otherwise. Also ignored when iscf = -2. refer to https://docs.abinit.org/variables/gstate/#occ")

    # magnetic related parameters
    gp = subparser.add_argument_group(title="magnetic:",
        description="magnetic related parameters")

    gp.add_argument("--nsppol", type=int, default=None,
        choices=[1, 2],
        help="Give the number of INDEPENDENT spin polarisations, for which there are non- related wavefunctions. Can take the values 1 or 2. for more information, refer to https://docs.abinit.org/variables/basic/#nsppol")

    # vdw related parameters
    gp = subparser.add_argument_group(title="vdw:",
        description="setting of van der valls exchange-correlation functional")

    gp.add_argument("--vdw-xc", type=int, default=None,
        choices=[0, 1, 2, 5, 6, 7, 10, 11, 14],
        help="Van Der Waals exchange-correlation functional. 0: no correction, 1: vdW-DF1, 2: vdW-DF2, 5: DFT-D2, 6: DFT-D3, 7: DFT-D3(BJ). for more information, refer to https://docs.abinit.org/variables/vdw/#vdw_xc")

    gp.add_argument("--vdw-tol", type=float,
        default=None,
        help="Van Der Waals tolerance, only work when vdw_xc == 5 or 6 or 7. to be included in the potential a pair of atom must have contribution to the energy larger than vdw_tol. default value is 1.0e-10. fore more information, refer to https://docs.abinit.org/variables/vdw/#vdw_tol")

    gp = subparser.add_argument_group(title="write parameters:",
        description="control writing behavior of abinit")

    gp.add_argument("--prtden", type=int ,default=1,
        choices=[0, 1],
        help="print the density. for more information, refer to https://docs.abinit.org/variables/files/#prtden")

    gp.add_argument("--prtdos", type=int, default=None,
        choices=[0, 1, 2, 3],
        help="can be 0, 1, 2, 3. for more information, refer to https://docs.abinit.org/variables/files/#prtdos")

    gp.add_argument("--properties", nargs="+", type=int,
        default=[],
        help="options for properties calculation")

    #                        ions moving related parameters
    # -----------------------------------------------------------
    gp = subparser.add_argument_group(title="ions:",
        description="setting of ions related parameters")

    gp.add_argument("--ionmov", type=int, default=3,
        choices=[2, 3, 4, 5],
        help="type of ionmov algorithm. fore more information, refer to https://docs.abinit.org/variables/rlx/#ionmov")

    gp.add_argument("--optcell", type=int,
        choices=[0, 1, 2, 3, 4, 5, 6, 7, 8, 9],
        default=0,
        help="whether to optimize the cell shape and dimension. fore more information, refer to https://docs.abinit.org/variables/rlx/#optcell")

    gp.add_argument("--chkdilatmx", type=int, default=None,
        choices=[0, 1],
        help="check dilatmx. fore more information, refer to https://docs.abinit.org/variables/rlx/#chkdilatmx")

    gp.add_argument("--dilatmx", type=float, default=None,
        help="lattice dilation maximal value. fore more information, refer to https://docs.abinit.org/variables/rlx/#dilatmx")

    gp.add_argument("--ecutsm", type=float, default=None,
        help="when optcell != 0, must specify encutsm larser than zero. for more information refer to https://docs.abinit.org/variables/rlx/#ecutsm")

    gp.add_argument("--ntime", type=int, default=None,
        help="Default value: 0 if ionmvov == 0, set to 1000 if ionvmov != 0 and imgmov != 0 and the variable is not specified., for more information refer to https://docs.abinit.org/variables/rlx/#ntime")

    # na stepa nc stepc
    # ------------------------------------------------
    gp = subparser.add_argument_group(title="cell optimization",
        description="setting of parameters needed by cubic, hexagonal, tetragonal cell parameters optimization needed by matflow")

    gp.add_argument("--na", type=int, default=10,
        help="number of a to run")
    
    gp.add_argument("--stepa", type=float, default=0.05,
        help="step of a in unit of Angstrom")
    
    gp.add_argument("--nc", type=int, default=10,
        help="number of c to run")
            
    gp.add_argument("--stepc", type=float, default=0.05,
        help="step of c in unit of Angstrom")

    gp = subparser.add_argument_group(title="phonopy:",
        description="setting of parameters needed by phonopy")

    gp.add_argument("--supercell-n", type=int, nargs="+",
        default=[1, 1, 1],
        help="supercell build for phonopy.")

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

            

def abinitDriver(args):
    from pymatflow.cmd.matflow import getXyzFile
    xyzfile, images = getXyzFile(args)
    # server
    # xxx.set_run can only deal with pbs, llhpc, lsf_sz server now 
    # however both guangzhou chaosuan llhpc are build on tianhe2, so they can use the same job system(yhbatch...)
    # we add tianhe2 option to args.server which cannot be handled by xxx.set_run. so we convert it to llhpc if tianhe2 is chosen
    server = args.server if args.server != "tianhe2" else "llhpc"
          
    params = {}
    kpoints = {}

    params["chkprim"] = args.chkprim
    params["ecut"] = args.ecut
    params["nstep"] = args.nstep
    params["ixc"] = args.ixc
    params["vdw_xc"] = args.vdw_xc
    params["vdw_tol"] = args.vdw_tol

    kpoints["kptopt"] = args.kptopt
    kpoints["ngkpt"] = args.ngkpt

    params["occopt"] = args.occopt
    params["nband"] = args.nband
    params["occ"] = args.occ

    params["autoparal"] = args.autoparal


    if args.runtype == 0:
        # static
        params["nsppol"] = args.nsppol
        params["prtden"] = args.prtden
        params["prtdos"] = args.prtdos

        params["ecutsm"] = args.ecutsm

        from pymatflow.abinit.static import StaticRun
        from pymatflow.cmd.matflow import get_kpath
        task = StaticRun()
        if get_kpath(args.kpath_manual, args.kpath_file) == None:
            print("================================================\n")
            print("Warning: matflow abinit\n")
            print("in abinit static runing you must provide kpath\n")
            sys.exit(1)
        task.dataset[3].electrons.kpoints.set_band(kpath=get_kpath(args.kpath_manual, args.kpath_file))
        task.get_xyz(xyzfile)
        task.set_params(params=params)
        task.set_kpoints(kpoints=kpoints, ndtset=0)
        task.set_kpoints(kpoints=kpoints, ndtset=1)
        task.set_kpoints(kpoints=kpoints, ndtset=2)
        task.set_run(mpi=args.mpi, server=server, jobname=args.jobname, nodes=args.nodes, ppn=args.ppn, queue=args.queue)
        task.set_llhpc(partition=args.partition, nodes=args.nodes, ntask=args.ntask, jobname=args.jobname, stdout=args.stdout, stderr=args.stderr)
        task.run(directory=args.directory, runopt=args.runopt, auto=args.auto)
    elif args.runtype == 1:
        # optimization
        params["optcell"] = args.optcell
        params["chkdilatmx"] = args.chkdilatmx
        params["dilatmx"] = args.dilatmx
        params["ionmov"] = args.ionmov
        params["ecutsm"] = args.ecutsm
        params["ntime"] = args.ntime
        from pymatflow.abinit.opt import OptRun
        task = OptRun()
        task.get_xyz(xyzfile)
        task.set_params(params=params)
        task.set_kpoints(kpoints=kpoints)
        task.set_run(mpi=args.mpi, server=server, jobname=args.jobname, nodes=args.nodes, ppn=args.ppn, queue=args.queue)
        task.set_llhpc(partition=args.partition, nodes=args.nodes, ntask=args.ntask, jobname=args.jobname, stdout=args.stdout, stderr=args.stderr)
        task.run(directory=args.directory, runopt=args.runopt, auto=args.auto)
    elif args.runtype == 2:
        ## cubic optimization
        params["optcell"] = 0 # must be 0
        params["ionmov"] = args.ionmov
        params["ntime"] = args.ntime

        params["ecutsm"] = args.ecutsm

        from pymatflow.abinit.opt import OptRun
        task = OptRun()
        task.get_xyz(xyzfile)
        task.set_params(params=params)
        task.set_kpoints(kpoints=kpoints)
        task.set_run(mpi=args.mpi, server=server, jobname=args.jobname, nodes=args.nodes, ppn=args.ppn, queue=args.queue)
        task.set_llhpc(partition=args.partition, nodes=args.nodes, ntask=args.ntask, jobname=args.jobname, stdout=args.stdout, stderr=args.stderr)
        task.cubic(directory=args.directory, runopt=args.runopt, auto=args.auto, na=args.na, stepa=args.stepa)
    elif args.runtype == 3:
        # hexagonal optimization
        params["optcell"] = 0 # must be 0
        params["ionmov"] = args.ionmov
        params["ntime"] = args.ntime

        params["ecutsm"] = args.ecutsm

        from pymatflow.abinit.opt import OptRun
        task = OptRun()
        task.get_xyz(xyzfile)
        task.set_params(params=params)
        task.set_kpoints(kpoints=kpoints)
        task.set_run(mpi=args.mpi, server=server, jobname=args.jobname, nodes=args.nodes, ppn=args.ppn, queue=args.queue)
        task.set_llhpc(partition=args.partition, nodes=args.nodes, ntask=args.ntask, jobname=args.jobname, stdout=args.stdout, stderr=args.stderr)
        task.hexagonal(directory=args.directory, runopt=args.runopt, auto=args.auto, na=args.na, stepa=args.stepa, nc=args.nc, stepc=args.stepc)
    elif args.runtype == 4:
        # tetragonal optimization
        params["optcell"] = 0 # must be 0
        params["ionmov"] = args.ionmov
        params["ntime"] = args.ntime
        
        params["ecutsm"] = args.ecutsm

        from pymatflow.abinit.opt import OptRun
        task = OptRun()
        task.get_xyz(xyzfile)
        task.set_params(params=params)
        task.set_kpoints(kpoints=kpoints)
        task.set_run(mpi=args.mpi, server=server, jobname=args.jobname, nodes=args.nodes, ppn=args.ppn, queue=args.queue)
        task.set_llhpc(partition=args.partition, nodes=args.nodes, ntask=args.ntask, jobname=args.jobname, stdout=args.stdout, stderr=args.stderr)
        task.tetragonal(directory=args.directory, runopt=args.runopt, auto=args.auto, na=args.na, stepa=args.stepa, nc=args.nc, stepc=args.stepc)
    elif args.runtype == 5:
        # dfpt-elastic-piezo-dielec
        params["ecutsm"] = args.ecutsm

        from pymatflow.abinit.dfpt import DfptElasticPiezoDielec
        task = DfptElasticPiezoDielec()
        task.get_xyz(xyzfile)
        task.set_params(params=params)
        task.set_kpoints(kpoints=kpoints)
        task.set_run(mpi=args.mpi, server=server, jobname=args.jobname, nodes=args.nodes, ppn=args.ppn, queue=args.queue)
        task.set_llhpc(partition=args.partition, nodes=args.nodes, ntask=args.ntask, jobname=args.jobname, stdout=args.stdout, stderr=args.stderr)
        task.run(directory=args.directory, runopt=args.runopt, auto=args.auto)
    elif args.runtype == 6:
        # dfpt-phonon
        params["ecutsm"] = args.ecutsm

        from pymatflow.abinit.dfpt import DfptPhonon
        from pymatflow.cmd.matflow import get_kpath
        task = DfptPhonon()
        task.get_qpath(get_kpath(args.kpath_manual, args.kpath_file))

        task.get_xyz(xyzfile)
        task.set_params(params=params)
        task.set_kpoints(kpoints=kpoints)
        task.set_properties(properties=args.properties)
        task.set_run(mpi=args.mpi, server=server, jobname=args.jobname, nodes=args.nodes, ppn=args.ppn, queue=args.queue)
        task.set_llhpc(partition=args.partition, nodes=args.nodes, ntask=args.ntask, jobname=args.jobname, stdout=args.stdout, stderr=args.stderr)
        task.run(directory=args.directory, runopt=args.runopt, auto=args.auto)
    elif args.runtype == 7:
        # phonopy phonon

        params["ecutsm"] = args.ecutsm

        from pymatflow.abinit.phonopy import PhonopyRun
        task = PhonopyRun()
        task.get_xyz(xyzfile)
        task.set_params(params=params)
        task.set_kpoints(kpoints=kpoints)
        task.supercell_n = args.supercell_n
        task.set_run(mpi=args.mpi, server=server, jobname=args.jobname, nodes=args.nodes, ppn=args.ppn, queue=args.queue)
        task.set_llhpc(partition=args.partition, nodes=args.nodes, ntask=args.ntask, jobname=args.jobname, stdout=args.stdout, stderr=args.stderr)
        task.phonopy(directory=args.directory, runopt=args.runopt, auto=args.auto)
    elif args.runtype == 8:
        # abc opt
        params["optcell"] = 0 # must be 0
        params["ionmov"] = args.ionmov

        params["ecutsm"] = args.ecutsm

        from pymatflow.abinit.opt import OptRun
        task = OptRun()
        task.get_xyz(xyzfile)
        task.set_params(params=params)
        task.set_kpoints(kpoints=kpoints)
        task.set_run(mpi=args.mpi, server=server, jobname=args.jobname, nodes=args.nodes, ppn=args.ppn, queue=args.queue)
        task.set_llhpc(partition=args.partition, nodes=args.nodes, ntask=args.ntask, jobname=args.jobname, stdout=args.stdout, stderr=args.stderr)
        task.batch_a = args.batch_a     
        task.batch_b = args.batch_b
        task.batch_c = args.batch_c     
        task.abc(directory=args.directory, runopt=args.runopt, auto=args.auto, range_a=args.range_a, range_b=args.range_b, range_c=args.range_c)
    elif args.runtype == 9:
        # optic
        params["ecutsm"] = args.ecutsm

        from pymatflow.abinit.optic import OpticRun
        task = OpticRun()
        task.get_xyz(xyzfile)
        task.set_params(params=params)
        task.set_kpoints(kpoints=kpoints)
        task.set_run(mpi=args.mpi, server=server, jobname=args.jobname, nodes=args.nodes, ppn=args.ppn, queue=args.queue)
        task.set_llhpc(partition=args.partition, nodes=args.nodes, ntask=args.ntask, jobname=args.jobname, stdout=args.stdout, stderr=args.stderr)
        task.run(directory=args.directory, runopt=args.runopt, auto=args.auto)
    elif args.runtype == 10:
        # optic
        params["ecutsm"] = args.ecutsm

        from pymatflow.abinit.bse import BseRun
        task = BseRun()
        task.get_xyz(xyzfile)
        task.set_params(params=params)
        task.set_kpoints(kpoints=kpoints)
        task.set_run(mpi=args.mpi, server=server, jobname=args.jobname, nodes=args.nodes, ppn=args.ppn, queue=args.queue)
        task.set_llhpc(partition=args.partition, nodes=args.nodes, ntask=args.ntask, jobname=args.jobname, stdout=args.stdout, stderr=args.stderr)
        task.run(directory=args.directory, runopt=args.runopt, auto=args.auto)    
    else:
        pass