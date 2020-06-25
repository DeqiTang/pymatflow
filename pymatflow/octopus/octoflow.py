#!/usr/bin/env python

import os
import sys
import argparse

from pymatflow.cmd.matflow import get_kpath



def main():
    parser = argparse.ArgumentParser()
    # --------------------------------------------------------------------------
    # Octopus
    # --------------------------------------------------------------------------

    gp = parser.add_argument_group(title="overall running control",
            description="control the overall running parameters")

    #                      run params
    # -----------------------------------------------------------------
    gp = parser.add_argument_group(title="overall running control", description="control the task running parameters")

    gp.add_argument("-d", "--directory", type=str, default="matflow-running",
            help="directory to generate all the files, do not specify the current directory")

    gp.add_argument("-r", "--runtype", type=int, default=0,
            choices=[0, 1, 2, 3, 4, 5, 6, 7, 8, 9],
            help="choices of runtype. 0->static_run; 1->optimization; 2->cubic-cell; 3->hexagonal-cell; 4->tetragonal-cell; 5->neb; 6->vasp-phonon; 7->phonopy; 8->surf pes; 9->abc")

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
            choices=["pbs", "llhpc", "lsf_sz"],
            help="type of remote server, can be pbs or llhpc or lsf_sz")

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
    gp = parser.add_mutually_exclusive_group(required=True) # at leaset one of cif and xyz is provided
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
    #gp = parser.add_argument_group(title="pseudopotential")

    #gp.add_argument("--pot", type=str, default="./",
    #        help="specify the path to the potential file, default is ./. if you pass 'auto' to it, matflow will build the POTCAR foryou(need simple configuration, see manual)")
        
    #gp.add_argument("--pot-type", type=str, default="PAW_PBE",
    #        choices=["PAW_PBE", "PAW_LDA", "PAW_PW91", "paw_pbe", "paw_lda", "paw_pw91"],
    #        help="choose type of POT for POTCAR")

    # --------------------------------------------------------
    #                   inp parameters
    # --------------------------------------------------------
    
    gp = parser.add_argument_group(title="Special",
            description="Special parameters")
            
    gp.add_argument("--kpoints-mp", type=int, nargs=6,
            default=[1, 1, 1, 0, 0, 0],
            help="set kpoints like -k 1 1 1 0 0 0")

    gp.add_argument("--kpoints-mp-nscf", type=int, nargs="+",
            default=None,
            help="set kpoints like -k 1 1 1 0 0 0")

    gp.add_argument("--kpath-manual", type=str, nargs="+", default=None,
            help="set kpoints for band structure calculation manually")

    gp.add_argument("--kpath-file", type=str, default="kpath.txt",
            help="set kpoints for band structure calculation manually from file")
    
    # Calculation Modes
    gp = parser.add_argument_group(title="Calculation Modes",
            description="Calculation Modes related parameters")

    gp.add_argument("--calculationmode", type=str, default=None,
            help="Decides what kind of calculation is to be performed.")
            
    gp.add_argument("--geocenter", type=str, default=None,
            help="If set to yes, Octopus centers the geometry at every optimization step.")

    gp.add_argument("--gofireintegrator", type=str, default=None,
            choices=["euler", "verlet"],
            help="The Fire algorithm (GOMethod = fire) uses a molecular dynamics integrator to compute new geometries and velocities.")

    gp.add_argument("--golinetol", type=float, default=None,
            help="Tolerance for line-minimization. Applies only to GSL methods that use the forces.")

    gp.add_argument("--gomaxiter", type=int, default=None,
            help="Even if the convergence criterion is not satisfied, the minimization will stop after this number of iterations.")
            
    gp.add_argument("--gomethod", type=str, default=None,
            choices=["steep_native", "steep", "cg_fr", "cg_pr", "cg_bfgs", "cg_bfgs2", "simplex", "fire"],
            help="Method by which the minimization is performed.")
    
    gp.add_argument("--gominimummove", type=float, default=None,
            help="Convergence criterion, for stopping the minimization. In units of length; minimization is stopped when the coordinates of all species change less than GOMinimumMove, or the GOTolerance criterion is satisfied.")
            
    gp.add_argument("--goobjective", type=str, default=None,
            choices=["minimize_energy", "minimize_forces"],
            help="This rather esoteric option allows one to choose which objective function to minimize during a geometry minimization.")
            
    gp.add_argument("--gostep", type=float, default=None,
            help="Initial step for the geometry optimizer. The default is 0.5.")
            
    gp.add_argument("--gotolerance", type=float, default=None,
            help="Convergence criterion, for stopping the minimization. In units of force; minimization is stopped when all forces on ions are smaller than this criterion, or the GOMinimumMove is satisfied. If GOTolerance < 0, this criterion is ignored.")
            
    gp.add_argument("--invertksmethod", type=str, default=None,
            choices=["two_particle", "iterative", "iter_stella", "iter_godby"],
            help="Selects whether the exact two-particle method or the iterative scheme is used to invert the density to get the KS potential.")
            
    # Execution
    gp = parser.add_argument_group(title="Execution",
            description="Execution related parameters")

    gp.add_argument("--splines", type=str, default=None,
            choices=["gsl", "native"],
            help="Selects the implementation of the spline interpolation.")

    # Hamiltonian
    gp = parser.add_argument_group(title="Hamiltonian",
            description="Hamiltonian related parameters")

    gp.add_argument("--theorylevel", type=str, default=None,
            choices=["hartree", "independent_particles", "hartree_fock", "dft", "rdmft"],
            help="The calculations can be run with different \"theory levels\" that control how electrons are simulated. The default is dft. When hybrid functionals are requested, through the XCFunctional variable, the default is hartree_fock.")
  
    # Linear Response
    gp = parser.add_argument_group(title="Linear Response",
            description="Linear Response related parameters")

    gp.add_argument("--responsemethod", type=str, default=None,
            choices=["sternheimer", "finite_differences"],
            help="Some response properties can be calculated either via Sternheimer linear response or by using finite differences. You can use this variable to select how you want them to be calculated, it applies to em_resp and vib_modes calculation modes. By default, the Sternheimer linear-response technique is used.")  

    # Math
    gp = parser.add_argument_group(title="Math",
            description="Math related parameters")

    gp.add_argument("--rootsolver", type=str, default=None,
            choices=["root_newton", "root_watterstrom"],
            help="Specifies what kind of root solver will be used.")

    # Mesh
    gp = parser.add_argument_group(title="Mesh",
            description="Mesh related parameters")

    gp.add_argument("--spacing", type=float, default=None,
            help="The spacing between the points in the mesh. This controls the quality of the discretization: smaller spacing gives more precise results but increased computational cost.")

    # Output
    gp = parser.add_argument_group(title="Output",
            description="Output related parameters")

    gp.add_argument("--outputformat", type=str, default=None,
            choices=["axis_x", "axis_y", "matlab", "meshgrid", "boundary_points", "binary", "etsf",
                    "xyz", "cube", "bild", "axis_z", "vtk", "integrate_xy", "integrate_xz", "integrate_yz",
                    "plane_x", "plane_y", "plane_z", "dx", "netcdf", "mesh_index", "xcrysden"],
            help="Specifies what kind of root solver will be used.")

    # SCF
    gp = parser.add_argument_group(title="SCF",
            description="SCF related parameters")

    gp.add_argument("--scfcalculatedipole", type=str, default=None,
            choices=["yes", "no"],
            help="This variable controls whether the dipole is calculated at the end of a self-consistent iteration. For finite systems the default is yes. For periodic systems the default is no, unless an electric field is being applied in a periodic direction. The single-point Berry`s phase approximation is used for periodic directions.")


    # States
    gp = parser.add_argument_group(title="States",
            description="States related parameters")

    gp.add_argument("--smearing", type=float, default=None,
            help="If Occupations is not set, Smearing is the smearing width used in the SmearingFunction to distribute the electrons among the existing states.")
            
    gp.add_argument("--smearingfunction", type=str, default=None,
            choices=["semiconducting", "fermi_dirac", "cold_smearing", "methfessel_paxton", "spline_smearing"],
            help="This is the function used to smear the electronic occupations. It is ignored if the Occupations block is set.")

    # System
    gp = parser.add_argument_group(title="System",
            description="System related parameters")

    gp.add_argument("--dimensions", type=int, default=None,
            choices=[0, 1, 2, 3],
            help="Octopus can run in 1, 2 or 3 dimensions, depending on the value of this variable (or more, if configured with --with-max-dim=4 or higher). ")

    # Time-Dependent
    gp = parser.add_argument_group(title="Time-Dependent",
            description="Time-Dependent related parameters")

    gp.add_argument("--tdfunctions", type=str, default=None,
            choices=["tdf_cw", "tdf_gaussian", "tdf_cosinoidal", "tdf_trapezoidal", "tdf_from_file", "tdf_from_expr"],
            help="This block specifies the shape of a \"time-dependent function\", such as the envelope needed when using the TDExternalFields block.")

    # Utilities
    gp = parser.add_argument_group(title="Utilities",
            description="Utilities related parameters")

    gp.add_argument("--volume", type=str, default=None,
            choices=["vol_sphere", "vol_slab"],
            help="Describes a volume in space defined through the addition and substraction of spheres. ")


    # range_a range_b range_c
    # ----------------------------------------------
    gp = parser.add_argument_group(title="cell optimization",
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

    # inp template
    gp = parser.add_argument_group(title="template", 
            description="read in inp template")

    gp.add_argument("--inp", type=str, default=None,
            help="specify inp template to set parameters")

    # surf pes
    gp = parser.add_argument_group(title="surf pes",
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
            help="0 -> do not fix any z of the atoms, 1: only fix z of the buttom atoms, 2: fix z of both the buttom and the moving atoms. note x y are all fixed")

    gp.add_argument("--batch-x-y", type=int, nargs=2, default=None,
            help="number of structures to calculate each batch x and y, default is all in one batch")

    # fix atoms
    gp = parser.add_argument_group(title="fix atoms",
            description="specify atoms to fix in optimization, only used when --runtype=1")

    gp.add_argument("--fix", help="list of fixed atoms, index start from 1", nargs='+', type=int)


    # static calc related setting
    gp = parser.add_argument_group(title="static calc",
            description="setting type of static calculation when -r 0")

    gp.add_argument("--static", type=str, default="band",
            choices=["scf", "band", "dos", "optics", "bse", "stm"],
            help="in case of band(default), run scf, nscf(bands) in a single run; in case of scf, run scf only, in case of optics, run scf and optics calc in a single run")

    gp.add_argument("--hse-in-scf", type=str, default="false",
            choices=["true", "false", "True", "False"],
            help="choose whether to use HSE in both scf and nscf or only in nscf, when calc band structure")

    gp.add_argument("--bse-level", type=int, default=0,
            choices=[0, 1, 2],
            help="0 -> bse on standard DFT; 1 -> bse on hybrid functional; 2 -> bse on GW")

    gp.add_argument("--algo-gw", type=str, default="EVGW",
            choices=["EVGW", "GW", "GW0", "QPGW0", "QPGW"],
            help="ALGO used for GW")

            
    # ==========================================================
    # transfer parameters from the arg subparser to static_run setting
    # ==========================================================

    args = parser.parse_args()

    # if no argument passed to matflow
    if len(sys.argv) == 1:
        # display help message when no args provided
        parser.print_help()
        sys.exit(1)

    # check directory
    if os.getcwd() == os.path.abspath(args.directory):
        print("================================================================\n")
        print("                    Warning from matflow\n")
        print("----------------------------------------------------------------\n")
        print("you are trying to specify current directory for task running\n")
        print("which is not allowed.\n")
        print("try specify it to a new directory like -d matflow-running\n")
        sys.exit(1)

    # dealing wich structure files
    if args.xyz != None:
        xyzfile = args.xyz
    elif args.cif != None:
        os.system("structflow convert -i %s -o %s.xyz" % (args.cif, args.cif))
        xyzfile = "%s.xyz" % args.cif
    elif args.xsd != None:
        os.system("structflow convert -i %s -o %s.xyz" % (args.xsd, args.xsd))
        xyzfile = "%s.xyz" % args.xsd
    elif args.xsf != None:
        os.sytem("structflow convert -i % -o %s.xyz" % (args.xsf, args.xsf))
        xyzfile = "%s.xyz" % args.xsf
    else:
        # neb caculattion with
        images = []
        for image in args.images:
            if image.split(".")[-1] == "xyz":
                images.append(image)
            else:
                os.system("structflow convert -i %s -o %s.xyz" % (image, image))
                images.append("%s.xyz" % image)
        xyzfile = images[0] # this set only for dealing with pseudo potential file
        #

    # octopus run type 
    params = {}
    # deal with octopus inp template specified by --inp
    if args.inp == None:
        pass
    else:
        if not os.path.exists(args.incar):
            print("====================================================\n")
            print("                  Warning !!!!\n")
            print("----------------------------------------------------\n")
            print("matflow vasp:\n")
            print("the specified incar file by --inp doesn't exist\n")
            print("go and check it\n")
            sys.exit(1)
        with open(args.incar, 'r') as fin:
            incar = fin.readlines()
        for line in incar:
            if len(line.split()) == 0:
                continue
            if line[0] == "#":
                continue
            if len(line.split("\n")[0].split("#")[0].split("=")) == 2:
                # in case of single value INCAR variable
                params[line.split("=")[0].split()[0].upper()] = line.split("\n")[0].split("#")[0].split("=")[1].split()[0]
            else:
                params[line.split("=")[0].split()[0].upper()] = line.split("\n")[0].split("#")[0].split("=")[1].split()
    #
    # if xxx is alraedy in params(set from --inp) and args.xxx is None
    # params[xxx] will not be control by args.xxx
    params["Calculation Modes/CalculationMode"] = args.calculationmode
    params["Calculation Modes/Geometry Optimization/GOCenter"] = args.geocenter
    params["Calculation Modes/Geometry Optimization/GOFireIntegrator"] = args.gofireintegrator
    params["Calculation Modes/Geometry Optimization/GOLineTol"] = args.golinetol
    params["Calculation Modes/Geometry Optimization/GOMaxIter"] = args.gomaxiter
    params["Calculation Modes/Geometry Optimization/GOMethod"] = args.gomethod
    params["Calculation Modes/Geometry Optimization/GOMinimumMove"] = args.gominimummove
    params["Calculation Modes/Geometry Optimization/GOObjective"] = args.goobjective
    params["Calculation Modes/Geometry Optimization/GOStep"] = args.gostep
    params["Calculation Modes/Geometry Optimization/GOTolerance"] = args.gotolerance
    params["Calculation Modes/Invert KS/InvertKSmethod"] = args.invertksmethod
    params["Execution/Splines"] = args.splines
    params["Hamiltonian/TheoryLevel"] = args.theorylevel
    params["Linear Response/ResponseMethod"] = args.responsemethod
    params["Math/RootSolver/RootSolver"] = args.rootsolver
    params["Mesh/Spacing"] = args.spacing
    params["Output/OutputFormat"] = args.outputformat
    params["SCF/SCFCalculateDipole"] = args.scfcalculatedipole
    params["States/Smearing"] = args.smearing
    params["States/SmearingFunction"] = args.smearingfunction
    params["System/Dimensions"] = args.dimensions
    params["Time-Dependent/TDFunctions"] = args.tdfunctions
    params["Utilities/Volume"] = args.volume

    
    if args.runtype == 0:
        # static
        from pymatflow.octopus.static import static_run
        task = static_run()
        task.get_xyz(xyzfile)
        task.set_params(params, runtype="static")
        task.set_kpoints(kpoints_mp=args.kpoints_mp)
        task.set_run(mpi=args.mpi, server=args.server, jobname=args.jobname, nodes=args.nodes, ppn=args.ppn, queue=args.queue)
        task.set_llhpc(partition=args.partition, nodes=args.nodes, ntask=args.ntask, jobname=args.jobname, stdout=args.stdout, stderr=args.stderr)
        
        if args.static == "scf":
            task.set_kpoints(kpoints_mp=args.kpoints_mp)
            task.scf(directory=args.directory, runopt=args.runopt, auto=args.auto)
        elif args.static == "band":
            if args.hse_in_scf.lower() == "true":
                hse_in_scf = True
            elif args.hse_in_scf.lower() == "false":
                hse_in_scf = False
            task.band(directory=args.directory, runopt=args.runopt, auto=args.auto, kpath=get_kpath(args.kpath_manual, args.kpath_file), hse_in_scf=hse_in_scf)                
        elif args.static == "dos":
            if args.hse_in_scf.lower() == "true":
                hse_in_scf = True
            elif args.hse_in_scf.lower() == "false":
                hse_in_scf = False
            if args.kpoints_mp_nscf == None:
                kpoints_mp_nscf = args.kpoints_mp #[2*k for k in args.kpoints_mp]
            else:
                kpoints_mp_nscf = args.kpoints_mp_nscf
            task.dos(directory=args.directory, runopt=args.runopt, auto=args.auto, hse_in_scf=hse_in_scf, kpoints_mp_nscf=kpoints_mp_nscf)
        elif args.static == "optics":
            task.set_kpoints(kpoints_mp=args.kpoints_mp)                    
            task.optics(directory=args.directory, runopt=args.runopt, auto=args.auto)
        elif args.static == "bse":
            task.set_kpoints(kpoints_mp=args.kpoints_mp)
            task.bse(directory=args.directory, runopt=args.runopt, auto=args.auto, bse_level=args.bse_level, algo_gw=args.algo_gw)
        elif args.static == "stm":
            if args.hse_in_scf.lower() == "true":
                hse_in_scf = True
            elif args.hse_in_scf.lower() == "false":
                hse_in_scf = False            
            task.set_kpoints(kpoints_mp=args.kpoints_mp)
            task.stm(directory=args.directory, runopt=args.runopt, auto=args.auto, hse_in_scf=hse_in_scf)                
    elif args.runtype == 1:
        # optimization
        from pymatflow.octopus.opt import opt_run
        #
        if args.fix != None:
            fix_str = ""
            for i in args.fix:
                    fix_str += "%d " % i
            os.system("xyz-fix-atoms.py -i %s -o %s --fix %s" % (xyzfile, xyzfile, fix_str))
            args.selective_dynamics = "T"
        #            
        task = opt_run()
        task.get_xyz(xyzfile)
        task.set_params(params=params, runtype="opt")
        task.set_kpoints(kpoints_mp=args.kpoints_mp)
        task.poscar.selective_dynamics = True if args.selective_dynamics.upper()[0] == "T" else False
        task.set_run(mpi=args.mpi, server=args.server, jobname=args.jobname, nodes=args.nodes, ppn=args.ppn, queue=args.queue)
        task.set_llhpc(partition=args.partition, nodes=args.nodes, ntask=args.ntask, jobname=args.jobname, stdout=args.stdout, stderr=args.stderr)
        task.optimize(directory=args.directory, runopt=args.runopt, auto=args.auto)
    elif args.runtype == 2:
        # cubic cell
        from pymatflow.octopus.opt import opt_run
        # some must set parameters 
        if params["IBRION"] == None:
            params["IBRION"] = 2
        params["ISIF"] = 2
        if params["NSW"] == None:
            params["NSW"] = 100
        task = opt_run()
        task.get_xyz(xyzfile)
        task.set_params(params=params, runtype="opt")
        task.set_kpoints(kpoints_mp=args.kpoints_mp)
        task.set_run(mpi=args.mpi, server=args.server, jobname=args.jobname, nodes=args.nodes, ppn=args.ppn, queue=args.queue)
        task.set_llhpc(partition=args.partition, nodes=args.nodes, ntask=args.ntask, jobname=args.jobname, stdout=args.stdout, stderr=args.stderr)
        task.batch_a = args.batch_a
        task.cubic(directory=args.directory, runopt=args.runopt, auto=args.auto, range_a=args.range_a)
    elif args.runtype == 3:
        # hexagonal cell
        from pymatflow.octopus.opt import opt_run
        task = opt_run()
        task.get_xyz(xyzfile)
        task.set_params(params=params, runtype="opt")
        task.set_kpoints(kpoints_mp=args.kpoints_mp)
        task.set_run(mpi=args.mpi, server=args.server, jobname=args.jobname, nodes=args.nodes, ppn=args.ppn, queue=args.queue)
        task.set_llhpc(partition=args.partition, nodes=args.nodes, ntask=args.ntask, jobname=args.jobname, stdout=args.stdout, stderr=args.stderr)
        task.batch_a = args.batch_a
        task.batch_c = args.batch_c            
        task.hexagonal(directory=args.directory, runopt=args.runopt, auto=args.auto, range_a=args.range_a, range_c=args.range_c)
    elif args.runtype == 4:
        # tetragonal cell
        from pymatflow.octopus.opt import opt_run
        task = opt_run()
        task.get_xyz(xyzfile)
        task.set_params(params=params, runtype="opt")
        task.set_kpoints(kpoints_mp=args.kpoints_mp)
        task.set_run(mpi=args.mpi, server=args.server, jobname=args.jobname, nodes=args.nodes, ppn=args.ppn, queue=args.queue)
        task.set_llhpc(partition=args.partition, nodes=args.nodes, ntask=args.ntask, jobname=args.jobname, stdout=args.stdout, stderr=args.stderr)
        task.batch_a = args.batch_a     
        task.batch_c = args.batch_c            
        task.tetragonal(directory=args.directory, runopt=args.runopt, auto=args.auto, range_a=args.range_a, range_c=args.range_c)
    elif args.runtype == 5:
        # neb
        # we better set NSW manually in VTST neb calc. 
        # if not set, pymatflow.vasp.neb will set it to 100 automatically
        from pymatflow.octopus.neb import neb_run
        task = neb_run()
        task.get_images(images)
        task.set_params(params=params, runtype="neb")
        task.set_kpoints(kpoints_mp=args.kpoints_mp)
        task.nimage = args.nimage
        if args.nebmake == 1 and args.moving_atom == None:
            print("============================================\n")
            print("when using nebmake.py to generate inter image\n")
            print("you have to specify the moving atoms.\n")
            print("index start from 0\n")
            sys.exit(1)
        task.nebmake = "nebmake.pl" if args.nebmake == 0 else "nebmake.py"
        task.moving_atom = args.moving_atom
        task.set_run(mpi=args.mpi, server=args.server, jobname=args.jobname, nodes=args.nodes, ppn=args.ppn, queue=args.queue)
        task.set_llhpc(partition=args.partition, nodes=args.nodes, ntask=args.ntask, jobname=args.jobname, stdout=args.stdout, stderr=args.stderr)
        task.neb(directory=args.directory, runopt=args.runopt, auto=args.auto)
        # move the OUTCAR for initial stucture and final structure to the corresponding dir
        # if they are specified
        if args.outcars != None and len(args.outcars) > 0:
            os.system("cp %s %s" % (args.outcars[0], os.path.join(args.directory, "00/")))
            os.system("cp %s %s" % (args.outcars[-1], os.path.join(args.directory, "%.2d/" % (args.nimage+1))))                
    elif args.runtype == 6:
        # vasp phonon
        from pymatflow.octopus.phonon import phonon_run
        task = phonon_run() 
        task.get_xyz(xyzfile)
        task.set_params(params=params, runtype="phonon")
        task.set_kpoints(kpoints_mp=args.kpoints_mp)
        task.supercell_n = args.supercell_n
        task.set_run(mpi=args.mpi, server=args.server, jobname=args.jobname, nodes=args.nodes, ppn=args.ppn, queue=args.queue)
        task.set_llhpc(partition=args.partition, nodes=args.nodes, ntask=args.ntask, jobname=args.jobname, stdout=args.stdout, stderr=args.stderr)
        task.phonon(directory=args.directory, runopt=args.runopt, auto=args.auto)
    elif args.runtype == 7:
        # phonopy
        from pymatflow.octopus.phonopy import phonopy_run
        task = phonopy_run()
        task.get_xyz(xyzfile)
        task.set_params(params=params, runtype="phonopy")
        task.set_kpoints(kpoints_mp=args.kpoints_mp)
        task.supercell_n = args.supercell_n
        task.set_run(mpi=args.mpi, server=args.server, jobname=args.jobname, nodes=args.nodes, ppn=args.ppn, queue=args.queue)
        task.set_llhpc(partition=args.partition, nodes=args.nodes, ntask=args.ntask, jobname=args.jobname, stdout=args.stdout, stderr=args.stderr)
        task.phonopy(directory=args.directory, runopt=args.runopt, auto=args.auto)
    elif args.runtype == 8:
        # sur pes
        properties_resp_slab_sampling_surf_direction
    elif args.runtype == 9:
        # abc cell
        from pymatflow.octopus.opt import opt_run
        task = opt_run()
        task.get_xyz(xyzfile)
        task.set_params(params=params, runtype="opt")
        task.set_kpoints(kpoints_mp=args.kpoints_mp)
        task.set_run(mpi=args.mpi, server=args.server, jobname=args.jobname, nodes=args.nodes, ppn=args.ppn, queue=args.queue)
        task.set_llhpc(partition=args.partition, nodes=args.nodes, ntask=args.ntask, jobname=args.jobname, stdout=args.stdout, stderr=args.stderr)
        task.batch_a = args.batch_a     
        task.batch_b = args.batch_b
        task.batch_c = args.batch_c            
        task.abc(directory=args.directory, runopt=args.runopt, auto=args.auto, range_a=args.range_a, range_b=args.range_b, range_c=args.range_c)
# --------------------------------------------------------------------------



if __name__ == "__main__":
    main()
