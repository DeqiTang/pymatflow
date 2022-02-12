
import os

def siestaSubparser(subparsers):
    # --------------------------------------------------------------------------
    # SIESTA
    # --------------------------------------------------------------------------
    subparser = subparsers.add_parser("siesta", help="using siesta as calculator")

    gp = subparser.add_argument_group(title="overall running control:")

    gp.add_argument("-r", "--runtype", type=int, default=0,
        choices=[0, 1, 2, 3, 4, 5, 6, 7],
        help="choices of runtype. 0->static_run; 1->optimization; 2->cubic-cell; 3->hexagonal-cell; 4->tetragonal-cell; 5->phonopy; 6->molecular dynamics; 7->abc")

    gp.add_argument("-d", "--directory", type=str, default="matflow-running",
        help="Directory for the running.")

    gp.add_argument("--runopt", type=str, default="gen",
        choices=["gen", "run", "genrun"],
        help="Generate or run or both at the same time.")

    gp.add_argument("--auto", type=int, default=3,
        choices=[0, 1, 2, 3],
        help="auto:0 nothing, 1: copying files to server, 2: copying and executing, 3: pymatflow run inserver with direct submit,  in order use auto=1, 2, you must make sure there is a working ~/.pymatflow/server_[pbs|llhpc].conf")

    # -----------------------------------------------------------------
    #                      run params
    # -----------------------------------------------------------------

    gp.add_argument("--mpi", type=str, default="",
        help="MPI command: like 'mpirun -np 4'")

    gp.add_argument("--server", type=str, default="pbs",
        choices=["pbs", "llhpc", "tianhe2"],
        help="type of remote server, can be pbs or llhpc")

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
    structfile = subparser.add_mutually_exclusive_group(required=True) # at leaset one of cif and xyz is provided
    # argparse will make sure only one of argument in structfile(xyz, cif) appear on command line
    structfile.add_argument("--xyz", type=str, default=None,
        help="The xyz structure file with the second line specifying the cell parameter")

    structfile.add_argument("--cif", type=str, default=None,
        help="The cif structure file")

    structfile.add_argument("--xsd", type=str, default=None,
        help="The xsd structure file")

    structfile.add_argument("--xsf", type=str, default=None,
        help="The xsf structure file")

    # potential file
    gp = subparser.add_argument_group(title="pseudopotential")

    gp.add_argument("--pot", type=str, default="./",
        help="specify the path to dir containing all the needed pseudopotential, default behavior is find them in the current directory automatically. if you pass 'auto' to it, matflow will get the pots automatically(need simple configuration, see manual)")

    # --------------------------------------------------------------------------

    gp = subparser.add_argument_group(title="electronic")

    gp.add_argument("--meshcutoff", type=int, default=200,
        help="MeshCutoff (Ry)")

    gp.add_argument("--solution-method", type=str, default="diagon",
        choices=["diagon", "OMM", "OrderN", "PEXSI"],
        help="SolutionMethod(diagon, OMM, OrderN, PEXSI)")

    gp.add_argument("--functional", type=str, default="GGA",
        help="XC.functional")

    gp.add_argument("--authors", type=str, default="PBE",
        help="XC.authors")

    gp.add_argument("--tolerance", type=float, default=1.0e-6,
        help="DM.Tolerance")

    gp.add_argument("--numberpulay", type=int, default=8,
        help="DM.NumberPulay")

    gp.add_argument("--mixing", type=float, default=0.1,
        help="DM.MixingWeight")

    gp.add_argument("--kpoints-mp", type=int, nargs="+",
        default=[3, 3, 3],
        help="set kpoints like '3 3 3'")

    gp.add_argument("--kpath-manual", type=str, nargs="+", default=None,
        help="manual input kpath for band structure calculation")

    gp.add_argument("--kpath-file", type=str,
        help="file to read the kpath for band structure calculation")

    gp.add_argument("--occupation", type=str, default="FD",
        choices=["FD", "MP"],
        help="OccupationFunction(FD or MP)")

    gp.add_argument("--electronic-temperature", type=int, default=300,
        help="Electronic Temperature")

    gp.add_argument("--pao-fix-split-table", type=str, default=None,
        choices=["T", "F"],
        help="can fix problem with small split_norm WARNING")
    # properties related parameter
    # ------------------------------
    gp = subparser.add_argument_group(title="properties")

    gp.add_argument("-p", "--properties" ,nargs="+", type=int, default=[],
        help="Option for properties calculation. 1->PDOS; 2->LDOS; 3->Bands; 4->Charge Density; 5->Chemical analysis; 6->Macro Polarization; 7->Net Charge Dipole Electric Field; 8->Optical; 9->Wannier90 ")

    gp.add_argument("--pdos-block", type=float, nargs="+",
        default=[-20, 10, 0.2, 500])

    gp.add_argument("--polarization-grids", nargs="+", type=str,
        default=["10 3 3 no", "2 20 2 no", "4 4 15 no"],
        help="PolarizationGrids")

    gp.add_argument("--external-electric-field", nargs="+", type=float,
        default=[0.0, 0.0, 0.5],
        help="External Electric field")

    gp.add_argument("--optical-energy-minimum", type=float,
        default=0.0,
        help="Optical.Energy.Minimum")

    gp.add_argument("--optical-energy-maximum", type=float,
        default=10.0,
        help="Optical.Energy.Maximum")

    gp.add_argument("--optical-broaden", type=float,
        default=0.0,
        help="Optical.Broaden")

    gp.add_argument("--optical-scissor", type=float,
        default=0.0,
        help="Optical.Scissor")

    gp.add_argument("--optical-mesh", nargs="+", type=int,
        default=[5, 5, 5],
        help="Optical.Mesh")

    gp.add_argument("--optical-polarization-type", type=str,
        default="unpolarized",
        help="Optical.PolarizationType")

    gp.add_argument("--optical-vector", nargs="+", type=float,
        default=[1.0, 0.0, 0.5],
        help="Optical.Vector")

    gp.add_argument("--wannier90-unkgrid", nargs="+", type=int,
        default=[10, 10, 10],
        help="Siesta2Wannier90.UnkGrid[1-3]")

    #           ions relaed parameter
    # ==================================================
    gp = subparser.add_argument_group(title="ions")

    gp.add_argument("--variablecell", type=str, default="false",
        choices=["true", "false"],
        help="MD.VariableCell")

    gp.add_argument("--forcetol", type=float, default=0.04,
        help="Force tolerance in coordinate optimization. default=0.04 eV/Ang")

    gp.add_argument("--stresstol", type=float, default=1,
        help="Stress tolerance in variable-cell CG optimization. default=1 GPa")

    gp.add_argument("--targetpressure", type=float, default=0,
        help="Target pressure for Parrinello-Rahman method, variable cell optimizations, and annealing options.")

    gp.add_argument("--mdstep", type=int, default=1000,
        help="Final time step of the MD simulation.")

    gp.add_argument("--timestep", type=float, default=1.0,
        help="Length of the time step of the MD simulation.")

    gp.add_argument("--initial-temp", type=float, default=0,
        help="Initial temperature for the MD run.")

    gp.add_argument("--target-temp", type=float, default=0,
        help="arget temperature for Nose thermostat and annealing options.")

    # na nc stepa stepc
    # --------------------------------------------------------------------------
    gp = subparser.add_argument_group(title="cell optimization",
        description="setting of parameters needed in matflow cubic, hexagonal, tetragonal cell parameters optimization")

    gp.add_argument("--na", type=int, default=10,
        help="number of a used")

    gp.add_argument("--nc", type=int, default=10,
        help="number of c used")

    gp.add_argument("--stepa", type=float, default=0.05,
        help="a step")

    gp.add_argument("--stepc", type=float, default=0.05,
        help="c step")

    #      Phonopy
    # -------------------------------
    gp = subparser.add_argument_group(title="phonopy")

    gp.add_argument("-n", "--supercell-n", type=int, nargs="+",
        default=[1, 1,1],
        help="supercell option for phonopy, like '2 2 2'")

    # range_a range_c
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


def siestaDriver(args):
    # ==============================================================================
    # SIESTA SIESTA SIESTA SIESTA SIESTA SIESTA SIESTA SIESTA SIESTA SIESTA SIESTA
    # ==============================================================================
    from pymatflow.cmd.matflow import getXyzFile
    xyzfile, images = getXyzFile(args)
    # server
    # xxx.set_run can only deal with pbs, llhpc, lsf_sz server now 
    # however both guangzhou chaosuan llhpc are build on tianhe2, so they can use the same job system(yhbatch...)
    # we add tianhe2 option to args.server which cannot be handled by xxx.set_run. so we convert it to llhpc if tianhe2 is chosen
    server = args.server if args.server != "tianhe2" else "llhpc"
          
    from pymatflow.siesta.base import default_units
    params = {}

    params["MeshCutoff"] = args.meshcutoff
    params["SolutionMethod"] = args.solution_method
    params["XC.funtional"] = args.functional
    params["XC.authors"] = args.authors
    params["DM.Tolerance"] = args.tolerance
    params["DM.NumberPulay"] = args.numberpulay
    params["DM.MixingWeight"] = args.mixing
    params["OccupationFunction"] = args.occupation
    params["ElectronicTemperature"] = args.electronic_temperature

    params["PAO.FixSplitTable"] = args.pao_fix_split_table
    
    
    if args.runtype == 0:
        # static
        from pymatflow.cmd.matflow import get_kpath
        from pymatflow.siesta.static import StaticRun
        task = StaticRun()
        task.get_xyz(xyzfile)

        task.properties.set_params(
            #bandlines = args.bandlines,
            #bandpoints = args.bandpoints,
            polarization_grids = args.polarization_grids,
            external_electric_field = args.external_electric_field,
            optical_energy_minimum = args.optical_energy_minimum,
            optical_energy_maximum = args.optical_energy_maximum,
            optical_broaden = args.optical_broaden,
            optical_scissor = args.optical_scissor,
            optical_mesh = args.optical_mesh,
            optical_polarization_type = args.optical_polarization_type,
            optical_vector = args.optical_vector,
            wannier90_unkgrid = args.wannier90_unkgrid,
            )

        if 3 in args.properties:
            task.properties.bandlines = get_kpath(args.kpath_manual, args.kpath_file)

        task.set_params(params=params, units=default_units)
        task.set_kpoints(kpoints_mp=args.kpoints_mp)
        task.set_run(mpi=args.mpi, server=server, jobname=args.jobname, nodes=args.nodes, ppn=args.ppn, queue=args.queue)
        task.set_llhpc(partition=args.partition, nodes=args.nodes, ntask=args.ntask, jobname=args.jobname, stdout=args.stdout, stderr=args.stderr)
        task.scf(directory=args.directory, runopt=args.runopt, auto=args.auto, properties=args.properties)
    elif args.runtype == 1:
        # optimization
        from pymatflow.siesta.opt import OptRun
        params["MD.VariableCell"] = args.variablecell
        params["MD.MaxForceTol"] = args.forcetol
        params["MD.MaxStressTol"] = args.stresstol
        params["MD.TargetPressure"] = args.targetpressure
        task = OptRun()
        task.get_xyz(xyzfile)
        task.set_params(params=params, units=default_units)
        task.set_kpoints(kpoints_mp=args.kpoints_mp)
        task.set_run(mpi=args.mpi, server=server, jobname=args.jobname, nodes=args.nodes, ppn=args.ppn, queue=args.queue)
        task.set_llhpc(partition=args.partition, nodes=args.nodes, ntask=args.ntask, jobname=args.jobname, stdout=args.stdout, stderr=args.stderr)
        task.opt(directory=args.directory, runopt=args.runopt, auto=args.auto)
    elif args.runtype == 2:
        # cubic cell
        from pymatflow.siesta.opt import OptRun
        params["MD.VariableCell"] = "false"
        params["MD.MaxForceTol"] = args.forcetol
        params["MD.MaxStressTol"] = args.stresstol
        params["MD.TargetPressure"] = args.targetpressure
        task = OptRun()
        task.get_xyz(xyzfile)
        task.set_params(params=params, units=default_units)
        task.set_kpoints(kpoints_mp=args.kpoints_mp)
        task.set_run(mpi=args.mpi, server=server, jobname=args.jobname, nodes=args.nodes, ppn=args.ppn, queue=args.queue)
        task.set_llhpc(partition=args.partition, nodes=args.nodes, ntask=args.ntask, jobname=args.jobname, stdout=args.stdout, stderr=args.stderr)
        task.cubic(directory=args.directory, runopt=args.runopt, auto=args.auto, na=args.na, stepa=args.stepa)
    elif args.runtype == 3:
        # hexagonal cell
        from pymatflow.siesta.opt import OptRun
        params["MD.VariableCell"] = "false"
        params["MD.MaxForceTol"] = args.forcetol
        params["MD.MaxStressTol"] = args.stresstol
        params["MD.TargetPressure"] = args.targetpressure
        task = OptRun()
        task.get_xyz(xyzfile)
        task.set_params(params=params, units=default_units)
        task.set_kpoints(kpoints_mp=args.kpoints_mp)
        task.set_run(mpi=args.mpi, server=server, jobname=args.jobname, nodes=args.nodes, ppn=args.ppn, queue=args.queue)
        task.set_llhpc(partition=args.partition, nodes=args.nodes, ntask=args.ntask, jobname=args.jobname, stdout=args.stdout, stderr=args.stderr)
        task.hexagonal(directory=args.directory, runopt=args.runopt, auto=args.auto, na=args.na, nc=args.nc, stepa=args.stepa, stepc=args.stepc)
    elif args.runtype == 4:
        # tetragonal cell
        from pymatflow.siesta.opt import OptRun
        params["MD.VariableCell"] = "false"
        params["MD.MaxForceTol"] = args.forcetol
        params["MD.MaxStressTol"] = args.stresstol
        params["MD.TargetPressure"] = args.targetpressure
        task = OptRun()
        task.get_xyz(xyzfile)
        task.set_params(params=params, units=default_units)
        task.set_kpoints(kpoints_mp=args.kpoints_mp)
        task.set_run(mpi=args.mpi, server=server, jobname=args.jobname, nodes=args.nodes, ppn=args.ppn, queue=args.queue)
        task.set_llhpc(partition=args.partition, nodes=args.nodes, ntask=args.ntask, jobname=args.jobname, stdout=args.stdout, stderr=args.stderr)
        task.tetragonal(directory=args.directory, runopt=args.runopt, auto=args.auto, na=args.na, nc=args.nc, stepa=args.stepa, stepc=args.stepc)
    elif args.runtype == 5:
        # phonopy
        from pymatflow.siesta.phonopy import PhonopyRun
        task = PhonopyRun()
        task.get_xyz(xyzfile)
        task.set_params(params=params, units=default_units)
        task.set_kpoints(kpoints_mp=args.kpoints_mp)
        task.supercell_n = args.supercell_n
        task.set_run(mpi=args.mpi, server=server, jobname=args.jobname, nodes=args.nodes, ppn=args.ppn, queue=args.queue)
        task.set_llhpc(partition=args.partition, nodes=args.nodes, ntask=args.ntask, jobname=args.jobname, stdout=args.stdout, stderr=args.stderr)
        task.phonopy(directory=args.directory, runopt=args.runopt, auto=args.auto)
    elif args.runtype == 6:
        # molecular dynamics
        from pymatflow.siesta.md import DdRun
        params["MD.FinalTimeStep"] = args.mdstep
        params["MD.LengthTimeStep"] = args.timestep
        params["MD.InitialTemperature"] = args.initial_temp
        params["MD.TargetTemperature"] = args.target_temp         
        task = MdRun()
        task.get_xyz(xyzfile)
        task.set_params(params=params, units=default_units)
        task.set_kpoints(kpoints_mp=args.kpoints_mp)
        task.set_run(mpi=args.mpi, server=server, jobname=args.jobname, nodes=args.nodes, ppn=args.ppn, queue=args.queue)
        task.set_llhpc(partition=args.partition, nodes=args.nodes, ntask=args.ntask, jobname=args.jobname, stdout=args.stdout, stderr=args.stderr)
        task.md(directory=args.directory, runopt=args.runopt, auto=args.auto)
    elif args.runtype == 7:
        # abc cell opt
        from pymatflow.siesta.opt import OptRun
        task = OptRun()
        task.get_xyz(xyzfile)
        task.set_relax()
        task.set_kpoints(kpoints_mp=args.kpoints_mp)
        task.set_params(params=params, units=default_units)
        task.set_run(mpi=args.mpi, server=server, jobname=args.jobname, nodes=args.nodes, ppn=args.ppn, queue=args.queue)
        task.set_llhpc(partition=args.partition, nodes=args.nodes, ntask=args.ntask, jobname=args.jobname, stdout=args.stdout, stderr=args.stderr)
        task.batch_a = args.batch_a     
        task.batch_b = args.batch_b
        task.batch_c = args.batch_c     
        task.abc(directory=args.directory, runopt=args.runopt, auto=args.auto, range_a=args.range_a, range_b=args.range_b, range_c=args.range_c)
    else:
        pass        