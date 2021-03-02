#!/usr/bin/env python

import os
import sys
import argparse



def get_kpath(kpath_manual=None, kpath_file=None):
    """
    :param kpath_manual: manual input kpath like --kpath-manual '0.000000 0.000000 0.000000 GAMMA 5' '0.500000 0.000000 0.000000 X 5' '0.0000 0.000 0.50000 A |' '0.5 0.5 0.5 R '
    :param kpath_file: manual input kpath read from the file
    :return kpath or None(when kpath_manual and kpath_file are both None)
    """
    # dealing with standard kpath
    kpath = None
    if kpath_manual != None:
        # kpath from script argument args.kpath
        kpath = []
        for kpoint in kpath_manual:
            if kpoint.split()[4] != "|":
                kpath.append([
                    float(kpoint.split()[0]),
                    float(kpoint.split()[1]),
                    float(kpoint.split()[2]),
                    kpoint.split()[3].upper(),
                    int(kpoint.split()[4]),
                    ])
            elif kpoint.split()[4] == "|":
                kpath.append([
                    float(kpoint.split()[0]),
                    float(kpoint.split()[1]),
                    float(kpoint.split()[2]),
                    kpoint.split()[3].upper(),
                    "|",
                    ])
        return kpath
    elif kpath_file != None:
        # kpath read from file specified by kpath_file
        # file is in format like this
        """
        5
        0.0 0.0 0.0 #GAMMA 15
        x.x x.x x.x #XXX |
        x.x x.x x.x #XXX 10
        x.x x.x x.x #XXX 15
        x.x x.x x.x #XXX 20
        """
        # if there is a '|' behind the label it means the path is
        # broken after that point!!!
        kpath = []
        with open(kpath_file, 'r') as fin:
            lines = fin.readlines()
        nk = int(lines[0])
        for i in range(nk):
            if lines[i+1].split("\n")[0].split()[4] != "|":
                kpath.append([
                    float(lines[i+1].split()[0]),
                    float(lines[i+1].split()[1]),
                    float(lines[i+1].split()[2]),
                    lines[i+1].split()[3].split("#")[1].upper(),
                    int(lines[i+1].split()[4]),
                    ])
            elif lines[i+1].split("\n")[0].split()[4] == "|":
                kpath.append([
                    float(lines[i+1].split()[0]),
                    float(lines[i+1].split()[1]),
                    float(lines[i+1].split()[2]),
                    lines[i+1].split()[3].split("#")[1].upper(),
                    '|',
                    ])
        return kpath
    else:
        pass
    # -------------------------------------------------------------------



def main():
    parser = argparse.ArgumentParser()

    # --------------------------------------------------------------------------
    # Octopus
    # --------------------------------------------------------------------------
    #parser = subparsers.add_parser("octopus", help="using abinit as calculator")

    gp = parser.add_argument_group(title="overall running control",
            description="control the overall running parameters")

    gp.add_argument("-r", "--runtype", type=int, default=0,
            choices=[0, 1, 2, 3, 4, 5, 6, 7],
            help="choices of runtype. 0->static_run; 1->optimization; 2->cubic-opt; 3->hexagonal-opt; 4->tetragonal-opt;")

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
            choices=["pbs", "llhpc"],
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
    # actually this can be put in the main parser, but it will make the command not like git sub-cmmand
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
    gp = parser.add_argument_group(title="pseudopotential")

    gp.add_argument("--pot", type=str, default="./",
            help="specify the path to the directory containing all the needed pseudopotential, default behavior is find them in the current directory automatically. if you pass 'auto' to it, matflow will get the pots automatically(need simple configuration, see manual)")

    # calculation mode
    gp = parser.add_argument_group(title="Calculation Mode",
            description="setting of calculation mode related parameters")

    # execution
    gp = parser.add_argument_group(title="Execution",
            description="setting of execution related parameters")

    # hamiltonian
    gp = parser.add_argument_group(title="Hamiltonian",
            description="setting of hamiltonian related parameters")

    # linear_response
    gp = parser.add_argument_group(title="Linear Response",
            description="setting of linear response related parameters")

    # math
    gp = parser.add_argument_group(title="Math",
            description="setting of math related parameters")

    # mesh
    gp = parser.add_argument_group(title="Mesh",
            description="setting of mesh related parameters")

    gp.add_argument("--kpath-manual", type=str, nargs="+", default=None,
            help="manual input kpath for band structure calculation")

    gp.add_argument("--kpath-file", type=str,
            help="file to read the kpath for band structure calculation")

    gp.add_argument("--kpoints-mp", type=int, nargs=6,
            default=[1, 1, 1, 0, 0, 0],
            help="monkhorst-pack type k mesh generation using KPointsGrid")

    # output
    gp = parser.add_argument_group(title="Output",
            description="setting of output related parameters")
    
    # scf
    gp = parser.add_argument_group(title="Scf",
            description="setting of scf related parameters")

    # states
    gp = parser.add_argument_group(title="States",
            description="setting of states related parameters")

    # system
    gp = parser.add_argument_group(title="System",
            description="setting of system related parameters")

    # time dependent
    gp = parser.add_argument_group(title="Time dependent",
            description="setting of time dependent related parameters")

    # utilities
    gp = parser.add_argument_group(title="Utilities",
            description="setting of utilities related parameters")

    # na stepa nc stepc
    # ------------------------------------------------
    gp = parser.add_argument_group(title="cell optimization",
            description="setting of parameters needed by cubic, hexagonal, tetragonal cell parameters optimization needed by matflow")

    gp.add_argument("--na", type=int, default=10,
            help="number of a to run")
    gp.add_argument("--stepa", type=float, default=0.05,
            help="step of a in unit of Angstrom")
    gp.add_argument("--nc", type=int, default=10,
            help="number of c to run")
    gp.add_argument("--stepc", type=float, default=0.05,
            help="step of c in unit of Angstrom")


    # ==========================================================
    # transfer parameters from the arg parser to static_run setting
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




    # dealing with pseudo potential file
    if args.pot == "./":
        #TODO make a simple check, whether there exists the potential file
        pass
    elif args.pot == "auto":
        pass
    else:
        os.system("cp %s/* ./" % args.pot)



    params = {}

    params[""] = None


    if args.runtype == 0:
        # static
        from pymatflow.octopus.static import static_run
        task = static_run()
        if get_kpath(args.kpath_manual, args.kpath_file) == None:
            print("================================================\n")
            print("Warning: octopus\n")
            print("in abinit static runing you must provide kpath\n")
            sys.exit(1)
        task.dataset[3].electrons.kpoints.set_band(kpath=get_kpath(args.kpath_manual, args.kpath_file))
        task.get_xyz(xyzfile)
        task.set_params(params=params)
        task.set_kpoints(kpoints=kpoints)
        task.set_run(mpi=args.mpi, server=args.server, jobname=args.jobname, nodes=args.nodes, ppn=args.ppn, queue=args.queue)
        task.set_llhpc(partition=args.partition, nodes=args.nodes, ntask=args.ntask, jobname=args.jobname, stdout=args.stdout, stderr=args.stderr)
        task.run(directory=args.directory, runopt=args.runopt, auto=args.auto)
    elif args.runtype == 1:
        # optimization
        params["optcell"] = args.optcell
        params["chkdilatmx"] = args.chkdilatmx
        params["dilatmx"] = args.dilatmx
        params["ionmov"] = args.ionmov
        params["ecutsm"] = args.ecutsm
        from pymatflow.octopus.opt import opt_run
        task = opt_run()
        task.get_xyz(xyzfile)
        task.set_params(params=params)
        task.set_kpoints(kpoints=kpoints)
        task.set_run(mpi=args.mpi, server=args.server, jobname=args.jobname, nodes=args.nodes, ppn=args.ppn, queue=args.queue)
        task.set_llhpc(partition=args.partition, nodes=args.nodes, ntask=args.ntask, jobname=args.jobname, stdout=args.stdout, stderr=args.stderr)
        task.run(directory=args.directory, runopt=args.runopt, auto=args.auto)
    elif args.runtype == 2:
        ## cubic optimization
        params["optcell"] = 0 # must be 0
        params["ionmov"] = args.ionmov
        from pymatflow.octopus.opt import opt_run
        task = opt_run()
        task.get_xyz(xyzfile)
        task.set_params(params=params)
        task.set_kpoints(kpoints=kpoints)
        task.set_run(mpi=args.mpi, server=args.server, jobname=args.jobname, nodes=args.nodes, ppn=args.ppn, queue=args.queue)
        task.set_llhpc(partition=args.partition, nodes=args.nodes, ntask=args.ntask, jobname=args.jobname, stdout=args.stdout, stderr=args.stderr)
        task.cubic(directory=args.directory, runopt=args.runopt, auto=args.auto, na=args.na, stepa=args.stepa)
    elif args.runtype == 3:
        # hexagonal optimization
        params["optcell"] = 0 # must be 0
        params["ionmov"] = args.ionmov
        from pymatflow.octopus.opt import opt_run
        task = opt_run()
        task.get_xyz(xyzfile)
        task.set_params(params=params)
        task.set_kpoints(kpoints=kpoints)
        task.set_run(mpi=args.mpi, server=args.server, jobname=args.jobname, nodes=args.nodes, ppn=args.ppn, queue=args.queue)
        task.set_llhpc(partition=args.partition, nodes=args.nodes, ntask=args.ntask, jobname=args.jobname, stdout=args.stdout, stderr=args.stderr)
        task.hexagonal(directory=args.directory, runopt=args.runopt, auto=args.auto, na=args.na, stepa=args.stepa, nc=args.nc, stepc=args.stepc)
    elif args.runtype == 4:
        # tetragonal optimization
        params["optcell"] = 0 # must be 0
        params["ionmov"] = args.ionmov
        from pymatflow.octopus.opt import opt_run
        task = opt_run()
        task.get_xyz(xyzfile)
        task.set_params(params=params)
        task.set_kpoints(kpoints=kpoints)
        task.set_run(mpi=args.mpi, server=args.server, jobname=args.jobname, nodes=args.nodes, ppn=args.ppn, queue=args.queue)
        task.set_llhpc(partition=args.partition, nodes=args.nodes, ntask=args.ntask, jobname=args.jobname, stdout=args.stdout, stderr=args.stderr)
        task.tetragonal(directory=args.directory, runopt=args.runopt, auto=args.auto, na=args.na, stepa=args.stepa, nc=args.nc, stepc=args.stepc)
    else:
        pass


if __name__ == "__main__":
    main()
