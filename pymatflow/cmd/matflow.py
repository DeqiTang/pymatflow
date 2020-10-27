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
    subparsers = parser.add_subparsers(dest="driver", title="subcommands", description="choose one and only one calculator")

    # --------------------------------------------------------------------------
    # Abinit
    # --------------------------------------------------------------------------
    subparser = subparsers.add_parser("abinit", help="using abinit as calculator")

    gp = subparser.add_argument_group(title="overall running control",
            description="control the overall running parameters")

    gp.add_argument("-r", "--runtype", type=int, default=0,
            choices=[0, 1, 2, 3, 4, 5, 6, 7, 8],
            help="choices of runtype. 0->static_run; 1->optimization; 2->cubic-opt; 3->hexagonal-opt; 4->tetragonal-opt; 5->dfpt-elastic-piezo-dielec; 6->dfpt-phonon; 7->phonopy; 8->abc")

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
    #        default=[1, 1, 1, 0, 0, 0],
    #        help="monkhorst-pack type k mesh generation using ngkpt.")

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
            

    # --------------------------------------------------------------------------
    # CP2K
    # --------------------------------------------------------------------------
    subparser = subparsers.add_parser("cp2k", help="using cp2k as calculator")

    gp = subparser.add_argument_group(title="overall running control")

    gp.add_argument("-r", "--runtype", type=int, default=0,
            choices=[0, 1, 2, 3, 4 ,5, 6, 7, 8, 9, 10, 11, 12],
            help="choices of runtype. 0->static_run; 1->geo-opt; 2->cell-opt; 3->cubic-cell; 4->hexagonal-cell; 5->tetragonal-cell; 6-neb; 7->phonopy; 8->vibrational_analysis; 9->converge test; 10->aimd; 11->abc; 12->metadynamics")

    gp.add_argument("-d", "--directory", type=str, default="matflow-running",
            help="Directory to do the calculation")

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
    gp = subparser.add_mutually_exclusive_group(required=True) # at leaset one of cif and xyz is provided
    # argparse will make sure only one of argument in structfile(xyz, cif) appear on command line
    gp.add_argument("--xyz", type=str, default=None,
            help="The xyz structure file with the second line specifying the cell parameter")

    gp.add_argument("--cif", type=str, default=None,
            help="The cif structure file")

    gp.add_argument("--xsf", type=str, default=None,
            help="The xsf structure file")

    gp.add_argument("--xsd", type=str, default=None,
            help="The xsd structure file")

    gp.add_argument("--images", type=str, nargs="+",
            help="the image stucture file(--images first.cif final.xsd), can only be cif, xsd, xsd, or xyz(second line is cell parameter) format")

    # potential file and basis set
    gp = subparser.add_argument_group(title="pseudopotential and basis set")

    gp.add_argument("--pot", type=str, default="auto",
            choices=["auto"],
            help="setting pseudopotential file, in cp2k can only be auto(no need to set)")

    gp.add_argument("--kind-pot", type=str, nargs="+", default=None,
            help="set like this: --kind-pot H GTH-PBE Li GTH-PBE ...[they have default setting]")

    gp.add_argument("--kind-basis", type=str, nargs="+", default=None,
            help="set like this: --kind-basis H DZVP-MOLOPT-SR-GTH Li DZVP-MOLOPT-SR-GTH ...[they have default setting]")

    gp.add_argument("--basis-file", type=str, default="BASIS_MOLOPT",
            help="set the BASIS_SET_FILE_NAME")

    gp.add_argument("--pot-file", type=str, default="GTH_POTENTIALS",
            help="set the POTENTIAL_FILE_NAME")

    # GLOBAL
    gp = subparser.add_argument_group(title="GLOBAL")

    gp.add_argument("--print-level", type=str, default=None,
            choices=["DEBUG", "HIGH", "LOW", "MEDIUM", "SILENT", "debug", "high", "low", "medium", "silent"],
            help="How much output is written out.")

    # FORCE_EVAL/SUBSYS
    gp = subparser.add_argument_group(title="FORCE_EVAL/SUBSYS")

    gp.add_argument("--cell-periodic", type=str, default=None,
            choices=["NONE", "X", "XY", "XYZ", "XZ", "Y", "YZ", "Z", "none", "x", "xy", "xyz", "xz", "y", "yz", "z"],
            help="Specify the directions for which periodic boundary conditions (PBC) will be applied. Important notice: This applies to the generation of the pair lists as well as to the application of the PBCs to positions. See the POISSON section to specify the periodicity used for the electrostatics. Typically the settings should be the same. DEFAULT is XYZ")

    gp.add_argument("--cell-symmetry", type=str, default=None,
            help="Imposes an initial cell symmetry. must be set when you do KEEP_SYMMETRY cell_opt")

    # FORCE_EVAL/DFT
    gp = subparser.add_argument_group(title="FORCE_EVAL/DFT")

    gp.add_argument("--qs-method", type=str, default=None, #"gpw",
            choices=["am1", "dftb", "gapw", "gapw_xc", "gpw", "lrigpw", "mndo", "mndod",
                "ofgpw", "pdg", "pm3", "pm6", "pm6-fm", "pnnl", "rigpw", "rm1"],
            help="specify the electronic structure method that should be employed, default is gpw")
            
    gp.add_argument("--lsd", type=str, default=None,
            choices=["TRUE", "FALSE", "true", "false"],
            help="Requests a spin-polarized calculation using alpha and beta orbitals, i.e. no spin restriction is applied.Alias names for this keyword: UNRESTRICTED_KOHN_SHAM, UKS, SPIN_POLARIZED")

    gp.add_argument("--charge", type=int, default=None,
            help="The total charge of the system")

    gp.add_argument("--surface-dipole-correction", type=str, default=None,
            choices=["TRUE", "FALSE", "true", "false"],
            help="For slab calculations with asymmetric geometries, activate the correction of the electrostatic potential with by compensating for the surface dipole. Implemented only for slabs with normal parallel to one Cartesian axis. The normal direction is given by the keyword SURF_DIP_DIR")

    gp.add_argument("--surf-dip-dir", type=str, default=None,
            choices=["X", "Y", "Z", "x", "y", "z"],
            help="Cartesian axis parallel to surface normal")

    # FORCE_EVAL/DFT/POISSON
    gp = subparser.add_argument_group(title="FORCE_EVAL/DFT/POISSON")

    gp.add_argument("--poisson-periodic", type=str, default=None,
            choices=["NONE", "X", "XY", "XYZ", "XZ", "Y", "YZ", "Z", "none", "x", "xy", "xyz", "xz", "y", "yz", "z"],
            help="Specify the directions in which PBC apply. Important notice, this only applies to the electrostatics. See the CELL section to specify the periodicity used for e.g. the pair lists. Typically the settings should be the same. default is XYZ")
            
    gp.add_argument("--poisson-solver", type=str, default=None,
            choices=["ANALYTIC", "IMPLICIT", "MT", "MULTIPOLE", "PERIODIC", "WAVELET", "analytic", "implicit", "mt", "multipole", "periodic", "wavelet"],
            help="Specify which kind of solver to use to solve the Poisson equation.")

    # FORCE_EVAL/DFT/SCF
    gp = subparser.add_argument_group(title="FORCE_EVAL/DFT/SCF")

    gp.add_argument("--max-scf", type=int, default=50,
            help="Maximum number of SCF iteration to be performed for one optimization.")

    gp.add_argument("--eps-default", type=float, default=1.0e-14,
            help="Try setting all EPS_xxx to values leading to an energy correct up to EPS_DEFAULT")

    gp.add_argument("--eps-scf", type=float, default=1.0e-6,
            help="target accuracy for the scf convergence, default is 1.0e-6")

    gp.add_argument("--level-shift", type=float, default=None,
            help="Use level shifting to improve convergence, default is 0.0")

    gp.add_argument("--max-scf-history", type=int, default=None, #0
            help="Maximum number of SCF iterations after the history pipeline is filled")
            
    gp.add_argument("--max-diis", type=int, default=None, #4
            help="Maximum number of DIIS vectors to be used")

    gp.add_argument("--xc-functional", type=str, default="pbe",
            choices=["B3LYP", "BEEFVDW", "BLYP", "BP", "LDA", "PBE", "PADE", "PBE0", "TPSS",
                "b3lyp", "beefvdw", "blyp", "bp", "lda", "pbe", "pade", "pbe0", "tpss"],
            help="shortcut for the most common functional combinations, default is PBE")

    gp.add_argument("--cutoff", type=int, default=None, #100,
            help="The cutoff of the finest grid level, default value: 100 Ry")

    gp.add_argument("--rel-cutoff", type=int, default=None, #60,
            help="determines the grid at which a Gaussian is mapped, giving the cutoff used for a gaussian with alpha=1. A value 50+-10Ry might be required for highly accurate results, Or for simulations with a variable cell, default value: 60 Ry")

    gp.add_argument("--ngrids", type=int, default=4,
            help="The number of multigrids to use, default is 4")

    gp.add_argument("--kpoints-scheme", type=str,
            default=None, #"GAMMA",
            help="kpoint scheme to be used, can be NONE, GAMMA, MONKHORST-PACK, MACDONALD, GENERAL. when you set MONKHORST-PACK, you should also add the three integers like 'monkhorst-pack 3 3 3'")

    gp.add_argument("--kpath-manual", type=str, nargs="+", default=None,
            help="manual input kpath for band structure calculation")

    gp.add_argument("--kpath-file", type=str,
            help="file to read the kpath for band structure calculation")

    gp.add_argument("--diag", type=str, default="TRUE",
            choices=["TRUE", "FALSE", "true", "false"],
            help="whether choosing tranditional diagonalization for SCF")

    gp.add_argument("--diag-algo", type=str, default="STANDARD",
            choices=["DAVIDSON", "FILTER_MATRIX", "OT", "STANDARD",
            "davidson", "filter_matrix", "ot", "standard"],
            help="Algorithm to be used for diagonalization")

    gp.add_argument("--ot", type=str, default="FALSE",
            choices=["TRUE", "FALSE", "true", "false"],
            help="whether choosing orbital transformation for SCF")

    gp.add_argument("--mixing-method", type=str, default=None,
            choices=["BROYDEN_MIXING", "BROYDEN_MIXING_NEW", "DIRECT_P_MIXING", "KERKER_MIXING", "MULTISECANT_MIXING",
            "NONE", "PULAY_MIXING", "broyden_mixing", "broyden_mixing_new", "direct_p_mixing", "kerker_mixing", "multisecant_mixing",
            "none", "pulay_mixing"],
            help="Mixing method to be applied")

    gp.add_argument("--mixing-alpha", type=float, default=0.4,
            help="fraction of new density to be included, default is 0.4")

    gp.add_argument("--mixing-beta", type=float, default=0.5, 
            help="Denominator parameter in Kerker damping introduced to suppress charge sloshing: rho_mix(g) =rho_in(g) + alpha*g^2/(g^2 + beta^2)*(rho_out(g)-rho_in(g))")

    gp.add_argument("--mixing-nbuffer", type=int, default=None, # default is 4
            help="Number of previous steps stored for the actual mixing scheme. default is 4")

    # outer scf
    gp.add_argument("--outer-scf", type=str, default=None, # "FALSE",
            choices=["TRUE", "FALSE", "true", "false"],
            help="controls the activation of the outer SCF loop")

    gp.add_argument("--outer-scf-bisect-trust-count", type=int, default=None,
            help="Maximum number of times the same point will be used in bisection, a small number guards against the effect of wrongly converged states.")

    gp.add_argument("--outer-scf-diis-buffer-length", type=int, default=None,
            help="Maximum number of DIIS vectors used")

    gp.add_argument("--outer-scf-extrapolation-order", type=int, default=None,
            choices=[1, 2, 3, 4],
            help="Number of past states used in the extrapolation of the variables during e.g. MD")

    gp.add_argument("--outer-scf-eps-scf", type=str, default=None,
            help="The target gradient of the outer SCF variables. Notice that the EPS_SCF of the inner loop also determines the value that can be reached in the outer loop, typically EPS_SCF of the outer loop must be smaller than or equal to EPS_SCF of the inner loop.")

    gp.add_argument("--outer-scf-max-scf", type=int, default=None,
            help="The maximum number of outer loops")

    gp.add_argument("--outer-scf-optimizer", type=str, default=None,
            choices=["BISECT", "BROYDEN", "DIIS", "NEWTON", "NEWTON_LS", "NONE", "SD", "SECANT", "bisect", "broyden", "diis", "newton", "newton_ls", "none", "sd", "secant"],
            help="Method used to bring the outer loop to a stationary point")

    gp.add_argument("--outer-scf-type", type=str, default=None,
            choices=["BASIS_CENTER_OPT", "CDFT_CONSTRAINT", "DDAPC_CONSTRAINT", "NONE", "S2_CONSTRAINT", "basis_center_opt", "cdft_constraint", "ddapc_constraint", "none"],
            help="Specifies which kind of outer SCF should be employed ")
    
    # smear
    gp.add_argument("--smear", type=str, default=None, #"FALSE",
            choices=["TRUE", "FALSE", "true", "false"],
            help="controls the activation of smearing")

    gp.add_argument("--smear-method", type=str, default="FERMI_DIRAC",
            help="Smearing method to be applied, can be fermi_dirac(default) or energy_window")

    gp.add_argument("--added-mos", type=int, default=None, #0,
            help="Number of additional MOS added for each spin")

    gp.add_argument("--electronic-temp", type=float, default=None, #300,
            help="Electronic temperature in the case of Fermi-Dirac smearing in unit of [K], default is 300")

    gp.add_argument("--eps-fermi-dirac", type=float, default=None,
            help="Accuracy checks on occupation numbers use this as a tolerance. default is 1.0E-10")

    gp.add_argument("--window-size", type=float, default=None, #0,
            help="size of the energy window centred at the Fermi level for energy_window type smearing")

    gp.add_argument("--ls-scf", type=str, default=None, #"false",
            choices=["true", "false", "true", "false"],
            help="use linear scaling scf method")

    gp.add_argument("--ot-preconditioner", type=str, default=None,
            choices=["FULL_ALL", "FULL_KINETIC", "FULL_SINLGE", "FULL_SINGLE_INVERSE", "FULL_S_INVERSE", "NONE", "full_all", "full_kinetic", "full_single", "full_single_inverse", "full_s_inverse", "none"],
            help="Type of preconditioner to be used with all minimization schemes. They differ in effectiveness, cost of construction, cost of application. Properly preconditioned minimization can be orders of magnitude faster than doing nothing.")

    gp.add_argument("--ot-energy-gap", type=float, default=None, 
            help="Should be an estimate for the energy gap [a.u.] (HOMO-LUMO) and is used in preconditioning, especially effective with the FULL_ALL preconditioner, in which case it should be an underestimate of the gap (can be a small number, e.g. 0.002). FULL_SINGLE_INVERSE takes it as lower bound (values below 0.05 can cause stability issues). In general, heigher values will tame the preconditioner in case of poor initial guesses. A negative value will leave the choice to CP2K depending on type of preconditioner. ")            

    gp.add_argument("--ot-minimizer", type=str, default=None,
            choices=["BROYDEN", "CG", "DIIS", "SD", "broyden", "cg", "diis", "sd"],
            help="Minimizer to be used with the OT method")
            
    gp.add_argument("--ot-algorithm", type=str, default=None,
            choices=["IRAC", "STRICT", "irac", "strict"],
            help="Algorithm to be used for OT")
            
    gp.add_argument("--ot-linesearch", type=str, default=None,
            choices=["2PNT", "3PNT", "GOLD", "NONE"],
            help="1D line search algorithm to be used with the OT minimizer, in increasing order of robustness and cost. MINIMIZER CG combined with LINESEARCH GOLD should always find an electronic minimum. Whereas the 2PNT minimizer is almost always OK, 3PNT might be needed for systems in which successive OT CG steps do not decrease the total energy.")

    # vdw correction related
    gp = subparser.add_argument_group(title="FORCE_EVAL/DFT/XC")

    gp.add_argument("--vdw-potential-type", type=str, default="NONE",
            choices=["PAIR_POTENTIAL", "NON_LOCAL", "NONE", "pair_potential", "non_local", "none"],
            help="Type of dispersion/vdW functional or potential to use")

    gp.add_argument("--pair-type", type=str, default=None, #"DFTD3",
            choices=["DFTD2", "DFTD3", "DFTD3(BJ)", "dftd2", "dftd3", "dftd3(bj)"],
            help="Type of potential(VDW)")

    gp.add_argument("--reference-functional", type=str, default="PBE",
            choices=["PBE", "pbe"],
            help="Use parameters for this specific density functional. For available D3 and D3(BJ) parameters. Must be set if using VDW")
  
    gp.add_argument("--r-cutoff", type=float, default=None, #1.05835442E+001,
            help="Range of potential. The cutoff will be 2 times this value")

    # printout option
    gp = subparser.add_argument_group(title="printout option")

    gp.add_argument("-p", "--printout-option", nargs="+", type=int,
            default=[],
            choices = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14],
            help=
            """
            Properties printout option, you can also activate multiple prinout-option at the same time.
            1: printout pdos
            2: printout band
            3: printout electron densities
            4: printout electron local function(ELF)
            5: printout molecular orbitals
            6: printout molecular orbital cube files
            7: printout mulliken populaltion analysis
            8: printout cubes for generation of STM images
            9: printout cube file with total density(electrons+atomic core)
           10: printout v_hartree_cube
           11: printout v_xc_cube
           12: printout xray_diffraction_spectrum
           13: request a RESP fit of charges.
           14: printout MOMENTS
           default is no printout of these properties.
           """)

    gp.add_argument("--dft-print-elf-cube-stride", type=int, nargs="+",
            default=None, #[1, 1, 1],
            help="The stride (X,Y,Z) used to write the cube file (larger values result in smaller cube files). You can provide 3 numbers (for X,Y,Z) or 1 number valid for all components.")

    gp.add_argument("--dft-print-e-density-cube-stride", type=int, nargs="+",
            default=None, #[1, 1, 1],
            help="The stride (X,Y,Z) used to write the cube file (larger values result in smaller cube files). You can provide 3 numbers (for X,Y,Z) or 1 number valid for all components.")

    gp.add_argument("--dft-print-moments-periodic", type=str, default=None,
            choices=["TRUE", "FALSE", "true", "false"],
            help="Use Berry phase formula (PERIODIC=T) or simple operator (PERIODIC=F). The latter normally requires that the CELL is periodic NONE.")

    gp.add_argument("--dft-print-moments-reference", type=str, default=None,
            choices=["COAC", "COM", "USER_DEFINED", "ZERO", "coac", "com", "user_defined", "zero"],
            help="Define the reference point for the calculation of the electrostatic moment. default is ZERO")

    # FORCE_EVAL/PROPERTIES
    gp = subparser.add_argument_group(title="FORCE_EVAL/PROPERTIES")

    gp.add_argument("--properties-resp-slab-sampling-range", type=float, nargs="+",
            default=None, #[0.3, 3.0],
            help="Range where the fitting points are sampled. A range of 3 to 5 Angstroms means that the fitting points are sampled in the region of 3 to 5 Angstroms above the surface which is defined by atom indexes given in ATOM_LIST.")

    gp.add_argument("--properties-resp-slab-sampling-surf-direction", type=str, default=None, #"Z",
            choices=["X", "Y", "Z", "x", "y", "z", "-X", "-Y", "-Z", "-x", "-y", "-z"],
            help="Specifies what above the surface means. Defines the direction")

    gp.add_argument("--properties-resp-slab-sampling-atom-list", type=int, nargs="+",
            default=None, #[1],
            help="Specifies the list of indexes of atoms used to define the region for the RESP fitting. The list should contain indexes of atoms of the first surface layer")

    # MOTION/CELL_OPT related parameters
    # --------------------------------------------------------------------------
    gp = subparser.add_argument_group(title="MOTION/CELL_OPT")

    gp.add_argument("--cell-opt-optimizer", type=str, default=None, #"BFGS",
            choices=["BFGS", "CG", "LBFGS"],
            help="Specify which method to use to perform a geometry optimization, can be: BFGS(default), CG, LBFGS")

    gp.add_argument("--cell-opt-max-iter", type=int, default=None, #200,
            help="	Specifies the maximum number of geometry optimization steps. One step might imply several force evaluations for the CG and LBFGS optimizers")

    gp.add_argument("--cell-opt-type", type=str, default=None, #"DIRECT_CELL_OPT",
            choices=["DIRECT_CELL_OPT", "GEO_OPT", "MD"],
            help="Specify which kind of method to use for the optimization of the simulation cell, can be: DIRECT_CELL_OPT(default), GEO_OPT, MD")

    gp.add_argument("--cell-opt-max-dr", type=float, default=None, #3e-3,
            help="Convergence criterion for the maximum geometry change between the current and the last optimizer iteration in unit of bohr, default is 3.e-3")

    gp.add_argument("--cell-opt-max-force", type=float, default=None, #4.50000000E-004,
            help="Convergence criterion for the maximum force component of the current configuration in unit of bohr^-1*hartree")

    gp.add_argument("--cell-opt-rms-dr", type=float, default=None, #1.50000000E-003,
            help="Convergence criterion for the root mean square (RMS) geometry change between the current and the last optimizer iteration in unit of bohr, default is 1.5e-3")

    gp.add_argument("--cell-opt-rms-force", type=float, default=None, #3.00000000E-004,
            help="Convergence criterion for the root mean square (RMS) force of the current configuration, default is 3.e-4")

    gp.add_argument("--cell-opt-pressure-tolerance", type=float, default=None, #1.00000000E+002,
            help="Specifies the Pressure tolerance (compared to the external pressure) to achieve during the cell optimization, default is 1.0e2")
        
    gp.add_argument("--cell-opt-keep-angles", type=str, default=None,
            help="Keep angles between the cell vectors constant, but allow the lenghts of the cell vectors to change independently.")

    gp.add_argument("--cell-opt-keep-symmetry", type=str, default=None,
            help="Keep the requested initial cell symmetry (e.g. during a cell optimisation). The initial symmetry must be specified in the &CELL section.")


    # MOTION/GEO_OPT related parameters
    # --------------------------------------------------------------------------
    gp = subparser.add_argument_group(title="MOTION/GEO_OPT")

    gp.add_argument("--geo-opt-optimizer", type=str, default=None, #"BFGS",
            choices=["BFGS", "CG", "LBFGS", "bfgs", "cg", "lbfgs", "cg"],
            help="Specify which method to use to perform a geometry optimization, can be: BFGS, CG, LBFGS")

    gp.add_argument("--geo-opt-max-iter", type=int, default=None, #200,
            help="Specifies the maximum number of geometry optimization steps. One step might imply several force evaluations for the CG and LBFGS optimizers, default is 200")

    gp.add_argument("--geo-opt-type", type=str, default=None, #"MINIMIZATION",
            choices=["MINIMIZATION", "TRANSITION_STATE", "minimization", "transition_state"],
            help="Specify which kind of geometry optimization to perform, default is MINIMIZATION")

    gp.add_argument("--geo-opt-max-dr", type=float, default=None, #3e-3,
            help="Convergence criterion for the maximum geometry change between the current and the last optimizer iteration, default is 3e-3")

    gp.add_argument("--geo-opt-max-force", type=float, default=None, #4.50000000E-004,
            help="Convergence criterion for the maximum force component of the current configuration, default is 4.0e-4")

    gp.add_argument("--geo-opt-rms-dr", type=float, default=None, #1.50000000E-003,
            help="Convergence criterion for the root mean square (RMS) geometry change between the current and the last optimizer iteration, default is 1.5e-3")

    gp.add_argument("--geo-opt-rms-force", type=float, default=None, #3.00000000E-004,
            help="Convergence criterion for the root mean square (RMS) force of the current configuration.")

    # MOTION/MD 
    # --------------------------------------------------------------------------
    gp = subparser.add_argument_group(title="MOTION/MD")

    gp.add_argument("--md-steps", type=int, default=None, #1000,
            help="The number of MD steps to perform")

    gp.add_argument("--timestep", type=float, default=None, #5.0e-1,
            help="The length of an integration step (in case RESPA the large TIMESTEP), default and also recommended is 0.5 fs.")

    gp.add_argument("--ensemble", type=str, default=None, #"NVE",
            choices=["NVE",  "NVT","HYDROSTATICSHOCK", "ISOKIN", "LANGEVIN", "MSST", "MSST_DAMPED"],
            help="The ensemble/integrator that you want to use for MD propagation")

    gp.add_argument("--temperature", type=str, default=None, #300,
            help="The temperature in K used to initialize the velocities with init and pos restart, and in the NPT/NVT simulations")

    gp.add_argument("--temp-tol", type=float, default=None, #0.0,
            help="The maximum accepted deviation of the (global) temperaturefrom the desired target temperature before a rescaling of the velocites is performed. If it is 0 no rescaling is performed. NOTE: This keyword is obsolescent; Using a CSVR thermostat with a short timeconstant is recommended as a better alternative")

    gp.add_argument("--traj-format", type=str, default=None, #"XMOL",
            help="type of output trajectory for MOTION, note: DCD format can be visualized by vmd")

    # MOTION/BAND
    # --------------------------------
    gp = subparser.add_argument_group(title="MOTION/BAND")

    gp.add_argument("--band-type", type=str, default=None, #"CI-NEB",
            choices=["IT-NEB", "CI-NEB", "B-NEB", "D-NEB", "SM", "it-neb", "ci-neb", "b-neb", "d-neb", "sm"],
            help="specify the type of band calculation")

    gp.add_argument("--number-of-replica", type=int, default=None, #5,
            help="number of replicas")

    gp.add_argument("--k-spring", type=float, default=None, #2.0e-2,
            help="value of the spring constant")

    gp.add_argument("--nproc-rep", type=int, default=None,
            help="Specify the number of processors to be used per replica environment (for parallel runs) default is 1.")

    gp.add_argument("--align-frames", type=str, default=None, #"TRUE",
            choices=["TRUE", "FALSE", "true", "false"],
            help="Enables the alignment of the frames at the beginning of a BAND calculation. This keyword does not affect the rotation of the replicas during a BAND calculation.")

    gp.add_argument("--rotate-frames", type=str, default=None, #"TRUE",
            choices=["TRUE", "FALSE", "true", "false"],
            help="Compute at each BAND step the RMSD and rotate the frames in order to minimize it.")

    gp.add_argument("--optimize-end-points", type=str, default=None,
            choices=["TRUE", "FALSE", "true", "false"],
            help="Performs also an optimization of the end points of the band. default is FALSE")

    gp.add_argument("--convergence-control-max-dr", type=float, default=None,
            help="Tolerance on the maximum value of the displacement on the BAND.")

    gp.add_argument("--convergence-control-max-force", type=float, default=None,
            help="Tolerance on the maximum value of Forces on the BAND. ")

    gp.add_argument("--convergence-control-rms-dr", type=float, default=None,
            help="Tolerance on RMS displacements on the BAND.")

    gp.add_argument("--convergence-control-rms-force", type=float, default=None,
            help="Tolerance on RMS Forces on the BAND.")

    gp.add_argument('--optimize-band-opt-type', type=str, default=None,
            choices=["DIIS", "MD", "diis", "md"],
            help="Specifies the type optimizer used for the band")

    gp.add_argument("--optimize-band-diis-max-steps", type=int, default=None,
            help="Specify the maximum number of optimization steps for neb run")

    gp.add_argument("--optimize-band-md-max-steps", type=int, default=None,
            help="Specify the maximum number of optimization steps for neb run")

    # MOTION/FREE_ENERGY
    gp = subparser.add_argument_group(title="MOTION/FREE_ENERGY")
    
    gp.add_argument("--metadyn-delta-t", type=float, default=None,
            help="If Well-tempered metaD is used, the temperature parameter must be specified.")

    gp.add_argument("--metadyn-do-hills", type=str, default=None,
            choices=["TRUE", "FALSE", "true", "false"],
            help="This keyword enables the spawning of the hills. Default FALSE")

    gp.add_argument("--metadyn-nt-hills", type=int, default=None,
            help="Specify the maximum MD step interval between spawning two hills.")
            
    #             vibrational_analysis related parameters
    # ---------------------------------------------------------------
    gp = subparser.add_argument_group(title="VIBRATIONAL_ANALYSIS")

    gp.add_argument("--dx", type=float, default=None, #1.0e-2,
            help="specify the increment to be used to construct the HESSIAN with finite difference method")

    gp.add_argument("--fully-periodic", type=str, default=None, #"FALSE",
            choices=["TRUE", "FALSE", "true", "false"],
            help="avoids to clean rotations from the Hessian matrix")

    gp.add_argument("--intensities", type=str, default=None, #"FALSE",
            choices=["TRUE", "FALSE", "true", "false"],
            help="Calculation of the IR-Intensities. Calculation of dipoles has to be specified explicitly"
            )

    gp.add_argument("--tc-pressure", type=float, default=None, #1.01325000E+005,
            help="Pressure for the calculation of the thermochemical data in unit of [Pa]")

    gp.add_argument("--tc-temperature", type=float, default=None, #2.73150000E+002,
            help="Temperature for the calculation of the thermochemical data in unit of [K]")

    gp.add_argument("--thermochemistry", type=str, default=None, #"FALSE",
            help="Calculation of the thermochemical data. Valid for molecules in the gas phase.")

    gp.add_argument("--vib-nproc-rep", type=int, default=None,
            help="Specify the number of processors to be used per replica environment (for parallel runs). default is 1")

    #                   PHONOPY related parameters
    # ------------------------------------------------------------------
    gp = subparser.add_argument_group(title="phonopy:",
            description="setting parameters needed by phonopy")

    gp.add_argument("--supercell-n", nargs="+", type=int, default=[1, 1, 1],
            help="Supercell for Phonopy calculation.")

    # na nc stepa stepc
    # -----------------------------------------
    gp = subparser.add_argument_group(title="cell optimization",
            description="setting of parameters needed in cubic, hexagonal, tetragonal cell optimizaiton inf matflow")

    gp.add_argument("--na", type=int, default=10,
            help="number of a used")

    gp.add_argument("--nc", type=int, default=10,
            help="number of c used")

    gp.add_argument("--stepa", type=float, default=0.05,
            help="a step")

    gp.add_argument("--stepc", type=float, default=0.05,
            help="c step")

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

    # converge test
    gp = subparser.add_argument_group(title="converge test",
            description="setting of parameters needed in converge test")

    gp.add_argument("--converge", type=str, default="cutoff",
            choices=["cutoff", "rel_cutoff", "kpoints_auto", "kpoints_manual"],
            help="choose type of converge test")

    gp.add_argument("--cutoff-range", type=int, nargs=3,
            help="cutoff converge test range: emin emax step")

    gp.add_argument("--rel-cutoff-range", type=int, nargs=3,
            help="rel_cutoff converge test range: emin emax step")

    gp.add_argument("--krange", type=int, nargs=3,
            help="kpoints(auto) converge test range: kmin kmax step")

    gp.add_argument("--klist", type=int, nargs="+",
            help="kpoints(manual) converge test, list of kpoints")
            
    # inp template
    gp = subparser.add_argument_group(title="template inp",
            description="read inp template")
            
    gp.add_argument("--inp", type=str, default=None,
            help="read parameters from cp2k inp template")

    # --------------------------------------------------------------------------
    # Quantum ESPRESSO
    # --------------------------------------------------------------------------
    subparser = subparsers.add_parser("qe", help="using quantum espresso as calculator")

    gp = subparser.add_argument_group(title="overall running control")

    gp.add_argument("-r", "--runtype", type=int, default=0,
            choices=[0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11],
            help="choices of runtype. 0->static_run; 1->relax; 2->vc-relax; 3->cubic-cell; 4->hexagonal-cell; 5->tetragonal-cell; 6->neb; 7->dfpt; 8->phonopy; 9->pp.x; 10->abc; 11->converge")

    gp.add_argument("--static", type=str, default="all",
            choices=["all", "scf"],
            help="in case of all(default), run scf, nscf, bands in a single run; in case of scf, run scf only")


    gp.add_argument("-d", "--directory", type=str, default="matflow-running",
            help="Directory for the running.")

    # run option
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
            help="specify the path to the dir containing all the needed pseudopotential, default behavior is find them in the current directory automatically. if you pass 'auto' to it, matflow will get the pots automatically(need simple configuration, see manual)")

    gp.add_argument("--pot-type", type=str, default="sssp_efficiency",
            choices=["paw_pbe", "PAW_PBE", "sssp_efficiency", "SSSP_EFFICIENCY", "sssp_precision", "SSSP_PRECISION"],
            help="specify type of pseudo potentials to prepare, when --pot auto")

    # -------------------------------------------------------------------
    #                       scf related parameters
    # -------------------------------------------------------------------
    # &control
    gp = subparser.add_argument_group(title="pw.x->control")

    gp.add_argument("--etot-conv-thr", type=float, default=None, #1.0e-4,
            help="convergence threshold of energy for geometric optimization")

    gp.add_argument("--forc-conv-thr", type=float, default=None, #1.0e-3,
            help="convergence threshold for force in optimization,(usually it is more important than energy)")

    gp.add_argument("--tstress", type=str, default=".false.",
            choices=[".true.", ".false."],
            help="calculate stress. default=.false.")

    # &system
    gp = subparser.add_argument_group(title="pw.x->system")

    gp.add_argument("--occupations", type=str, default="smearing",
            choices=["smearing", "tetrahedra", "tetrahedra_lin", "tetrahedra_opt", "fixed", "from_input"],
            help="Occupation method for the calculation.")

    gp.add_argument("--smearing", type=str, default="gaussian",
            choices=["gaussian", "methfessel-paxton", "marzari-vanderbilt", "fermi-dirac"],
            help="Smearing type for occupations by smearing, default is gaussian in this script")

    gp.add_argument("--degauss", type=float, default=0.001,
            help="Value of the gaussian spreading (Ry) for brillouin-zone integration in metals.(defualt: 0.001 Ry)")

    gp.add_argument("--nbnd", type=int, default=None,
            help="Number of electronic states (bands) to be calculated")

    gp.add_argument("--ecutwfc",
            type=int, default=100)

    gp.add_argument("--ecutrho", type=int, default=None,
            help="Kinetic energy cutoff for charge density and potential in unit of Rydberg, default value: None")

    gp.add_argument("--vdw-corr", help="vdw_corr = dft-d, dft-d3, ts-vdw, xdm", type=str, default="none")

    gp.add_argument("--tot-charge", type=int, default=None,
            help="Total charge of the system. Useful for simulations with charged cells. tot_charge=+1 means one electron missing from the system, tot_charge=-1 means one additional electron, and so on.")

    gp.add_argument("--nosym", type=str, default=None,
            choices=[".true.", ".false."],
            help="Do not use this option unless you know exactly what you want and what you get. May be useful in the following cases: - in low-symmetry large cells, if you cannot afford a k-point grid with the correct symmetry - in MD simulations - in calculations for isolated atoms")
    
    gp.add_argument("--nosym-evc", type=str, default=None,
            choices=[".true.", ".false."],
            help="if (.TRUE.) symmetry is not used, and k points are forced to have the symmetry of the Bravais lattice;")

    gp.add_argument("--noinv", type=str, default=None,
            choices=[".true.", ".false"],
            help="if (.TRUE.) disable the usage of k => -k symmetry(time reversal) in k-point generation")
    

    # magnetic related parameters
    gp.add_argument("--nspin", type=int, default=None,
            choices=[1, 2],
            help="choose either 1 or 2, and 4 should not be used as suggested by pwscf official documentation.")

    gp.add_argument("--starting-magnetization", type=float, nargs="+", default=None,
            help="starting_magnetization(i), i=1,ntyp -> Starting spin polarization on atomic type i in a spin polarized calculation. Values range between -1 (all spins down for the valence electrons of atom type i) to 1 (all spins up).")

    gp.add_argument("--noncolin", type=str, default=None,
            choices=[".true.", ".false."],
            help="if .true. the program will perform a noncollinear calculation.")
    # DFT+U
    gp.add_argument("--lda-plus-u", type=str, default=None,
            choices=[".true.", ".false."],
            help="DFT+U (formerly known as LDA+U) currently works only for a few selected elements. Modify Modules/set_hubbard_l.f90 and PW/src/tabd.f90 if you plan to use DFT+U with an element that is not configured there.")

    gp.add_argument("--lda-plus-u-kind", type=int, default=None,
            choices=[0, 1, 2],
            help="0   DFT+U simplified version of Cococcioni and de Gironcoli, using Hubbard_U;  1   DFT+U rotationally invariant scheme of Liechtenstein et al.,using Hubbard U and Hubbard J; 2   DFT+U+V simplified version of Campo Jr and Cococcioni, using Hubbard V")

    gp.add_argument("--hubbard-u", type=float, nargs="+", default=None,
            help="Hubbard_U(i): U parameter (eV) for species i, DFT+U calculation")

    gp.add_argument("--hubbard-j0", type=float, nargs="+", default=None,
            help="Hubbard_J0(i): J0 parameter (eV) for species i, DFT+U+J calculation")

    gp.add_argument("--hubbard-alpha", type=float, default=None,
            help="Hubbard_alpha(i) is the perturbation (on atom i, in eV) used to compute U with the linear-response method of Cococcioni and de Gironcoli(only for lda_plus_u_kind=0)")

    gp.add_argument("--hubbard-beta", type=float, default=None,
            help="Hubbard_beta(i) is the perturbation (on atom i, in eV) used to compute J0 with the linear-response method of Cococcioni and de Gironcoli(only for lda_plus_u_kind=0)")

    gp.add_argument("--u-projection-type", type=str, default=None,
            choices=["atomic", "ortho-atomic", "norm-atomic", "file", "pseudo"],
            help="Only active when lda_plus_U is .true., specifies the type of projector on localized orbital to be used in the DFT+U scheme. default is atomic")

    # Hybrid functional
    gp.add_argument("--input-dft", type=str, default="pbe",
            choices=["pbe", "pbe0", "b3lyp", "hse", "vdw-DF"], 
            help="Exchange-correlation functional: eg 'PBE', 'BLYP' etc")

    gp.add_argument("--ace", type=str, default=None,
            choices=[".true.", ".false."],
            help="Use Adaptively Compressed Exchange operator as in Lin Lin, J. Chem. Theory Comput. 2016, 12, 2242--2249, Set to false to use standard Exchange (much slower)")

    gp.add_argument("--exx-fraction", type=float, default=None,
            help="Fraction of EXX for hybrid functional calculations. In the case of input_dft='PBE0', the default value is 0.25, while for input_dft='B3LYP' the exx_fraction default value is 0.20.")

    gp.add_argument("--screening-parameter", type=float, default=None,
            help="screening_parameter for HSE like hybrid functionals., default is 0.106")

    gp.add_argument("--exxdiv-treatment", type=str, default=None,
            choices=["gygi-baldereschi", "vcut_spherical", "vcut_ws", "none"],
            help="Specific for EXX. It selects the kind of approach to be used for treating the Coulomb potential divergencies at small q vectors.")
    
    gp.add_argument("--x-gamma-extrapolation", type=str, default=None,
            choices=[".true.", ".false."],
            help="Specific for EXX. If .true., extrapolate the G=0 term of the potential ")

    gp.add_argument("--ecutvcut", type=float, default=None,
            help="Reciprocal space cutoff for correcting Coulomb potential divergencies at small q vectors.")

    gp.add_argument("--nqx", type=float, nargs=3, default=[None, None, None], 
            help="Three-dimensional mesh for q (k1-k2) sampling of the Fock operator (EXX). Can be smaller than the number of k-points. Currently this defaults to the size of the k-point mesh used. In QE =< 5.0.2 it defaulted to nqx1=nqx2=nqx3=1.")

    # &electrons
    gp = subparser.add_argument_group(title="pw.x->electrons")

    gp.add_argument("--electron-maxstep", type=int, default=None,
            help="maximum number of iterations in a scf step")

    gp.add_argument("--conv-thr", type=float, default=1.0e-6,
            help="the conv_thr for scf, when doing geometric optimization better use a strict covnergec for scf")

    gp.add_argument("--mixing-beta", type=float, default=None,
            help="mixing factor for self-consistency, default is 0.7")

    gp.add_argument("--mixing-ndim", type=float, default=None,
            help="number of iterations used in mixing scheme. If you are tight with memory, you may reduce it to 4 or so.")

    gp.add_argument("--diagonalization", type=float, default=None,
            help="Available options are: david cg ppcg paro")

    gp.add_argument("--scf-must-converge", type=str, default=None,
            choices=[".true.", ".false."],
            help="If .false. do not stop molecular dynamics or ionic relaxation when electron_maxstep is reached. Use with care.")

    # &ions
    gp = subparser.add_argument_group(title="pw.x->ions")

    gp.add_argument("--ion-dynamics", type=str, default=None,
            choices=["bfgs", "damp", "verlet", "langevin", "langevin-smc", "beeman"],
            help="For different type of calculation different possibilities are allowed and different default values apply.")

    gp.add_argument("--nstep",
            type=int, default=50,
            help="maximum ion steps for geometric optimization")

    gp.add_argument("--pot-extrapolation", type=str, default=None,
            choices=["none", "atomic", "first_order", "second_order"],
            help="Used to extrapolate the potential from preceding ionic steps. first_order and second-order extrapolation make sense only for molecular dynamics calculations")
   
    gp.add_argument("--wfc-extrapolation", type=str, default=None,
            choices=["none", "first_order", "second_order"],
            help="Used to extrapolate the wavefunctions from preceding ionic steps. first_order and second-order extrapolation make sense only for molecular dynamics calculations")

    gp.add_argument("--ion-temperature", type=str, default=None,
            choices=["rescaling", "rescaling-v", "rescaling-T", "reduce-T", "berendsen", "andersen", "svr", "initial", "not_controlled"],
            help="for molecular dynamics")

    gp.add_argument("--tempw", type=float, default=None,
            help="Starting temperature (Kelvin) in MD runs target temperature for most thermostats.")

    # &cell
    gp = subparser.add_argument_group(title="pw.x->cells")

    gp.add_argument("--cell-dofree", type=str, default=None,
            choices=['all', 'ibrav', 'x', 'y', 'z', 'xy', 'xz', 'yz', 'xyz', 'shape', 'volume', '2Dxy', '2Dshape', 'epitaxial_ab', 'epitaxial_ac', 'epitaxial_bc'],
            help="cell_dofree for &cell/")

    gp.add_argument("--press-conv-thr", type=float, default=None,
            help="convergence threshold on the pressure for variable cell relaxation")

    # KPOINTS
    gp = subparser.add_argument_group(title="pw.x-K_POINTS")

    gp.add_argument("--kpoints-option", type=str, default="automatic",
            choices=["automatic", "gamma", "crystal_b"],
            help="Kpoints generation scheme option for the SCF or non-SCF calculation")

    gp.add_argument("--kpoints-mp", type=int, nargs=6,
            default=[1, 1, 1, 0, 0, 0],
            help="Monkhorst-Pack kpoint grid, in format like --kpoints-mp 1 1 1 0 0 0")

    gp.add_argument("--kpoints-mp-nscf", type=int, nargs=6,
            default=[3, 3, 3, 0, 0, 0],
            help="Monkhorst-Pack kpoint grid, in format like --kpoints-mp 3 3 3 0 0 0")

    gp.add_argument("--kpath-manual", type=str, nargs="+", default=None,
            help="manual input kpath in crystal_b, like --kpath-manual '0.000000 0.000000 0.000000 GAMMA 5' '0.500000 0.000000 0.000000 X 5' '0.0000 0.000 0.50000 A |' '0.5 0.5 0.5 R '")

    gp.add_argument("--kpath-file", type=str,
            help="manual input kpath in crystal_b read from the file")


    # ATOMIC_FORCES
    gp = subparser.add_argument_group(title="ATOMIC_FORCES")

    gp.add_argument("--pressure", type=float, default=None,
            help="specify pressure acting on system in unit of Pa")

    gp.add_argument("--pressuredir", type=str, default=None,
            choices=["x", "y", "z"],
            help="specify direction of pressure acting on system.")

    # projwfc
    gp = subparser.add_argument_group(title="projwfc")

    gp.add_argument("--projwfc-filpdos", type=str, default="projwfc",
            help="output projected dos file name")

    gp.add_argument("--projwfc-ngauss", type=str, default="default",
            help="gaussian broadening type")

    gp.add_argument("--projwfc-degauss", type=str, default='default',
            help="gaussian broadening")

    gp.add_argument("--projwfc-emin", type=str, default='default',
            help="min energy for DOS")

    gp.add_argument("--projwfc-emax", type=str, default='default',
            help="max energy for DOS")

    gp.add_argument("--projwfc-deltae", type=str, default='default',
            help="DeltaE: energy grid step (eV)")

    #         bands.x related parameters
    # -----------------------------------------
    gp = subparser.add_argument_group(title="bands.x")

    gp.add_argument("--lsym", type=str, default=".true.",
            choices=[".true.", ".false."],
            help="set lsym variable in bands.x input.")


    # na nc stepa stepc
    # -----------------------------------------
    gp = subparser.add_argument_group(title="cell optimization",
            description="setting parameters needed by matflow cubic, hexagonal, tetragonal cell optimization")

    gp.add_argument("--na", type=int, default=10,
            help="number of a used")

    gp.add_argument("--nc", type=int, default=10,
            help="number of c used")

    gp.add_argument("--stepa", type=float, default=0.05,
            help="a step")

    gp.add_argument("--stepc", type=float, default=0.05,
            help="c step")

    # neb
    gp = subparser.add_argument_group(title="neb")

    gp.add_argument("--string-method", type=str, default="neb",
            help="string_method")

    gp.add_argument("--nstep-path", type=int, default=100,
            help="nstep_path")

    gp.add_argument("--opt-scheme", type=str, default="broyden",
            help="Specify the type of optimization scheme(sd, broyden, broyden2, quick-min, langevin)")

    gp.add_argument("--num-of-images", type=int, default=5,
            help="number of total images(including the initial and final image). about how to set proper number of images: usually the inter-image distance between 1~2Bohr is OK")

    gp.add_argument("--k-max", type=float, default=0.3e0,
            help="Set them to use a Variable Elastic Constants scheme elastic constants are in the range [ k_min, k_max  ], this is useful to rise the resolution around the saddle point")

    gp.add_argument("--k-min", type=float, default=0.2e0,
            help="Set them to use a Variable Elastic Constants scheme elastic constants are in the range [ k_min, k_max  ], this is useful to rise the resolution around the saddle point")

    gp.add_argument("--ci-scheme", type=str, default="auto",
            help="Specify the type of Climbing Image scheme(no-CI, auto, manual)")

    gp.add_argument("--path-thr", type=float, default=0.05,
            help="path_thr")

    gp.add_argument("--ds", type=float, default=1.e0, help="Optimisation step length ( Hartree atomic units )")

    gp.add_argument("--first-last-opt", type=str, default=None,
            choices=[".true.", ".false."],
            help="whether to optimize the first and last image")


    # for phx
    # --------------------------------------------------------------
    gp = subparser.add_argument_group(title="ph.x")

    gp.add_argument("--tr2-ph", type=float, default=1.0e-14,
            help="threshold for self-consistency.")

    gp.add_argument("--nq", type=int, nargs="+",
            default=[0, 0, 0],
            help="set value of nq1 nq2 nq3.")

    gp.add_argument("--epsil", type=str, default=None,
            choices=[".true.", ".false."],
            help="set epsil in inputph")

    gp.add_argument("--lraman", type=str, default=None,
            choices=[".true.", ".false."],
            help="set lraman, can be 'true' or 'false' only. default is None which means 'false' in real world.")

    gp.add_argument("--search-sym", type=str, default=None,
            choices=[".true.", ".false."],
            help="set it to .false. if you want to disable the mode symmetry analysis.")

    # Phonopy
    # ---------------------------------------------------------
    gp = subparser.add_argument_group(title="phonopy")

    gp.add_argument("--supercell-n", type=int, nargs="+",
            default=[1, 1, 1],
            help="supercell build for Phonopy.")

    gp = subparser.add_argument_group(title="pp.x")

    # pp.x
    gp.add_argument("--plot-num", type=int, nargs="+", default=[0],
            choices=[0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 17, 18, 19, 20, 21],
            help="""
                type of analysis stored in the filplot file for later plot, 0: electron-pseudo-charge-density,
                    1: total-potential,
                    2: local-ionic-potential,
                    3: ldos,
                    4: local-density-of-electronic-entropy,
                    5: stm,
                    6: spin-polar,
                    7: molecular-orbitals,
                    8: electron-local-function,
                    9: charge-density-minus-superposition-of-atomic-densities,
                    10: ILDOS,
                    11: v_bare+v_H-potential,
                    12: sawtooth-electric-field-potential,
                    13: nocollinear-magnetization,
                    17: all-electron-charge-density-paw-only,
                    18: exchage-correlation-magnetic-field-noncollinear-case,
                    19: reduced-density-gradient,
                    20: product-of-charge-density-with-hessian,
                    21: all-electron-density-paw-only,""")

    gp.add_argument("--iflag", type=int,
            default=3,
            choices=[0, 1, 2, 3, 4],
            help="dimension of the plot. 0: 1D plot of the spherical average, 1: 1D plot, 2: 2D plot, 3: 3D plot, 4: 2D polar plot on a sphere")

    gp.add_argument("--output-format", type=int, default=5,
            choices=[0, 1, 2, 3, 4, 5, 6, 7],
            help="output file format for visualization. 0: gnuplot(1D), 1: no longer supported, 2: plotrho(2D), 3: XCRYSDEN(2d), 4: no longer supported, 5: XCRYSDEN(3D), 6: gaussian cube(3D), 7: gnuplot(2D)")

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

    # fix atoms
    gp = subparser.add_argument_group(title="fix atoms",
            description="specify atoms to fix in optimization, only used when --runtype=1")

    gp.add_argument("--fix", help="list of fixed atoms, index start from 1, have privilege over --fix-around-z", nargs='+', type=int)

    gp.add_argument("--fix-around-z", type=float, nargs=3, default=None,
        help="select atoms around specified z in Angstrom with tolerance, like this --fix-around-z 10 -0.5 0.5")            
            
    gp.add_argument("--color-fix", type=str, default="white",
        choices=["red", "green", "blue", "white"],
        help="select color to color the fixed atoms in xsd file, can be: red green blue and white")
        
    # inp template
    gp = subparser.add_argument_group(title="template input",
            description="read template input")
            
    gp.add_argument("--pwin", type=str, default=None,
            help="read parameters from pwscf input template")

    gp.add_argument("--nebin", type=str, default=None,
            help="read parameters from neb.x input template")
            
    gp.add_argument("--phin", type=str, default=None,
            help="read parameters from ph.x input template")            

    # converge test
    gp = subparser.add_argument_group(title="converge test",
            description="converge test for ecutrho ecutwfc degauss kpoints")
    
    gp.add_argument("--converge", type=str, default="ecutwfc",
            choices=["ecutwfc", "ecutrho", "degauss", "kpoints"],
            help="choose what to do converge test, ecutwfc or ecutrho or degauss, or kpoints")
            
    gp.add_argument("--ecutwfc-range", type=int, nargs=3,
            help="test range for ecutwfc, like --ecutwfc-range 30 120 10")

    gp.add_argument("--ecutrho-range", type=int, nargs=3,
            help="test range for ecutrho, like --ecutrho-range 600 1000 50")

    gp.add_argument("--degauss-range", type=float, nargs=3,
            help="test range for degauss, like --degauss-range 0.001 0.02 0.001")

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


    # --------------------------------------------------------------------------
    # VASP
    # --------------------------------------------------------------------------
    subparser = subparsers.add_parser("vasp", help="using vasp as calculator")

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
            choices=["pbs", "llhpc", "lsf_sz", "tianhe2"],
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
            help="specify the path to the POTCAR, default is ./. if you pass 'auto' to it, matflow will build the POTCAR foryou(need simple configuration, see manual)")
        
    gp.add_argument("--pot-type", type=str, default="PAW_PBE",
            choices=["PAW_PBE", "PAW_LDA", "PAW_PW91", "paw_pbe", "paw_lda", "paw_pw91"],
            help="choose type of POT for POTCAR")

    # --------------------------------------------------------
    #                   INCAR PARAMETERS
    # --------------------------------------------------------
    # incar->start parameters
    gp = subparser.add_argument_group(title="incar->start parameters",
            description="start parameters to be set in INCAR")

    gp.add_argument("--nwrite", type=int, default=None,
            help=" This flag determines how much will be written to the file OUTCAR (verbosity flag)")

    gp.add_argument("--prec", type=str, default=None,
            choices=["Normal", "Accurate", "A", "N", "Low", "L","" "Single"],
            help="PREC, default value: Normal")

    gp.add_argument("--ncore", type=int, default=None,
            help="NCORE determines the number of compute cores that work on an individual orbital ")

    # incar->electrons
    gp = subparser.add_argument_group(title="incar->electron",
            description="electrons calculation related parameters")

    gp.add_argument("--encut", type=int, default=None,
            help="ENCUT, default value: 300 eV")

    gp.add_argument("--ediff", type=float, default=None,
            help="EDIFF, default value: 1.0e-4")

    gp.add_argument("--nelm", type=int, default=None,
            help="NELM sets the maximum number of electronic SC (selfconsistency) steps which may be performed")

    gp.add_argument("--nfree", type=int, default=None,
            help="NFREE specifies the number of remembered steps in the history of ionic convergence runs, or the number of ionic displacements in frozen phonon calculations")

    gp.add_argument("--kpoints-mp", type=int, nargs=6,
            default=[1, 1, 1, 0, 0, 0],
            help="set kpoints like -k 1 1 1 0 0 0")

    #gp.add_argument("--kpoints-mp-scf", type=int, nargs="+",
    #        default=[1, 1, 1, 0, 0, 0],
    #        help="set kpoints like -k 1 1 1 0 0 0")

    gp.add_argument("--kpoints-mp-nscf", type=int, nargs="+",
            default=None,
            help="set kpoints like -k 1 1 1 0 0 0")

    gp.add_argument("--kpath-manual", type=str, nargs="+", default=None,
            help="set kpoints for band structure calculation manually")

    gp.add_argument("--kpath-file", type=str, default="kpath.txt",
            help="set kpoints for band structure calculation manually from file")

    gp.add_argument("--kpath-intersections", type=int, default=15,
            help="intersection of the line mode kpoint for band calculation")

    gp.add_argument("--ismear", type=int, default=None,
            help="smearing type(methfessel-paxton(>0), gaussian(0), fermi-dirac(-1), tetra(-4), tetra-bloch-dorrected(-5)), default: 0")

    gp.add_argument("--sigma", type=float, default=None,
            help="determines the width of the smearing in eV.")

    gp.add_argument("--ivdw", type=int, default=None,
            choices=[0, 11, 12, 21, 202, 4],
            help="IVDW = 0(no correction), 1(dft-d2), 11(dft-d3 Grimme), 12(dft-d3 Becke-Jonson), default: None which means 0, no correction")
    # -----------------------------------------------------------------

    gp.add_argument("--lorbit", type=int, default=None,
            choices=[0, 1, 2, 5, 10, 11, 12],
            help="together with an appropriate RWIGS, determines whether the PROCAR or PROOUT files are written")

    # optics related
    gp.add_argument("--loptics", type=str, default=None,
            choices=["TRUE", "FALSE"],
            help="calculates the frequency dependent dielectric matrix after the electronic ground state has been determined.")

    gp.add_argument("--cshift", type=float, default=None,
            help="CSHIFT sets the (small) complex shift η in the Kramers-Kronig transformation")

    gp.add_argument("--nedos", type=int, default=None,
            help="NEDOS specifies number of gridpoints on which the DOS is evaluated")
            
    # magnetic related
    gp.add_argument("--ispin", type=int, default=None,
            choices=[1, 2],
            help="specifies spin polarization: 1->no spin polarized, 2->spin polarized(collinear). combine SIPIN with MAGMOM to study collinear magnetism.")

    gp.add_argument("--magmom", type=float, nargs="+", default=None,
            help="Specifies the initial magnetic moment for each atom, if and only if ICHARG=2, or if ICHARG=1 and the CHGCAR file contains no magnetisation density.")

    gp.add_argument("--lnoncollinear", type=str, default=None,
            choices=["T", "F", ".TRUE.", ".FALSE."],
            help="specifies whether fully non-collinear magnetic calculations are performed")

    gp.add_argument("--lsorbit", type=str, default=None,
            choices=["T", "F", ".TRUE.", ".FALSE."],
            help="specifies whether spin-orbit coupling is taken into account.")

    gp.add_argument("--saxis", type=float, nargs=3, default=None,
            help="SAXIS specifies the quantisation axis for noncollinear spins")

    gp.add_argument("--lmaxmix", type=int, default=None,
            help="LMAXMIX controls up to which l-quantum number the one-center PAW charge densities are passed through the charge density mixer and written to the CHGCAR file.")

    # hybrid functional
    gp = subparser.add_argument_group(title="incar->Exchange correlation")
    gp.add_argument("--lhfcalc", type=str, default=None,
            choices=["T", "F", ".TRUE.", ".FALSE."],
            help=" specifies whether Hartree-Fock/DFT hybrid functional type calculations are performed")

    gp.add_argument("--hfscreen", type=float, default=None,
            choices=[0.3, 0.2],
            help=" specifies the range-separation parameter in range separated hybrid functionals: HSE03->0.3, HSE06->0.2, must also set LHFCALC=.TRUE.")

    gp.add_argument("--aexx", type=float, default=None,
            help="AEXX specifies the fraction of exact exchange in a Hartree-Fock/DFT hybrid functional type calculation")

    gp.add_argument("--lsubrot", type=str, default=None,
            choices=["T", "F", ".TRUE.", ".FALSE."],
            help="This flag can be set for hybrid functionals (HF-type calculations).")

    gp.add_argument("--nsw", type=int, default=None,
            help="NSW sets the maximum number of ionic steps")

    gp.add_argument("--ediffg", type=float, default=None,
            help="EDIFFG, default value: 10*EDIFF")

    gp = subparser.add_argument_group(title="incar->ions",
            description="setting ions related parameters")

    gp.add_argument("--ibrion", type=int, default=None,
            choices=[-1, 0, 1, 2, 3, 5, 6, 7, 8, 44],
            help="IBRION = refer to https://cms.mpi.univie.ac.at/wiki/index.php/IBRION for how to set the algorithm of optimization you need!")

    gp.add_argument("--isif", type=int, default=None,
            choices=[0, 1, 2, 3, 4, 5, 6, 7],
            help="ISIF = 0-7: refer to https://cms.mpi.univie.ac.at/wiki/index.php/ISIF for how to set the type of Geometri Optimization you need!")

    gp.add_argument("--potim", type=float, default=None,
            help="step width scaling (ionic relaxations), default: None = 0.015 in phonon calculation")

    gp.add_argument("--selective-dynamics", type=str, default="False",
            choices=["True", "False", "T", "F"],
            help="whether use selective dyanmics")

    gp = subparser.add_argument_group(title="molecular dynamics",
            description="molecular dynamics related")
            
    gp.add_argument("--smass", type=int, default=None,
            help="controls the velocities during an ab-initio molecular dynamics run.")
            
    gp.add_argument("--mdalgo", type=int, default=None,
            choices=[0, 1, 2, 3, 11, 21, 13],
            help="specifies the molecular dynamics simulation protocol (in case IBRION=0 and VASP was compiled with -Dtbdyn).")

    gp.add_argument("--anderson-prob", type=float, default=None,
            help=" sets the collision probability for the Anderson thermostat (in case VASP was compiled with -Dtbdyn).")
            
    gp.add_argument("--tebeg", type=float, default=None,
            help=" sets the start temperature for an ab-initio molecular dynamics run (IBRION=0) and other routines (e.g. Electron-phonon interactions from Monte-Carlo sampling).")
            
    gp.add_argument("--teend", type=float, default=None,
            help="sets the final temperature for an ab-initio molecular dynamics run (IBRION=0; SMASS=−1).")            

    # incar-miscellaneous
    gp = subparser.add_argument_group(title="incar-miscellaneous",
            description="miscellaneous input parameters")

    gp.add_argument("--algo", type=str, default=None,
            choices=["N", "D", "V", "F", "VeryFast"],  #"Exact", "G0W0", "GW0", "GW"],
            help=" a convenient option to specify the electronic minimisation algorithm (as of VASP.4.5) and/or to select the type of GW calculations")

    gp.add_argument("--ialgo", type=int, default=None,
            choices=[5, 6, 7, 8, 38, 44, 46, 48],
            help="IALGO selects the algorithm used to optimize the orbitals.Mind: We strongly urge the users to select the algorithms via ALGO. Algorithms other than those available via ALGO are subject to instabilities.")

    gp.add_argument("--addgrid", type=str, default=None,
            choices=[".TRUE.", ".FALSE.", "T", "F"],
            help="ADDGRID determines whether an additional support grid is used for the evaluation of the augmentation charges.")

    gp.add_argument("--isym", type=int, default=None,
            choices=[-1, 0, 1, 2, 3],
            help=" ISYM determines the way VASP treats symmetry.")

    gp.add_argument('--lreal', type=str, default=None,
            choices=["T", "F", ".TRUE.", ".FALSE.", "O", "On", "A", "Auto"],
            help="LREAL determines whether the projection operators are evaluated in real-space or in reciprocal space.")

    gp.add_argument("--pstress", type=float, default=None,
            help="controls whether Pulay corrections are added to the stress tensor or not.")

    # properties parameters
    gp.add_argument("--lelf", type=str, default=None,
            choices=["T", "F", ".TRUE.", ".FALSE."],
            help="LELF determines whether to create an ELFCAR file or not.")

    # write PARAMETERS
    gp = subparser.add_argument_group(title="incar->write parameters",
            description="set writing parameters")

    gp.add_argument("--lwave", type=str, default=None,
            choices=['T', 'F', ".TRUE.", '.FALSE.'],
            help="LWAVE determines whether the wavefunctions are written to the WAVECAR file at the end of a run.")

    gp.add_argument("--lcharg", type=str, default=None,
            choices=['T', 'F', ".TRUE.", '.FALSE.'],
            help="LCHARG determines whether the charge densities (files CHGCAR and CHG) are written.")


    #                     neb related PARAMETERS
    # --------------------------------------------------------------------------
    gp = subparser.add_argument_group(title="incar->neb",
            description="nudged elastic band related setting")

    gp.add_argument("--iopt", type=int, default=None,
            choices=[0, 1, 2, 3],
            help="chioce for optimizer: 0->vasp, 1, 2->vtst")

    gp.add_argument("--lclimb", type=str, default=None,
            choices=["T", "F"],
            help="whether use climbing image")

    gp.add_argument("--lnebcell", type=str, default=None,
            choices=["T", "F"],
            help="flag to turn on SS-NEB, used with ISIF=3 and IOPT=3")

    gp.add_argument("--spring", type=int, default=None,
            help="gives the spring constant between the images as used in the elastic band method")

    gp.add_argument("--maxmove", type=float, default=None,
            help="maximum allowed step size for translation, default is None which means 0.2")

    gp.add_argument("--lglobal", type=str, default=None,
            choices=["T", "F"],
            help="optimize the NEB globally instead of image-by-image, default is None which means .TRUE.")

    gp.add_argument("--lautoscale", type=str, default=None,
            choices=["T", "F"],
            help="automatically determines INVCURV, default is T")

    gp.add_argument("--invcurv", type=float, default=None,
            help="initial inverse curvature, used to construct the inverse Hessian matrix. default is None which means 0.01")

    gp.add_argument("--llineopt", type=str, default=None,
            choices=["T", "F"],
            help="use a force based line minimizer for translation. default is None(means F)")

    gp.add_argument("--fdstep", type=float, default=None,
            help="finite difference step size for line optimizer, default is None(5E-3)")

    gp.add_argument("--nimage", type=int, default=None,
            help="number of image to interpolate. total image will be nimage+2.")

    gp.add_argument("--outcars", type=str, nargs="+",
            help="OUTCAR for the initial and final structure in neb calc")

    gp.add_argument("--nebmake", type=int, default=0,
            choices=[0, 1],
            help="0(default): use nebmake.pl, 1: use nebmake.py")
    
    gp.add_argument("--moving-atom", type=int, nargs="+",
            help="spedifying moving atoms, only needed when --nebmake=1 using nebmake.py. index start from 0")

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

    # incar template
    gp = subparser.add_argument_group(title="template", 
            description="read in INCAR template")

    gp.add_argument("--incar", type=str, default=None,
            help="specify incar template to set parameters")

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

    # static calc related setting
    gp = subparser.add_argument_group(title="static calc",
            description="setting type of static calculation when -r 0")

    gp.add_argument("--static", type=str, default="band",
            choices=["scf", "band", "dos", "optics", "bse", "parchg","stm"],
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

    # VASP PARCHG STM
    gp = subparser.add_argument_group(title="PARCHG(STM) related",
            description="PARCHG(STM) calc related parameters")
            
    gp.add_argument("--lpard", type=str, default=None,
            choices=[".TRUE.", "T", "F", ".FALSE."],
            help="Determines whether partial (band or k-point decomposed) charge densities are evaluated.")
            
    gp.add_argument("--lsepk", type=str, default=None,
            choices=[".TRUE.", "T", "F", ".FALSE."],
            help="Specifies whether the charge density of every k-point is write to the files PARCHG.*.nk (LSEPK=.TRUE.) or whether it is merged to a single file.")
            
    gp.add_argument("--lsepb", type=str, default=None,
            choices=[".TRUE.", "T", "F", ".FALSE."],
            help="Specifies whether the charge density is calculated for every band separately and written to a file PARCHG.nb.* (LSEPB=.TRUE.) or whether charge density is merged for all selected bands and written to the files PARCHG.ALLB.* or PARCHG.")

    gp.add_argument("--nbmod", type=int, default=None,
            help="Controls which bands are used in the calculation of Band decomposed charge densities.")
            
    gp.add_argument("--eint", type=float, nargs=2,
            help="Specifies the energy range of the bands that are used for the evaluation of the partial charge density needed in Band decomposed charge densities.")
            
    # miscellaneous
    gp = subparser.add_argument_group(title="miscellaneous",
            description="miscallaneous setting")
    
    gp.add_argument("--symprec", type=float, default=None,
            help="determines how accurately the positions in the POSCAR file must be specified. The default, SYMPREC=10-5, is usually large enough, even if the POSCAR file has been generated with single precision accuracy. Increasing SYMPREC means that the positions in the POSCAR file can be specified with less accuracy (increasing fuzziness).")

    gp.add_argument("--amix", type=float, default=None,
            help="specifies the linear mixing parameter.")

    gp.add_argument("--bmix", type=float, default=None,
            help="sets the cutoff wave vector for Kerker mixing scheme")            

    gp.add_argument("--nelect", type=int, default=None,
            help="sets the number of electrons")

    gp.add_argument("--laechg", type=str, default=None,
            choices=[".TRUE.", ".FALSE.", "T", "F"],
            help="for bader analysis. when LAECHG=.TRUE. the all-electron charge density will be reconstructed explicitly and written out to file.")
            
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




    # dealing with pseudo potential file
    if args.pot == "./":
        #TODO make a simple check, whether there exists the potential file
        pass
    elif args.pot == "auto":
        if args.driver == "abinit":
            os.system("pot-from-xyz-modified.py -i %s -d ./ -p abinit --abinit-type=ncpp" % xyzfile)
        elif args.driver == "qe":
            if args.runtype == 6:
                os.system("pot-from-xyz-modified.py -i %s -d ./ -p qe --qe-type=%s" % (images[0]), args.pot_type)
            else:
                os.system("pot-from-xyz-modified.py -i %s -d ./ -p qe --qe-type=%s" % (xyzfile, args.pot_type))
        elif args.driver == "siesta":
            print("=============================================================\n")
            print("                     WARNING\n")
            print("-------------------------------------------------------------\n")
            print("support for auto preparation of pseudopotential file for siesta\n")
            print("is not fully implemented now!\n")
            print("please prepare it yourself\n")
            sys.exit(1)
        elif args.driver == "vasp":
            os.system("vasp-potcar-from-xyz.py --type %s -i %s -o ./POTCAR" % (args.pot_type, xyzfile))
    else:
        os.system("cp %s/* ./" % args.pot)


    # server
    # xxx.set_run can only deal with pbs, llhpc, lsf_sz server now 
    # however both guangzhou chaosuan llhpc are build on tianhe2, so they can use the same job system(yhbatch...)
    # we add tianhe2 option to args.server which cannot be handled by xxx.set_run. so we convert it to llhpc if tianhe2 is chosen
    server = args.server if args.server != "tianhe2" else "llhpc"

    if args.driver == "abinit":
        params = {}
        kpoints = {}

        params["chkprim"] = args.chkprim
        params["ecut"] = args.ecut
        params["ixc"] = args.ixc
        params["vdw_xc"] = args.vdw_xc
        params["vdw_tol"] = args.vdw_tol

        kpoints["kptopt"] = args.kptopt
        kpoints["ngkpt"] = args.ngkpt

        params["occopt"] = args.occopt
        params["nband"] = args.nband
        params["occ"] = args.occ


        if args.runtype == 0:
            # static
            params["nsppol"] = args.nsppol
            params["prtden"] = args.prtden
            params["prtdos"] = args.prtdos
            from pymatflow.abinit.static import static_run
            task = static_run()
            if get_kpath(args.kpath_manual, args.kpath_file) == None:
                print("================================================\n")
                print("Warning: matflow abinit\n")
                print("in abinit static runing you must provide kpath\n")
                sys.exit(1)
            task.dataset[3].electrons.kpoints.set_band(kpath=get_kpath(args.kpath_manual, args.kpath_file))
            task.get_xyz(xyzfile)
            task.set_params(params=params)
            task.set_kpoints(kpoints=kpoints)
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
            from pymatflow.abinit.opt import opt_run
            task = opt_run()
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
            from pymatflow.abinit.opt import opt_run
            task = opt_run()
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
            from pymatflow.abinit.opt import opt_run
            task = opt_run()
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
            from pymatflow.abinit.opt import opt_run
            task = opt_run()
            task.get_xyz(xyzfile)
            task.set_params(params=params)
            task.set_kpoints(kpoints=kpoints)
            task.set_run(mpi=args.mpi, server=server, jobname=args.jobname, nodes=args.nodes, ppn=args.ppn, queue=args.queue)
            task.set_llhpc(partition=args.partition, nodes=args.nodes, ntask=args.ntask, jobname=args.jobname, stdout=args.stdout, stderr=args.stderr)
            task.tetragonal(directory=args.directory, runopt=args.runopt, auto=args.auto, na=args.na, stepa=args.stepa, nc=args.nc, stepc=args.stepc)
        elif args.runtype == 5:
            # dfpt-elastic-piezo-dielec
            from pymatflow.abinit.dfpt import dfpt_elastic_piezo_dielec
            task = dfpt_elastic_piezo_dielec()
            task.get_xyz(xyzfile)
            task.set_params(params=params)
            task.set_kpoints(kpoints=kpoints)
            task.set_run(mpi=args.mpi, server=server, jobname=args.jobname, nodes=args.nodes, ppn=args.ppn, queue=args.queue)
            task.set_llhpc(partition=args.partition, nodes=args.nodes, ntask=args.ntask, jobname=args.jobname, stdout=args.stdout, stderr=args.stderr)
            task.run(directory=args.directory, runopt=args.runopt, auto=args.auto)
        elif args.runtype == 6:
            # dfpt-phonon
            from pymatflow.abinit.dfpt import dfpt_phonon
            task = dfpt_phonon()
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
            from pymatflow.abinit.phonopy import phonopy_run
            task = phonopy_run()
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
            from pymatflow.abinit.opt import opt_run
            task = opt_run()
            task.get_xyz(xyzfile)
            task.set_params(params=params)
            task.set_kpoints(kpoints=kpoints)
            task.set_run(mpi=args.mpi, server=server, jobname=args.jobname, nodes=args.nodes, ppn=args.ppn, queue=args.queue)
            task.set_llhpc(partition=args.partition, nodes=args.nodes, ntask=args.ntask, jobname=args.jobname, stdout=args.stdout, stderr=args.stderr)
            task.batch_a = args.batch_a     
            task.batch_b = args.batch_b
            task.batch_c = args.batch_c     
            task.abc(directory=args.directory, runopt=args.runopt, auto=args.auto, range_a=args.range_a, range_b=args.range_b, range_c=args.range_c)
        else:
            pass
# ==============================================================================
# CP2K CP2K CP2K CP2K CP2K CP2K CP2K CP2K CP2K CP2K CP2K CP2K CP2K CP2K C2PK CP2K
# ==============================================================================
    elif args.driver == "cp2k":
        from pymatflow.cmd.cp2k_parser import read_inp
        if args.inp == None:
            params = {}
        else:
            params = read_inp(args.inp)
                
        #print(params)
        params["GLOBAL-PRINT_LEVEL"] = args.print_level if "GLOBAL-PRINT_LEVEL" not in params or args.print_level != None else params["GLOBAL-PRINT_LEVEL"]
        
        params["FORCE_EVAL-SUBSYS-CELL-PERIODIC"] = args.cell_periodic if "FORCE_EVAL-SUBSYS-CELL-PERIODIC" not in params or args.cell_periodic != None else params["FORCE_EVAL-SUBSYS-CELL-PERIODIC"]
        params["FORCE_EVAL-SUBSYS-CELL-SYMMETRY"] = args.cell_symmetry if "FORCE_EVAL-SUBSYS-CELL-SYMMETRY" not in params or args.cell_symmetry != None else params["FORCE_EVAL-SUBSYS-CELL-SYMMETRY"]

        params["FORCE_EVAL-DFT-LS_SCF"] = args.ls_scf if "FORCE_EVAL-DFT-LS_SCF" not in params or args.ls_scf != None else params["FORCE_EVAL-DFT-LS_SCF"]
        params["FORCE_EVAL-DFT-LSD"] = args.lsd if "FORCE_EVAL-DFT-LSD" not in params or args.lsd != None else params["FORCE_EVAL-DFT-LSD"]
        params["FORCE_EVAL-DFT-CHARGE"] = args.charge if "FORCE_EVAL-DFT-CHARGE" not in params or args.charge != None else params["FORCE_EVAL-DFT-CHARGE"]
        params["FORCE_EVAL-DFT-SURFACE_DIPOLE_CORRECTION"] = args.surface_dipole_correction if "FORCE_EVAL-DFT-SURFACE_DIPOLE_CORRECTION" not in params or args.surface_dipole_correction != None else params["FORCE_EVAL-DFT-SURFACE_DIPOLE_CORRECTION"]
        params["FORCE_EVAL-DFT-SURF_DIP_DIR"] = args.surf_dip_dir if "FORCE_EVAL-DFT-SURF_DIP_DIR" not in params or args.surf_dip_dir != None else params["FORCE_EVAL-DFT-SURF_DIP_DIR"]        
        params["FORCE_EVAL-DFT-POISSON-PERIODIC"] = args.poisson_periodic if "FORCE_EVAL-DFT-POISSON-PERIODIC" not in params or args.poisson_periodic != None else params["FORCE_EVAL-DFT-POISSON-PERIODIC"]
        params["FORCE_EVAL-DFT-POISSON-POISSON_SOLVER"] = args.poisson_solver if "FORCE_EVAL-DFT-POISSON-POISSON_SOLVER" not in params or args.poisson_solver != None else params["FORCE_EVAL-DFT-POISSON-POISSON_SOLVER"]
        params["FORCE_EVAL-DFT-QS-METHOD"] = args.qs_method if "FORCE_EVAL-DFT-QS-METHOD" not in params or args.qs_method != None else params["FORCE_EVAL-DFT-QS-METHOD"]
        params["FORCE_EVAL-DFT-MGRID-CUTOFF"] = args.cutoff if "FORCE_EVAL-DFT-MGRID-CUTOFF" not in params or args.cutoff != None else params["FORCE_EVAL-DFT-MGRID-CUTOFF"]
        params["FORCE_EVAL-DFT-MGRID-REL_CUTOFF"] = args.rel_cutoff if "FORCE_EVAL-DFT-MGRID-REL_CUTOFF" not in params or args.rel_cutoff != None else params["FORCE_EVAL-DFT-MGRID-REL_CUTOFF"]
        params["FORCE_EVAL-DFT-MGRID-NGRIDS"] = args.ngrids if "FORCE_EVAL-DFT-MGRID-NGRIDS" not in params or  args.ngrids != None else params["FORCE_EVAL-DFT-MGRID-NGRIDS"]
        params["FORCE_EVAL-DFT-XC-XC_FUNCTIONAL"] = args.xc_functional if "FORCE_EVAL-DFT-XC-XC_FUNCTIONAL" not in params or  args.xc_functional != None else params["FORCE_EVAL-DFT-XC-XC_FUNCTIONAL"]
        params["FORCE_EVAL-DFT-SCF-MAX_SCF"] = args.max_scf if "FORCE_EVAL-DFT-SCF-MAX_SCF" not in params or args.max_scf != None else params["FORCE_EVAL-DFT-SCF-MAX_SCF"]
        params["FORCE_EVAL-DFT-QS-EPS_DEFAULT"] = args.eps_default if "FORCE_EVAL-DFT-QS-EPS_DEFAULT" not in params or args.eps_default != None else params["FORCE_EVAL-DFT-QS-EPS_DEFAULT"]
        params["FORCE_EVAL-DFT-SCF-EPS_SCF"] = args.eps_scf if "FORCE_EVAL-DFT-SCF-EPS_SCF" not in params or  args.eps_scf != None else params["FORCE_EVAL-DFT-SCF-EPS_SCF"]
        params["FORCE_EVAL-DFT-SCF-LEVEL_SHIFT"] = args.level_shift if "FORCE_EVAL-DFT-SCF-LEVEL_SHIFT" not in params or  args.level_shift != None else params["FORCE_EVAL-DFT-SCF-LEVEL_SHIFT"]
        params["FORCE_EVAL-DFT-SCF-MAX_SCF_HISTORY"] = args.max_scf_history if "FORCE_EVAL-DFT-SCF-MAX_SCF_HISTORY" not in params or args.max_scf_history != None else params["FORCE_EVAL-DFT-SCF-MAX_SCF_HISTORY"]
        params["FORCE_EVAL-DFT-SCF-MAX_DIIS"] = args.max_diis if "FORCE_EVAL-DFT-SCF-MAX_DIIS" not in params or args.max_diis != None else params["FORCE_EVAL-DFT-SCF-MAX_DIIS"]
        params["FORCE_EVAL-DFT-SCF-ADDED_MOS"] = args.added_mos if "FORCE_EVAL-DFT-SCF-ADDED_MOS" not in params or  args.added_mos != None else params["FORCE_EVAL-DFT-SCF-ADDED_MOS"]
       

        params["FORCE_EVAL-DFT-SCF-OUTER_SCF"] = args.outer_scf if "FORCE_EVAL-DFT-SCF-OUTER_SCF" not in params or  args.outer_scf != None else params["FORCE_EVAL-DFT-SCF-OUTER_SCF"]
        params["FORCE_EVAL-DFT-SCF-OUTER_SCF-BISECT_TRUST_COUNT"] = args.outer_scf_bisect_trust_count if "FORCE_EVAL-DFT-SCF-OUTER_SCF-BISECT_TRUST_COUNT" not in params or  args.outer_scf_bisect_trust_count != None else params["FORCE_EVAL-DFT-SCF-OUTER_SCF-BISECT_TRUST_COUNT"]
        params["FORCE_EVAL-DFT-SCF-OUTER_SCF-DIIS_BUFFER_LENGTH"] = args.outer_scf_diis_buffer_length if "FORCE_EVAL-DFT-SCF-OUTER_SCF-DIIS_BUFFER_LENGTH" not in params or  args.outer_scf_diis_buffer_length != None else params["FORCE_EVAL-DFT-SCF-OUTER_SCF-DIIS_BUFFER_LENGTH"]
        params["FORCE_EVAL-DFT-SCF-OUTER_SCF-EXTRAPOLATION_ORDER"] = args.outer_scf_extrapolation_order if "FORCE_EVAL-DFT-SCF-OUTER_SCF-EXTRAPOLATION_ORDER" not in params or  args.outer_scf_extrapolation_order != None else params["FORCE_EVAL-DFT-SCF-OUTER_SCF-EXTRAPOLATION_ORDER"]
        params["FORCE_EVAL-DFT-SCF-OUTER_SCF-EPS_SCF"] = args.outer_scf_eps_scf if "FORCE_EVAL-DFT-SCF-OUTER_SCF-EPS_SCF" not in params or  args.outer_scf_eps_scf != None else params["FORCE_EVAL-DFT-SCF-OUTER_SCF-EPS_SCF"]
        params["FORCE_EVAL-DFT-SCF-OUTER_SCF-MAX_SCF"] = args.outer_scf_max_scf if "FORCE_EVAL-DFT-SCF-OUTER_SCF-MAX_SCF" not in params or  args.outer_scf_max_scf != None else params["FORCE_EVAL-DFT-SCF-OUTER_SCF-MAX_SCF"]
        params["FORCE_EVAL-DFT-SCF-OUTER_SCF-OPTIMIZER"] = args.outer_scf_optimizer if "FORCE_EVAL-DFT-SCF-OUTER_SCF-OPTIMIZER" not in params or  args.outer_scf_optimizer != None else params["FORCE_EVAL-DFT-SCF-OUTER_SCF-OPTIMIZER"]
        params["FORCE_EVAL-DFT-SCF-OUTER_SCF-TYPE"] = args.outer_scf_type if "FORCE_EVAL-DFT-SCF-OUTER_SCF-TYPE" not in params or  args.outer_scf_type != None else params["FORCE_EVAL-DFT-SCF-OUTER_SCF-TYPE"]

        params["FORCE_EVAL-DFT-SCF-SMEAR"] = args.smear if "FORCE_EVAL-DFT-SCF-SMEAR" not in params or  args.smear != None else params["FORCE_EVAL-DFT-SCF-SMEAR"]
        params["FORCE_EVAL-DFT-SCF-SMEAR-METHOD"] = args.smear_method if "FORCE_EVAL-DFT-SCF-SMEAR-METHOD" not in params or args.smear_method != None else params["FORCE_EVAL-DFT-SCF-SMEAR-METHOD"]
        params["FORCE_EVAL-DFT-SCF-SMEAR-ELECTRONIC_TEMPERATURE"] = args.electronic_temp if "FORCE_EVAL-DFT-SCF-SMEAR-ELECTRONIC_TEMPERATURE" not in params or  args.electronic_temp != None else params["FORCE_EVAL-DFT-SCF-SMEAR-ELECTRONIC_TEMPERATURE"]
        params["FORCE_EVAL-DFT-SCF-SMEAR-EPS_FERMI_DIRAC"] = args.eps_fermi_dirac if "FORCE_EVAL-DFT-SCF-SMEAR-EPS_FERMI_DIRAC" not in params or  args.esp_fermi_dirac != None else params["FORCE_EVAL-DFT-SCF-SMEAR-EPS_FERMI_DIRAC"]
        params["FORCE_EVAL-DFT-SCF-SMEAR-WINDOW_SIZE"] = args.window_size if "FORCE_EVAL-DFT-SCF-SMEAR-WINDOW_SIZE" not in params or  args.window_size != None else params["FORCE_EVAL-DFT-SCF-SMEAR-WINDOW_SIZE"]
        
        params["FORCE_EVAL-DFT-SCF-DIAGONALIZATION"] = args.diag if "FORCE_EVAL-DFT-SCF-DIAGONALIZATION" not in params or  args.diag != None else params["FORCE_EVAL-DFT-SCF-DIAGONALIZATION"]
        params["FORCE_EVAL-DFT-SCF-DIAGONALIZATION-ALGORITHM"] = args.diag_algo if "FORCE_EVAL-DFT-SCF-DIAGONALIZATION-ALGORITHM" not in params or  args.diag_algo != None else params["FORCE_EVAL-DFT-SCF-DIAGONALIZATION-ALGORITHM"]
        params["FORCE_EVAL-DFT-SCF-OT"] = args.ot if "FORCE_EVAL-DFT-SCF-OT" not in params or  args.ot != None else params["FORCE_EVAL-DFT-SCF-OT"]
        params["FORCE_EVAL-DFT-SCF-MIXING-METHOD"] = args.mixing_method if "FORCE_EVAL-DFT-SCF-MIXING-METHOD" not in params or  args.mixing_method != None else params["FORCE_EVAL-DFT-SCF-MIXING-METHOD"]
        params["FORCE_EVAL-DFT-SCF-MIXING-ALPHA"] = args.mixing_alpha if "FORCE_EVAL-DFT-SCF-MIXING-ALPHA" not in params or  args.mixing_alpha != None else params["FORCE_EVAL-DFT-SCF-MIXING-ALPHA"]
        params["FORCE_EVAL-DFT-SCF-MIXING-BETA"] = args.mixing_beta if "FORCE_EVAL-DFT-SCF-MIXING-BETA" not in params or  args.mixing_beta != None else params["FORCE_EVAL-DFT-SCF-MIXING-BETA"]
        params["FORCE_EVAL-DFT-SCF-MIXING-NBUFFER"] = args.mixing_nbuffer if "FORCE_EVAL-DFT-SCF-MIXING-NBUFFER" not in params or  args.mixing_nbuffer != None else params["FORCE_EVAL-DFT-SCF-MIXING-NBUFFER"]        
        
        params["FORCE_EVAL-DFT-KPOINTS-SCHEME"] = args.kpoints_scheme if "FORCE_EVAL-DFT-KPOINTS-SCHEME" not in params or  args.kpoints_scheme != None else params["FORCE_EVAL-DFT-KPOINTS-SCHEME"]
        params["FORCE_EVAL-DFT-XC-VDW_POTENTIAL-POTENTIAL_TYPE"] = args.vdw_potential_type if "FORCE_EVAL-DFT-XC-VDW_POTENTIAL-POTENTIAL_TYPE" not in params or  args.vdw_potential_type != None else params["FORCE_EVAL-DFT-XC-VDW_POTENTIAL-POTENTIAL_TYPE"]
        params["FORCE_EVAL-DFT-SCF-OT-PRECONDITIONER"] = args.ot_preconditioner if "FORCE_EVAL-DFT-SCF-OT-PRECONDITIONER" not in params or args.ot_preconditioner != None else params["FORCE_EVAL-DFT-SCF-OT-PRECONDITIONER"]
        params["FORCE_EVAL-DFT-SCF-OT-ENERGY_GAP"] = args.ot_energy_gap if "FORCE_EVAL-DFT-SCF-OT-ENERGY_GAP" not in params or args.ot_energy_gap != None else params["FORCE_EVAL-DFT-SCF-OT-ENERGY_GAP"]
        params["FORCE_EVAL-DFT-SCF-OT-MINIMIZER"] = args.ot_minimizer if "FORCE_EVAL-DFT-SCF-OT-MINIMIZER" not in params or args.ot_minimizer != None else params["FORCE_EVAL-DFT-SCF-OT-MINIMIZER"]
        params["FORCE_EVAL-DFT-SCF-OT-ALGORITHM"] = args.ot_algorithm if "FORCE_EVAL-DFT-SCF-OT-ALGORITHM" not in params or args.ot_algorithm != None else params["FORCE_EVAL-DFT-SCF-OT-ALGORITHM"]
        params["FORCE_EVAL-DFT-SCF-OT-LINESEARCH"] = args.ot_linesearch if "FORCE_EVAL-DFT-SCF-OT-LINESEARCH" not in params or args.ot_linesearch != None else params["FORCE_EVAL-DFT-SCF-OT-LINESEARCH"]
        params["FORCE_EVAL-DFT-XC-VDW_POTENTIAL-PAIR_POTENTIAL-TYPE"] = args.pair_type if "FORCE_EVAL-DFT-XC-VDW_POTENTIAL-PAIR_POTENTIAL-TYPE" not in params or  args.pair_type != None else params["FORCE_EVAL-DFT-XC-VDW_POTENTIAL-PAIR_POTENTIAL-TYPE"]
        params["FORCE_EVAL-DFT-XC-VDW_POTENTIAL-PAIR_POTENTIAL-REFERENCE_FUNCTIONAL"] = args.reference_functional if "FORCE_EVAL-DFT-XC-VDW_POTENTIAL-PAIR_POTENTIAL-REFERENCE_FUNCTIONAL" not in params or  args.reference_functional != None else params["FORCE_EVAL-DFT-XC-VDW_POTENTIAL-PAIR_POTENTIAL-REFERENCE_FUNCTIONAL"]
        params["FORCE_EVAL-DFT-XC-VDW_POTENTIAL-PAIR_POTENTIAL-R_CUTOFF"] = args.r_cutoff if "FORCE_EVAL-DFT-XC-VDW_POTENTIAL-PAIR_POTENTIAL-R_CUTOFF" not in params or  args.r_cutoff != None else params["FORCE_EVAL-DFT-XC-VDW_POTENTIAL-PAIR_POTENTIAL-R_CUTOFF"]
        params["FORCE_EVAL-DFT-PRINT-ELF_CUBE-STRIDE"] = args.dft_print_elf_cube_stride if "FORCE_EVAL-DFT-PRINT-ELF_CUBE-STRIDE" not in params or  args.dft_print_elf_cube_stride != None else params["FORCE_EVAL-DFT-PRINT-ELF_CUBE-STRIDE"]
        params["FORCE_EVAL-DFT-PRINT-E_DENSITY_CUBE-STRIDE"] = args.dft_print_e_density_cube_stride if "FORCE_EVAL-DFT-PRINT-E_DENSITY_CUBE-STRIDE" not in params or  args.dft_print_e_density_cube_stride != None else params["FORCE_EVAL-DFT-PRINT-E_DENSITY_CUBE-STRIDE"]
        params["FORCE_EVAL-PROPERTIES-RESP-SLAB_SAMPLING-RANGE"] = args.properties_resp_slab_sampling_range if "FORCE_EVAL-PROPERTIES-RESP-SLAB_SAMPLING-RANGE" not in params or  args.properties_resp_slab_sampling_range != None else params["FORCE_EVAL-PROPERTIES-RESP-SLAB_SAMPLING-RANGE"]
        params["FORCE_EVAL-PROPERTIES-RESP-SLAB_SAMPLING-SURF_DIRECTION"] = args.properties_resp_slab_sampling_surf_direction if "FORCE_EVAL-PROPERTIES-RESP-SLAB_SAMPLING-SURF_DIRECTION" not in params or  args.properties_resp_slab_sampling_surf_direction != None else params["FORCE_EVAL-PROPERTIES-RESP-SLAB_SAMPLING-SURF_DIRECTION"]
        params["FORCE_EVAL-PROPERTIES-RESP-SLAB_SAMPLING-ATOM_LIST"] = args.properties_resp_slab_sampling_atom_list if "FORCE_EVAL-PROPERTIES-RESP-SLAB_SAMPLING-ATOM_LIST" not in params or  args.properties_resp_slab_sampling_atom_list != None else params["FORCE_EVAL-PROPERTIES-RESP-SLAB_SAMPLING-ATOM_LIST"]
        params["FORCE_EVAL-DFT-PRINT-MOMENTS-PERIODIC"] = args.dft_print_moments_periodic if "FORCE_EVAL-DFT-PRINT-MOMENTS-PERIODIC" not in params or args.dft_print_moments_periodic != None else params["FORCE_EVAL-DFT-PRINT-MOMENTS-PERIODIC"]
        params["FORCE_EVAL-DFT-PRINT-MOMENTS-REFERENCE"] = args.dft_print_moments_reference if "FORCE_EVAL-DFT-PRINT-MOMENTS-REFERENCE" not in params or args.dft_print_moments_reference != None else params["FORCE_EVAL-DFT-PRINT-MOMENTS-REFERENCE"]

        params["MOTION-GEO_OPT-MAX_ITER"] = args.geo_opt_max_iter if "MOTION-GEO_OPT-MAX_ITER" not in params or  args.geo_opt_max_iter != None else params["MOTION-GEO_OPT-MAX_ITER"]
        params["MOTION-GEO_OPT-OPTIMIZER"] = args.geo_opt_optimizer if "MOTION-GEO_OPT-OPTIMIZER" not in params or  args.geo_opt_optimizer != None else params["MOTION-GEO_OPT-OPTIMIZER"]
        params["MOTION-GEO_OPT-TYPE"] = args.geo_opt_type if "MOTION-GEO_OPT-TYPE" not in params or  args.geo_opt_type != None else params["MOTION-GEO_OPT-TYPE"]
        params["MOTION-GEO_OPT-MAX_DR"] = args.geo_opt_max_dr if "MOTION-GEO_OPT-MAX_DR" not in params or  args.geo_opt_max_dr != None else params["MOTION-GEO_OPT-MAX_DR"]
        params["MOTION-GEO_OPT-MAX_FORCE"] = args.geo_opt_max_force if "MOTION-GEO_OPT-MAX_FORCE" not in params or  args.geo_opt_max_force != None else params["MOTION-GEO_OPT-MAX_FORCE"]
        params["MOTION-GEO_OPT-RMS_DR"] = args.geo_opt_rms_dr if "MOTION-GEO_OPT-RMS_DR" not in params or  args.geo_opt_rms_dr != None else params["MOTION-GEO_OPT-RMS_DR"]
        params["MOTION-GEO_OPT-RMS_FORCE"] = args.geo_opt_rms_force if "MOTION-GEO_OPT-RMS_FORCE" not in params or  args.geo_opt_rms_force != None else params["MOTION-GEO_OPT-RMS_FORCE"]

        params["MOTION-CELL_OPT-MAX_ITER"] = args.cell_opt_max_iter if "MOTION-CELL_OPT-MAX_ITER" not in params or  args.cell_opt_max_iter != None else params["MOTION-CELL_OPT-MAX_ITER"]
        params["MOTION-CELL_OPT-OPTIMIZER"] = args.cell_opt_optimizer if "MOTION-CELL_OPT-OPTIMIZER" not in params or  args.cell_opt_optimizer != None else params["MOTION-CELL_OPT-OPTIMIZER"]
        params["MOTION-CELL_OPT-TYPE"] = args.cell_opt_type if "MOTION-CELL_OPT-TYPE" not in params or  args.cell_opt_type != None else params["MOTION-CELL_OPT-TYPE"]
        params["MOTION-CELL_OPT-MAX_DR"] = args.cell_opt_max_dr if "MOTION-CELL_OPT-MAX_DR" not in params or  args.cell_opt_max_dr != None else params["MOTION-CELL_OPT-MAX_DR"]
        params["MOTION-CELL_OPT-MAX_FORCE"] = args.cell_opt_max_force if "MOTION-CELL_OPT-MAX_FORCE" not in params or  args.cell_opt_max_force != None else params["MOTION-CELL_OPT-MAX_FORCE"]
        params["MOTION-CELL_OPT-RMS_DR"] = args.cell_opt_rms_dr if "MOTION-CELL_OPT-RMS_DR" not in params or  args.cell_opt_rms_dr != None else params["MOTION-CELL_OPT-RMS_DR"]
        params["MOTION-CELL_OPT-RMS_FORCE"] = args.cell_opt_rms_force if "MOTION-CELL_OPT-RMS_FORCE" not in params or  args.cell_opt_rms_force != None else params["MOTION-CELL_OPT-RMS_FORCE"]
        params["MOTION-CELL_OPT-PRESSURE_TOLERANCE"] = args.cell_opt_pressure_tolerance if "MOTION-CELL_OPT-PRESSURE_TOLERANCE" not in params or  args.cell_opt_pressure_tolerance != None else params["MOTION-CELL_OPT-PRESSURE_TOLERANCE"]
        params["MOTION-CELL_OPT-KEEP_ANGLES"] = args.cell_opt_keep_angles if "MOTION-CELL_OPT-KEEP_ANGLES" not in params or  args.cell_opt_keep_angles != None else params["MOTION-CELL_OPT-KEEP_ANGLES"]
        params["MOTION-CELL_OPT-KEEP_SYMMETRY"] = args.cell_opt_keep_symmetry if "MOTION-CELL_OPT-KEEP_SYMMETRY" not in params or  args.cell_opt_keep_symmetry != None else params["MOTION-CELL_OPT-KEEP_SYMMETRY"]

        params["MOTION-BAND-BAND_TYPE"] = args.band_type if "MOTION-BAND-BAND_TYPE" not in params or  args.band_type != None else params["MOTION-BAND-BAND_TYPE"]
        params["MOTION-BAND-NUMBER_OF_REPLICA"] = args.number_of_replica if "MOTION-BAND-NUMBER_OF_REPLICA" not in params or  args.number_of_replica != None else params["MOTION-BAND-NUMBER_OF_REPLICA"]
        params["MOTION-BAND-ALIGN_FRAMES"] = args.align_frames if "MOTION-BAND-ALIGN_FRAMES" not in params or  args.align_frames != None else params["MOTION-BAND-ALIGN_FRAMES"]
        params["MOTION-BAND-ROTATE-FRAMES"] = args.rotate_frames if "MOTION-BAND-ROTATE-FRAMES" not in params or  args.rotate_frames != None else params["MOTION-BAND-ROTATE-FRAMES"]
        params["MOTION-BAND-K_SPRING"] = args.k_spring if "MOTION-BAND-K_SPRING" not in params or  args.k_spring != None else params["MOTION-BAND-K_SPRING"]
        params["MOTION-BAND-NPROC_REP"] = args.nproc_rep  if "MOTION-BAND-NPROC_REP" not in params or  args.nproc_rep != None else params["MOTION-BAND-NPROC_REP"]
        params["MOTION-BAND-OPTIMIZE_BAND-OPTIMIZE_END_POINTS"] = args.optimize_end_points if "MOTION-BAND-OPTIMIZE_BAND-OPTIMIZE_END_POINTS" not in params or  args.optimize_end_points != None else params["MOTION-BAND-OPTIMIZE_BAND-OPTIMIZE_END_POINTS"]
        params["MOTION-BAND-CONVERGENCE_CONTROL-MAX_DR"] = args.convergence_control_max_dr if "MOTION-BAND-OPTIMIZE_BAND-CONVERGENCE_CONTROL-MAX_DR" not in params or  args.convergence_control_max_dr != None else params["MOTION-BAND-OPTIMIZE_BAND-CONVERGENCE_CONTROL-MAX_DR"]
        params["MOTION-BAND-CONVERGENCE_CONTROL-MAX_FORCE"] = args.convergence_control_max_force if "MOTION-BAND-OPTIMIZE_BAND-CONVERGENCE_CONTROL-MAX_FORCE" not in params or  args.convergence_control_max_force != None else params["MOTION-BAND-OPTIMIZE_BAND-CONVERGENCE_CONTROL-MAX_FORCE"]
        params["MOTION-BAND-CONVERGENCE_CONTROL-RMS_DR"] = args.convergence_control_rms_dr if "MOTION-BAND-OPTIMIZE_BAND-CONVERGENCE_CONTROL-RMS_DR" not in params or  args.convergence_control_rms_dr != None else params["MOTION-BAND-OPTIMIZE_BAND-CONVERGENCE_CONTROL-RMS_DR"]
        params["MOTION-BAND-CONVERGENCE_CONTROL-RMS_FORCE"] = args.convergence_control_rms_force if "MOTION-BAND-OPTIMIZE_BAND-CONVERGENCE_CONTROL-RMS_FORCE" not in params or  args.convergence_control_rms_force != None else params["MOTION-BAND-OPTIMIZE_BAND-CONVERGENCE_CONTROL-RMS_FORCE"]
        params["MOTION-BAND-OPTIMIZE_BAND-OPT_TYPE"] = args.optimize_band_opt_type if "MOTION-BAND-OPTIMIZE_BAND-OPTIMIZE_BAND-OPT_TYPE" not in params or  args.optimize_band_opt_type != None else params["MOTION-BAND-OPTIMIZE_BAND-OPTIMIZE_BAND-OPT_TYPE"]
        params["MOTION-BAND-OPTIMIZE_BAND-DIIS-MAX_STEPS"] = args.optimize_band_diis_max_steps if "MOTION-BAND-OPTIMIZE_BAND-OPTIMIZE_BAND-DIIS-MAX_STEPS" not in params or  args.optimize_band_diis_max_steps != None else params["MOTION-BAND-OPTIMIZE_BAND-OPTIMIZE_BAND-DIIS-MAX_STEPS"]
        params["MOTION-BAND-OPTIMIZE_BAND-MD-MAX_STEPS"] = args.optimize_band_md_max_steps if "MOTION-BAND-OPTIMIZE_BAND-OPTIMIZE_BAND-MD-MAX_STEPS" not in params or  args.optimize_band_md_max_steps != None else params["MOTION-BAND-OPTIMIZE_BAND-OPTIMIZE_BAND-MD-MAX_STEPS"]


        params["MOTION-MD-STEPS"] = args.md_steps if "MOTION-MD-STEPS" not in params or  args.md_steps != None else params["MOTION-MD-STEPS"]
        params["MOTION-MD-TIMESTEP"] = args.timestep if "MOTION-MD-TIMESTEP" not in params or  args.timestep != None else params["MOTION-MD-TIMESTEP"]
        params["MOTION-MD-ENSEMBLE"] = args.ensemble if "MOTION-MD-ENSEMBLE" not in params or  args.ensemble != None else params["MOTION-MD-ENSEMBLE"]
        params["MOTION-MD-TEMPERATURE"] = args.temperature if "MOTION-MD-TEMPERATURE" not in params or  args.temperature != None else params["MOTION-MD-TEMPERATURE"]
        params["MOTION-MD-TEMP_TOL"] = args.temp_tol if "MOTION-MD-TEMP_TOL" not in params or  args.temp_tol != None else params["MOTION-MD-TEMP_TOL"]
        params["MOTION-PRINT-TRAJECTORY-FORMAT"] = args.traj_format if "MOTION-PRINT-TRAJECTORY-FORMAT" not in params or  args.traj_format != None else params["MOTION-PRINT-TRAJECTORY-FORMAT"]

        params["MOTION-FREE_ENERGY-METADYN-DELTA_T"] = args.metadyn_delta_t if "MOTION-FREE_ENERGY-METADYN-DELTA_T" not in params or args.metadyn_delta_t != None else params["MOTION-FREE_ENERGY-METADYN-DELTA_T"]
        params["MOTION-FREE_ENERGY-METADYN-DO_HILLS"] =  args.metadyn_do_hills if "MOTION-FREE_ENERGY-METADYN-DO_HILLS" not in params or args.metadyn_do_hills != None else params["MOTION-FREE_ENERGY-METADYN-DO_HILLS"]
        params["MOTION-FREE_ENERGY-METADYN-NT_HILLS"] =  args.metadyn_nt_hills if "MOTION-FREE_ENERGY-METADYN-NT_HILLS" not in params or args.metadyn_nt_hills != None else params["MOTION-FREE_ENERGY-METADYN-NT_HILLS"]

        params["VIBRATIONAL_ANALYSIS-DX"] = args.dx if "VIBRATIONAL_ANALYSIS-DX" not in params or  args.dx != None else params["VIBRATIONAL_ANALYSIS-DX"]
        params["VIBRATIONAL_ANALYSIS-FULLY_PERIODIC"] = args.fully_periodic if "VIBRATIONAL_ANALYSIS-FULLY_PERIODIC" not in params or  args.fully_periodic != None else params["VIBRATIONAL_ANALYSIS-FULLY_PERIODIC"]
        params["VIBRATIONAL_ANALYSIS-INTENSITIES"] = args.intensities if "VIBRATIONAL_ANALYSIS-INTENSITIES" not in params or  args.intensities != None else params["VIBRATIONAL_ANALYSIS-INTENSITIES"]
        params["VIBRATIONAL_ANALYSIS-TC_PRESSURE"] = args.tc_pressure if "VIBRATIONAL_ANALYSIS-TC_PRESSURE" not in params or  args.tc_pressure != None else params["VIBRATIONAL_ANALYSIS-TC_PRESSURE"]
        params["VIBRATIONAL_ANALYSIS-TC_TEMPERATURE"] = args.tc_temperature if "VIBRATIONAL_ANALYSIS-TC_TEMPERATURE" not in params or  args.tc_temperature != None else params["VIBRATIONAL_ANALYSIS-TC_TEMPERATURE"]
        params["VIBRATIONAL_ANALYSIS-THERMOCHEMISTRY"] = args.thermochemistry if "VIBRATIONAL_ANALYSIS-THERMOCHEMISTRY" not in params or  args.thermochemistry != None else params["VIBRATIONAL_ANALYSIS-THERMOCHEMISTRY"]
        params["VIBRATIONAL_ANALYSIS-NPROC_REP"] = args.vib_nproc_rep if "VIBRATIONAL_ANALYSIS-NPROC_REP" not in params or  args.vib_nproc_rep != None else params["VIBRATIONAL_ANALYSIS-NPROC_REP"]

        
        # deal with POTENTIAL and BASIS SET
        kind_pot = {}
        if args.kind_pot == None:
            pass
        else:
            n_ele = int(len(args.kind_pot.split()) / 2)
            for i in range(n_ele):
                kind_pot[args.kind_pot.split()[2*i]] = args.kind_pot.split()[2*i+1]
        kind_basis = {}
        if args.kind_basis == None:
            pass
        else:
            n_ele = int(len(args.kind_basis.split()) / 2)
            for i in range(n_ele):
                kind_basis[args.kind_basis.split()[2*i]] = args.kind_basis.split()[2*i+1]
        #
        basis_file = os.path.abspath(args.basis_file) if os.path.exists(args.basis_file) else os.path.basename(args.basis_file)
        pot_file = os.path.abspath(args.pot_file) if os.path.exists(args.pot_file) else os.path.basename(args.pot_file)


        # do some check
        if params["MOTION-CELL_OPT-KEEP_SYMMETRY"] == None:
            pass
        elif params["MOTION-CELL_OPT-KEEP_SYMMETRY"].upper() == ".TRUE." or params["MOTION-CELL_OPT-KEEP_SYMMETRY"].upper() == "TURE":
            if params["FORCE_EVAL-SUBSYS-CELL-SYMMETRY"] == None:
                print("==============================================================\n")
                print("                    WARNING!!!\n")
                print("--------------------------------------------------------------")
                print("you are trying to use KEEP_SYMMETRY in CELL_OPT\n")
                print("in which case you must set an initial symmetry in &CELL\n")
                print("go and check it\n")
                sys.exit(1)
                
        if args.runtype == 0:
            from pymatflow.cp2k.static import static_run
            task = static_run()
            task.get_xyz(xyzfile)
            task.set_params(params=params)
            task.set_pot_basis(kind_basis=kind_basis, kind_pot=kind_pot, basis_set_file=basis_file, potential_file=pot_file)
            task.set_printout(option=args.printout_option)
            if 2 in args.printout_option:
                task.force_eval.dft.printout.band_structure.set_band(kpath=get_kpath(args.kpath_manual, args.kpath_file))
            task.set_vdw(usevdw=True if args.vdw_potential_type.lower() != "none" else False)
            task.set_run(mpi=args.mpi, server=server, jobname=args.jobname, nodes=args.nodes, ppn=args.ppn, queue=args.queue)
            task.set_llhpc(partition=args.partition, nodes=args.nodes, ntask=args.ntask, jobname=args.jobname, stdout=args.stdout, stderr=args.stderr)
            task.scf(directory=args.directory, runopt=args.runopt, auto=args.auto)
        elif args.runtype == 1:
            # geo opt
            from pymatflow.cp2k.opt import opt_run
            task = opt_run()
            task.get_xyz(xyzfile)
            task.set_geo_opt()
            task.set_params(params=params)
            task.set_pot_basis(kind_basis=kind_basis, kind_pot=kind_pot, basis_set_file=basis_file, potential_file=pot_file)
            task.set_vdw(usevdw=True if args.vdw_potential_type.lower() != "none" else False)
            task.set_run(mpi=args.mpi, server=server, jobname=args.jobname, nodes=args.nodes, ppn=args.ppn, queue=args.queue)
            task.set_llhpc(partition=args.partition, nodes=args.nodes, ntask=args.ntask, jobname=args.jobname, stdout=args.stdout, stderr=args.stderr)
            task.geo_opt(directory=args.directory, runopt=args.runopt, auto=args.auto)
        elif args.runtype == 2:
            # cell opt
            from pymatflow.cp2k.opt import opt_run
            task = opt_run()
            task.get_xyz(xyzfile)
            task.set_cell_opt()
            task.set_params(params=params)
            task.set_pot_basis(kind_basis=kind_basis, kind_pot=kind_pot, basis_set_file=basis_file, potential_file=pot_file)
            task.set_vdw(usevdw=True if args.vdw_potential_type.lower() != "none" else False)
            task.set_run(mpi=args.mpi, server=server, jobname=args.jobname, nodes=args.nodes, ppn=args.ppn, queue=args.queue)
            task.set_llhpc(partition=args.partition, nodes=args.nodes, ntask=args.ntask, jobname=args.jobname, stdout=args.stdout, stderr=args.stderr)
            task.cell_opt(directory=args.directory, runopt=args.runopt, auto=args.auto)
        elif args.runtype == 3:
            # cubic cell opt
            from pymatflow.cp2k.opt import opt_run
            task = opt_run()
            task.get_xyz(xyzfile)
            task.set_geo_opt()
            task.set_params(params=params)
            task.set_pot_basis(kind_basis=kind_basis, kind_pot=kind_pot, basis_set_file=basis_file, potential_file=pot_file)
            task.batch_a = args.batch_a     
            task.batch_b = args.batch_b
            task.batch_c = args.batch_c     
            task.set_vdw(usevdw=True if args.vdw_potential_type.lower() != "none" else False)
            task.set_run(mpi=args.mpi, server=server, jobname=args.jobname, nodes=args.nodes, ppn=args.ppn, queue=args.queue)
            task.set_llhpc(partition=args.partition, nodes=args.nodes, ntask=args.ntask, jobname=args.jobname, stdout=args.stdout, stderr=args.stderr)
            #task.cubic(directory=args.directory, runopt=args.runopt, auto=args.auto, na=args.nc, stepa=args.stepa)
            task.cubic(directory=args.directory, runopt=args.runopt, auto=args.auto, range_a=args.range_a)
        elif args.runtype == 4:
            # hexagonal cell opt
            from pymatflow.cp2k.opt import opt_run
            task = opt_run()
            task.get_xyz(xyzfile)
            task.set_geo_opt()
            task.set_params(params=params)
            task.set_pot_basis(kind_basis=kind_basis, kind_pot=kind_pot, basis_set_file=basis_file, potential_file=pot_file)
            task.batch_a = args.batch_a     
            task.batch_b = args.batch_b
            task.batch_c = args.batch_c     
            task.set_vdw(usevdw=True if args.vdw_potential_type.lower() != "none" else False)
            task.set_run(mpi=args.mpi, server=server, jobname=args.jobname, nodes=args.nodes, ppn=args.ppn, queue=args.queue)
            task.set_llhpc(partition=args.partition, nodes=args.nodes, ntask=args.ntask, jobname=args.jobname, stdout=args.stdout, stderr=args.stderr)
            #task.hexagonal(directory=args.directory, runopt=args.runopt, auto=args.auto, na=args.na, nc=args.nc, stepa=args.stepa, stepc=args.stepc)
            task.hexagonal(directory=args.directory, runopt=args.runopt, auto=args.auto, range_a=args.range_a, range_c=args.range_c)
        elif args.runtype == 5:
            # tetragonal cell opt
            from pymatflow.cp2k.opt import opt_run
            task = opt_run()
            task.get_xyz(xyzfile)
            task.set_geo_opt()
            task.set_params(params=params)
            task.set_pot_basis(kind_basis=kind_basis, kind_pot=kind_pot, basis_set_file=basis_file, potential_file=pot_file)
            task.batch_a = args.batch_a     
            task.batch_b = args.batch_b
            task.batch_c = args.batch_c     
            task.set_vdw(usevdw=True if args.vdw_potential_type.lower() != "none" else False)
            task.set_run(mpi=args.mpi, server=server, jobname=args.jobname, nodes=args.nodes, ppn=args.ppn, queue=args.queue)
            task.set_llhpc(partition=args.partition, nodes=args.nodes, ntask=args.ntask, jobname=args.jobname, stdout=args.stdout, stderr=args.stderr)
            #task.tetragonal(directory=args.directory, runopt=args.runopt, auto=args.auto, na=args.na, nc=args.nc, stepa=args.stepa, stepc=args.stepc)
            task.tetragonal(directory=args.directory, runopt=args.runopt, auto=args.auto, range_a=args.range_a, range_c=args.range_c)
        elif args.runtype == 6:
            # neb
            from pymatflow.cp2k.neb import neb_run
            task = neb_run()
            task.get_images(images=images)
            task.set_params(params=params)
            task.check_neb()
            task.set_pot_basis(kind_basis=kind_basis, kind_pot=kind_pot, basis_set_file=basis_file, potential_file=pot_file)
            task.set_vdw(usevdw=True if args.vdw_potential_type.lower() != "none" else False)
            task.set_run(mpi=args.mpi, server=server, jobname=args.jobname, nodes=args.nodes, ppn=args.ppn, queue=args.queue)
            task.set_llhpc(partition=args.partition, nodes=args.nodes, ntask=args.ntask, jobname=args.jobname, stdout=args.stdout, stderr=args.stderr)
            task.neb(directory=args.directory, runopt=args.runopt, auto=args.auto)
        elif args.runtype == 7:
            # phonopy
            from pymatflow.cp2k.phonopy import phonopy_run
            task = phonopy_run()
            task.get_xyz(xyzfile)
            task.supercell_n = args.supercell_n
            task.set_params(params=params)
            task.set_pot_basis(kind_basis=kind_basis, kind_pot=kind_pot, basis_set_file=basis_file, potential_file=pot_file)
            task.set_vdw(usevdw=True if args.vdw_potential_type.lower() != "none" else False)
            task.set_run(mpi=args.mpi, server=server, jobname=args.jobname, nodes=args.nodes, ppn=args.ppn, queue=args.queue)
            task.set_llhpc(partition=args.partition, nodes=args.nodes, ntask=args.ntask, jobname=args.jobname, stdout=args.stdout, stderr=args.stderr)
            task.phonopy(directory=args.directory, runopt=args.runopt, auto=args.auto)
        elif args.runtype == 8:
            # vibrational analysis
            from pymatflow.cp2k.vib import vib_run
            task = vib_run()
            task.get_xyz(xyzfile)
            task.set_printout(option=args.printout_option)
            task.set_params(params=params)
            task.set_pot_basis(kind_basis=kind_basis, kind_pot=kind_pot, basis_set_file=basis_file, potential_file=pot_file)
            task.set_vdw(usevdw=True if args.vdw_potential_type.lower() != "none" else False)
            task.set_run(mpi=args.mpi, server=server, jobname=args.jobname, nodes=args.nodes, ppn=args.ppn, queue=args.queue)
            task.set_llhpc(partition=args.partition, nodes=args.nodes, ntask=args.ntask, jobname=args.jobname, stdout=args.stdout, stderr=args.stderr)
            task.vib(directory=args.directory, runopt=args.runopt, auto=args.auto)
        elif args.runtype == 9:
            # converge test
            from pymatflow.cp2k.static import static_run
            task = static_run()
            task.get_xyz(xyzfile)
            task.set_params(params=params)
            task.set_pot_basis(kind_basis=kind_basis, kind_pot=kind_pot, basis_set_file=basis_file, potential_file=pot_file)
            task.set_vdw(usevdw=True if args.vdw_potential_type.lower() != "none" else False)
            task.set_run(mpi=args.mpi, server=server, jobname=args.jobname, nodes=args.nodes, ppn=args.ppn, queue=args.queue)
            task.set_llhpc(partition=args.partition, nodes=args.nodes, ntask=args.ntask, jobname=args.jobname, stdout=args.stdout, stderr=args.stderr)
            if args.converge.lower() == "cutoff":
                task.converge_cutoff(directory=args.directory, runopt=args.runopt, auto=args.auto, emin=args.cutoff_range[0], emax=args.cutoff_range[1], step=args.cutoff_range[2])
            elif args.converge.lower() == "rel_cutoff":
                task.converge_rel_cutoff(directory=args.directory, runopt=args.runopt, auto=args.auto, emin=args.rel_cutoff_range[0], emax=args.rel_cutoff_range[1], step=args.rel_cutoff_range[2])
            elif args.converge.lower() == "kpoints_auto":
                task.converge_kpoints_auto(directory=args.directory, runopt=args.runopt, auto=args.auto, kmin=args.krange[0], kmax=args.krange[1], step=args.krange[2])
            elif args.converge.lower() == "kpoints_manual":
                kpoints = []
                i = 0
                while i + 2 <= len(args.klist) - 1:
                        kpoints.append([args.klist[i], args.klist[i+1], args.klist[i+2]])
                        i = i + 3
                task.converge_kpoints_manual(directory=args.directory, runopt=args.runopt, auto=args.auto, kpoints_list=kpoints)
        elif args.runtype == 10:
            # aimd
            from pymatflow.cp2k.md import md_run
            task = md_run()
            task.get_xyz(xyzfile)
            task.set_printout(option=args.printout_option)
            task.set_params(params=params)
            task.set_pot_basis(kind_basis=kind_basis, kind_pot=kind_pot, basis_set_file=basis_file, potential_file=pot_file)
            task.set_vdw(usevdw=True if args.vdw_potential_type.lower() != "none" else False)
            task.set_run(mpi=args.mpi, server=server, jobname=args.jobname, nodes=args.nodes, ppn=args.ppn, queue=args.queue)
            task.set_llhpc(partition=args.partition, nodes=args.nodes, ntask=args.ntask, jobname=args.jobname, stdout=args.stdout, stderr=args.stderr)
            task.aimd(directory=args.directory, runopt=args.runopt, auto=args.auto)
        elif args.runtype == 11:
            # abc cell opt
            from pymatflow.cp2k.opt import opt_run
            task = opt_run()
            task.get_xyz(xyzfile)
            task.set_geo_opt()
            task.set_params(params=params)
            task.set_pot_basis(kind_basis=kind_basis, kind_pot=kind_pot, basis_set_file=basis_file, potential_file=pot_file)
            task.set_vdw(usevdw=True if args.vdw_potential_type.lower() != "none" else False)
            task.set_run(mpi=args.mpi, server=server, jobname=args.jobname, nodes=args.nodes, ppn=args.ppn, queue=args.queue)
            task.set_llhpc(partition=args.partition, nodes=args.nodes, ntask=args.ntask, jobname=args.jobname, stdout=args.stdout, stderr=args.stderr)
            task.batch_a = args.batch_a     
            task.batch_b = args.batch_b
            task.batch_c = args.batch_c     
            task.abc(directory=args.directory, runopt=args.runopt, auto=args.auto, range_a=args.range_a, range_b=args.range_b, range_c=args.range_c)
        elif args.runtype == 12:
            # metadynamics
            from pymatflow.cp2k.md import md_run
            task = md_run()
            task.get_xyz(xyzfile)
            task.set_params(params=params)
            task.set_pot_basis(kind_basis=kind_basis, kind_pot=kind_pot, basis_set_file=basis_file, potential_file=pot_file)
            task.set_vdw(usevdw=True if args.vdw_potential_type.lower() != "none" else False)
            task.set_run(mpi=args.mpi, server=server, jobname=args.jobname, nodes=args.nodes, ppn=args.ppn, queue=args.queue)
            task.set_llhpc(partition=args.partition, nodes=args.nodes, ntask=args.ntask, jobname=args.jobname, stdout=args.stdout, stderr=args.stderr)
            task.metadynamics(directory=args.directory, runopt=args.runopt, auto=args.auto)
        else:                
            pass
# ====================================================================================
# Quantum ESPERSSO Quantum ESPERSSO Quantum ESPERSSO Quantum ESPERSSO Quantum ESPERSSO
# ====================================================================================
    elif args.driver == "qe":
        from pymatflow.cmd.qe_parser import read_pwscf_in, read_neb_in, read_ph_in
        control = {}
        electrons = {}
        system = {}
        ions = {}

        if args.pwin != None:
            control, electrons, system, ions, cell = read_pwscf_in(args.pwin)
        
        control["tstress"] = args.tstress if "tstress" not in control or args.tstress != None else control["tstress"]
        control["nstep"] = args.nstep if "nstep" not in control or args.nstep != None else control["nstep"]
        control["etot_conv_thr"] = args.etot_conv_thr if "etot_conv_thr" not in control or args.etot_conv_thr != None else control["etot_conv_thr"]
        control["forc_conv_thr"] = args.forc_conv_thr if "forc_conv_thr" not in control or args.forc_conv_thr != None else control["forc_conv_thr"]
        system["ecutwfc"] = args.ecutwfc if "ecutwfc" not in system or args.ecutwfc != None else system["ecutwfc"]
        system["ecutrho"] = args.ecutrho if "ecutrho" not in system or args.ecutrho != None else system["ecutrho"]
        system["occupations"] = args.occupations if "occupations" not in system or args.occupations != None else system["occupations"]
        system["smearing"] = args.smearing if "smearing" not in system or args.smearing != None else system["smearing"]
        system["degauss"] = args.degauss if "degauss" not in system or args.degauss != None else system["degauss"]
        system["vdw_corr"] = args.vdw_corr if "vde_corr" not in system or args.vdw_corr != None  else system["vdw_corr"]
        system["nbnd"] = args.nbnd if "nbnd" not in system or args.nbnd != None else system["nbnd"]
        system["tot_charge"] = args.tot_charge if "tot_charge" not in system or args.tot_charge != None else system["tot_charge"]
        system["nosym"] = args.nosym if "nosym" not in system or args.nosym != None else system["nosym"]
        system["nosym_evc"] = args.nosym_evc if "nosym_evc" not in system or args.nosym_evc != None else system["nosym_evc"]
        system["noinv"] = args.noinv if "noinv" not in system or args.noinv != None else system["noinv"]

        system["nspin"] = args.nspin if "nspin" not in system or args.nspin != None else system["nspin"]
        system["starting_magnetization"] = args.starting_magnetization if "starting_magnetization" not in system or args.starting_magnetization != None else system["starting_magnetization"]
        system["noncolin"] = args.noncolin if "noncolin" not in system or args.nnoncolin != None else system["noncolin"]
        
        system["lda_plus_u"] = args.lda_plus_u if "lda_plus_u" not in system or args.lda_plus_u != None else system["lda_plus_u"]
        system["lda_plus_u_kind"] = args.lda_plus_u_kind if "lda_plus_u_kind" not in system or args.lda_plus_u_Kind != None else system["lda_plus_u_kind"]
        system["Hubbard_U"] = args.hubbard_u if "Hubbard_U" not in system or args.hubbard_u != None else system["Hubbard_U"]
        system["Hubbard_J0"] = args.hubbard_j0 if "Hubbard_J0" not in system or args.hubbard_j0 != None else system["Hubbard_J0"]
        system["Hubbard_alpha"] = args.hubbard_alpha if "Hubbard_alpha" not in system or args.hubbard_alpha != None else system["Hubbard_alpha"]
        system["Hubbard_beta"] = args.hubbard_beta if "Hubbard_beta" not in system or args.hubbard_beta != None else system["Hubbard_beta"]
        system["U_projection_type"] = args.u_projection_type if "U_projection_type" not in system or args.u_projection_type != None else system["U_projection_type"]
        
        system["input_dft"] = args.input_dft if "input_dft" not in system or args.input_dft != None else params["input_dft"]
        system["ace"] = args.ace if "ace" not in system or args.ace != None else params["ace"]
        system["exx_fraction"]  = args.exx_fraction if "exx_fraction" not in system or args.exx_fraction != None else params["exx_fraction"]
        system["screening_parameter"] = args.screening_parameter if "screening_parameter" not in system or args.screening_parameter != None else params["screening_parameter"]
        system["exxdiv_treatment"] = args.exxdiv_treatment if "exxdiv_treatment" not in system or args.exxdiv_treatment != None else params["exxdiv_treatment"]
        system["x_gamma_extrapolation"] = args.x_gamma_extrapolation if "x_gamma_extrapolation" not in system or args.x_gamma_extrapolation != None else params["x_gamma_extrapolation"]
        system["ecutvcut"] = args.ecutvcut if "ecutvcut" not in system or args.ecutvcut != None else params["ecutvcut"]
        system["nqx1"] = args.nqx[0] if "nqx1" not in system or args.nqx[0] != None else params["nqx1"]
        system["nqx2"] = args.nqx[1] if "nqx2" not in system or args.nqx[1] != None else params["nqx2"]
        system["nqx3"] = args.nqx[2] if "nqx3" not in system or args.nqx[2] != None else params["nqx3"]
        

        electrons["electron_maxstep"] = args.electron_maxstep if "electron_maxstep" not in electrons or args.electron_maxstep != None else electrons["electron_maxstep"]
        electrons["conv_thr"] = args.conv_thr if "conv_thr" not in electrons or args.conv_thr != None else electrons["conv_thr"]
        electrons["mixing_beta"] = args.mixing_beta if "mixing_beta" not in electrons or args.mixing_beta != None else electrons["mixing_beta"]
        electrons["mixing_ndim"] = args.mixing_ndim if "mixing_ndim" not in electrons or args.mixing_ndim != None else electrons["mixing_ndim"]
        electrons["diagonalization"] = args.diagonalization if "diagonalization" not in electrons or args.diagonalization != None else electrons["diagonalization"]
        electrons["scf_must_converge"] = args.scf_must_converge if "scf_must_converge" not in electrons or args.scf_must_converge != None else electrons["scf_must_converge"]

        ions["ion_dynamics"] = args.ion_dynamics if "ion_dynamics" not in ions or args.ion_cynamics != None else ions["ion_dynamics"]
        ions["pot_extrapolation"] = args.pot_extrapolation if "pot_extrapolation" not in ions or args.pot_extrapolation != None else ions["pot_extrapolation"]
        ions["wfc_extrapolation"] = args.wfc_extrapolation if "wfc_extrapolation" not in ions or args.wfc_extrapolation != None else ions["wfc_extrapolation"]
        ions["ion_temperature"] = args.ion_temperature if "ion_temperature" not in ions or args.ion_temperature != None else ions["ion_temperature"]
        ions["tempw"] = args.tempw if "tempw" not in ions or args.tempw != None else ions["tempw"]
        

        path = {}
        if args.nebin != None:
            path = read_neb_in(args.nebin)
        path["string_method"] = args.string_method if "string_method" not in path or args.string_method != None else path["string_method"]
        path["nstep_path"] = args.nstep_path if "nstep_path" not in path or args.nstep_path != None else path["nstep_path"]
        path["opt_scheme"] = args.opt_scheme if "opt_scheme" not in path or args.opt_scheme != None else path["opt_scheme"]
        path["num_of_images"] = args.num_of_images if "num_of_images" not in path or args.num_of_images != None else path["num_of_images"]
        path["k_max"] = args.k_max if "k_max" not in path or args.k_max != None else path["k_max"]
        path["k_min"] = args.k_min if "k_min" not in path or args.k_min != None else path["k_min"]
        path["CI_scheme"] = args.ci_scheme if "CI_scheme" not in path or args.ci_scheme != None else path["CI_scheme"]
        path["path_thr"] = args.path_thr if "path_thr" not in path or args.path_thr != None else path["path_thr"]
        path["ds"] = args.ds if "ds" not in path or args.ds != None else path["ds"]
        path["first_last_opt"] = args.first_last_opt if "first_last_opt" not in path or args.first_last_opt != None else path["first_last_opt"]

        # for ph.x
        inputph = {}
        if args.phin != None:
            inputph = read_ph_in(args.phin)
        inputph["tr2_ph"] = args.tr2_ph if "tr2_ph" not in inputph or args.tr2_ph != None else inputph["tr2_ph"]
        inputph["lraman"] = args.lraman if "lraman" not in inputph or args.lraman != None else inputph["lraman"]
        inputph["epsil"] = args.epsil if "epsil" not in inputph or args.epsil != None else inputph["epsil"]
        inputph["nq1"] = args.nq[0] if "nq1" not in inputph or args.nq[0] != None else inputph["nq1"]
        inputph["nq2"] = args.nq[1] if "nq2" not in inputph or args.nq[1] != None else inputph["nq2"]
        inputph["nq3"] = args.nq[2] if "nq3" not in inputph or args.nq[2] != None else inputph["nq3"]
        inputph["search_sym"] = args.search_sym if "search_sym" not in inputph or args.search_sysm != None else inputph["search_sym"]


        if args.runtype == 0:
            # static scf nscf projwfc bands pp.x in a single run
            from pymatflow.qe.static import static_run
            projwfc_input = {}
            if args.projwfc_ngauss == 'default':
                ngauss = args.projwfc_ngauss
            else:
                ngauss = int(args.projwfc_ngauss)
            if args.projwfc_degauss == 'default':
                degauss = args.projwfc_degauss
            else:
                degauss = float(args.projwfc_degauss)
            if args.projwfc_emin == 'default':
                emin = args.projwfc_emin
            else:
                emin = float(args.projwfc_emin)
            if args.projwfc_emax == 'default':
                emax = args.projwfc_emax
            else:
                emax = float(args.projwfc_emax)
            if args.projwfc_deltae == 'default':
                deltae = args.projwfc_deltae
            else:
                deltae = float(args.projwfc_deltae)

            projwfc_input["filpdos"] = args.projwfc_filpdos
            projwfc_input["ngauss"] = ngauss
            projwfc_input["degauss"] = degauss
            projwfc_input["emin"] = emin
            projwfc_input["emax"] = emax
            projwfc_input["deltae"] = deltae
            bands = {}
            bands["lsym"] = args.lsym
            inputpp = {}
            plotpp = {}
            inputpp["plot_num"] = args.plot_num
            plotpp["iflag"] = args.iflag
            plotpp["output_format"] = args.output_format
            task = static_run()
            task.get_xyz(xyzfile)
            task.set_kpoints(kpoints_option=args.kpoints_option, kpoints_mp=args.kpoints_mp)
            task.set_params(control=control, system=system, electrons=electrons)
            task.set_atomic_forces(pressure=args.pressure, pressuredir=args.pressuredir)
            task.set_projwfc(projwfc_input=projwfc_input)
            task.set_bands(bands_input=bands)
            task.set_pp(inputpp=inputpp, plotpp=plotpp)
            task.set_run(mpi=args.mpi, server=server, jobname=args.jobname, nodes=args.nodes, ppn=args.ppn, queue=args.queue)
            task.set_llhpc(partition=args.partition, nodes=args.nodes, ntask=args.ntask, jobname=args.jobname, stdout=args.stdout, stderr=args.stderr)
            if args.static == "all":
                task.run(directory=args.directory, runopt=args.runopt, auto=args.auto, kpath=get_kpath(args.kpath_manual, args.kpath_file), kpoints_mp_scf=args.kpoints_mp, kpoints_mp_nscf=args.kpoints_mp_nscf)
            elif args.static == "scf":
                task.scf(directory=args.directory, runopt=args.runopt, auto=args.auto)
        elif args.runtype == 1:
            # relax
            from pymatflow.qe.opt import opt_run
            
            if args.fix != None or args.fix_around_z != None:
                # can only write xyz and poscar file
                from pymatflow.cmd.structflow import read_structure                
                a = read_structure(filepath=xyzfile)
                if args.fix != None:
                    fix = args.fix
                elif args.fix_around_z != None:
                    atoms_index_from_1 = []
                    for i in range(len(a.atoms)):
                        if a.atoms[i].z > (args.fix_around_z[0] + args.fix_around_z[1]) and a.atoms[i].z < (args.fix_around_z[0] + args.fix_around_z[2]):
                            atoms_index_from_1.append(i+1)
                    fix = atoms_index_from_1
                else:
                    fix = []                      
                    
                fix_str = ""
                for i in fix:
                        fix_str += "%d " % i
                os.system("xyz-fix-atoms.py -i %s -o %s --fix %s" % (xyzfile, xyzfile, fix_str))
                #args.selective_dynamics = "T"
                
                # output an xsd file with fixed atoms colored specifically so that user can check the atoms fixed
                from xml.etree.ElementTree import parse
                from pymatflow.cmd.structflow import write_structure
                os.system("mkdir -p /tmp/structflow/fix")
                write_structure(a, filepath="/tmp/structflow/fix/tmp.xsd")
                # read xsd file
                xsd = parse("/tmp/structflow/fix/tmp.xsd")
    
                # ID of Atom3D in xsd file start from 4
                imap = xsd.getroot().find("AtomisticTreeRoot").find("SymmetrySystem").find("MappingSet").find("MappingFamily").find("IdentityMapping")
                atoms = imap.findall("Atom3d")
                if args.color_fix == "white":
                    RGB = [255, 255, 255]
                elif args.color_fix == "red":
                    RGB = [255, 0, 0]
                elif args.color_fix == "green":
                    RGB = [0, 255, 0]
                elif args.color_fix == "blue":
                    RGB = [0, 0, 255]
                else:
                    RGB = [255, 255, 255] # default
            
                for i in fix:
                    atoms[i-1].set("Color", "%f, %f, %f, %f" % (RGB[0], RGB[1], RGB[2], 1))
                    
                # write xsd file
                xsd.write(xyzfile+".coloring.atoms.fixed.xsd")

            task = opt_run()
            task.get_xyz(xyzfile)
            task.set_relax()
            task.set_kpoints(kpoints_option=args.kpoints_option, kpoints_mp=args.kpoints_mp)
            task.set_params(control=control, system=system, electrons=electrons, ions=ions)
            task.set_run(mpi=args.mpi, server=server, jobname=args.jobname, nodes=args.nodes, ppn=args.ppn, queue=args.queue)
            task.set_llhpc(partition=args.partition, nodes=args.nodes, ntask=args.ntask, jobname=args.jobname, stdout=args.stdout, stderr=args.stderr)
            task.relax(directory=args.directory, runopt=args.runopt, auto=args.auto)
        elif args.runtype == 2:
            # vc-relax
            from pymatflow.qe.opt import opt_run
            task = opt_run()
            task.get_xyz(xyzfile)
            task.set_vc_relax()
            task.set_kpoints(kpoints_option=args.kpoints_option, kpoints_mp=args.kpoints_mp)
            task.set_params(control=control, system=system, electrons=electrons, ions=ions)
            task.set_run(mpi=args.mpi, server=server, jobname=args.jobname, nodes=args.nodes, ppn=args.ppn, queue=args.queue)
            task.set_llhpc(partition=args.partition, nodes=args.nodes, ntask=args.ntask, jobname=args.jobname, stdout=args.stdout, stderr=args.stderr)
            task.vc_relax(directory=args.directory, runopt=args.runopt, auto=args.auto)
        elif args.runtype == 3:
            # cubic cell opt
            from pymatflow.qe.opt import opt_run
            task = opt_run()
            task.get_xyz(xyzfile)
            task.set_relax()
            task.batch_a = args.batch_a     
            task.batch_b = args.batch_b
            task.batch_c = args.batch_c                
            task.set_kpoints(kpoints_option=args.kpoints_option, kpoints_mp=args.kpoints_mp)
            task.set_params(control=control, system=system, electrons=electrons, ions=ions)
            task.set_run(mpi=args.mpi, server=server, jobname=args.jobname, nodes=args.nodes, ppn=args.ppn, queue=args.queue)
            task.set_llhpc(partition=args.partition, nodes=args.nodes, ntask=args.ntask, jobname=args.jobname, stdout=args.stdout, stderr=args.stderr)
            #task.cubic(directory=args.directory, runopt=args.runopt, auto=args.auto, na=args.na, stepa=args.stepa)
            task.cubic(directory=args.directory, runopt=args.runopt, auto=args.auto, range_a=args.range_a)
        elif args.runtype == 4:
            # hexagonal cell opt
            from pymatflow.qe.opt import opt_run
            task = opt_run()
            task.get_xyz(xyzfile)
            task.set_relax()
            task.batch_a = args.batch_a     
            task.batch_b = args.batch_b
            task.batch_c = args.batch_c  
            task.set_kpoints(kpoints_option=args.kpoints_option, kpoints_mp=args.kpoints_mp)
            task.set_params(control=control, system=system, electrons=electrons, ions=ions)
            task.set_run(mpi=args.mpi, server=server, jobname=args.jobname, nodes=args.nodes, ppn=args.ppn, queue=args.queue)
            task.set_llhpc(partition=args.partition, nodes=args.nodes, ntask=args.ntask, jobname=args.jobname, stdout=args.stdout, stderr=args.stderr)
            #task.hexagonal(directory=args.directory, runopt=args.runopt, auto=args.auto, na=args.na, nc=args.nc, stepa=args.stepa, stepc=args.stepc)
            task.hexagonal(directory=args.directory, runopt=args.runopt, auto=args.auto, range_a=args.range_a, range_c=args.range_c)
        elif args.runtype == 5:
            # tetragonal cell opt
            from pymatflow.qe.opt import opt_run
            task = opt_run()
            task.get_xyz(xyzfile)
            task.set_relax()
            task.batch_a = args.batch_a     
            task.batch_b = args.batch_b
            task.batch_c = args.batch_c              
            task.set_kpoints(kpoints_option=args.kpoints_option, kpoints_mp=args.kpoints_mp)
            task.set_params(control=control, system=system, electrons=electrons, ions=ions)
            task.set_run(mpi=args.mpi, server=server, jobname=args.jobname, nodes=args.nodes, ppn=args.ppn, queue=args.queue)
            task.set_llhpc(partition=args.partition, nodes=args.nodes, ntask=args.ntask, jobname=args.jobname, stdout=args.stdout, stderr=args.stderr)
            #task.tetragonal(directory=args.directory, runopt=args.runopt, auto=args.auto, na=args.na, nc=args.nc, stepa=args.stepa, stepc=args.stepc)
            task.tetragonal(directory=args.directory, runopt=args.runopt, auto=args.auto, range_a=args.range_a, range_c=args.range_c)
        elif args.runtype == 6:
            from pymatflow.qe.neb import neb_run
            task = neb_run()
            task.get_images(images=images)
            task.set_kpoints(kpoints_option=args.kpoints_option, kpoints_mp=args.kpoints_mp)
            task.set_path(path=path)
            task.set_run(mpi=args.mpi, server=server, jobname=args.jobname, nodes=args.nodes, ppn=args.ppn, queue=args.queue)
            task.set_llhpc(partition=args.partition, nodes=args.nodes, ntask=args.ntask, jobname=args.jobname, stdout=args.stdout, stderr=args.stderr)
            task.neb(directory=args.directory, runopt=args.runopt, auto=args.auto)
        elif args.runtype == 7:
            from pymatflow.qe.dfpt import dfpt_run
            task = dfpt_run()
            task.get_xyz(xyzfile)
            task.set_inputph(inputph=inputph)
            task.set_run(mpi=args.mpi, server=server, jobname=args.jobname, nodes=args.nodes, ppn=args.ppn, queue=args.queue)
            task.set_llhpc(partition=args.partition, nodes=args.nodes, ntask=args.ntask, jobname=args.jobname, stdout=args.stdout, stderr=args.stderr)
            task.phx(directory=args.directory, runopt=args.runopt, auto=args.auto)
        elif args.runtype == 8:
            # phonopy
            from pymatflow.qe.phonopy import phonopy_run
            task = phonopy_run()
            task.get_xyz(xyzfile)
            task.set_kpoints(kpoints_option=args.kpoints_option, kpoints_mp=args.kpoints_mp)
            task.set_params(control=control, system=system, electrons=electrons)
            task.supercell_n = args.supercell_n
            task.set_run(mpi=args.mpi, server=server, jobname=args.jobname, nodes=args.nodes, ppn=args.ppn, queue=args.queue)
            task.set_llhpc(partition=args.partition, nodes=args.nodes, ntask=args.ntask, jobname=args.jobname, stdout=args.stdout, stderr=args.stderr)
            task.phonopy(directory=args.directory, runopt=args.runopt, auto=args.auto)
        elif args.runtype == 9:
            # pp.x
            from pymatflow.qe.static import static_run

            inputpp = {}
            plotpp = {}
            
            inputpp["plot_num"] = args.plot_num
            plotpp["iflag"] = args.iflag
            plotpp["output_format"] = args.output_format

            task = static_run()
            task.get_xyz(xyzfile)
            task.set_pp(inputpp=inputpp, plotpp=plotpp)
            task.set_run(mpi=args.mpi, server=server, jobname=args.jobname, nodes=args.nodes, ppn=args.ppn, queue=args.queue)
            task.set_llhpc(partition=args.partition, nodes=args.nodes, ntask=args.ntask, jobname=args.jobname, stdout=args.stdout, stderr=args.stderr)
            task.pp(directory=args.directory, runopt=args.runopt, auto=args.auto)
        elif args.runtype == 10:
            # abc cell opt
            from pymatflow.qe.opt import opt_run
            task = opt_run()
            task.get_xyz(xyzfile)
            task.set_relax()
            task.set_kpoints(kpoints_option=args.kpoints_option, kpoints_mp=args.kpoints_mp)
            task.set_params(control=control, system=system, electrons=electrons, ions=ions)
            task.set_run(mpi=args.mpi, server=server, jobname=args.jobname, nodes=args.nodes, ppn=args.ppn, queue=args.queue)
            task.set_llhpc(partition=args.partition, nodes=args.nodes, ntask=args.ntask, jobname=args.jobname, stdout=args.stdout, stderr=args.stderr)
            task.batch_a = args.batch_a     
            task.batch_b = args.batch_b
            task.batch_c = args.batch_c     
            task.abc(directory=args.directory, runopt=args.runopt, auto=args.auto, range_a=args.range_a, range_b=args.range_b, range_c=args.range_c)
        elif args.runtype == 11:
            # converge test
            from pymatflow.qe.static import static_run
            
            task = static_run()
            task.get_xyz(xyzfile)
            task.set_kpoints(kpoints_option=args.kpoints_option, kpoints_mp=args.kpoints_mp)
            task.set_params(control=control, system=system, electrons=electrons, ions=ions)
            task.set_run(mpi=args.mpi, server=server, jobname=args.jobname, nodes=args.nodes, ppn=args.ppn, queue=args.queue)
            task.set_llhpc(partition=args.partition, nodes=args.nodes, ntask=args.ntask, jobname=args.jobname, stdout=args.stdout, stderr=args.stderr)
            if args.converge == "ecutwfc":            
                task.converge_ecutwfc(args.ecutwfc_range[0], args.ecutwfc_range[1], args.ecutwfc_range[2], directory=args.directory, runopt=args.runopt, auto=args.auto)
            elif args.converge == "ecutrho":
                task.converge_ecutrho(args.ecutrho_range[0], args.ecutrho_range[1], args.ecutrho_range[2], args.ecutwfc, directory=args.directory, runopt=args.runopt, auto=args.auto)
            elif args.converge == "degauss":
                task.converge_degauss(round(args.degauss_range[0], 6), round(args.degauss_range[1], 6), round(args.degauss_range[2], 6), directory=args.directory, runopt=args.runopt, auto=args.auto)
            elif args.converge == "kpoints":
                pass
        else:
            pass

# ==============================================================================
# SIESTA SIESTA SIESTA SIESTA SIESTA SIESTA SIESTA SIESTA SIESTA SIESTA SIESTA
# ==============================================================================
    elif args.driver == "siesta":
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
            from pymatflow.siesta.static import static_run
            task = static_run()
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

            task.set_params(params=params)
            task.set_kpoints(kpoints_mp=args.kpoints_mp)
            task.set_run(mpi=args.mpi, server=server, jobname=args.jobname, nodes=args.nodes, ppn=args.ppn, queue=args.queue)
            task.set_llhpc(partition=args.partition, nodes=args.nodes, ntask=args.ntask, jobname=args.jobname, stdout=args.stdout, stderr=args.stderr)
            task.scf(directory=args.directory, runopt=args.runopt, auto=args.auto, properties=args.properties)
        elif args.runtype == 1:
            # optimization
            from pymatflow.siesta.opt import opt_run
            params["MD.VariableCell"] = args.variablecell
            params["MD.MaxForceTol"] = args.forcetol
            params["MD.MaxStressTol"] = args.stresstol
            params["MD.TargetPressure"] = args.targetpressure
            task = opt_run()
            task.get_xyz(xyzfile)
            task.set_params(params=params)
            task.set_kpoints(kpoints_mp=args.kpoints_mp)
            task.set_run(mpi=args.mpi, server=server, jobname=args.jobname, nodes=args.nodes, ppn=args.ppn, queue=args.queue)
            task.set_llhpc(partition=args.partition, nodes=args.nodes, ntask=args.ntask, jobname=args.jobname, stdout=args.stdout, stderr=args.stderr)
            task.opt(directory=args.directory, runopt=args.runopt, auto=args.auto)
        elif args.runtype == 2:
            # cubic cell
            from pymatflow.siesta.opt import opt_run
            params["MD.VariableCell"] = "false"
            params["MD.MaxForceTol"] = args.forcetol
            params["MD.MaxStressTol"] = args.stresstol
            params["MD.TargetPressure"] = args.targetpressure
            task = opt_run()
            task.get_xyz(xyzfile)
            task.set_params(params=params)
            task.set_kpoints(kpoints_mp=args.kpoints_mp)
            task.set_run(mpi=args.mpi, server=server, jobname=args.jobname, nodes=args.nodes, ppn=args.ppn, queue=args.queue)
            task.set_llhpc(partition=args.partition, nodes=args.nodes, ntask=args.ntask, jobname=args.jobname, stdout=args.stdout, stderr=args.stderr)
            task.cubic(directory=args.directory, runopt=args.runopt, auto=args.auto, na=args.na, stepa=args.stepa)
        elif args.runtype == 3:
            # hexagonal cell
            from pymatflow.siesta.opt import opt_run
            params["MD.VariableCell"] = "false"
            params["MD.MaxForceTol"] = args.forcetol
            params["MD.MaxStressTol"] = args.stresstol
            params["MD.TargetPressure"] = args.targetpressure
            task = opt_run()
            task.get_xyz(xyzfile)
            task.set_params(params=params)
            task.set_kpoints(kpoints_mp=args.kpoints_mp)
            task.set_run(mpi=args.mpi, server=server, jobname=args.jobname, nodes=args.nodes, ppn=args.ppn, queue=args.queue)
            task.set_llhpc(partition=args.partition, nodes=args.nodes, ntask=args.ntask, jobname=args.jobname, stdout=args.stdout, stderr=args.stderr)
            task.hexagonal(directory=args.directory, runopt=args.runopt, auto=args.auto, na=args.na, nc=args.nc, stepa=args.stepa, stepc=args.stepc)
        elif args.runtype == 4:
            # tetragonal cell
            from pymatflow.siesta.opt import opt_run
            params["MD.VariableCell"] = "false"
            params["MD.MaxForceTol"] = args.forcetol
            params["MD.MaxStressTol"] = args.stresstol
            params["MD.TargetPressure"] = args.targetpressure
            task = opt_run()
            task.get_xyz(xyzfile)
            task.set_params(params=params)
            task.set_kpoints(kpoints_mp=args.kpoints_mp)
            task.set_run(mpi=args.mpi, server=server, jobname=args.jobname, nodes=args.nodes, ppn=args.ppn, queue=args.queue)
            task.set_llhpc(partition=args.partition, nodes=args.nodes, ntask=args.ntask, jobname=args.jobname, stdout=args.stdout, stderr=args.stderr)
            task.tetragonal(directory=args.directory, runopt=args.runopt, auto=args.auto, na=args.na, nc=args.nc, stepa=args.stepa, stepc=args.stepc)
        elif args.runtype == 5:
            # phonopy
            from pymatflow.siesta.phonopy import phonopy_run
            task = phonopy_run()
            task.get_xyz(xyzfile)
            task.set_params(params=params)
            task.set_kpoints(kpoints_mp=args.kpoints_mp)
            task.supercell_n = args.supercell_n
            task.set_run(mpi=args.mpi, server=server, jobname=args.jobname, nodes=args.nodes, ppn=args.ppn, queue=args.queue)
            task.set_llhpc(partition=args.partition, nodes=args.nodes, ntask=args.ntask, jobname=args.jobname, stdout=args.stdout, stderr=args.stderr)
            task.phonopy(directory=args.directory, runopt=args.runopt, auto=args.auto)
        elif args.runtype == 6:
            # molecular dynamics
            from pymatflow.siesta.md import md_run
            params["MD.FinalTimeStep"] = args.mdstep
            params["MD.LengthTimeStep"] = args.timestep
            params["MD.InitialTemperature"] = args.initial_temp
            params["MD.TargetTemperature"] = args.target_temp         
            task = md_run()
            task.get_xyz(xyzfile)
            task.set_params(params=params)
            task.set_kpoints(kpoints_mp=args.kpoints_mp)
            task.set_run(mpi=args.mpi, server=server, jobname=args.jobname, nodes=args.nodes, ppn=args.ppn, queue=args.queue)
            task.set_llhpc(partition=args.partition, nodes=args.nodes, ntask=args.ntask, jobname=args.jobname, stdout=args.stdout, stderr=args.stderr)
            task.md(directory=args.directory, runopt=args.runopt, auto=args.auto)
        elif args.runtype == 7:
            # abc cell opt
            from pymatflow.siesta.opt import opt_run
            task = opt_run()
            task.get_xyz(xyzfile)
            task.set_relax()
            task.set_kpoints(kpoints_option=args.kpoints_option, kpoints_mp=args.kpoints_mp)
            task.set_params(control=control, system=system, electrons=electrons, ions=ions)
            task.set_run(mpi=args.mpi, server=server, jobname=args.jobname, nodes=args.nodes, ppn=args.ppn, queue=args.queue)
            task.set_llhpc(partition=args.partition, nodes=args.nodes, ntask=args.ntask, jobname=args.jobname, stdout=args.stdout, stderr=args.stderr)
            task.batch_a = args.batch_a     
            task.batch_b = args.batch_b
            task.batch_c = args.batch_c     
            task.abc(directory=args.directory, runopt=args.runopt, auto=args.auto, range_a=args.range_a, range_b=args.range_b, range_c=args.range_c)
        else:
            pass
    elif args.driver == "vasp":
        params = {}
        # deal with INCAR template specified by --incar
        if args.incar == None:
            pass
        else:
            if not os.path.exists(args.incar):
                print("====================================================\n")
                print("                  Warning !!!!\n")
                print("----------------------------------------------------\n")
                print("matflow vasp:\n")
                print("the specified incar file by --incar doesn't exist\n")
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
        # if xxx is alraedy in params(set from --incar) and args.xxx is None
        # params[xxx] will not be control by args.xxx
        params["NWRITE"] = args.nwrite if "NWRITE" not in params or args.nwrite != None else params["NWRITE"]
        params["PREC"] = args.prec if "PREC" not in params or args.prec != None else params["PREC"]
        params["NCORE"] = args.ncore if "NCORE" not in params or args.ncore != None else params["NCORE"]
        params["ENCUT"] = args.encut if "ENCUT" not in params or args.encut != None else params["ENCUT"]
        params["EDIFF"] = args.ediff if "EDIFF" not in params or args.ediff != None else params["EDIFF"]
        params["NELM"] = args.nelm if "NELM" not in params or args.nelm != None else params["NELM"]
        params["NFREE"] = args.nfree if "NFREE" not in params or args.nfree != None else params["NFREE"]
        params["ISMEAR"] = args.ismear if "ISMEAR" not in params or args.ismear != None else params["ISMEAR"]
        params["SIGMA"] = args.sigma if "SIGMA" not in params or args.sigma != None else params["SIGMA"]
        params["IVDW"] = args.ivdw if "IVDW" not in params or args.ivdw != None else params["IVDW"]
        params["EDIFFG"] = args.ediffg if "EDIFFG" not in params or args.ediffg != None else params["EDIFFG"]
        params["NSW"] = args.nsw if "NSW" not in params or args.nsw != None else params["NSW"]
        params["IBRION"] = args.ibrion if "IBRION" not in params or args.ibrion != None else params["IBRION"]
        params["ISIF"] = args.isif if "ISIF" not in params or args.isif != None else params["ISIF"]
        params["POTIM"] = args.potim if "POTIM" not in params or args.potim != None else params["POTIM"]
        params["LORBIT"] = args.lorbit if "LORBIT" not in params or args.potim != None else params["LORBIT"]
        params["LOPTICS"] = args.loptics if "LOPTICS" not in params or args.loptics != None else params["LOPTICS"]
        params["CSHIFT"] = args.cshift if "CSHIFT" not in params or args.cshift != None else params["CSHIFT"]
        params["NEDOS"] = args.nedos if "NEDOS"  not in params or args.nedos != None else params["NEDOS"]
        params["LSUBROT"] = args.lsubrot if "LSUBROT" not in params or args.lsubrot != None else params["LSUBROT"]
        params["SAXIS"] = args.saxis if "SAXIS" not in params or args.saxis != None else params["SAXIS"]
        params["LMAXMIX"] = args.lmaxmix if "LMAXMIX" not in params or args.lmaxmix != None else params["LMAXMIX"]
        params["MDALGO"] = args.mdalgo if "MDALGO" not in params or args.mdalgo != None else params["MDALGO"]
        params["SMASS"] = args.smass if "SMASS" not in params or args.smass != None else params["SMASS"]
        params["ANDERSON_PROB"] = args.anderson_prob if "ANDERSON_PROB" not in params or args.anderson_prob != None else params["ANDERSON_PROB"]
        params["TEBEG"] = args.tebeg if "TEBEG" not in params or args.tebeg != None else params["TEBEG"]
        params["TEEND"] = args.teend if "TEEND" not in params or args.teend != None else params["TEEND"]
        params["ALGO"] = args.algo if "ALGO" not in params or args.algo != None else params["ALGO"]
        params["IALGO"] = args.ialgo if "IALGO" not in params or args.ialgo != None else params["IALGO"]
        params["ADDGRID"] = args.addgrid if "ADDGRID" not in params or args.addgrid != None else params["ADDGRID"]
        params["ISYM"] = args.isym if "ISYM" not in params or args.isym != None else params["ISYM"]
        params["LREAL"] = args.lreal if "LREAL" not in params or args.lreal != None else params["LREAL"]
        params["PSTRESS"] = args.pstress if "PSTRESS" not in params or args.pstress != None else params["PSTRESS"]
        params["LWAVE"] = args.lwave if "LWAVE" not in params or args.lwave != None else params["LWAVE"]
        params["LCHARG"] = args.lcharg if "LCHARG" not in params or args.lcharg != None else params["LCHARG"]
        params["ISPIN"] = args.ispin if "ISPIN" not in params or args.ispin != None else params["ISPIN"]
        params["MAGMOM"] = args.magmom if "MAGMOM" not in params or args.magmom != None else params["MAGMOM"] # magmom can be a list that can be automatically dealt with by base.incar.to_incar()
        params["LNONCOLLINEAR"] = args.lnoncollinear if "LNONCOLLINEAR" not in params or args.lnoncollinear != None else params["LNONCOLLINEAR"]
        params["LSORBIT"] = args.lsorbit if "LSORBIT" not in params or args.lsorbit != None else params["LSORBIT"]
        params["ALGO"] = args.algo if "ALGO" not in params or args.algo != None else params["ALGO"]
        params["LHFCALC"] = args.lhfcalc if "LHFCALC" not in params or args.lhfcalc != None else params["LHFCALC"]
        params["HFSCREEN"] = args.hfscreen if "HFSCREEN" not in params or args.hfscreen != None else params["HFSCREEN"]
        params["AEXX"] = args.aexx if "AEXX" not in params or args.aexx != None else params["AEXX"]
        params["LELF"] = args.lelf if "LELF" not in params or args.lelf != None else params["LELF"]
        params["IOPT"] = args.iopt if "IOPT" not in params or args.iopt != None else params["IOPT"]
        params["LCLIMB"] = args.lclimb if "LCLIMB" not in params or args.lclimb != None else params["LCLIMB"]
        params["LNEBCELL"] = args.lnebcell if "LNEBCELL" not in params or args.lnebcell != None else params["LNEBCELL"]
        params["SPRING"] = args.spring if "SPRING" not in params or args.spring != None else params["SPRING"]
        params["MAXMOVE"] = args.maxmove if "MAXMOVE" not in params or args.maxmove != None else params["MAXMOVE"]
        params["LGLOBAL"] = args.lglobal if "LGLOBAL" not in params or args.lglobal != None else params["LGLOBAL"]
        params["LAUTOSCALE"] = args.lautoscale if "LAUTOSCALE" not in params or args.lautoscale != None else params["LAUTOSCALE"]
        params["INVCURV"] = args.invcurv if "INVCURV" not in params or args.invcurv != None else params["INVCURV"]
        params["IMAGES"] = args.nimage if "IMAGES" not in params or args.images != None else params["IMAGES"]
        params["LLINEOPT"] = args.llineopt if "LLINEOPT" not in params or args.llineopt != None else params["LLINEOPT"]
        params["FDSTEP"] = args.fdstep if "FDSTEP" not in params or args.fdstep != None else params["FDSTEP"]
        params["SYMPREC"] = args.symprec if "SYMPREC" not in params or args.symprec != None else params["SYMPREC"]
        params["AMIX"] = args.amix if "AMIX" not in params or args.amix != None else params["AMIX"]
        params["BMIX"] = args.bmix if "BMIX" not in params or args.bmix != None else params["BMIX"]
        params["NELECT"] = args.nelect if "NELECT" not in params or args.nelect != None else params["NELECT"]
        params["LAECHG"] = args.laechg if "LAECHG" not in params or args.laechg != None else params["LAECHG"]
        params["LPARD"] = args.lpard if "LPARD" not in params or args.lpard != None else params["LPARD"]
        params["LSEPK"] = args.lsepk if "LSEPK" not in params or args.lsepk != None else params["LSEPK"]
        params["LSEPB"] = args.lsepb if "LSEPB" not in params or args.lsepb != None else params["LSEPB"]
        params["NBMOD"] = args.nbmod if "NBMOD" not in params or args.nbmod != None else params["NBMOD"]
        params["EINT"] = args.eint if "EINT" not in params or args.eint != None else params["EINT"]
        
        
        if args.runtype == 0:
            # static
            from pymatflow.vasp.static import static_run
            task = static_run()
            task.get_xyz(xyzfile)
            task.set_params(params, runtype="static")
            task.set_kpoints(kpoints_mp=args.kpoints_mp)
            task.set_run(mpi=args.mpi, server=server, jobname=args.jobname, nodes=args.nodes, ppn=args.ppn, queue=args.queue)
            task.set_llhpc(partition=args.partition, nodes=args.nodes, ntask=args.ntask, jobname=args.jobname, stdout=args.stdout, stderr=args.stderr)
            if params["LNONCOLLINEAR"] != None:
                if params["LNONCOLLINEAR"].upper() == ".TRUE." or params["LNONCOLLINEAR"].upper() == "T":
                    task.magnetic_status = "non-collinear"

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
            elif args.static == "parchg" or args.static == "stm":
                if args.hse_in_scf.lower() == "true":
                    hse_in_scf = True
                elif args.hse_in_scf.lower() == "false":
                    hse_in_scf = False
                task.set_kpoints(kpoints_mp=args.kpoints_mp)
                task.parchg_stm(directory=args.directory, runopt=args.runopt, auto=args.auto, hse_in_scf=hse_in_scf)                
        elif args.runtype == 1:
            # optimization
            from pymatflow.vasp.opt import opt_run
            #  
            
            if args.fix != None or args.fix_around_z != None:
                # can only write xyz and poscar file
                from pymatflow.cmd.structflow import read_structure                
                a = read_structure(filepath=xyzfile)
                if args.fix != None:
                    fix = args.fix
                elif args.fix_around_z != None:
                    atoms_index_from_1 = []
                    for i in range(len(a.atoms)):
                        if a.atoms[i].z > (args.fix_around_z[0] + args.fix_around_z[1]) and a.atoms[i].z < (args.fix_around_z[0] + args.fix_around_z[2]):
                            atoms_index_from_1.append(i+1)
                    fix = atoms_index_from_1
                else:
                    fix = []                      
                    
                fix_str = ""
                for i in fix:
                        fix_str += "%d " % i
                os.system("xyz-fix-atoms.py -i %s -o %s --fix %s" % (xyzfile, xyzfile, fix_str))
                args.selective_dynamics = "T"
                
                # output an xsd file with fixed atoms colored specifically so that user can check the atoms fixed
                from xml.etree.ElementTree import parse
                from pymatflow.cmd.structflow import write_structure
                os.system("mkdir -p /tmp/structflow/fix")
                write_structure(a, filepath="/tmp/structflow/fix/tmp.xsd")
                # read xsd file
                xsd = parse("/tmp/structflow/fix/tmp.xsd")
    
                # ID of Atom3D in xsd file start from 4
                imap = xsd.getroot().find("AtomisticTreeRoot").find("SymmetrySystem").find("MappingSet").find("MappingFamily").find("IdentityMapping")
                atoms = imap.findall("Atom3d")
                if args.color_fix == "white":
                    RGB = [255, 255, 255]
                elif args.color_fix == "red":
                    RGB = [255, 0, 0]
                elif args.color_fix == "green":
                    RGB = [0, 255, 0]
                elif args.color_fix == "blue":
                    RGB = [0, 0, 255]
                else:
                    RGB = [255, 255, 255] # default
            
                for i in fix:
                    atoms[i-1].set("Color", "%f, %f, %f, %f" % (RGB[0], RGB[1], RGB[2], 1))
                    
                # write xsd file
                xsd.write(xyzfile+".coloring.atoms.fixed.xsd")
                        
            #            
            task = opt_run()
            task.get_xyz(xyzfile)
            task.set_params(params=params, runtype="opt")
            task.set_kpoints(kpoints_mp=args.kpoints_mp)
            task.poscar.selective_dynamics = True if args.selective_dynamics.upper()[0] == "T" else False
            task.set_run(mpi=args.mpi, server=server, jobname=args.jobname, nodes=args.nodes, ppn=args.ppn, queue=args.queue)
            task.set_llhpc(partition=args.partition, nodes=args.nodes, ntask=args.ntask, jobname=args.jobname, stdout=args.stdout, stderr=args.stderr)
            task.optimize(directory=args.directory, runopt=args.runopt, auto=args.auto)
        elif args.runtype == 2:
            # cubic cell
            from pymatflow.vasp.opt import opt_run
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
            task.set_run(mpi=args.mpi, server=server, jobname=args.jobname, nodes=args.nodes, ppn=args.ppn, queue=args.queue)
            task.set_llhpc(partition=args.partition, nodes=args.nodes, ntask=args.ntask, jobname=args.jobname, stdout=args.stdout, stderr=args.stderr)
            task.batch_a = args.batch_a
            task.cubic(directory=args.directory, runopt=args.runopt, auto=args.auto, range_a=args.range_a)
        elif args.runtype == 3:
            # hexagonal cell
            from pymatflow.vasp.opt import opt_run
            task = opt_run()
            task.get_xyz(xyzfile)
            task.set_params(params=params, runtype="opt")
            task.set_kpoints(kpoints_mp=args.kpoints_mp)
            task.set_run(mpi=args.mpi, server=server, jobname=args.jobname, nodes=args.nodes, ppn=args.ppn, queue=args.queue)
            task.set_llhpc(partition=args.partition, nodes=args.nodes, ntask=args.ntask, jobname=args.jobname, stdout=args.stdout, stderr=args.stderr)
            task.batch_a = args.batch_a
            task.batch_c = args.batch_c            
            task.hexagonal(directory=args.directory, runopt=args.runopt, auto=args.auto, range_a=args.range_a, range_c=args.range_c)
        elif args.runtype == 4:
            # tetragonal cell
            from pymatflow.vasp.opt import opt_run
            task = opt_run()
            task.get_xyz(xyzfile)
            task.set_params(params=params, runtype="opt")
            task.set_kpoints(kpoints_mp=args.kpoints_mp)
            task.set_run(mpi=args.mpi, server=server, jobname=args.jobname, nodes=args.nodes, ppn=args.ppn, queue=args.queue)
            task.set_llhpc(partition=args.partition, nodes=args.nodes, ntask=args.ntask, jobname=args.jobname, stdout=args.stdout, stderr=args.stderr)
            task.batch_a = args.batch_a     
            task.batch_c = args.batch_c            
            task.tetragonal(directory=args.directory, runopt=args.runopt, auto=args.auto, range_a=args.range_a, range_c=args.range_c)
        elif args.runtype == 5:
            # neb
            # we better set NSW manually in VTST neb calc. 
            # if not set, pymatflow.vasp.neb will set it to 100 automatically
            from pymatflow.vasp.neb import neb_run
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
            task.set_run(mpi=args.mpi, server=server, jobname=args.jobname, nodes=args.nodes, ppn=args.ppn, queue=args.queue)
            task.set_llhpc(partition=args.partition, nodes=args.nodes, ntask=args.ntask, jobname=args.jobname, stdout=args.stdout, stderr=args.stderr)
            task.neb(directory=args.directory, runopt=args.runopt, auto=args.auto)
            # move the OUTCAR for initial stucture and final structure to the corresponding dir
            # if they are specified
            if args.outcars != None and len(args.outcars) > 0:
                os.system("cp %s %s" % (args.outcars[0], os.path.join(args.directory, "00/")))
                os.system("cp %s %s" % (args.outcars[-1], os.path.join(args.directory, "%.2d/" % (args.nimage+1))))                
        elif args.runtype == 6:
            # vasp phonon
            from pymatflow.vasp.phonon import phonon_run
            task = phonon_run() 
            task.get_xyz(xyzfile)
            task.set_params(params=params, runtype="phonon")
            task.set_kpoints(kpoints_mp=args.kpoints_mp)
            task.supercell_n = args.supercell_n
            task.set_run(mpi=args.mpi, server=server, jobname=args.jobname, nodes=args.nodes, ppn=args.ppn, queue=args.queue)
            task.set_llhpc(partition=args.partition, nodes=args.nodes, ntask=args.ntask, jobname=args.jobname, stdout=args.stdout, stderr=args.stderr)
            task.phonon(directory=args.directory, runopt=args.runopt, auto=args.auto)
        elif args.runtype == 7:
            # phonopy
            from pymatflow.vasp.phonopy import phonopy_run
            task = phonopy_run()
            task.get_xyz(xyzfile)
            task.set_params(params=params, runtype="phonopy")
            task.set_kpoints(kpoints_mp=args.kpoints_mp)
            task.supercell_n = args.supercell_n
            task.set_run(mpi=args.mpi, server=server, jobname=args.jobname, nodes=args.nodes, ppn=args.ppn, queue=args.queue)
            task.set_llhpc(partition=args.partition, nodes=args.nodes, ntask=args.ntask, jobname=args.jobname, stdout=args.stdout, stderr=args.stderr)
            task.phonopy(directory=args.directory, runopt=args.runopt, auto=args.auto)
        elif args.runtype == 8:
            # sur pes
            from pymatflow.flow.surface_pes import vasp_run
            task = vasp_run()
            task.get_xyz(xyzfile)
            task.set_params(params=params, runtype="opt")
            task.set_kpoints(kpoints_mp=args.kpoints_mp)
            #task.poscar.selective_dynamics = True if args.selective_dynamics.upper()[0] == "T" else False
            task.poscar.selective_dynamics = True # always use selective_dynamics            
            task.set_run(mpi=args.mpi, server=server, jobname=args.jobname, nodes=args.nodes, ppn=args.ppn, queue=args.queue)
            #task.set_llhpc(partition=args.partition, nodes=args.nodes, ntask=args.ntask, jobname=args.jobname, stdout=args.stdout, stderr=args.stderr)
            task.batch_x_y = args.batch_x_y
            task.set_pes(move_atom=args.move_atom, xrange=args.xrange, yrange=args.yrange, zshift=args.zshift, fix_z=args.fix_z, fix_y=args.fix_y, fix_x=args.fix_x)
            task.run(directory=args.directory, runopt=args.runopt, auto=args.auto)
        elif args.runtype == 9:
            # abc cell
            from pymatflow.vasp.opt import opt_run
            task = opt_run()
            task.get_xyz(xyzfile)
            task.set_params(params=params, runtype="opt")
            task.set_kpoints(kpoints_mp=args.kpoints_mp)
            task.set_run(mpi=args.mpi, server=server, jobname=args.jobname, nodes=args.nodes, ppn=args.ppn, queue=args.queue)
            task.set_llhpc(partition=args.partition, nodes=args.nodes, ntask=args.ntask, jobname=args.jobname, stdout=args.stdout, stderr=args.stderr)
            task.batch_a = args.batch_a     
            task.batch_b = args.batch_b
            task.batch_c = args.batch_c            
            task.abc(directory=args.directory, runopt=args.runopt, auto=args.auto, range_a=args.range_a, range_b=args.range_b, range_c=args.range_c)
        elif args.runtype == 10:
            # AIMD
            from pymatflow.vasp.md import md_run
            task = md_run()
            task.get_xyz(xyzfile)
            task.set_params(params=params, runtype="md")
            task.set_kpoints(kpoints_mp=args.kpoints_mp)
            """
            ensemble:
                    0: NVE
                    1: NVT
                    2: NPT
                    3: NPH
            thermostat:
                    0: Anderson
                    1: Nose-Hoover
                    2: Langevin
                    3: Multiple Anderson
            """
            #task.incar.set_md(ensemble=ensemble, thermostat=thermostat)
            task.set_run(mpi=args.mpi, server=server, jobname=args.jobname, nodes=args.nodes, ppn=args.ppn, queue=args.queue)
            task.set_llhpc(partition=args.partition, nodes=args.nodes, ntask=args.ntask, jobname=args.jobname, stdout=args.stdout, stderr=args.stderr)
            task.md(directory=args.directory, runopt=args.runopt, auto=args.auto)
    # --------------------------------------------------------------------------



if __name__ == "__main__":
    main()
