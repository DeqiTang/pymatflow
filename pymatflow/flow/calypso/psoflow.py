#!/usr/bin/env python

import os
import sys
import argparse



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
            choices=[0, 1, 2, 3, 4, 5, 6, 7],
            help="choices of runtype. 0->static_run; 1->optimization; 2->cubic-opt; 3->hexagonal-opt; 4->tetragonal-opt; 5->dfpt-elastic-piezo-dielec; 6->dfpt-phonon; 7->phonopy")

    gp.add_argument("-d", "--directory", type=str, default="matflow-running",
            help="Directory to do the calculation")
            
    # calypso input.dat template
    gp = subparser.add_argument_group(title="template", 
            description="read in Calypso input.dat template")

    gp.add_argument("--input-dat", type=str, default=None,
            help="specify Calypso input.dat template to set parameters")
    
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

    # --------------------------------------------------------------------------
    # CP2K
    # --------------------------------------------------------------------------
    subparser = subparsers.add_parser("cp2k", help="using cp2k as calculator")

    gp = subparser.add_argument_group(title="overall running control")

    gp.add_argument("-r", "--runtype", type=int, default=0,
            choices=[0, 1, 2, 3, 4 ,5, 6, 7, 8, 9, 10],
            help="choices of runtype. 0->static_run; 1->geo-opt; 2->cell-opt; 3->cubic-cell; 4->hexagonal-cell; 5->tetragonal-cell; 6-neb; 7->phonopy; 8->vibrational_analysis; 9->converge test; 10->aimd")

    gp.add_argument("-d", "--directory", type=str, default="matflow-running",
            help="Directory to do the calculation")
    
    # calypso input.dat template
    gp = subparser.add_argument_group(title="template", 
            description="read in Calypso input.dat template")

    gp.add_argument("--input-dat", type=str, default=None,
            help="specify Calypso input.dat template to set parameters")


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

    # potential file
    gp = subparser.add_argument_group(title="pseudopotential")

    gp.add_argument("--pot", type=str, default="auto",
            choices=["auto"],
            help="setting pseudopotential file, in cp2k can only be auto(no need to set)")

    # GLOBAL
    gp = subparser.add_argument_group(title="GLOBAL")

    gp.add_argument("--print-level", type=str, default=None,
            choices=["DEBUG", "HIGH", "LOW", "MEDIUM", "SILENT", "debug", "high", "low", "medium", "silent"],
            help="How much output is written out.")

    # FORCE_EVAL/SUBSYS
    gp = subparser.add_argument_group(title="FORCE_EVAL/SUBSYS")

    gp.add_argument("--cell-symmetry", type=str, default=None,
            help="Imposes an initial cell symmetry. must be set when you do KEEP_SYMMETRY cell_opt")

    # FORCE_EVAL/DFT
    gp = subparser.add_argument_group(title="FORCE_EVAL/DFT")

    gp.add_argument("--qs-method", type=str, default="gpw",
            choices=["am1", "dftb", "gapw", "gapw_xc", "gpw", "lrigpw", "mndo", "mndod",
                "ofgpw", "pdg", "pm3", "pm6", "pm6-fm", "pnnl", "rigpw", "rm1"],
            help="specify the electronic structure method that should be employed, default is gpw")

    # FORCE_EVAL/DFT/SCF
    gp = subparser.add_argument_group(title="FORCE_EVAL/DFT/SCF")

    gp.add_argument("--eps-scf", type=float, default=1.0e-6,
            help="target accuracy for the scf convergence, default is 1.0e-6")

    gp.add_argument("--xc-functional", type=str, default="pbe",
            choices=["B3LYP", "BEEFVDW", "BLYP", "BP", "LDA", "PBE", "PADE", "PBE0", "TPSS",
                "b3lyp", "beefvdw", "blyp", "bp", "lda", "pbe", "pade", "pbe0", "tpss"],
            help="shortcut for the most common functional combinations, default is PBE")

    gp.add_argument("--cutoff", type=int, default=100,
            help="The cutoff of the finest grid level, default value: 100 Ry")

    gp.add_argument("--rel-cutoff", type=int, default=60,
            help="determines the grid at which a Gaussian is mapped, giving the cutoff used for a gaussian with alpha=1. A value 50+-10Ry might be required for highly accurate results, Or for simulations with a variable cell, default value: 60 Ry")

    gp.add_argument("--ngrids", type=int, default=4,
            help="The number of multigrids to use, default is 4")

    gp.add_argument("--kpoints-scheme", type=str,
            default="GAMMA",
            help="kpoint scheme to be used, can be NONE, GAMMA, MONKHORST-PACK, MACDONALD, GENERAL. when you set MONKHORST-PACK, you should also add the three integers like 'monkhorst-pack 3 3 3'")

    gp.add_argument("--kpath-manual", type=str, nargs="+", default=None,
            help="manual input kpath for band structure calculation")

    gp.add_argument("--kpath-file", type=str,
            help="file to read the kpath for band structure calculation")

    gp.add_argument("--diag", type=str, default="TRUE",
            choices=["TRUE", "FALSE", "true", "false"],
            help="whether choosing tranditional diagonalization for SCF")

    gp.add_argument("--ot", type=str, default="FALSE",
            choices=["TRUE", "FALSE", "true", "false"],
            help="whether choosing orbital transformation for SCF")

    gp.add_argument("--mixing-alpha", type=float, default=0.4,
            help="fraction of new density to be included, default is 0.4")

    gp.add_argument("--smear", type=str, default="FALSE",
            choices=["TRUE", "FALSE", "true", "false"],
            help="controls the activation of smearing")

    gp.add_argument("--smear-method", type=str, default="FERMI_DIRAC",
            help="Smearing method to be applied, can be fermi_dirac(default) or energy_window")

    gp.add_argument("--added-mos", type=int, default=0,
            help="Number of additional MOS added for each spin")

    gp.add_argument("--electronic-temp", type=float, default=300,
            help="Electronic temperature in the case of Fermi-Dirac smearing in unit of [K], default is 300")

    gp.add_argument("--window-size", type=float, default=0,
            help="size of the energy window centred at the Fermi level for energy_window type smearing")

    gp.add_argument("--ls-scf", type=str, default="false",
            choices=["true", "false", "true", "false"],
            help="use linear scaling scf method")

    # vdw correction related
    gp = subparser.add_argument_group(title="FORCE_EVAL/DFT/XC")

    gp.add_argument("--vdw-potential-type", type=str, default="NONE",
            choices=["PAIR_POTENTIAL", "NON_LOCAL", "NONE", "pair_potential", "non_local", "none"],
            help="Type of dispersion/vdW functional or potential to use")

    gp.add_argument("--pair-type", type=str, default="DFTD3",
            choices=["DFTD2", "DFTD3", "DFTD3(BJ)"],
            help="Type of potential(VDW)")

    gp.add_argument("--r-cutoff", type=float, default=1.05835442E+001,
            help="Range of potential. The cutoff will be 2 times this value")

    # printout option
    gp = subparser.add_argument_group(title="printout option")

    gp.add_argument("-p", "--printout-option", nargs="+", type=int,
            default=[],
            choices = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13],
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
           default is no printout of these properties.
           """)

    gp.add_argument("--dft-print-elf-cube-stride", type=int, nargs="+",
            default=[1, 1, 1],
            help="The stride (X,Y,Z) used to write the cube file (larger values result in smaller cube files). You can provide 3 numbers (for X,Y,Z) or 1 number valid for all components.")

    gp.add_argument("--dft-print-e-density-cube-stride", type=int, nargs="+",
            default=[1, 1, 1],
            help="The stride (X,Y,Z) used to write the cube file (larger values result in smaller cube files). You can provide 3 numbers (for X,Y,Z) or 1 number valid for all components.")

    # FORCE_EVAL/PROPERTIES
    gp = subparser.add_argument_group(title="FORCE_EVAL/PROPERTIES")

    gp.add_argument("--properties-resp-slab-sampling-range", type=float, nargs="+",
            default=[0.3, 3.0],
            help="Range where the fitting points are sampled. A range of 3 to 5 Angstroms means that the fitting points are sampled in the region of 3 to 5 Angstroms above the surface which is defined by atom indexes given in ATOM_LIST.")

    gp.add_argument("--properties-resp-slab-sampling-surf-direction", type=str, default="Z",
            choices=["X", "Y", "Z", "x", "y", "z", "-X", "-Y", "-Z", "-x", "-y", "-z"],
            help="Specifies what above the surface means. Defines the direction")

    gp.add_argument("--properties-resp-slab-sampling-atom-list", type=int, nargs="+",
            default=[1],
            help="Specifies the list of indexes of atoms used to define the region for the RESP fitting. The list should contain indexes of atoms of the first surface layer")

    # MOTION/CELL_OPT related parameters
    # --------------------------------------------------------------------------
    gp = subparser.add_argument_group(title="MOTION/CELL_OPT")

    gp.add_argument("--cell-opt-optimizer", type=str, default="BFGS",
            choices=["BFGS", "CG", "LBFGS"],
            help="Specify which method to use to perform a geometry optimization, can be: BFGS(default), CG, LBFGS")

    gp.add_argument("--cell-opt-max-iter", type=int, default=200,
            help="	Specifies the maximum number of geometry optimization steps. One step might imply several force evaluations for the CG and LBFGS optimizers")

    gp.add_argument("--cell-opt-type", type=str, default="DIRECT_CELL_OPT",
            choices=["DIRECT_CELL_OPT", "GEO_OPT", "MD"],
            help="Specify which kind of method to use for the optimization of the simulation cell, can be: DIRECT_CELL_OPT(default), GEO_OPT, MD")

    gp.add_argument("--cell-opt-max-dr", type=float, default=3e-3,
            help="Convergence criterion for the maximum geometry change between the current and the last optimizer iteration in unit of bohr, default is 3.e-3")

    gp.add_argument("--cell-opt-max-force", type=float, default=4.50000000E-004,
            help="Convergence criterion for the maximum force component of the current configuration in unit of bohr^-1*hartree")

    gp.add_argument("--cell-opt-rms-dr", type=float, default=1.50000000E-003,
            help="Convergence criterion for the root mean square (RMS) geometry change between the current and the last optimizer iteration in unit of bohr, default is 1.5e-3")

    gp.add_argument("--cell-opt-rms-force", type=float, default=3.00000000E-004,
            help="Convergence criterion for the root mean square (RMS) force of the current configuration, default is 3.e-4")

    gp.add_argument("--cell-opt-pressure-tolerance", type=float, default=1.00000000E+002,
            help="Specifies the Pressure tolerance (compared to the external pressure) to achieve during the cell optimization, default is 1.0e2")
        
    gp.add_argument("--cell-opt-keep-angles", type=str, default=None,
            help="Keep angles between the cell vectors constant, but allow the lenghts of the cell vectors to change independently.")

    gp.add_argument("--cell-opt-keep-symmetry", type=str, default=None,
            help="Keep the requested initial cell symmetry (e.g. during a cell optimisation). The initial symmetry must be specified in the &CELL section.")


    # MOTION/GEO_OPT related parameters
    # --------------------------------------------------------------------------
    gp = subparser.add_argument_group(title="MOTION/GEO_OPT")

    gp.add_argument("--geo-opt-optimizer", type=str, default="BFGS",
            choices=["BFGS", "CG", "LBFGS", "bfgs", "cg", "lbfgs", "cg"],
            help="Specify which method to use to perform a geometry optimization, can be: BFGS, CG, LBFGS")

    gp.add_argument("--geo-opt-max-iter", type=int, default=200,
            help="Specifies the maximum number of geometry optimization steps. One step might imply several force evaluations for the CG and LBFGS optimizers, default is 200")

    gp.add_argument("--geo-opt-type", type=str, default="MINIMIZATION",
            choices=["MINIMIZATION", "TRANSITION_STATE", "minimization", "transition_state"],
            help="Specify which kind of geometry optimization to perform, default is MINIMIZATION")

    gp.add_argument("--geo-opt-max-dr", type=float, default=3e-3,
            help="Convergence criterion for the maximum geometry change between the current and the last optimizer iteration, default is 3e-3")

    gp.add_argument("--geo-opt-max-force", type=float, default=4.50000000E-004,
            help="Convergence criterion for the maximum force component of the current configuration, default is 4.0e-4")

    gp.add_argument("--geo-opt-rms-dr", type=float, default=1.50000000E-003,
            help="Convergence criterion for the root mean square (RMS) geometry change between the current and the last optimizer iteration, default is 1.5e-3")

    gp.add_argument("--geo-opt-rms-force", type=float, default=3.00000000E-004,
            help="Convergence criterion for the root mean square (RMS) force of the current configuration.")

    # MOTION/MD 
    # --------------------------------------------------------------------------
    gp = subparser.add_argument_group(title="MOTION/MD")

    gp.add_argument("--md-steps", type=int, default=1000,
            help="The number of MD steps to perform")

    gp.add_argument("--timestep", type=float, default=5.0e-1,
            help="The length of an integration step (in case RESPA the large TIMESTEP), default and also recommended is 0.5 fs.")

    gp.add_argument("--ensemble", type=str, default="NVE",
            choices=["NVE",  "NVT","HYDROSTATICSHOCK", "ISOKIN", "LANGEVIN", "MSST", "MSST_DAMPED"],
            help="The ensemble/integrator that you want to use for MD propagation")

    gp.add_argument("--temperature", type=str, default=300,
            help="The temperature in K used to initialize the velocities with init and pos restart, and in the NPT/NVT simulations")

    gp.add_argument("--temp-tol", type=float, default=0.0,
            help="The maximum accepted deviation of the (global) temperaturefrom the desired target temperature before a rescaling of the velocites is performed. If it is 0 no rescaling is performed. NOTE: This keyword is obsolescent; Using a CSVR thermostat with a short timeconstant is recommended as a better alternative")

    gp.add_argument("--traj-format", type=str, default="XMOL",
            help="type of output trajectory for MOTION, note: DCD format can be visualized by vmd")

    # MOTION/BAND
    # --------------------------------
    gp = subparser.add_argument_group(title="MOTION/BAND")

    gp.add_argument("--band-type", type=str, default="CI-NEB",
            help="specify the type of band calculation")

    gp.add_argument("--number-of-replica", type=int, default=5,
            help="number of replicas")

    gp.add_argument("--k-spring", type=float, default=2.0e-2,
            help="value of the spring constant")

    gp.add_argument("--align-frames", type=str, default="TRUE",
            choices=["TRUE", "FALSE", "true", "false"],
            help="Enables the alignment of the frames at the beginning of a BAND calculation. This keyword does not affect the rotation of the replicas during a BAND calculation.")

    gp.add_argument("--rotate-frames", type=str, default="TRUE",
            choices=["TRUE", "FALSE", "true", "false"],
            help="Compute at each BAND step the RMSD and rotate the frames in order to minimize it.")

    #             vibrational_analysis related parameters
    # ---------------------------------------------------------------
    gp = subparser.add_argument_group(title="VIBRATIONAL_ANALYSIS")

    gp.add_argument("--dx", type=float, default=1.0e-2,
            help="specify the increment to be used to construct the HESSIAN with finite difference method")

    gp.add_argument("--fully-periodic", type=str, default="FALSE",
            choices=["TRUE", "FALSE", "true", "false"],
            help="avoids to clean rotations from the Hessian matrix")

    gp.add_argument("--intensities", type=str, default="FALSE",
            choices=["TRUE", "FALSE", "true", "false"],
            help="Calculation of the IR-Intensities. Calculation of dipoles has to be specified explicitly"
            )

    gp.add_argument("--tc-pressure", type=float, default=1.01325000E+005,
            help="Pressure for the calculation of the thermochemical data in unit of [Pa]")

    gp.add_argument("--tc-temperature", type=float, default=2.73150000E+002,
            help="Temperature for the calculation of the thermochemical data in unit of [K]")

    gp.add_argument("--thermochemistry", type=str, default="FALSE",
            help="Calculation of the thermochemical data. Valid for molecules in the gas phase.")

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

    # --------------------------------------------------------------------------
    # Quantum ESPRESSO
    # --------------------------------------------------------------------------
    subparser = subparsers.add_parser("qe", help="using quantum espresso as calculator")

    gp = subparser.add_argument_group(title="overall running control")

    gp.add_argument("-r", "--runtype", type=int, default=0,
            choices=[0, 1, 2, 3, 4, 5, 6, 7, 8, 9],
            help="choices of runtype. 0->static_run; 1->relax; 2->vc-relax; 3->cubic-cell; 4->hexagonal-cell; 5->tetragonal-cell; 6->neb; 7->dfpt; 8->phonopy; 9->pp.x")

    gp.add_argument("--static", type=str, default="all",
            choices=["all", "scf"],
            help="in case of all(default), run scf, nscf, bands in a single run; in case of scf, run scf only")


    gp.add_argument("-d", "--directory", type=str, default="matflow-running",
            help="Directory for the running.")

    # calypso input.dat template
    gp = subparser.add_argument_group(title="template", 
            description="read in Calypso input.dat template")

    gp.add_argument("--input-dat", type=str, default=None,
            help="specify Calypso input.dat template to set parameters")

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
            help="specify type of pseudo potentials to prepare, when --pot auto")

    # -------------------------------------------------------------------
    #                       scf related parameters
    # -------------------------------------------------------------------
    # &control
    gp = subparser.add_argument_group(title="pw.x->control")

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


    # magnetic related parameters
    gp.add_argument("--nspin", type=int, default=None,
            choices=[1, 2],
            help="choose either 1 or 2, and 4 should not be used as suggested by pwscf official documentation.")

    gp.add_argument("--starting-magnetization", type=float, nargs="+", default=None,
            help="starting_magnetization(i), i=1,ntyp -> Starting spin polarization on atomic type i in a spin polarized calculation. Values range between -1 (all spins down for the valence electrons of atom type i) to 1 (all spins up).")

    gp.add_argument("--noncolin", type=str, default=None,
            choices=[".true.", ".false."],
            help="if .true. the program will perform a noncollinear calculation.")
    # &electrons
    gp = subparser.add_argument_group(title="pw.x->electrons")

    gp.add_argument("--conv-thr", type=float, default=1.0e-6,
            help="the conv_thr for scf, when doing geometric optimization better use a strict covnergec for scf")

    # &ions
    gp = subparser.add_argument_group(title="pw.x->ions")

    gp.add_argument("--etot-conv-thr",
            type=float, default=1.0e-4,
            help="convergence threshold of energy for geometric optimization")

    gp.add_argument("--forc-conv-thr",
            type=float, default=1.0e-3,
            help="convergence threshold for force in optimization,(usually it is more important than energy)")

    gp.add_argument("--nstep",
            type=int, default=50,
            help="maximum ion steps for geometric optimization")
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

    gp.add_argument("--kpoints-mp-scf", type=int, nargs=6,
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

    # --------------------------------------------------------------------------
    # SIESTA
    # --------------------------------------------------------------------------
    subparser = subparsers.add_parser("siesta", help="using siesta as calculator")

    gp = subparser.add_argument_group(title="overall running control:")

    gp.add_argument("-r", "--runtype", type=int, default=0,
            choices=[0, 1, 2, 3, 4, 5],
            help="choices of runtype. 0->static_run; 1->optimization; 2->cubic-cell; 3->hexagonal-cell; 4->tetragonal-cell; 5->phonopy; 6->molecular dynamics")

    gp.add_argument("-d", "--directory", type=str, default="matflow-running",
            help="Directory for the running.")

    gp.add_argument("--runopt", type=str, default="gen",
            choices=["gen", "run", "genrun"],
            help="Generate or run or both at the same time.")

    gp.add_argument("--auto", type=int, default=3,
            choices=[0, 1, 2, 3],
            help="auto:0 nothing, 1: copying files to server, 2: copying and executing, 3: pymatflow run inserver with direct submit,  in order use auto=1, 2, you must make sure there is a working ~/.pymatflow/server_[pbs|llhpc].conf")

    # calypso input.dat template
    gp = subparser.add_argument_group(title="template", 
            description="read in Calypso input.dat template")

    gp.add_argument("--input-dat", type=str, default=None,
            help="specify Calypso input.dat template to set parameters")

    # -----------------------------------------------------------------
    #                      run params
    # -----------------------------------------------------------------

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
            choices=[0, 1, 2, 3, 4, 5, 6, 7, 8],
            help="choices of runtype. 0->static_run; 1->optimization; 2->cubic-cell; 3->hexagonal-cell; 4->tetragonal-cell; 5->neb; 6->vasp-phonon; 7->phonopy; 8->surf pes")


    # calypso input.dat template
    gp = subparser.add_argument_group(title="template", 
            description="read in Calypso input.dat template")

    gp.add_argument("--input-dat", type=str, default=None,
            help="specify Calypso input.dat template to set parameters")
            
            
            
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

    gp.add_argument("--nwrite", type=int, nargs="+", default=None,
            help=" This flag determines how much will be written to the file OUTCAR (verbosity flag)")

    gp.add_argument("--prec", type=str, nargs="+", default=None,
            choices=["Normal", "Accurate", "A", "N"],
            help="PREC, default value: Normal")

    gp.add_argument("--ncore", type=int, nargs="+", default=None,
            help="NCORE determines the number of compute cores that work on an individual orbital ")

    # incar->electrons
    gp = subparser.add_argument_group(title="incar->electron",
            description="electrons calculation related parameters")

    gp.add_argument("--encut", type=int, nargs="+", default=None,
            help="ENCUT, default value: 300 eV")

    gp.add_argument("--ediff", type=float, nargs="+", default=None,
            help="EDIFF, default value: 1.0e-4")

    gp.add_argument("--nelm", type=int, nargs="+", default=None,
            help="NELM sets the maximum number of electronic SC (selfconsistency) steps which may be performed")

    gp.add_argument("--nfree", type=int, nargs="+", default=None,
            help="NFREE specifies the number of remembered steps in the history of ionic convergence runs, or the number of ionic displacements in frozen phonon calculations")

    gp.add_argument("--kpoints-mp", type=int, nargs="+",
            default=[1, 1, 1, 0, 0, 0],
            help="set kpoints like -k 1 1 1 0 0 0")

    gp.add_argument("--kspacing", type=float, nargs="+", default=None,
            help="determines the number of k-points if the KPOINTS file is not present. default is 0.5")
            
    #gp.add_argument("--kpoints-mp-scf", type=int, nargs="+",
    #        default=[1, 1, 1, 0, 0, 0],
    #        help="set kpoints like -k 1 1 1 0 0 0")

    #gp.add_argument("--kpoints-mp-nscf", type=int, nargs="+",
    #        default=[3, 3, 3, 0, 0, 0],
    #        help="set kpoints like -k 1 1 1 0 0 0")

    gp.add_argument("--kpath-manual", type=str, nargs="+", default=None,
            help="set kpoints for band structure calculation manually")

    gp.add_argument("--kpath-file", type=str, nargs="+", default="kpath.txt",
            help="set kpoints for band structure calculation manually from file")

    gp.add_argument("--kpath-intersections", type=int, nargs="+", default=15,
            help="intersection of the line mode kpoint for band calculation")

    gp.add_argument("--ismear", type=int, nargs="+", default=None,
            help="smearing type(methfessel-paxton(>0), gaussian(0), fermi-dirac(-1), tetra(-4), tetra-bloch-dorrected(-5)), default: 0")

    gp.add_argument("--sigma", type=float, nargs="+", default=None,
            help="determines the width of the smearing in eV.")

    gp.add_argument("--ivdw", type=int, nargs="+", default=None,
            choices=[0, 11, 12, 21, 202, 4],
            help="IVDW = 0(no correction), 1(dft-d2), 11(dft-d3 Grimme), 12(dft-d3 Becke-Jonson), default: None which means 0, no correction")
    # -----------------------------------------------------------------

    gp.add_argument("--lorbit", type=int, nargs="+", default=None,
            choices=[0, 1, 2, 5, 10, 11, 12],
            help="together with an appropriate RWIGS, determines whether the PROCAR or PROOUT files are written")

    # optics related
    gp.add_argument("--loptics", type=str, nargs="+", default=None,
            choices=["TRUE", "FALSE"],
            help="calculates the frequency dependent dielectric matrix after the electronic ground state has been determined.")

    gp.add_argument("--cshift", type=float, nargs="+", default=None,
            help="CSHIFT sets the (small) complex shift η in the Kramers-Kronig transformation")

    gp.add_argument("--nedos", type=int, nargs="+", default=None,
            help="NEDOS specifies number of gridpoints on which the DOS is evaluated")
            
    # magnetic related
    gp.add_argument("--ispin", type=int, nargs="+", default=None,
            choices=[1, 2],
            help="specifies spin polarization: 1->no spin polarized, 2->spin polarized(collinear). combine SIPIN with MAGMOM to study collinear magnetism.")

    gp.add_argument("--magmom", type=float, nargs="+", default=None,
            help="Specifies the initial magnetic moment for each atom, if and only if ICHARG=2, or if ICHARG=1 and the CHGCAR file contains no magnetisation density.")

    gp.add_argument("--lnoncollinear", type=str, nargs="+", default=None,
            choices=["T", "F", ".TRUE.", ".FALSE."],
            help="specifies whether fully non-collinear magnetic calculations are performed")

    gp.add_argument("--lsorbit", type=str, nargs="+", default=None,
            choices=["T", "F", ".TRUE.", ".FALSE."],
            help="specifies whether spin-orbit coupling is taken into account.")

    gp.add_argument("--saxis", type=float, nargs="+", default=None,
            help="SAXIS specifies the quantisation axis for noncollinear spins")

    gp.add_argument("--lmaxmix", type=int, nargs="+", default=None,
            help="LMAXMIX controls up to which l-quantum number the one-center PAW charge densities are passed through the charge density mixer and written to the CHGCAR file.")

    # hybrid functional
    gp = subparser.add_argument_group(title="incar->Exchange correlation")
    gp.add_argument("--lhfcalc", type=str, nargs="+", default=None,
            choices=["T", "F", ".TRUE.", ".FALSE."],
            help=" specifies whether Hartree-Fock/DFT hybrid functional type calculations are performed")

    gp.add_argument("--hfscreen", type=float, nargs="+", default=None,
            choices=[0.3, 0.2],
            help=" specifies the range-separation parameter in range separated hybrid functionals: HSE03->0.3, HSE06->0.2, must also set LHFCALC=.TRUE.")

    gp.add_argument("--aexx", type=float, nargs="+", default=None,
            help="AEXX specifies the fraction of exact exchange in a Hartree-Fock/DFT hybrid functional type calculation")

    gp.add_argument("--lsubrot", type=str, nargs="+", default=None,
            choices=["T", "F", ".TRUE.", ".FALSE."],
            help="This flag can be set for hybrid functionals (HF-type calculations).")

    gp.add_argument("--nsw", type=int, nargs="+", default=None,
            help="NSW sets the maximum number of ionic steps")

    gp.add_argument("--ediffg", type=float, nargs="+", default=None,
            help="EDIFFG, default value: 10*EDIFF")

    gp = subparser.add_argument_group(title="incar->ions",
            description="setting ions related parameters")

    gp.add_argument("--ibrion", type=int, nargs="+", default=None,
            choices=[-1, 0, 1, 2, 3, 5, 6, 7, 8, 44],
            help="IBRION = refer to https://cms.mpi.univie.ac.at/wiki/index.php/IBRION for how to set the algorithm of optimization you need!")

    gp.add_argument("--isif", type=int, nargs="+", default=None,
            choices=[0, 1, 2, 3, 4, 5, 6, 7],
            help="ISIF = 0-7: refer to https://cms.mpi.univie.ac.at/wiki/index.php/ISIF for how to set the type of Geometri Optimization you need!")

    gp.add_argument("--potim", type=float, nargs="+", default=None,
            help="step width scaling (ionic relaxations), default: None = 0.015 in phonon calculation")

    gp.add_argument("--selective-dynamics", type=str, nargs="+", default="False",
            choices=["True", "False", "T", "F"],
            help="whether use selective dyanmics")

    # incar-miscellaneous
    gp = subparser.add_argument_group(title="incar-miscellaneous",
            description="miscellaneous input parameters")

    gp.add_argument("--algo", type=str, nargs="+", default=None,
            choices=["N", "D", "V", "F"],  #"Exact", "G0W0", "GW0", "GW"],
            help=" a convenient option to specify the electronic minimisation algorithm (as of VASP.4.5) and/or to select the type of GW calculations")

    gp.add_argument("--ialgo", type=int, nargs="+", default=None,
            choices=[5, 6, 7, 8, 38, 44, 46, 48],
            help="IALGO selects the algorithm used to optimize the orbitals.Mind: We strongly urge the users to select the algorithms via ALGO. Algorithms other than those available via ALGO are subject to instabilities.")

    gp.add_argument("--addgrid", type=str, nargs="+", default=None,
            choices=[".TRUE.", ".FALSE.", "T", "F"],
            help="ADDGRID determines whether an additional support grid is used for the evaluation of the augmentation charges.")

    gp.add_argument("--isym", type=int, nargs="+", default=None,
            choices=[-1, 0, 1, 2, 3],
            help=" ISYM determines the way VASP treats symmetry.")

    gp.add_argument('--lreal', type=str, nargs="+", default=None,
            choices=["T", "F", ".TRUE.", ".FALSE.", "O", "On", "A", "Auto"],
            help="LREAL determines whether the projection operators are evaluated in real-space or in reciprocal space.")

    gp.add_argument("--pstress", type=float, nargs="+", default=None,
            help="controls whether Pulay corrections are added to the stress tensor or not.")

    # properties parameters
    gp.add_argument("--lelf", type=str, nargs="+", default=None,
            choices=["T", "F", ".TRUE.", ".FALSE."],
            help="LELF determines whether to create an ELFCAR file or not.")

    # write PARAMETERS
    gp = subparser.add_argument_group(title="incar->write parameters",
            description="set writing parameters")

    gp.add_argument("--lwave", type=str, nargs="+", default=None,
            choices=['T', 'F', ".TRUE.", '.FALSE.'],
            help="LWAVE determines whether the wavefunctions are written to the WAVECAR file at the end of a run.")

    gp.add_argument("--lcharg", type=str, nargs="+", default=None,
            choices=['T', 'F', ".TRUE.", '.FALSE.'],
            help="LCHARG determines whether the charge densities (files CHGCAR and CHG) are written.")


    #                     neb related PARAMETERS
    # --------------------------------------------------------------------------


    # PHONOPY parameters
    # ----------------------------------------



    # range_a range_c
    # ----------------------------------------------

    # incar template
    gp = subparser.add_argument_group(title="template", 
            description="read in INCAR template")

    gp.add_argument("--incar", type=str, nargs="+", default=None,
            help="specify incar template to set parameters")

    # fix atoms
    gp = subparser.add_argument_group(title="fix atoms",
            description="specify atoms to fix in optimization, only used when --runtype=1")

    gp.add_argument("--fix", help="list of fixed atoms", nargs='+', type=int)


    # static calc related setting
    gp = subparser.add_argument_group(title="static calc",
            description="setting type of static calculation when -r 0")

    gp.add_argument("--static", type=str, default="band",
            choices=["scf", "band", "optics", "bse"],
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

    # miscellaneous
    gp = subparser.add_argument_group(title="miscellaneous",
            description="miscallaneous setting")
    
    gp.add_argument("--symprec", type=float, default=None,
            help="determines how accurately the positions in the POSCAR file must be specified. The default, SYMPREC=10-5, is usually large enough, even if the POSCAR file has been generated with single precision accuracy. Increasing SYMPREC means that the positions in the POSCAR file can be specified with less accuracy (increasing fuzziness).")

    gp.add_argument("--amix", type=float, default=None,
            help="specifies the linear mixing parameter.")

    gp.add_argument("--bmix", type=float, default=None,
            help="sets the cutoff wave vector for Kerker mixing scheme")            

    # Control of calypso
    gp = subparser.add_argument_group(title="psoflow control",
            description="psoflow controlled setting")

    gp.add_argument("--gen-incar-n", type=int, default=1,
            help="control number of INCAR_* to generate")
            
    gp.add_argument("--split-batch", type=int, default=None,
            help="number of structure to calculate in one submit script, (several submit script are generated and submit at the same time)")

    # Calypso General
    gp = subparser.add_argument_group(title="Calypso General", 
            description="Claypso input.dat parameters")
            
    gp.add_argument("--systemname", type=str, default="PsoFlow",
            help="a descriptive name of the system")

    gp.add_argument("--numberofspecies", type=int, default=None, required=True)
    
    gp.add_argument("--nameofatoms", type=str, nargs="+", default=None, required=True)
    
    gp.add_argument("--numberofatoms", type=int, nargs="+", default=None, required=True)
    
    gp.add_argument("--numberofformula", type=int, nargs=2,
            default=None,
            help="define number of formulas used in simulation, default is 0 4")
            
    gp.add_argument("--volume", type=float, default=None,
            help="volume, calypso can set it so we do not need to set it")    
            
    gp.add_argument("--distanceofion", type=float, nargs="+", default=None,
            help="limit of the minimum distance between atoms, default is 0.7")
            
    gp.add_argument("--ialgo-pso", type=int, default=None, 
            help="optimization algorithm used in structure evolution, default is 2->local swarm particle group")
            
    gp.add_argument("--icode", type=int, default=1, 
            help="software used to relax the structures, default is 1 vasp")
            
    gp.add_argument("--numberoflocaloptim", type=int, default=None,
            help="time of optimization for every structure in the search, default is 3")
            
            
    gp.add_argument("--psoratio", type=float, default=None,
            help="ratio of number of structures generated by swarm algorithm against total structures, default is 0.6")
            
    gp.add_argument("--popsize", type=int, default=None,
            help="size of every group in structure evolution, default is 30")
            
    gp.add_argument("--kgrid", type=float, nargs=2,
            default=None,
            help="only used in vasp(version <5.3) and Quantum-Espresso and DFTB, default is 0.12 0.06")
            
    gp.add_argument("--maxstep", type=int, default=None,
            help="maximum generations for structure evolution")
            
    gp.add_argument("--pickup", type=str, default=None,
            choices=["True", "False"],
            help="whether to pick up the calculation from generation defined by PickStep")
            
    gp.add_argument("--pickstep", type=str, default=None,
            help="prediction from accidentally interrrupted generation")
            
    gp.add_argument("--maxtime", type=int, default=None,
            help="maximum time for structure relaxation, only apply to vasp now")
            
    gp.add_argument("--lmc", type=str, default=None, 
            help="whether to add Metropolis energy criteria in structure evolution")
            
    # Calypso 2D layer materials
    gp = subparser.add_argument_group(title="Calypso 2D layer materials",
            description="parameters needed in 2D layer structure prediction")
            
    gp.add_argument("--pso-2d", type=str, default=None, 
            help="whether do 2D layer structure prediction")
            
    gp.add_argument("--area", type=float, default=None,
            help="area of every formula, no default value, only has effect on the optimization effctiveness in first generation evolution")
            
    gp.add_argument("--deltaz", type=float, default=None,
            help="deformation along c, used to design wrinkle laer structure")
            
            
    gp.add_argument("--multilayer", type=int, default=None,
            help="number of layers in the system, default is 1")
            
            
    gp.add_argument("--layergap", type=float, default=None,
            help="distance between different layer")
            
            
    gp.add_argument("--vacuumgap", type=float, default=None, 
            help="vacuum gap, default is 10")
            
            
    gp.add_argument("--layertype", type=int, nargs="+", 
            help="a matrice of m x n. m is the number of layers in the system, n is the number of elements.")
    
    gp.add_argument("--latom-dis", type=float, default=None,
            help="minimum distance between atoms in layer structure")
            
    # Calypso cluster
    gp = subparser.add_argument_group(title="Calypso cluster",
            description="parameters needed in cluster prediction")
                       
    gp.add_argument("--cluster", type=str, default=None,
            choices=["True", "False"],
            help="whether to do cluster structure evolution")
   
    gp.add_argument("--vacancy", type=float, nargs=3,
            default=None,
            help="cell to put the cluster, default is 10.0 10.0 10.0")
            
    # Calypso rigid molecule
    gp = subparser.add_argument_group(title="Calypso rigid",
            description="Rigid Molecule structure search")
            
    gp.add_argument("--mol", type=str, default=None,
            choices=["True", "False"],
            help="whether to do rigid molecule structure prediction")
    
    # Calypso special constraint
    gp = subparser.add_argument_group(title="Calypso constraint", 
            description="Special Constraint in Calypso")
            
    # Calypso variation composition
    gp = subparser.add_argument_group(title="Calypso variation composition", 
            description="Variation Composition structure search")
            
    # Calypso surface stucture 
    gp = subparser.add_argument_group(title="Calypso surface structure", 
            description="Surface structure search")
            
    # Calypso hardness desing 
    gp = subparser.add_argument_group(title="Calypso hardness", 
            description="Hardness design structure search")

    # Calypso gap design
    gp = subparser.add_argument_group(title="Calypso gap", 
            description="Gap design structure search")
            
    # Calypso 2D adsorption 
    gp = subparser.add_argument_group(title="Calypso 2D adsorption", 
            description="2D adsorption structure prediction")

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


    # we use args.nameofatoms to construct a virtual xyz file in /tmp/matflow and use it to generate POTCAR with pot-from-xyz-modified.py
    os.system("mkdir -p /tmp/psoflow")
    with open("/tmp/psoflow/tmp.xyz", "w") as fout:
        fout.write("%d\n" % len(args.nameofatoms))
        fout.write("cell: 10 0 0 | 0 10 0 | 0 0 10\n")
        for name in args.nameofatoms:
            fout.write("%s 0.0 0.0 0.0\n" % name)
    fake_xyzfile = os.path.join("/tmp/psoflow", "tmp.xyz")
    
    # dealing with pseudo potential file
    if args.pot == "./":
        #TODO make a simple check, whether there exists the potential file
        pass
    elif args.pot == "auto":
        if args.driver == "abinit":
            os.system("pot-from-xyz-modified.py -i %s -d ./ -p abinit --abinit-type=ncpp" % fake_xyzfile)
        elif args.driver == "qe":
            if args.runtype == 6:
                os.system("pot-from-xyz-modified.py -i %s -d ./ -p qe --qe-type=%s" % (images[0]), args.pot_type)
            else:
                os.system("pot-from-xyz-modified.py -i %s -d ./ -p qe --qe-type=%s" % (fake_xyzfile, args.pot_type))
        elif args.driver == "siesta":
            print("=============================================================\n")
            print("                     WARNING\n")
            print("-------------------------------------------------------------\n")
            print("support for auto preparation of pseudopotential file for siesta\n")
            print("is not fully implemented now!\n")
            print("please prepare it yourself\n")
            sys.exit(1)
        elif args.driver == "vasp":
            os.system("vasp-potcar-from-xyz.py --type %s -i %s -o ./POTCAR" % (args.pot_type, fake_xyzfile))
    else:
        os.system("cp %s/* ./" % args.pot)



    # deal with Calypso input.dat template specified by --input-dat
    input_dat_params = {}
    if args.input_dat == None:
        pass
    else:
        if not os.path.exists(args.input_dat):
            print("====================================================\n")
            print("                  Warning !!!!\n")
            print("----------------------------------------------------\n")
            print("psoflow:\n")
            print("the specified input.dat file by --input-dat doesn't exist\n")
            print("go and check it\n")
            sys.exit(1)
        with open(args.input_dat, 'r') as fin:
            lines = fin.readlines()
        for i in range(len(lines)):
            if len(lines[i].split()) == 0:
                continue
            if lines[i][0] == "#":
                continue
            if lines[i].split()[0][0] == "@" and lines[i].split()[0] != "@End":
                key = lines[i].split("\n")[0].split("#")[0].split("@")[1]
                value = []
                nrow = 1
                while lines[i+nrow+1].split()[0] != "@End":
                    nrow += 1
                ncol = len(lines[i+1].split("\n")[0].split())
                for row in range(nrow):
                    vec = []
                    for col in range(ncol):
                        vec.append(lines[i+row+1].split("\n")[0].split()[col])
                    value.append(vec)
                input_dat_params[key] = value
            elif "@" not in lines[i] and "=" not in lines[i]:
                continue
            elif len(lines[i].split("\n")[0].split("#")[0].split("=")) == 2 and len(lines[i].split("\n")[0].split("#")[0].split("=")[1].split()) == 1:
                # in case of single value input parameters
                key = lines[i].split("=")[0].split()[0]
                value = lines[i].split("\n")[0].split("#")[0].split("=")[1].split()[0]
                input_dat_params[key] = value
            elif len(lines[i].split("\n")[0].split("#")[0].split("=")) == 2 and len(lines[i].split("\n")[0].split("#")[0].split("=")[1].split()) > 1:
                key = lines[i].split("=")[0].split()[0]
                value = lines[i].split("\n")[0].split("#")[0].split("=")[1].split()
                input_dat_params[key] = value
            else:
                pass
        #
        # if xxx is alraedy in input_dat_params(set from --input_dat_params) and args.xxx is None
        # input_dat_params[xxx] will not be control by args.xxx
        input_dat_params["SystemName"] = args.systemname if "SystemName" not in input_dat_params or args.systemname != None else input_dat_params["SystemName"]
        input_dat_params["NameOfAtoms"] = args.nameofatoms if "nameofatoms" not in input_dat_params or args.nameofatoms != None else input_dat_params["NameOfAtoms"]
        input_dat_params["NumberOfSpecies"] = args.numberofspecies if "NumberOfSpcies" not in input_dat_params or args.numberofspecies != None else input_dat_params["NumberOfSpecies"]
        input_dat_params["NumberOfAtoms"] = args.numberofatoms if "NumberOfAtoms" not in input_dat_params or args.numberofatoms != None else input_dat_params["NumberOfAtoms"]
        input_dat_params["NumberOfFormula"] = args.numberofformula if "NumberOfFormula" not in input_dat_params or args.numberofformula != None else input_dat_params["NumberOfFormula"]
        input_dat_params["Volume"] = args.volume if "Volume" not in input_dat_params or args.volume != None else input_dat_params["Volume"]
        input_dat_params["DistanceOfIon"] = args.distanceofion if "DistanceOfIon" not in input_dat_params or args.distanceofion != None else input_dat_params["DistanceOfIon"]
        input_dat_params["IAlgo"] = args.ialgo_pso if "IAlgo" not in input_dat_params or args.ialgo_pso != None else input_dat_params["IAlgo"]
        input_dat_params["ICode"] = args.icode if "ICode" not in input_dat_params or args.icode != None else input_dat_params["ICode"]
        input_dat_params["NumberOfLocalOptim"] = args.numberoflocaloptim if "NumberOfLocalOptim" not in input_dat_params or args.numberoflocaloptim != None else input_dat_params["NumberOfLocalOptim"]
        input_dat_params["PsoRatio"] = args.psoratio if "PsoRatio" not in input_dat_params or args.psoratio != None else input_dat_params["PsoRatio"]
        input_dat_params["PopSize"] = args.popsize if "PopSize" not in input_dat_params or args.popsize != None else input_dat_params["PopSize"]
        input_dat_params["Kgrid"] = args.kgrid if "Kgrid" not in input_dat_params or args.kgrid != None else input_dat_params["Kgrid"]
        input_dat_params["MaxStep"] = args.maxstep if "MaxStep" not in input_dat_params or args.maxstep != None else input_dat_params["MaxStep"]
        input_dat_params["PickUp"] = args.pickup if "PickUp" not in input_dat_params or args.pickup != None else input_dat_params["PickUp"]
        input_dat_params["PickStep"] = args.pickstep if "PickStep" not in input_dat_params or args.pickstep != None else input_dat_params["PickStep"]
        input_dat_params["MaxTime"] = args.maxtime if "MaxTime" not in input_dat_params or args.maxtime != None else input_dat_params["MaxTime"]
        input_dat_params["LMC"] = args.lmc if "LMC" not in input_dat_params or args.lmc != None else input_dat_params["LMC"]
        input_dat_params["2D"] = args.pso_2d if "2D" not in input_dat_params or args.pso_2d != None else input_dat_params["2D"]
        input_dat_params["Area"] = args.area if "Area" not in input_dat_params or args.area != None else input_dat_params["Area"]
        input_dat_params["DeltaZ"] = args.deltaz if "DeltaZ" not in input_dat_params or args.deltaz != None else input_dat_params["DeltaZ"]
        input_dat_params["MultiLayer"] = args.multilayer if "MultiLayer" not in input_dat_params or args.multilayer != None else input_dat_params["MultiLayer"]
        input_dat_params["LayerGap"] = args.layergap if "LayerGap" not in input_dat_params or args.layergap != None else input_dat_params["LayerGap"]
        input_dat_params["VacuumGap"] = args.vacuumgap if "VacuumGap" not in input_dat_params or args.vacuumgap != None else input_dat_params["VacuumGap"]
        input_dat_params["LayerType"] = args.layertype if "LayerType" not in input_dat_params or args.layertype != None else input_dat_params["LayerType"]
        input_dat_params["LAtom_Dis"] = args.latom_dis if "LAtom_Dis" not in input_dat_params or args.latom_dis != None else input_dat_params["LAtom_Dis"]
        input_dat_params["Cluster"] = args.cluster if "SystemName" not in input_dat_params or args.systemname != None else input_dat_params["SystemName"]
        input_dat_params["Vacancy"] = args.vacancy if "Vacancy" not in input_dat_params or args.vacancy != None else input_dat_params["Vacancy"]
        input_dat_params["Mol"] = args.mol if "Mol" not in input_dat_params or args.mol != None else input_dat_params["Mol"]
        

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
            from pymatflow.abinit.opt import opt_run
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
            from pymatflow.abinit.opt import opt_run
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
            from pymatflow.abinit.opt import opt_run
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
            from pymatflow.abinit.opt import opt_run
            task = opt_run()
            task.get_xyz(xyzfile)
            task.set_params(params=params)
            task.set_kpoints(kpoints=kpoints)
            task.set_run(mpi=args.mpi, server=args.server, jobname=args.jobname, nodes=args.nodes, ppn=args.ppn, queue=args.queue)
            task.set_llhpc(partition=args.partition, nodes=args.nodes, ntask=args.ntask, jobname=args.jobname, stdout=args.stdout, stderr=args.stderr)
            task.tetragonal(directory=args.directory, runopt=args.runopt, auto=args.auto, na=args.na, stepa=args.stepa, nc=args.nc, stepc=args.stepc)
        elif args.runtype == 5:
            # dfpt-elastic-piezo-dielec
            from pymatflow.abinit.dfpt import dfpt_elastic_piezo_dielec
            task = dfpt_elastic_piezo_dielec()
            task.get_xyz(xyzfile)
            task.set_params(params=params)
            task.set_kpoints(kpoints=kpoints)
            task.set_run(mpi=args.mpi, server=args.server, jobname=args.jobname, nodes=args.nodes, ppn=args.ppn, queue=args.queue)
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
            task.set_run(mpi=args.mpi, server=args.server, jobname=args.jobname, nodes=args.nodes, ppn=args.ppn, queue=args.queue)
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
            task.set_run(mpi=args.mpi, server=args.server, jobname=args.jobname, nodes=args.nodes, ppn=args.ppn, queue=args.queue)
            task.set_llhpc(partition=args.partition, nodes=args.nodes, ntask=args.ntask, jobname=args.jobname, stdout=args.stdout, stderr=args.stderr)
            task.phonopy(directory=args.directory, runopt=args.runopt, auto=args.auto)
        else:
            pass
# ==============================================================================
# CP2K CP2K CP2K CP2K CP2K CP2K CP2K CP2K CP2K CP2K CP2K CP2K CP2K CP2K C2PK CP2K
# ==============================================================================
    elif args.driver == "cp2k":
        params = {}

        params["GLOBAL-PRINT_LEVEL"] = args.print_level

        params["FORCE_EVAL-SUBSYS-CELL-SYMMETRY"] = args.cell_symmetry

        params["FORCE_EVAL-DFT-LS_SCF"] = args.ls_scf
        params["FORCE_EVAL-DFT-QS-METHOD"] = args.qs_method
        params["FORCE_EVAL-DFT-MGRID-CUTOFF"] = args.cutoff
        params["FORCE_EVAL-DFT-MGRID-REL_CUTOFF"] = args.rel_cutoff
        params["FORCE_EVAL-DFT-MGRID-NGRIDS"] = args.ngrids
        params["FORCE_EVAL-DFT-XC-XC_FUNCTIONAL"] = args.xc_functional
        params["FORCE_EVAL-DFT-SCF-EPS_SCF"] = args.eps_scf
        params["FORCE_EVAL-DFT-SCF-ADDED_MOS"] = args.added_mos
        params["FORCE_EVAL-DFT-SCF-SMEAR"] = args.smear
        params["FORCE_EVAL-DFT-SCF-SMEAR-METHOD"] = args.smear_method
        params["FORCE_EVAL-DFT-SCF-SMEAR-ELECTRONIC_TEMPERATURE"] = args.electronic_temp
        params["FORCE_EVAL-DFT-SCF-SMEAR-WINDOW_SIZE"] = args.window_size
        params["FORCE_EVAL-DFT-SCF-DIAGONALIZATION"] = args.diag
        params["FORCE_EVAL-DFT-SCF-OT"] = args.ot
        params["FORCE_EVAL-DFT-SCF-MIXING-ALPHA"] = args.mixing_alpha
        params["FORCE_EVAL-DFT-KPOINTS-SCHEME"] = args.kpoints_scheme
        params["FORCE_EVAL-DFT-XC-VDW_POTENTIAL-POTENTIAL_TYPE"] = args.vdw_potential_type
        params["FORCE_EVAL-DFT-XC-VDW_POTENTIAL-PAIR_POTENTIAL-TYPE"] = args.pair_type
        params["FORCE_EVAL-DFT-XC-VDW_POTENTIAL-PAIR_POTENTIAL-R_CUTOFF"] = args.r_cutoff
        params["FORCE_EVAL-DFT-PRINT-ELF_CUBE-STRIDE"] = args.dft_print_elf_cube_stride
        params["FORCE_EVAL-DFT-PRINT-E_DENSITY_CUBE-STRIDE"] = args.dft_print_e_density_cube_stride
        params["FORCE_EVAL-PROPERTIES-RESP-SLAB_SAMPLING-RANGE"] = args.properties_resp_slab_sampling_range
        params["FORCE_EVAL-PROPERTIES-RESP-SLAB_SAMPLING-SURF_DIRECTION"] = args.properties_resp_slab_sampling_surf_direction
        params["FORCE_EVAL-PROPERTIES-RESP-SLAB_SAMPLING-ATOM_LIST"] = args.properties_resp_slab_sampling_atom_list

        params["MOTION-GEO_OPT-MAX_ITER"] = args.geo_opt_max_iter
        params["MOTION-GEO_OPT-OPTIMIZER"] = args.geo_opt_optimizer
        params["MOTION-GEO_OPT-TYPE"] = args.geo_opt_type
        params["MOTION-GEO_OPT-MAX_DR"] = args.geo_opt_max_dr
        params["MOTION-GEO_OPT-MAX_FORCE"] = args.geo_opt_max_force
        params["MOTION-GEO_OPT-RMS_DR"] = args.geo_opt_rms_dr
        params["MOTION-GEO_OPT-RMS_FORCE"] = args.geo_opt_rms_force

        params["MOTION-CELL_OPT-MAX_ITER"] = args.cell_opt_max_iter
        params["MOTION-CELL_OPT-OPTIMIZER"] = args.cell_opt_optimizer
        params["MOTION-CELL_OPT-TYPE"] = args.cell_opt_type
        params["MOTION-CELL_OPT-MAX_DR"] = args.cell_opt_max_dr
        params["MOTION-CELL_OPT-MAX_FORCE"] = args.cell_opt_max_force
        params["MOTION-CELL_OPT-RMS_DR"] = args.cell_opt_rms_dr
        params["MOTION-CELL_OPT-RMS_FORCE"] = args.cell_opt_rms_force
        params["MOTION-CELL_OPT-PRESSURE_TOLERANCE"] = args.cell_opt_pressure_tolerance
        params["MOTION-CELL_OPT-KEEP_ANGLES"] = args.cell_opt_keep_angles
        params["MOTION-CELL_OPT-KEEP_SYMMETRY"] = args.cell_opt_keep_symmetry

        params["MOTION-BAND-BAND_TYPE"] = args.band_type
        params["MOTION-BAND-NUMBER_OF_REPLICA"] = args.number_of_replica
        params["MOTION-BAND-ALIGN_FRAMES"] = args.align_frames
        params["MOTION-BAND-ROTATE-FRAMES"] = args.rotate_frames
        params["MOTION-BAND-K_SPRING"] = args.k_spring

        params["MOTION-MD-STEPS"] = args.md_steps
        params["MOTION-MD-TIMESTEP"] = args.timestep
        params["MOTION-MD-ENSEMBLE"] = args.ensemble
        params["MOTION-MD-TEMPERATURE"] = args.temperature
        params["MOTION-MD-TEMP_TOL"] = args.temp_tol
        params["MOTION-PRINT-TRAJECTORY-FORMAT"] = args.traj_format

        params["VIBRATIONAL_ANALYSIS-DX"] = args.dx
        params["VIBRATIONAL_ANALYSIS-FULLY_PERIODIC"] = args.fully_periodic
        params["VIBRATIONAL_ANALYSIS-INTENSITIES"] = args.intensities
        params["VIBRATIONAL_ANALYSIS-TC_PRESSURE"] = args.tc_pressure
        params["VIBRATIONAL_ANALYSIS-TC_TEMPERATURE"] = args.tc_temperature
        params["VIBRATIONAL_ANALYSIS-THERMOCHEMISTRY"] = args.thermochemistry

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
            task.set_printout(option=args.printout_option)
            if 2 in args.printout_option and kpath != None:
                task.force_eval.dft.printout.band_structure.set_band(kpath=get_kpath(args.kpath_manual, args.kpath_file))
            task.set_vdw(usevdw=True if args.vdw_potential_type.lower() != "none" else False)
            task.set_run(mpi=args.mpi, server=args.server, jobname=args.jobname, nodes=args.nodes, ppn=args.ppn, queue=args.queue)
            task.set_llhpc(partition=args.partition, nodes=args.nodes, ntask=args.ntask, jobname=args.jobname, stdout=args.stdout, stderr=args.stderr)
            task.scf(directory=args.directory, runopt=args.runopt, auto=args.auto)
        elif args.runtype == 1:
            # geo opt
            from pymatflow.cp2k.opt import opt_run
            task = opt_run()
            task.get_xyz(xyzfile)
            task.set_geo_opt()
            task.set_params(params=params)
            task.set_run(mpi=args.mpi, server=args.server, jobname=args.jobname, nodes=args.nodes, ppn=args.ppn, queue=args.queue)
            task.set_llhpc(partition=args.partition, nodes=args.nodes, ntask=args.ntask, jobname=args.jobname, stdout=args.stdout, stderr=args.stderr)
            task.geo_opt(directory=args.directory, runopt=args.runopt, auto=args.auto)
        elif args.runtype == 2:
            # cell opt
            from pymatflow.cp2k.opt import opt_run
            task = opt_run()
            task.get_xyz(xyzfile)
            task.set_cell_opt()
            task.set_params(params=params)
            task.set_run(mpi=args.mpi, server=args.server, jobname=args.jobname, nodes=args.nodes, ppn=args.ppn, queue=args.queue)
            task.set_llhpc(partition=args.partition, nodes=args.nodes, ntask=args.ntask, jobname=args.jobname, stdout=args.stdout, stderr=args.stderr)
            task.cell_opt(directory=args.directory, runopt=args.runopt, auto=args.auto)
        elif args.runtype == 3:
            # cubic cell opt
            from pymatflow.cp2k.opt import opt_run
            task = opt_run()
            task.get_xyz(xyzfile)
            task.set_geo_opt()
            task.set_params(params=params)
            task.set_run(mpi=args.mpi, server=args.server, jobname=args.jobname, nodes=args.nodes, ppn=args.ppn, queue=args.queue)
            task.set_llhpc(partition=args.partition, nodes=args.nodes, ntask=args.ntask, jobname=args.jobname, stdout=args.stdout, stderr=args.stderr)
            task.cubic(directory=args.directory, runopt=args.runopt, auto=args.auto, na=args.nc, stepa=args.stepa)
        elif args.runtype == 4:
            # hexagonal cell opt
            from pymatflow.cp2k.opt import opt_run
            task = opt_run()
            task.get_xyz(xyzfile)
            task.set_geo_opt()
            task.set_params(params=params)
            task.set_run(mpi=args.mpi, server=args.server, jobname=args.jobname, nodes=args.nodes, ppn=args.ppn, queue=args.queue)
            task.set_llhpc(partition=args.partition, nodes=args.nodes, ntask=args.ntask, jobname=args.jobname, stdout=args.stdout, stderr=args.stderr)
            task.hexagonal(directory=args.directory, runopt=args.runopt, auto=args.auto, na=args.na, nc=args.nc, stepa=args.stepa, stepc=args.stepc)
        elif args.runtype == 5:
            # tetragonal cell opt
            from pymatflow.cp2k.opt import opt_run
            task = opt_run()
            task.get_xyz(xyzfile)
            task.set_geo_opt()
            task.set_params(params=params)
            task.set_run(mpi=args.mpi, server=args.server, jobname=args.jobname, nodes=args.nodes, ppn=args.ppn, queue=args.queue)
            task.set_llhpc(partition=args.partition, nodes=args.nodes, ntask=args.ntask, jobname=args.jobname, stdout=args.stdout, stderr=args.stderr)
            task.tetragonal(directory=args.directory, runopt=args.runopt, auto=args.auto, na=args.na, nc=args.nc, stepa=args.stepa, stepc=args.stepc)
        elif args.runtype == 6:
            # neb
            from pymatflow.cp2k.neb import neb_run
            task = neb_run()
            task.get_images(images=images)
            task.set_params(params=params)
            task.set_run(mpi=args.mpi, server=args.server, jobname=args.jobname, nodes=args.nodes, ppn=args.ppn, queue=args.queue)
            task.set_llhpc(partition=args.partition, nodes=args.nodes, ntask=args.ntask, jobname=args.jobname, stdout=args.stdout, stderr=args.stderr)
            task.neb(directory=args.directory, runopt=args.runopt, auto=args.auto)
        elif args.runtype == 7:
            # phonopy
            from pymatflow.cp2k.phonopy import phonopy_run
            task = phonopy_run()
            task.get_xyz(xyzfile)
            task.supercell_n = args.supercell_n
            task.set_params(params=params)
            task.set_run(mpi=args.mpi, server=args.server, jobname=args.jobname, nodes=args.nodes, ppn=args.ppn, queue=args.queue)
            task.set_llhpc(partition=args.partition, nodes=args.nodes, ntask=args.ntask, jobname=args.jobname, stdout=args.stdout, stderr=args.stderr)
            task.phonopy(directory=args.directory, runopt=args.runopt, auto=args.auto)
        elif args.runtype == 8:
            # vibrational analysis
            from pymatflow.cp2k.vib import vib_run
            task = vib_run()
            task.get_xyz(xyzfile)
            task.set_params(params=params)
            task.set_run(mpi=args.mpi, server=args.server, jobname=args.jobname, nodes=args.nodes, ppn=args.ppn, queue=args.queue)
            task.set_llhpc(partition=args.partition, nodes=args.nodes, ntask=args.ntask, jobname=args.jobname, stdout=args.stdout, stderr=args.stderr)
            task.vib(directory=args.directory, runopt=args.runopt, auto=args.auto)
        elif args.runtype == 9:
            # converge test
            from pymatflow.cp2k.static import static_run
            task = static_run()
            task.get_xyz(xyzfile)
            task.set_params(params=params)
            task.set_run(mpi=args.mpi, server=args.server, jobname=args.jobname, nodes=args.nodes, ppn=args.ppn, queue=args.queue)
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
            task.set_params(params=params)
            task.set_run(mpi=args.mpi, server=args.server, jobname=args.jobname, nodes=args.nodes, ppn=args.ppn, queue=args.queue)
            task.set_llhpc(partition=args.partition, nodes=args.nodes, ntask=args.ntask, jobname=args.jobname, stdout=args.stdout, stderr=args.stderr)
            task.aimd(directory=args.directory, runopt=args.runopt, auto=args.auto)
        else:
            pass
# ==============================================================================
# Quantum ESPERSSO Quantum ESPERSSO Quantum ESPERSSO Quantum ESPERSSO Quantum ESPERSSO
# ==============================================================================
    elif args.driver == "qe":
        control = {}
        electrons = {}
        system = {}
        ions = {}

        control["tstress"] = args.tstress
        control["nstep"] = args.nstep
        control["etot_conv_thr"] = args.etot_conv_thr
        control["forc_conv_thr"] = args.forc_conv_thr
        system["ecutwfc"] = args.ecutwfc
        system["ecutrho"] = args.ecutrho
        system["occupations"] = args.occupations
        system["smearing"] = args.smearing
        system["degauss"] = args.degauss
        system["vdw_corr"] = args.vdw_corr
        system["nbnd"] = args.nbnd
        electrons["conv_thr"] = args.conv_thr

        system["nspin"] = args.nspin
        system["starting_magnetization"] = args.starting_magnetization
        system["noncolin"] = args.noncolin

        path = {}
        path["string_method"] = args.string_method
        path["nstep_path"] = args.nstep_path
        path["opt_scheme"] = args.opt_scheme
        path["num_of_images"] = args.num_of_images
        path["k_max"] = args.k_max
        path["k_min"] = args.k_min
        path["CI_scheme"] = args.ci_scheme
        path["path_thr"] = args.path_thr
        path["ds"] = args.ds
        path["first_last_opt"] = args.first_last_opt

        # for ph.x
        inputph = {}
        inputph["tr2_ph"] = args.tr2_ph
        inputph["lraman"] = args.lraman
        inputph["epsil"] = args.epsil
        inputph["nq1"] = args.nq[0]
        inputph["nq2"] = args.nq[1]
        inputph["nq3"] = args.nq[2]
        inputph["search_sym"] = args.search_sym


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
            task.set_run(mpi=args.mpi, server=args.server, jobname=args.jobname, nodes=args.nodes, ppn=args.ppn, queue=args.queue)
            task.set_llhpc(partition=args.partition, nodes=args.nodes, ntask=args.ntask, jobname=args.jobname, stdout=args.stdout, stderr=args.stderr)
            if args.static == "all":
                task.run(directory=args.directory, runopt=args.runopt, auto=args.auto, kpath=get_kpath(args.kpath_manual, args.kpath_file), kpoints_mp_scf=args.kpoints_mp_scf, kpoints_mp_nscf=args.kpoints_mp_nscf)
            elif args.static == "scf":
                task.scf(directory=args.directory, runopt=args.runopt, auto=args.auto)
        elif args.runtype == 1:
            # relax
            from pymatflow.qe.opt import opt_run
            task = opt_run()
            task.get_xyz(xyzfile)
            task.set_relax()
            task.set_kpoints(kpoints_option=args.kpoints_option, kpoints_mp=args.kpoints_mp)
            task.set_params(control=control, system=system, electrons=electrons, ions=ions)
            task.set_run(mpi=args.mpi, server=args.server, jobname=args.jobname, nodes=args.nodes, ppn=args.ppn, queue=args.queue)
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
            task.set_run(mpi=args.mpi, server=args.server, jobname=args.jobname, nodes=args.nodes, ppn=args.ppn, queue=args.queue)
            task.set_llhpc(partition=args.partition, nodes=args.nodes, ntask=args.ntask, jobname=args.jobname, stdout=args.stdout, stderr=args.stderr)
            task.vc_relax(directory=args.directory, runopt=args.runopt, auto=args.auto)
        elif args.runtype == 3:
            # cubic cell opt
            from pymatflow.qe.opt import opt_run
            task = opt_run()
            task.get_xyz(xyzfile)
            task.set_relax()
            task.set_kpoints(kpoints_option=args.kpoints_option, kpoints_mp=args.kpoints_mp)
            task.set_params(control=control, system=system, electrons=electrons, ions=ions)
            task.set_run(mpi=args.mpi, server=args.server, jobname=args.jobname, nodes=args.nodes, ppn=args.ppn, queue=args.queue)
            task.set_llhpc(partition=args.partition, nodes=args.nodes, ntask=args.ntask, jobname=args.jobname, stdout=args.stdout, stderr=args.stderr)
            task.cubic(directory=args.directory, runopt=args.runopt, auto=args.auto, na=args.na, stepa=args.stepa)
        elif args.runtype == 4:
            # hexagonal cell opt
            from pymatflow.qe.opt import opt_run
            task = opt_run()
            task.get_xyz(xyzfile)
            task.set_relax()
            task.set_kpoints(kpoints_option=args.kpoints_option, kpoints_mp=args.kpoints_mp)
            task.set_params(control=control, system=system, electrons=electrons, ions=ions)
            task.set_run(mpi=args.mpi, server=args.server, jobname=args.jobname, nodes=args.nodes, ppn=args.ppn, queue=args.queue)
            task.set_llhpc(partition=args.partition, nodes=args.nodes, ntask=args.ntask, jobname=args.jobname, stdout=args.stdout, stderr=args.stderr)
            task.hexagonal(directory=args.directory, runopt=args.runopt, auto=args.auto, na=args.na, nc=args.nc, stepa=args.stepa, stepc=args.stepc)
        elif args.runtype == 5:
            # tetragonal cell opt
            from pymatflow.qe.opt import opt_run
            task = opt_run()
            task.get_xyz(xyzfile)
            task.set_relax()
            task.set_kpoints(kpoints_option=args.kpoints_option, kpoints_mp=args.kpoints_mp)
            task.set_params(control=control, system=system, electrons=electrons, ions=ions)
            task.set_run(mpi=args.mpi, server=args.server, jobname=args.jobname, nodes=args.nodes, ppn=args.ppn, queue=args.queue)
            task.set_llhpc(partition=args.partition, nodes=args.nodes, ntask=args.ntask, jobname=args.jobname, stdout=args.stdout, stderr=args.stderr)
            task.tetragonal(directory=args.directory, runopt=args.runopt, auto=args.auto, na=args.na, nc=args.nc, stepa=args.stepa, stepc=args.stepc)
        elif args.runtype == 6:
            from pymatflow.qe.neb import neb_run
            task = neb_run()
            task.get_images(images=images)
            task.set_kpoints(kpoints_option=args.kpoints_option, kpoints_mp=args.kpoints_mp)
            task.set_path(path=path)
            task.set_run(mpi=args.mpi, server=args.server, jobname=args.jobname, nodes=args.nodes, ppn=args.ppn, queue=args.queue)
            task.set_llhpc(partition=args.partition, nodes=args.nodes, ntask=args.ntask, jobname=args.jobname, stdout=args.stdout, stderr=args.stderr)
            task.neb(directory=args.directory, runopt=args.runopt, auto=args.auto)
        elif args.runtype == 7:
            from pymatflow.qe.dfpt import dfpt_run
            task = dfpt_run()
            task.get_xyz(xyzfile)
            task.set_inputph(inputph=inputph)
            task.set_run(mpi=args.mpi, server=args.server, jobname=args.jobname, nodes=args.nodes, ppn=args.ppn, queue=args.queue)
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
            task.set_run(mpi=args.mpi, server=args.server, jobname=args.jobname, nodes=args.nodes, ppn=args.ppn, queue=args.queue)
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
            task.set_run(mpi=args.mpi, server=args.server, jobname=args.jobname, nodes=args.nodes, ppn=args.ppn, queue=args.queue)
            task.set_llhpc(partition=args.partition, nodes=args.nodes, ntask=args.ntask, jobname=args.jobname, stdout=args.stdout, stderr=args.stderr)
            task.pp(directory=args.directory, runopt=args.runopt, auto=args.auto)

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
            task.set_run(mpi=args.mpi, server=args.server, jobname=args.jobname, nodes=args.nodes, ppn=args.ppn, queue=args.queue)
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
            task.set_run(mpi=args.mpi, server=args.server, jobname=args.jobname, nodes=args.nodes, ppn=args.ppn, queue=args.queue)
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
            task.set_run(mpi=args.mpi, server=args.server, jobname=args.jobname, nodes=args.nodes, ppn=args.ppn, queue=args.queue)
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
            task.set_run(mpi=args.mpi, server=args.server, jobname=args.jobname, nodes=args.nodes, ppn=args.ppn, queue=args.queue)
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
            task.set_run(mpi=args.mpi, server=args.server, jobname=args.jobname, nodes=args.nodes, ppn=args.ppn, queue=args.queue)
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
            task.set_run(mpi=args.mpi, server=args.server, jobname=args.jobname, nodes=args.nodes, ppn=args.ppn, queue=args.queue)
            task.set_llhpc(partition=args.partition, nodes=args.nodes, ntask=args.ntask, jobname=args.jobname, stdout=args.stdout, stderr=args.stderr)
            task.phonopy(directory=args.directory, runopt=args.runopt, auto=args.auto)
        elif args.runtype == 6:
            # molecular dynamics
            from pymatflow.siesta.md import md_run
            task = md_run()
            task.get_xyz(xyzfile)
            task.set_params(params=params)
            task.set_kpoints(kpoints_mp=args.kpoints_mp)
            task.set_run(mpi=args.mpi, server=args.server, jobname=args.jobname, nodes=args.nodes, ppn=args.ppn, queue=args.queue)
            task.set_llhpc(partition=args.partition, nodes=args.nodes, ntask=args.ntask, jobname=args.jobname, stdout=args.stdout, stderr=args.stderr)
            task.md(directory=args.directory, runopt=args.runopt, auto=args.auto)
        else:
            pass
    elif args.driver == "vasp":
        # deal with INCAR template specified by --incar
        if args.incar == None:
            pass
        else:
            params_all = []
            for i in range(len(args.incar)):
                if not os.path.exists(args.incar[i]):
                    print("====================================================\n")
                    print("                  Warning !!!!\n")
                    print("----------------------------------------------------\n")
                    print("matflow vasp:\n")
                    print("the specified incar file by --incar doesn't exist\n")
                    print("go and check it\n")
                    sys.exit(1)
                with open(args.incar[i], 'r') as fin:
                    incar = fin.readlines()
                params = {}                        
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
                params_all.append(params)
                #
                # if xxx is alraedy in params(set from --incar) and args.xxx is None
                # params[xxx] will not be control by args.xxx
        # if number of specified incar is less than args.gen_incar_n, make an supplement of params by copying the last params
        while len(args.incar) < args.gen_incar_n:
            params_all.append(params[-1])
        
        # if number of specified value is less than args.gen_incar_n, make an summplement of params by copying the last value
        args.nwrite = [args.nwrite[i] if i < len(args.nwrite) else args.nwrite[-1] for i in range(args.gen_incar_n)] if args.nwrite != None else [None for i in range(args.gen_incar_n)]
        args.prec = [args.prec[i] if i < len(args.prec) else args.prec[-1] for i in range(args.gen_incar_n)] if args.prec != None else [None for i in range(args.gen_incar_n)]
        args.ncore = [args.ncore[i] if i < len(args.ncore) else args.ncore[-1] for i in range(args.gen_incar_n)] if args.ncore != None else [None for i in range(args.gen_incar_n)]
        args.encut = [args.encut[i] if i < len(args.encut) else args.encut[-1] for i in range(args.gen_incar_n)] if args.encut != None else [None for i in range(args.gen_incar_n)]
        args.ediff = [args.ediff[i] if i < len(args.ediff) else args.ediff[-1] for i in range(args.gen_incar_n)] if args.ediff != None else [None for i in range(args.gen_incar_n)]
        args.nelm = [args.nelm[i] if i < len(args.nelm) else args.nelm[-1] for i in range(args.gen_incar_n)] if args. nelm != None else [None for i in range(args.gen_incar_n)]
        args.nfree = [args.nfree[i] if i < len(args.nfree) else args.nfree[-1] for i in range(args.gen_incar_n)] if args.nfree != None else [None for i in range(args.gen_incar_n)]
        args.kspacing = [args.kspacing[i] if i < len(args.kspacing) else args.kspacing[-1] for i in range(args.gen_incar_n)] if args.kspacing != None else [None for i in range(args.gen_incar_n)]
        args.ismear = [args.ismear[i] if i < len(args.ismear) else args.ismear[-1] for i in range(args.gen_incar_n)] if args.ismear != None else [None for i in range(args.gen_incar_n)]
        args.sigma = [args.sigma[i] if i < len(args.sigma) else args.sigma[-1] for i in range(args.gen_incar_n)] if args.sigma != None else [None for i in range(args.gen_incar_n)]
        args.ivdw = [args.ivdw[i] if i < len(args.ivdw) else args.ivdw[-1] for i in range(args.gen_incar_n)] if args.ivdw != None else [None for i in range(args.gen_incar_n)]
        args.ediffg = [args.ediffg[i] if i < len(args.ediffg) else args.ediffg[-1] for i in range(args.gen_incar_n)] if args.ediffg != None else [None for i in range(args.gen_incar_n)]
        args.nsw = [args.nsw[i] if i < len(args.nsw) else args.nsw[-1] for i in range(args.gen_incar_n)] if args.nsw != None else [None for i in range(args.gen_incar_n)]
        args.ibrion = [args.ibrion[i] if i < len(args.ibrion) else args.ibrion[-1] for i in range(args.gen_incar_n)] if args.ibrion != None else [None for i in range(args.gen_incar_n)]
        args.isif = [args.isif[i] if i < len(args.isif) else args.isif[-1] for i in range(args.gen_incar_n)] if args.isif != None else [None for i in range(args.gen_incar_n)]
        args.potim = [args.potim[i] if i < len(args.potim) else args.potim[-1] for i in range(args.gen_incar_n)] if args.potim != None else [None for i in range(args.gen_incar_n)]
        args.lorbit = [args.lorbit[i] if i < len(args.lorbit) else args.lorbit[-1] for i in range(args.gen_incar_n)] if args.lorbit != None else [None for i in range(args.gen_incar_n)]
        args.loptics = [args.loptics[i] if i < len(args.loptics) else args.loptics[-1] for i in range(args.gen_incar_n)] if args.loptics != None else [None for i in range(args.gen_incar_n)]
        args.cshift = [args.cshift[i] if i < len(args.cshift) else args.cshift[-1] for i in range(args.gen_incar_n)] if args.cshift != None else [None for i in range(args.gen_incar_n)]
        args.nedos = [args.nedos[i] if i < len(args.nedos) else args.nedos[-1] for i in range(args.gen_incar_n)] if args.nedos != None else [None for i in range(args.gen_incar_n)]
        args.lsubrot = [args.lsubrot[i] if i < len(args.lsubrot) else args.lsubrot[-1] for i in range(args.gen_incar_n)] if args.lsubrot != None else [None for i in range(args.gen_incar_n)]
        args.saxis = [args.saxis[i] if i < len(args.saxis) else args.saxis[-1] for i in range(args.gen_incar_n)] if args.saxis != None else [None for i in range(args.gen_incar_n)]
        args.lmaxmix = [args.lmaxmix[i] if i < len(args.lmaxmix) else args.lmaxmix[-1] for i in range(args.gen_incar_n)] if args.lmaxmix != None else [None for i in range(args.gen_incar_n)]
        args.algo = [args.algo[i] if i < len(args.algo) else args.algo[-1] for i in range(args.gen_incar_n)] if args.algo != None else [None for i in range(args.gen_incar_n)]
        args.ialgo = [args.ialgo[i] if i < len(args.ialgo) else args.ialgo[-1] for i in range(args.gen_incar_n)] if args.ialgo != None else [None for i in range(args.gen_incar_n)]
        args.addgrid = [args.addgrid[i] if i < len(args.addgrid) else args.addgrid[-1] for i in range(args.gen_incar_n)] if args.addgrid != None else [None for i in range(args.gen_incar_n)]
        args.isym = [args.isym[i] if i < len(args.isym) else args.isym[-1] for i in range(args.gen_incar_n)] if args.isym != None else [None for i in range(args.gen_incar_n)]
        args.lreal = [args.lreal[i] if i < len(args.lreal) else args.lreal[-1] for i in range(args.gen_incar_n)] if args.lreal != None else [None for i in range(args.gen_incar_n)]
        args.pstress = [args.pstress[i] if i < len(args.pstress) else args.pstress[-1] for i in range(args.gen_incar_n)] if args.pstress != None else [None for i in range(args.gen_incar_n)]
        args.lwave = [args.lwave[i] if i < len(args.lwave) else args.lwave[-1] for i in range(args.gen_incar_n)] if args.lwave != None else [None for i in range(args.gen_incar_n)]
        args.lcharg = [args.lcharg[i] if i < len(args.lcharg) else args.lcharg[-1] for i in range(args.gen_incar_n)] if args.lcharg != None else [None for i in range(args.gen_incar_n)]
        args.ispin = [args.ispin[i] if i < len(args.ispin) else args.ispin[-1] for i in range(args.gen_incar_n)] if args.ispin != None else [None for i in range(args.gen_incar_n)]
        args.magmom = [args.magmom[i] if i < len(args.magmom) else args.magmom[-1] for i in range(args.gen_incar_n)] if args.magmom != None else [None for i in range(args.gen_incar_n)]
        args.lnoncollinear = [args.lnoncollinear[i] if i < len(args.lnoncollinear) else args.nwrite[-1] for i in range(args.gen_incar_n)] if args.lnoncollinear != None else [None for i in range(args.gen_incar_n)]
        args.lsorbit = [args.lsorbit[i] if i < len(args.lsorbit) else args.lsorbit[-1] for i in range(args.gen_incar_n)] if args.lsorbit != None else [None for i in range(args.gen_incar_n)]
        args.algo = [args.algo[i] if i < len(args.algo) else args.algo[-1] for i in range(args.gen_incar_n)] if args.algo != None else [None for i in range(args.gen_incar_n)]
        args.lhfcalc = [args.lhfcalc[i] if i < len(args.lhfcalc) else args.lhfcalc[-1] for i in range(args.gen_incar_n)] if args.lhfcalc != None else [None for i in range(args.gen_incar_n)]
        args.hfscreen = [args.hfscreen[i] if i < len(args.hfscreen) else args.hfscreen[-1] for i in range(args.gen_incar_n)] if args.hfscreen != None else [None for i in range(args.gen_incar_n)]
        args.aexx = [args.aexx[i] if i < len(args.aexx) else args.aexx[-1] for i in range(args.gen_incar_n)] if args.aexx != None else [None for i in range(args.gen_incar_n)]
        args.lelf = [args.lelf[i] if i < len(args.lelf) else args.lelf[-1] for i in range(args.gen_incar_n)] if args.lelf != None else [None for i in range(args.gen_incar_n)]
        args.amix = [args.amix[i] if i < len(args.amix) else args.amix[-1] for i in range(args.gen_incar_n)] if args.amix != None else [None for i in range(args.gen_incar_n)]
        args.bmix = [args.bmix[i] if i < len(args.bmix) else args.bmix[-1] for i in range(args.gen_incar_n)] if args.bmix != None else [None for i in range(args.gen_incar_n)]
        
        for i in range(args.gen_incar_n):
            params_all[i]["NWRITE"] = args.nwrite[i] if "NWRITE" not in params_all[i] or args.nwrite[i] != None else params_all[i]["NWRITE"]
            params_all[i]["PREC"] = args.prec[i] if "PREC" not in params_all[i] or args.prec[i] != None else params_all[i]["PREC"]
            params_all[i]["NCORE"] = args.ncore[i] if "NCORE" not in params_all[i] or args.ncore[i] != None else params_all[i]["NCORE"]
            params_all[i]["ENCUT"] = args.encut[i] if "ENCUT" not in params_all[i] or args.encut[i] != None else params_all[i]["ENCUT"]
            params_all[i]["EDIFF"] = args.ediff[i] if "EDIFF" not in params_all[i] or args.ediff[i] != None else params_all[i]["EDIFF"]
            params_all[i]["NELM"] = args.nelm[i] if "NELM" not in params_all[i] or args.nelm[i] != None else params_all[i]["NELM"]
            params_all[i]["NFREE"] = args.nfree[i] if "NFREE" not in params_all[i] or args.nfree[i] != None else params_all[i]["NFREE"]
            params_all[i]["KSPACING"] = args.kspacing[i] if "KSPACING" not in params_all[i] or args.kspacing[i] != None else params_all[i]["KSPACING"]
            params_all[i]["ISMEAR"] = args.ismear[i] if "ISMEAR" not in params_all[i] or args.ismear[i] != None else params_all[i]["ISMEAR"]
            params_all[i]["SIGMA"] = args.sigma[i] if "SIGMA" not in params_all[i] or args.sigma[i] != None else params_all[i]["SIGMA"]
            params_all[i]["IVDW"] = args.ivdw[i] if "IVDW" not in params_all[i] or args.ivdw[i] != None else params_all[i]["IVDW"]
            params_all[i]["EDIFFG"] = args.ediffg[i] if "EDIFFG" not in params_all[i] or args.ediffg[i] != None else params_all[i]["EDIFFG"]
            params_all[i]["NSW"] = args.nsw[i] if "NSW" not in params_all[i] or args.nsw[i] != None else params_all[i]["NSW"]
            params_all[i]["IBRION"] = args.ibrion[i] if "IBRION" not in params_all[i] or args.ibrion[i] != None else params_all[i]["IBRION"]
            params_all[i]["ISIF"] = args.isif[i] if "ISIF" not in params_all[i] or args.isif[i] != None else params_all[i]["ISIF"]
            params_all[i]["POTIM"] = args.potim[i] if "POTIM" not in params_all[i] or args.potim[i] != None else params_all[i]["POTIM"]
            params_all[i]["LORBIT"] = args.lorbit[i] if "LORBIT" not in params_all[i] or args.potim[i] != None else params_all[i]["LORBIT"]
            params_all[i]["LOPTICS"] = args.loptics[i] if "LOPTICS" not in params_all[i] or args.loptics[i] != None else params_all[i]["LOPTICS"]
            params_all[i]["CSHIFT"] = args.cshift[i] if "CSHIFT" not in params_all[i] or args.cshift[i] != None else params_all[i]["CSHIFT"]
            params_all[i]["NEDOS"] = args.nedos[i] if "NEDOS"  not in params_all[i] or args.nedos[i] != None else params_all[i]["NEDOS"]
            params_all[i]["LSUBROT"] = args.lsubrot[i] if "LSUBROT" not in params_all[i] or args.lsubrot[i] != None else params_all[i]["LSUBROT"]
            params_all[i]["SAXIS"] = args.saxis[i] if "SAXIS" not in params_all[i] or args.saxis[i] != None else params_all[i]["SAXIS"]
            params_all[i]["LMAXMIX"] = args.lmaxmix[i] if "LMAXMIX" not in params_all[i] or args.lmaxmix[i] != None else params_all[i]["LMAXMIX"]
            params_all[i]["ALGO"] = args.algo[i] if "ALGO" not in params_all[i] or args.algo[i] != None else params_all[i]["ALGO"]
            params_all[i]["IALGO"] = args.ialgo[i] if "IALGO" not in params_all[i] or args.ialgo[i] != None else params_all[i]["IALGO"]
            params_all[i]["ADDGRID"] = args.addgrid[i] if "ADDGRID" not in params_all[i] or args.addgrid[i] != None else params_all[i]["ADDGRID"]
            params_all[i]["ISYM"] = args.isym[i] if "ISYM" not in params_all[i] or args.isym[i] != None else params_all[i]["ISYM"]
            params_all[i]["PSTRESS"] = args.pstress[i] if "PSTRESS" not in params_all[i] or args.pstress[i] != None else params_all[i]["PSTRESS"]
            params_all[i]["LREAL"] = args.lreal[i] if "LREAL" not in params_all[i] or args.lreal[i] != None else params_all[i]["LREAL"]
            params_all[i]["LWAVE"] = args.lwave[i] if "LWAVE" not in params_all[i] or args.lwave[i] != None else params_all[i]["LWAVE"]
            params_all[i]["LCHARG"] = args.lcharg[i] if "LCHARG" not in params_all[i] or args.lcharg[i] != None else params_all[i]["LCHARG"]
            params_all[i]["ISPIN"] = args.ispin[i] if "ISPIN" not in params_all[i] or args.ispin[i] != None else params_all[i]["ISPIN"]
            params_all[i]["MAGMOM"] = args.magmom[i] if "MAGMOM" not in params_all[i] or args.magmom[i] != None else params_all[i]["MAGMOM"] # magmom can be a list that can be automatically dealt with by base.incar.to_incar()
            params_all[i]["LNONCOLLINEAR"] = args.lnoncollinear[i] if "LNONCOLLINEAR" not in params_all[i] or args.lnoncollinear[i] != None else params_all[i]["LNONCOLLINEAR"]
            params_all[i]["LSORBIT"] = args.lsorbit[i] if "LSORBIT" not in params_all[i] or args.lsorbit[i] != None else params_all[i]["LSORBIT"]
            params_all[i]["ALGO"] = args.algo[i] if "ALGO" not in params_all[i] or args.algo[i] != None else params_all[i]["ALGO"]
            params_all[i]["LHFCALC"] = args.lhfcalc[i] if "LHFCALC" not in params_all[i] or args.lhfcalc[i] != None else params_all[i]["LHFCALC"]
            params_all[i]["HFSCREEN"] = args.hfscreen[i] if "HFSCREEN" not in params_all[i] or args.hfscreen[i] != None else params_all[i]["HFSCREEN"]
            params_all[i]["AEXX"] = args.aexx[i] if "AEXX" not in params_all[i] or args.aexx[i] != None else params_all[i]["AEXX"]
            params_all[i]["LELF"] = args.lelf[i] if "LELF" not in params_all[i] or args.lelf[i] != None else params_all[i]["LELF"]
            params_all[i]["AMIX"] = args.amix[i] if "AMIX" not in params_all[i] or args.amix[i] != None else params_all[i]["AMIX"]
            params_all[i]["BMIX"] = args.bmix[i] if "BMIX" not in params_all[i] or args.bmix[i] != None else params_all[i]["BMIX"]
        
        if args.runtype == 0:
            pass
        elif args.runtype == 1:
            # optimization
            from pymatflow.flow.calypso.calypso import calypso
            #
            if args.fix != None:
                fix_str = ""
                for i in args.fix:
                        fix_str += "%d " % i
                os.system("xyz-fix-atoms.py -i %s -o %s --fix %s" % (xyzfile, xyzfile, fix_str))
                args.selective_dynamics = "T"
            #            
            pso = calypso()
            pso.set_params(params=input_dat_params)
            #pso.vasp.get_xyz(xyzfile)
            #pso.vasp.set_params(params=params, runtype="opt")
            pso.gen_incar_n = args.gen_incar_n
            pso.multi_incar_params = params_all
            pso.vasp.set_kpoints(kpoints_mp=args.kpoints_mp)
            #pso.vasp.poscar.selective_dynamics = True if args.selective_dynamics.upper()[0] == "T" else False
            pso.set_run(mpi=args.mpi, server=args.server, jobname=args.jobname, nodes=args.nodes, ppn=args.ppn, queue=args.queue)
            pso.set_llhpc(partition=args.partition, nodes=args.nodes, ntask=args.ntask, jobname=args.jobname, stdout=args.stdout, stderr=args.stderr)
            if args.split_batch == None:
                pso.run_vasp(directory=args.directory, runopt=args.runopt, auto=args.auto)
            else:
                pso.run_vasp_split(directory=args.directory, runopt=args.runopt, auto=args.auto, split_batch=args.split_batch)
        elif args.runtype == 2:
            pass
        elif args.runtype == 3:
            pass
        elif args.runtype == 4:
            pass
        elif args.runtype == 5:
            pass        
        elif args.runtype == 6:
            pass
        elif args.runtype == 7:
            pass
        elif args.runtype == 8:
            pass
    # --------------------------------------------------------------------------



if __name__ == "__main__":
    main()