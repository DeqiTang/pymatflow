import os

def qeSubparser(subparsers):
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
        choices=["pbs", "llhpc", "tianhe2", "cdcloud"],
        help="type of remote server, can be pbs or llhpc, cdcloud")

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

    gp.add_argument("--lspinorb", type=str, default=None,
        choices=[".true.", ".false."],
        help="if .TRUE. the noncollinear code can use a pseudopotential with spin-orbit.")

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

    gp.add_argument("--diagonalization", type=str, default=None,
        choices=["david", "cg", "ppcg", "paro"],
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


def qeDriver(args):
    from pymatflow.cmd.matflow import getXyzFile
    xyzfile, images = getXyzFile(args)
    # server
    # xxx.set_run can only deal with pbs, llhpc, lsf_sz server now 
    # however both guangzhou chaosuan llhpc are build on tianhe2, so they can use the same job system(yhbatch...)
    # we add tianhe2 option to args.server which cannot be handled by xxx.set_run. so we convert it to llhpc if tianhe2 is chosen
    server = args.server if args.server != "tianhe2" else "llhpc"
          
    # ====================================================================================
    # Quantum ESPERSSO Quantum ESPERSSO Quantum ESPERSSO Quantum ESPERSSO Quantum ESPERSSO
    # ====================================================================================
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
    system["lspinorb"] = args.lspinorb if "lspinorb" not in system or args.lspinorb != None else system["lspinorb"]
    
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
        from pymatflow.qe.static import StaticRun
        from pymatflow.cmd.matflow import get_kpath
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
        task = StaticRun()
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
        from pymatflow.qe.opt import OptRun
        
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

        task = OptRun()
        task.get_xyz(xyzfile)
        task.set_relax()
        task.set_kpoints(kpoints_option=args.kpoints_option, kpoints_mp=args.kpoints_mp)
        task.set_params(control=control, system=system, electrons=electrons, ions=ions)
        task.set_run(mpi=args.mpi, server=server, jobname=args.jobname, nodes=args.nodes, ppn=args.ppn, queue=args.queue)
        task.set_llhpc(partition=args.partition, nodes=args.nodes, ntask=args.ntask, jobname=args.jobname, stdout=args.stdout, stderr=args.stderr)
        task.relax(directory=args.directory, runopt=args.runopt, auto=args.auto)
    elif args.runtype == 2:
        # vc-relax
        from pymatflow.qe.opt import OptRun
        task = OptRun()
        task.get_xyz(xyzfile)
        task.set_vc_relax()
        task.set_kpoints(kpoints_option=args.kpoints_option, kpoints_mp=args.kpoints_mp)
        task.set_params(control=control, system=system, electrons=electrons, ions=ions)
        task.set_run(mpi=args.mpi, server=server, jobname=args.jobname, nodes=args.nodes, ppn=args.ppn, queue=args.queue)
        task.set_llhpc(partition=args.partition, nodes=args.nodes, ntask=args.ntask, jobname=args.jobname, stdout=args.stdout, stderr=args.stderr)
        task.vc_relax(directory=args.directory, runopt=args.runopt, auto=args.auto)
    elif args.runtype == 3:
        # cubic cell opt
        from pymatflow.qe.opt import OptRun
        task = OptRun()
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
        from pymatflow.qe.opt import OptRun
        task = OptRun()
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
        from pymatflow.qe.opt import OptRun
        task = OptRun()
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
        from pymatflow.qe.neb import NebRun
        task = NebRun()
        task.get_images(images=images)
        task.set_kpoints(kpoints_option=args.kpoints_option, kpoints_mp=args.kpoints_mp)
        task.set_path(path=path)
        task.set_run(mpi=args.mpi, server=server, jobname=args.jobname, nodes=args.nodes, ppn=args.ppn, queue=args.queue)
        task.set_llhpc(partition=args.partition, nodes=args.nodes, ntask=args.ntask, jobname=args.jobname, stdout=args.stdout, stderr=args.stderr)
        task.neb(directory=args.directory, runopt=args.runopt, auto=args.auto)
    elif args.runtype == 7:
        from pymatflow.qe.dfpt import DfptRun
        task = DfptRun()
        task.get_xyz(xyzfile)
        task.set_inputph(inputph=inputph)
        task.set_run(mpi=args.mpi, server=server, jobname=args.jobname, nodes=args.nodes, ppn=args.ppn, queue=args.queue)
        task.set_llhpc(partition=args.partition, nodes=args.nodes, ntask=args.ntask, jobname=args.jobname, stdout=args.stdout, stderr=args.stderr)
        task.phx(directory=args.directory, runopt=args.runopt, auto=args.auto)
    elif args.runtype == 8:
        # phonopy
        from pymatflow.qe.phonopy import PhonopyRun
        task = PhonopyRun()
        task.get_xyz(xyzfile)
        task.set_kpoints(kpoints_option=args.kpoints_option, kpoints_mp=args.kpoints_mp)
        task.set_params(control=control, system=system, electrons=electrons)
        task.supercell_n = args.supercell_n
        task.set_run(mpi=args.mpi, server=server, jobname=args.jobname, nodes=args.nodes, ppn=args.ppn, queue=args.queue)
        task.set_llhpc(partition=args.partition, nodes=args.nodes, ntask=args.ntask, jobname=args.jobname, stdout=args.stdout, stderr=args.stderr)
        task.phonopy(directory=args.directory, runopt=args.runopt, auto=args.auto)
    elif args.runtype == 9:
        # pp.x
        from pymatflow.qe.static import StaticRun

        inputpp = {}
        plotpp = {}
        
        inputpp["plot_num"] = args.plot_num
        plotpp["iflag"] = args.iflag
        plotpp["output_format"] = args.output_format

        task = StaticRun()
        task.get_xyz(xyzfile)
        task.set_pp(inputpp=inputpp, plotpp=plotpp)
        task.set_run(mpi=args.mpi, server=server, jobname=args.jobname, nodes=args.nodes, ppn=args.ppn, queue=args.queue)
        task.set_llhpc(partition=args.partition, nodes=args.nodes, ntask=args.ntask, jobname=args.jobname, stdout=args.stdout, stderr=args.stderr)
        task.pp(directory=args.directory, runopt=args.runopt, auto=args.auto)
    elif args.runtype == 10:
        # abc cell opt
        from pymatflow.qe.opt import OptRun
        task = OptRun()
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
        from pymatflow.qe.static import StaticRun
        
        task = StaticRun()
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