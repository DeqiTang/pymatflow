

def cp2kSubparser(subparsers):
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

    gp.add_argument("--flow-cp2k-version", type=str, default="stable",
        choices=["stable", "dev"],
        help="choose pymatflow.cp2k version, eigher stable or dev.")

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

    gp.add_argument("--print-level", type=str, default="LOW", # None,
        choices=["DEBUG", "HIGH", "LOW", "MEDIUM", "SILENT", "debug", "high", "low", "medium", "silent"],
        help="How much output is written out.")

    gp.add_argument("--extended-fft-lengths", type=str, default=None,
        choices=["TRUE", "FALSE", "true", "false"],
        help="Use fft library specific values for the allows number of points in FFTs. The default is to use the internal FFT lengths. For external fft libraries this may create an error at the external library level, because the length provided by cp2k is not supported by the external library. In this case switch on this keyword to obtain, with certain fft libraries, lengths matching the external fft library lengths, or larger allowed grids, or grids that more precisely match a given cutoff. IMPORTANT NOTE: in this case, the actual grids used in CP2K depends on the FFT library. A change of FFT library must therefore be considered equivalent to a change of basis, which implies a change of total energy.")

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
           
    gp.add_argument("--qs-extrapolation", type=str, default=None,
        choices=["ASPC", "FROZEN", "LINEAR_P", "LINEAR_PS", "PS", "USE_GUESS", "USE_PREV_P", "USE_PREV_RHO_R", "USE_PREV_WF", "aspc", "frozen", "linear_p", "linear_ps", "ps", "use_guess", "use_prev_p", "ues_prev_rho_r", "use_prev_wf"],
        help="Extrapolation strategy for the wavefunction during e.g. MD. Not all options are available for all simulation methods. PS and ASPC are recommended, see also EXTRAPOLATION_ORDER. default is ASPC")

    gp.add_argument("--qs-extrapolation-order", type=int, default=None,
        help="Order for the PS or ASPC extrapolation (typically 2-4). Higher order might bring more accuracy, but comes, for large systems, also at some cost. In some cases, a high order extrapolation is not stable, and the order needs to be reduced. default is 3")
    
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

    gp.add_argument("--xc-grid-use-finer-grid", type=str, default=None,
        choices=["TRUE", "FALSE", "true", "false"],
        help="Uses a finer grid only to calculate the xc")

    gp.add_argument("--xc-grid-xc-deriv", type=str, default=None,
        choices=["COLLOCATE", "NN10_SMOOTH", "NN4_SMOOTH", "NN50_SMOOTH", "NN6_SMOOTH", "PW", "SPLINE2", "SPLINE2_SMOOTH", "SPLINE3", "SPLINE3_SMOOTH", "collocate", "nn10_smooth", "nn4_smooth", "nn50_smooth", "nn6_smooth", "pw", "spline2", "spline2_smooth", "spline3", "spline3_smooth"],
        help="The method used to compute the derivatives, default is PW")

    gp.add_argument("--xc-grid-xc-smooth-rho", type=str, default=None,
        choices=["NN10", "NN4", "NN50", "NN6", "NONE", "SPLINE2", "SPLINE3", "nn10", "nn4", "nn50", "nn6", "none", "spline2", "spline3"],
        help="The density smoothing used for the xc calculation")

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

    gp.add_argument("--mixing-broy-w0", type=float, default=None,
        help="w0 parameter used in Broyden mixing, default is 0.01")

    gp.add_argument("--mixing-broy-wmax", type=float, default=None,
        help="default is 30")

    gp.add_argument("--mixing-broy-wref", type=float, default=None,
        help="default is 100")

    gp.add_argument("--mixing-nbuffer", type=int, default=None, # default is 4
        help="Number of previous steps stored for the actual mixing scheme. default is 4")

    gp.add_argument("--mixing-gmix-p", type=str, default=None,
        choices=["TRUE", "FALSE", "true", "false"],
        help="Activate the mixing of the denisty matrix, using the same mixing coefficient applied for the g-space mixing. default is false")

    gp.add_argument("--mixing-max-gvec-exp", type=int, default=None,
        help="Restricts the G-space mixing to lower part of G-vector spectrum, up to a G0, by assigning the exponent of the Gaussian that can be represented by vectors smaller than G0 within a certain accuracy. default is -1")            

    gp.add_argument("--mixing-nmixing", type=int, default=None,
        help="Minimal number of density mixing (should be greater than 0),before starting DIIS, default is 2")

    gp.add_argument("--mixing-nskip", type=int, default=None,
        help="Number of initial iteration for which the mixing is skipped, default is 0")

    gp.add_argument("--mixing-n-simple-mix", type=int, default=None,
        help="Number of kerker damping iterations before starting other mixing procedures, default is 0")

    gp.add_argument("--mixing-pulay-alpha", type=float, default=None,
        help="Fraction of new density to be added to the Pulay expansion. default is 0")

    gp.add_argument("--mixing-pulay-beta", type=float, default=None,
        help="Fraction of residual contribution to be added to Pulay expansion. default is 1")

    gp.add_argument("--mixing-regularization", type=float, default=None,
        help="Regularization parameter to stabilize the inversion of the residual matrix {Yn^t Yn} in the multisecant mixing scheme (noise)")

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

    gp.add_argument("--scf-guess", type=str, default=None,
        choices=["ATOMIC", "CORE", "HISTORY_RESTART", "MOPAC", "NONE", "RANDOM", "RESTART", "SPARSE", "atomic", "core", "history_restart", "mopac", "none", "random", "restart", "sparse"],
        help="Change the initial guess for the wavefunction.")

    # smear
    gp.add_argument("--smear", type=str, default=None, #"FALSE",
        choices=["TRUE", "FALSE", "true", "false"],
        help="controls the activation of smearing")

    gp.add_argument("--smear-method", type=str, default="FERMI_DIRAC",
        help="Smearing method to be applied, can be fermi_dirac(default) or energy_window")

    gp.add_argument("--added-mos", type=int, default=None, #0,
        help="Number of additional MOS added for each spin")

    gp.add_argument("--dft-scf-cholesky", type=str, default=None,
        choices=["INVERSE", "INVERSE_DBCSR", "OFF", "REDUCE", "RESTORE", "inverse", "inverse_dbcsr", "off", "reduce", "restore"],
        help="f the cholesky method should be used for computing the inverse of S, and in this case calling which Lapack routines. default is RESTORE")

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

    # DFT/SCCS
    gp = subparser.add_argument_group(title="DFT/SCCS")

    gp.add_argument("--dft-sccs", type=str, default=None,
        choices=["true", "false", "TRUE", "FALSE"],
        help="Controls the activation of the SCCS section")
    
    gp.add_argument("--dft-sccs-alpha", type=float, default=None,
        help="Solvent specific tunable parameter for the calculation of the repulsion term.")

    gp.add_argument("--dft-sccs-beta", type=float, default=None,
        help="Solvent specific tunable parameter for the calculation of the dispersion term")

    gp.add_argument("--dft-sccs-delta-rho", type=float, default=None,
        help="Numerical increment for the calculation of the (quantum) surface of the solute cavity")

    gp.add_argument("--dft-sccs-dielectric-constant", type=float, default=None,
        help="Dielectric constant of the solvent, defaut is 80")

    gp.add_argument("--dft-sccs-eps-sccs", type=float, default=None,
        help="Tolerance for the convergence of the polarisation density, i.e. requested accuracy for the SCCS iteration cycle, default is 1.0E-6")

    gp.add_argument("--dft-sccs-eps-scf", type=float, default=None,
        help="The SCCS iteration cycle is activated only if the SCF iteration cycle is converged to this threshold value, default is 1.0E-1")

    gp.add_argument("--dft-sccs-gamma", type=float, default=None,
        help="Surface tension of the solvent used for the calculation of the cavitation term, default is 0")

    gp.add_argument("--dft-sccs-max-iter", type=int, default=None,
        help="Maximum number of SCCS iteration steps performed to converge within the given tolerance, default is 100")            

    gp.add_argument("--dft-sccs-method", type=str, default=None,
        choices=["ANDREUSSI", "FATTEBERT-GYGI", "andreussi", "fatterbert-gygi"],
        help="Method used for the smoothing of the dielectric function")

    gp.add_argument("--dft-sccs-mixing", type=float, default=None,
        help="Mixing parameter (Hartree damping) employed during the iteration procedure, default is 0.6")

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

    gp.add_argument("--cell-opt-constraint", type=str, default=None,
        choices=["NONE", "X", "XY", "XZ", "Y", "YZ", "Z", "none", "x", "xy", "xz", "y", "yz", "z"],
        help="Imposes a constraint on the pressure tensor by fixing the specified cell components.")

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



def cp2kDriver(args):
    # ==============================================================================
    # CP2K CP2K CP2K CP2K CP2K CP2K CP2K CP2K CP2K CP2K CP2K CP2K CP2K CP2K C2PK CP2K
    # ==============================================================================        
    from pymatflow.cmd.cp2k_parser import read_inp
    if args.inp == None:
        params = {}
    else:
        params = read_inp(args.inp)
            
    #print(params)
    params["GLOBAL-PRINT_LEVEL"] = args.print_level if "GLOBAL-PRINT_LEVEL" not in params or args.print_level != None else params["GLOBAL-PRINT_LEVEL"]
    params["GLOBAL-EXTENDED_FFT_LENGTHS"] = args.extended_fft_lengths if "GLOBAL-EXTENDED_FFT_LENGTHS" not in params or args.extended_fft_lengths != None else params["GLOBAL-EXTENDED_FFT_LENGTHS"]

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
    params["FORCE_EVAL-DFT-QS-EXTRAPOLATION"] = args.qs_extrapolation if "FORCE_EVAL-DFT-QS-EXTRAPOLATION" not in params or args.qs_extrapolation != None else params["FORCE_EVAL-DFT-QS-EXTRAPOLATION"]
    params["FORCE_EVAL-DFT-QS-EXTRAPOLATION_ORDER"] = args.qs_extrapolation_order if "FORCE_EVAL-DFT-QS-EXTRAPOLATION_ORDER" not in params or args.qs_extrapolation_order != None else params["FORCE_EVAL-DFT-QS-EXTRAPOLATION_ORDER"]
    params["FORCE_EVAL-DFT-MGRID-CUTOFF"] = args.cutoff if "FORCE_EVAL-DFT-MGRID-CUTOFF" not in params or args.cutoff != None else params["FORCE_EVAL-DFT-MGRID-CUTOFF"]
    params["FORCE_EVAL-DFT-MGRID-REL_CUTOFF"] = args.rel_cutoff if "FORCE_EVAL-DFT-MGRID-REL_CUTOFF" not in params or args.rel_cutoff != None else params["FORCE_EVAL-DFT-MGRID-REL_CUTOFF"]
    params["FORCE_EVAL-DFT-MGRID-NGRIDS"] = args.ngrids if "FORCE_EVAL-DFT-MGRID-NGRIDS" not in params or  args.ngrids != None else params["FORCE_EVAL-DFT-MGRID-NGRIDS"]
    params["FORCE_EVAL-DFT-XC-XC_FUNCTIONAL"] = args.xc_functional if "FORCE_EVAL-DFT-XC-XC_FUNCTIONAL" not in params or  args.xc_functional != None else params["FORCE_EVAL-DFT-XC-XC_FUNCTIONAL"]

    params["FORCE_EVAL-DFT-XC-XC_GRID-USE_FINER_GRID"] = args.xc_grid_use_finer_grid if "FORCE_EVAL-DFT-XC-XC_GRID-USE_FINER_GRID" not in params or args.xc_grid_use_finer_grid != None else params["FORCE_EVAL-DFT-XC-XC_GRID-USE_FINER_GRID"]

    params["FORCE_EVAL-DFT-XC-XC_GRID-XC_DERIV"] = args.xc_grid_xc_deriv if "FORCE_EVAL-DFT-XC-XC_GRID-XC_DERIV" not in params or args.xc_grid_xc_deriv != None else params["FORCE_EVAL-DFT-XC-XC_GRID-XC_DERIV"]

    params["FORCE_EVAL-DFT-XC-XC_GRID-XC_SMOOTH_RHO"] = args.xc_grid_xc_smooth_rho if "FORCE_EVAL-DFT-XC-XC_GRID-XC_SMOOTH_RHO" not in params or args.xc_grid_xc_smooth_rho != None else params["FORCE_EVAL-DFT-XC-XC_GRID-XC_SMOOTH_RHO"]

    params["FORCE_EVAL-DFT-SCF-MAX_SCF"] = args.max_scf if "FORCE_EVAL-DFT-SCF-MAX_SCF" not in params or args.max_scf != None else params["FORCE_EVAL-DFT-SCF-MAX_SCF"]
    params["FORCE_EVAL-DFT-QS-EPS_DEFAULT"] = args.eps_default if "FORCE_EVAL-DFT-QS-EPS_DEFAULT" not in params or args.eps_default != None else params["FORCE_EVAL-DFT-QS-EPS_DEFAULT"]
    params["FORCE_EVAL-DFT-SCF-EPS_SCF"] = args.eps_scf if "FORCE_EVAL-DFT-SCF-EPS_SCF" not in params or  args.eps_scf != None else params["FORCE_EVAL-DFT-SCF-EPS_SCF"]
    params["FORCE_EVAL-DFT-SCF-LEVEL_SHIFT"] = args.level_shift if "FORCE_EVAL-DFT-SCF-LEVEL_SHIFT" not in params or  args.level_shift != None else params["FORCE_EVAL-DFT-SCF-LEVEL_SHIFT"]
    params["FORCE_EVAL-DFT-SCF-MAX_SCF_HISTORY"] = args.max_scf_history if "FORCE_EVAL-DFT-SCF-MAX_SCF_HISTORY" not in params or args.max_scf_history != None else params["FORCE_EVAL-DFT-SCF-MAX_SCF_HISTORY"]
    params["FORCE_EVAL-DFT-SCF-MAX_DIIS"] = args.max_diis if "FORCE_EVAL-DFT-SCF-MAX_DIIS" not in params or args.max_diis != None else params["FORCE_EVAL-DFT-SCF-MAX_DIIS"]
    params["FORCE_EVAL-DFT-SCF-ADDED_MOS"] = args.added_mos if "FORCE_EVAL-DFT-SCF-ADDED_MOS" not in params or  args.added_mos != None else params["FORCE_EVAL-DFT-SCF-ADDED_MOS"]
    params["FORCE_EVAL-DFT-SCF-CHOLESKY"] = args.dft_scf_cholesky if "FORCE_EVAL-DFT-SCF-CHOLESKY" not in params or args.dft_scf_cholesky != None else params["FORCE_EVAL-DFT-SCF-CHOLESKY"]


    params["FORCE_EVAL-DFT-SCF-OUTER_SCF"] = args.outer_scf if "FORCE_EVAL-DFT-SCF-OUTER_SCF" not in params or  args.outer_scf != None else params["FORCE_EVAL-DFT-SCF-OUTER_SCF"]
    params["FORCE_EVAL-DFT-SCF-OUTER_SCF-BISECT_TRUST_COUNT"] = args.outer_scf_bisect_trust_count if "FORCE_EVAL-DFT-SCF-OUTER_SCF-BISECT_TRUST_COUNT" not in params or  args.outer_scf_bisect_trust_count != None else params["FORCE_EVAL-DFT-SCF-OUTER_SCF-BISECT_TRUST_COUNT"]
    params["FORCE_EVAL-DFT-SCF-OUTER_SCF-DIIS_BUFFER_LENGTH"] = args.outer_scf_diis_buffer_length if "FORCE_EVAL-DFT-SCF-OUTER_SCF-DIIS_BUFFER_LENGTH" not in params or  args.outer_scf_diis_buffer_length != None else params["FORCE_EVAL-DFT-SCF-OUTER_SCF-DIIS_BUFFER_LENGTH"]
    params["FORCE_EVAL-DFT-SCF-OUTER_SCF-EXTRAPOLATION_ORDER"] = args.outer_scf_extrapolation_order if "FORCE_EVAL-DFT-SCF-OUTER_SCF-EXTRAPOLATION_ORDER" not in params or  args.outer_scf_extrapolation_order != None else params["FORCE_EVAL-DFT-SCF-OUTER_SCF-EXTRAPOLATION_ORDER"]
    params["FORCE_EVAL-DFT-SCF-OUTER_SCF-EPS_SCF"] = args.outer_scf_eps_scf if "FORCE_EVAL-DFT-SCF-OUTER_SCF-EPS_SCF" not in params or  args.outer_scf_eps_scf != None else params["FORCE_EVAL-DFT-SCF-OUTER_SCF-EPS_SCF"]
    params["FORCE_EVAL-DFT-SCF-OUTER_SCF-MAX_SCF"] = args.outer_scf_max_scf if "FORCE_EVAL-DFT-SCF-OUTER_SCF-MAX_SCF" not in params or  args.outer_scf_max_scf != None else params["FORCE_EVAL-DFT-SCF-OUTER_SCF-MAX_SCF"]
    params["FORCE_EVAL-DFT-SCF-OUTER_SCF-OPTIMIZER"] = args.outer_scf_optimizer if "FORCE_EVAL-DFT-SCF-OUTER_SCF-OPTIMIZER" not in params or  args.outer_scf_optimizer != None else params["FORCE_EVAL-DFT-SCF-OUTER_SCF-OPTIMIZER"]
    params["FORCE_EVAL-DFT-SCF-OUTER_SCF-TYPE"] = args.outer_scf_type if "FORCE_EVAL-DFT-SCF-OUTER_SCF-TYPE" not in params or  args.outer_scf_type != None else params["FORCE_EVAL-DFT-SCF-OUTER_SCF-TYPE"]

    params["FORCE_EVAL-DFT-SCF-SCF_GUESS"] = args.scf_guess if "FORCE_EVAL-DFT-SCF-SCF_GUESS" not in params or args.scf_guess != None else params["FORCE_EVAL-DFT-SCF-SCF_GUESS"]

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
    params["FORCE_EVAL-DFT-SCF-MIXING-BROY_W0"] = args.mixing_broy_w0 if "FORCE_EVAL-DFT-SCF-MIXING-BROY_W0" not in params or  args.mixing_broy_w0 != None else params["FORCE_EVAL-DFT-SCF-MIXING-BROY_W0"]
    params["FORCE_EVAL-DFT-SCF-MIXING-BROY_WMAX"] = args.mixing_broy_wmax if "FORCE_EVAL-DFT-SCF-MIXING-BROY_WMAX" not in params or  args.mixing_broy_wmax != None else params["FORCE_EVAL-DFT-SCF-MIXING-BROY_WMAX"]
    params["FORCE_EVAL-DFT-SCF-MIXING-BROY_WREF"] = args.mixing_broy_wref if "FORCE_EVAL-DFT-SCF-MIXING-BROY_WREF" not in params or  args.mixing_broy_wref != None else params["FORCE_EVAL-DFT-SCF-MIXING-BROY_WREF"]
    params["FORCE_EVAL-DFT-SCF-MIXING-NBUFFER"] = args.mixing_nbuffer if "FORCE_EVAL-DFT-SCF-MIXING-NBUFFER" not in params or  args.mixing_nbuffer != None else params["FORCE_EVAL-DFT-SCF-MIXING-NBUFFER"]        
    params["FORCE_EVAL-DFT-SCF-MIXING-GMIX_P"] = args.mixing_gmix_p if "FORCE_EVAL-DFT-SCF-MIXING-GMIX_P" not in params or  args.mixing_gmix_p != None else params["FORCE_EVAL-DFT-SCF-MIXING-GMIX_P"]
    params["FORCE_EVAL-DFT-SCF-MIXING-MAX_GVEC_EXP"] = args.mixing_max_gvec_exp if "FORCE_EVAL-DFT-SCF-MIXING-MAX_GVEC_EXP" not in params or args.mixing_max_gvec_exp != None else params["FORCE_EVAL-DFT-SCF-MIXING-MAX_GVEC_EXP"]
    params["FORCE_EVAL-DFT-SCF-MIXING-N_SIMPLE_MIX"] = args.mixing_n_simple_mix if "FORCE_EVAL-DFT-SCF-MIXING-N_SIMPLE_MIX" not in params or args.mixing_n_simple_mix != None else params["FORCE_EVAL-DFT-SCF-MIXING-N_SIMPLE_MIX"]
    params["FORCE_EVAL-DFT-SCF-MIXING-NMIXING"] = args.mixing_nmixing if "FORCE_EVAL-DFT-SCF-MIXING-NMIXING" not in params or args.mixing_nmixing != None else params["FORCE_EVAL-DFT-SCF-MIXING-NMIXING"]
    params["FORCE_EVAL-DFT-SCF-MIXING-NSKIP"] = args.mixing_nskip if "FORCE_EVAL-DFT-SCF-MIXING-NSKIP" not in params or args.mixing_nskip != None else params["FORCE_EVAL-DFT-SCF-MIXING-NSKIP"]
    params["FORCE_EVAL-DFT-SCF-MIXING-PULAY_ALPHA"] = args.mixing_pulay_alpha if "FORCE_EVAL-DFT-SCF-MIXING-PULAY_ALPHA" not in params or args.mixing_pulay_alpha != None else params["FORCE_EVAL-DFT-SCF-MIXING-PULAY_ALPHA"]
    params["FORCE_EVAL-DFT-SCF-MIXING-PULAY_BETA"] = args.mixing_pulay_beta if "FORCE_EVAL-DFT-SCF-MIXING-PULAY_BETA" not in params or args.mixing_pulay_beta != None else params["FORCE_EVAL-DFT-SCF-MIXING-PULAY_BETA"]
    params["FORCE_EVAL-DFT-SCF-MIXING-REGULARIZATION"] = args.mixing_regularization if "FORCE_EVAL-DFT-SCF-MIXING-REGULARIZATION" not in params or args.mixing_regularization != None else params["FORCE_EVAL-DFT-SCF-MIXING-REGULARIZATION"]


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


    params["FORCE_EVAL-DFT-SCCS"] = args.dft_sccs if "FORCE_EVAL-DFT-SCCS" not in params or  args.dft_sccs != None else params["FORCE_EVAL-DFT-SCCS"]
    params["FORCE_EVAL-DFT-SCCS-ALPHA"] = args.dft_sccs_alpha if "FORCE_EVAL-DFT-SCCS-ALPHA" not in params or  args.dft_sccs_alpha != None else params["FORCE_EVAL-DFT-SCCS-ALPHA"]
    params["FORCE_EVAL-DFT-SCCS-BETA"] = args.dft_sccs_beta if "FORCE_EVAL-DFT-SCCS-BETA" not in params or  args.dft_sccs_beta != None else params["FORCE_EVAL-DFT-SCCS-BETA"]
    params["FORCE_EVAL-DFT-SCCS-DELTA_RHO"] = args.dft_sccs_delta_rho if "FORCE_EVAL-DFT-SCCS-DELTA_RHO" not in params or  args.dft_sccs_delta_rho != None else params["FORCE_EVAL-DFT-SCCS-DELTA_RHO"]
    params["FORCE_EVAL-DFT-SCCS-DIELECTRIC_CONSTANT"] = args.dft_sccs_dielectric_constant if "FORCE_EVAL-DFT-SCCS-DIELECTRIC_CONSTANT" not in params or  args.dft_sccs_dielectric_constant != None else params["FORCE_EVAL-DFT-SCCS-DIELECTRIC_CONSTANT"]
    params["FORCE_EVAL-DFT-SCCS-EPS_SCCS"] = args.dft_sccs_eps_sccs if "FORCE_EVAL-DFT-SCCS-EPS_SCCS" not in params or  args.dft_sccs_eps_sccs != None else params["FORCE_EVAL-DFT-SCCS-EPS_SCCS"]
    params["FORCE_EVAL-DFT-SCCS-EPS_SCF"] = args.dft_sccs_eps_scf if "FORCE_EVAL-DFT-SCCS-EPS_SCF" not in params or  args.dft_sccs_eps_scf != None else params["FORCE_EVAL-DFT-SCCS-EPS_SCF"]
    params["FORCE_EVAL-DFT-SCCS-GAMMA"] = args.dft_sccs_gamma if "FORCE_EVAL-DFT-SCCS-GAMMA" not in params or  args.dft_sccs_gamma != None else params["FORCE_EVAL-DFT-SCCS-GAMMA"]
    params["FORCE_EVAL-DFT-SCCS-MAX_ITER"] = args.dft_sccs_max_iter if "FORCE_EVAL-DFT-SCCS-MAX_ITER" not in params or  args.dft_sccs_max_iter != None else params["FORCE_EVAL-DFT-SCCS-MAX_ITER"]
    params["FORCE_EVAL-DFT-SCCS-METHOD"] = args.dft_sccs_method if "FORCE_EVAL-DFT-SCCS-METHOD" not in params or  args.dft_sccs_method != None else params["FORCE_EVAL-DFT-SCCS-METHOD"]
    params["FORCE_EVAL-DFT-SCCS-MIXING"] = args.dft_sccs_mixing if "FORCE_EVAL-DFT-SCCS-MIXING" not in params or  args.dft_sccs_mixing != None else params["FORCE_EVAL-DFT-SCCS-MIXING"]

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
    params["MOTION-CELL_OPT-CONSTRAINT"] = args.cell_opt_constraint if "MOTION-CELL_OPT-CONSTRAINT" not in params or args.cell_opt_constraint != None else params["MOTION-CELL_OPT-CONSTRAINT"]

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
        from pymatflow.cmd.matflow import get_kpath
        if "stable" == args.flow_cp2k_version:
            from pymatflow.cp2k.version1 import StaticRun
        elif "dev" == args.flow_cp2k_version:
            from pymatflow.cp2k.version2 import StaticRun
        
        task = StaticRun()
        task.get_xyz(xyzfile)
        task.set_params(params=params)
        task.set_pot_basis(kind_basis=kind_basis, kind_pot=kind_pot, basis_set_file=basis_file, potential_file=pot_file)
        task.set_printout(option=args.printout_option)
        if 2 in args.printout_option:
            task.set_band(kpath=get_kpath(args.kpath_manual, args.kpath_file))
            task.set_vdw(usevdw=True if args.vdw_potential_type.lower() != "none" else False)
            task.set_run(mpi=args.mpi, server=server, jobname=args.jobname, nodes=args.nodes, ppn=args.ppn, queue=args.queue)
            task.set_llhpc(partition=args.partition, nodes=args.nodes, ntask=args.ntask, jobname=args.jobname, stdout=args.stdout, stderr=args.stderr)
            task.scf(directory=args.directory, runopt=args.runopt, auto=args.auto)
    elif args.runtype == 1:
        # geo opt
        if "stable" == args.flow_cp2k_version:
            from pymatflow.cp2k.version1 import OptRun
        elif "dev" == args.flow_cp2k_version:
            from pymatflow.cp2k.version2 import OptRun
        task = OptRun()
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
        if "stable" == args.flow_cp2k_version:
            from pymatflow.cp2k.version1 import OptRun
        elif "dev" == args.flow_cp2k_version:
            from pymatflow.cp2k.version2 import OptRun
        task = OptRun()
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
        if "stable" == args.flow_cp2k_version:
            from pymatflow.cp2k.version1 import OptRun
        elif "dev" == args.flow_cp2k_version:
            from pymatflow.cp2k.version2 import OptRun
        task = OptRun()
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
        if "stable" == args.flow_cp2k_version:
            from pymatflow.cp2k.version1 import OptRun
        elif "dev" == args.flow_cp2k_version:
            from pymatflow.cp2k.version2 import OptRun
        task = OptRun()
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
        if "stable" == args.flow_cp2k_version:
            from pymatflow.cp2k.version1 import OptRun
        elif "dev" == args.flow_cp2k_version:
            from pymatflow.cp2k.version2 import OptRun
        task = OptRun()
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
        if "stable" == args.flow_cp2k_version:
            from pymatflow.cp2k.version1 import NebRun
        elif "dev" == args.flow_cp2k_version:
            from pymatflow.cp2k.version2 import NebRun
        task = NebRun()
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
        if "stable" == args.flow_cp2k_version:
            from pymatflow.cp2k.version1 import PhonopyRun
        elif "dev" == args.flow_cp2k_version:
            from pymatflow.cp2k.version2 import PhonopyRun
        task = PhonopyRun()
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
        if "stable" == args.flow_cp2k_version:
            from pymatflow.cp2k.version1 import VibRun
        elif "dev" == args.flow_cp2k_version:
            from pymatflow.cp2k.version2 import VibRun
        task = VibRun()
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
        if "stable" == args.flow_cp2k_version:
            from pymatflow.cp2k.version1 import StaticRun
        elif "dev" == args.flow_cp2k_version:
            from pymatflow.cp2k.version2 import StaticRun
        task = StaticRun()
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
        if "stable" == args.flow_cp2k_version:
            from pymatflow.cp2k.version1 import MdRun
        elif "dev" == args.flow_cp2k_version:
            from pymatflow.cp2k.version2 import MdRun
        task = MdRun()
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
        if "stable" == args.flow_cp2k_version:
            from pymatflow.cp2k.version1 import OptRun
        elif "dev" == args.flow_cp2k_version:
            from pymatflow.cp2k.version2 import OptRun
        task = OptRun()
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
        if "stable" == args.flow_cp2k_version:
            from pymatflow.cp2k.version1 import MdRun
        elif "dev" == args.flow_cp2k_version:
            from pymatflow.cp2k.version2 import MdRun
        task = MdRun()
        task.get_xyz(xyzfile)
        task.set_params(params=params)
        task.set_pot_basis(kind_basis=kind_basis, kind_pot=kind_pot, basis_set_file=basis_file, potential_file=pot_file)
        task.set_vdw(usevdw=True if args.vdw_potential_type.lower() != "none" else False)
        task.set_run(mpi=args.mpi, server=server, jobname=args.jobname, nodes=args.nodes, ppn=args.ppn, queue=args.queue)
        task.set_llhpc(partition=args.partition, nodes=args.nodes, ntask=args.ntask, jobname=args.jobname, stdout=args.stdout, stderr=args.stderr)
        task.metadynamics(directory=args.directory, runopt=args.runopt, auto=args.auto)
    else:                
        pass        