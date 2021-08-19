
def vaspSubparser(subparsers):
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
        choices=[0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11],
        help="choices of runtype. 0->static_run; 1->optimization; 2->cubic-cell; 3->hexagonal-cell; 4->tetragonal-cell; 5->neb; 6->vasp-phonon; 7->phonopy; 8->surf pes; 9->abc; 10->AIMD, 11->custom")

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
        choices=["pbs", "llhpc", "lsf_sz", "tianhe2", "lsf_sustc", "cdcloud"],
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
        help="NCORE determines the number of compute cores that work on an individual orbital.")

    gp.add_argument("--npar", type=int, default=None,
        help="NPAR determines the number of bands that are treated in parallel.")

    gp.add_argument("--kpar", type=int, default=None,
        help="KPAR determines the number of k-points that are to be treated in parallel")

    gp.add_argument("--lplane", type=str, default=None,
        choices=[".TRUE.", ".FALSE.", "T", "F"],
        help="LPLANE switches on the plane-wise data distribution in real space. default is .TRUE.")

    gp.add_argument("--nsim", type=int, default=None,
        help="NSIM sets the number of bands that are optimized simultaneously by the RMM-DIIS algorithm. default NSIM=4")

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

    gp.add_argument("--lepsilon", type=str, default=None,
        choices=["T", "F", ".TRUE.", ".FALSE."],
        help="determines the static dielectric matrix, ion-clamped piezoelectric tensor and the Born effective charges using density functional perturbation theory.")

    gp.add_argument("--lpead", type=str, default=None,
        choices=["T", "F", ".TRUE.", ".FALSE."],
        help="the derivative of the cell-periodic part of the orbitals w.r.t. k, |∇kunk〉, is calculated using finite differences.")

    gp.add_argument("--lrpa", type=str, default=None,
        choices=["T", "F", ".TRUE.", ".FALSE."],
        help="includes local field effect on the Hartree level only.")

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
           
    gp.add_argument("--lvhar", type=str, default=None,
        choices=[".TRUE.", ".FALSE.", "T", "F"],
        help="This tag determines whether the total local potential (saved in the file LOCPOT) contains the entire local potential (ionic + Hartree + exchange correlation) or the electrostatic contributions only (ionic + Hartree).")


def vaspDriver(args):

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
    params["NPAR"] = args.npar if "NPAR" not in params or args.npar != None else params["NPAR"]
    params["KPAR"] = args.kpar if "KPAR" not in params or args.kpar != None else params["KPAR"]
    params["LPLANE"] = args.lplane if "LPLANE" not in params or args.lplane != None else params["LPLANE"]
    params["NSIM"] = args.nsim if "NSIM" not in params or args.nsim != None else params["NSIM"]
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
    params["LVHAR"] = args.lvhar if "LVHAR" not in params or args.lvhar != None else params["LVHAR"] 
    params["LEPSILON"] = args.lepsilon if "LEPSILON" not in params or args.lepsilon != None else params["LEPSILON"]
    params["LPEAD"] = args.lpead if "LPEAD" not in params or args.lpead != None else params["LPEAD"]
    params["LRPA"] = args.lrpa if "LRPA" not in params or args.lrpa != None else params["LRPA"]

    if args.runtype == 0:
        # static
        from pymatflow.cmd.matflow import get_kpath
        from pymatflow.vasp.static import StaticRun
        task = StaticRun()
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
        from pymatflow.vasp.opt import OptRun
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
        task = OptRun()
        task.get_xyz(xyzfile)
        task.set_params(params=params, runtype="opt")
        task.set_kpoints(kpoints_mp=args.kpoints_mp)
        task.poscar.selective_dynamics = True if args.selective_dynamics.upper()[0] == "T" else False
        task.set_run(mpi=args.mpi, server=server, jobname=args.jobname, nodes=args.nodes, ppn=args.ppn, queue=args.queue)
        task.set_llhpc(partition=args.partition, nodes=args.nodes, ntask=args.ntask, jobname=args.jobname, stdout=args.stdout, stderr=args.stderr)
        task.set_cdcloud(partition=args.partition, nodes=args.nodes, ntask=args.ntask, jobname=args.jobname, stdout=args.stdout, stderr=args.stderr)
        task.optimize(directory=args.directory, runopt=args.runopt, auto=args.auto)
    elif args.runtype == 2:
        # cubic cell
        from pymatflow.vasp.opt import OptRun
        # some must set parameters 
        if params["IBRION"] == None:
            params["IBRION"] = 2
        params["ISIF"] = 2
        if params["NSW"] == None:
            params["NSW"] = 100
        task = OptRun()
        task.get_xyz(xyzfile)
        task.set_params(params=params, runtype="opt")
        task.set_kpoints(kpoints_mp=args.kpoints_mp)
        task.set_run(mpi=args.mpi, server=server, jobname=args.jobname, nodes=args.nodes, ppn=args.ppn, queue=args.queue)
        task.set_llhpc(partition=args.partition, nodes=args.nodes, ntask=args.ntask, jobname=args.jobname, stdout=args.stdout, stderr=args.stderr)
        task.set_cdcloud(partition=args.partition, nodes=args.nodes, ntask=args.ntask, jobname=args.jobname, stdout=args.stdout, stderr=args.stderr)
        task.batch_a = args.batch_a
        task.cubic(directory=args.directory, runopt=args.runopt, auto=args.auto, range_a=args.range_a)
    elif args.runtype == 3:
        # hexagonal cell
        from pymatflow.vasp.opt import OptRun
        task = OptRun()
        task.get_xyz(xyzfile)
        task.set_params(params=params, runtype="opt")
        task.set_kpoints(kpoints_mp=args.kpoints_mp)
        task.set_run(mpi=args.mpi, server=server, jobname=args.jobname, nodes=args.nodes, ppn=args.ppn, queue=args.queue)
        task.set_llhpc(partition=args.partition, nodes=args.nodes, ntask=args.ntask, jobname=args.jobname, stdout=args.stdout, stderr=args.stderr)
        task.set_cdcloud(partition=args.partition, nodes=args.nodes, ntask=args.ntask, jobname=args.jobname, stdout=args.stdout, stderr=args.stderr)
        task.batch_a = args.batch_a
        task.batch_c = args.batch_c            
        task.hexagonal(directory=args.directory, runopt=args.runopt, auto=args.auto, range_a=args.range_a, range_c=args.range_c)
    elif args.runtype == 4:
        # tetragonal cell
        from pymatflow.vasp.opt import OptRun
        task = OptRun()
        task.get_xyz(xyzfile)
        task.set_params(params=params, runtype="opt")
        task.set_kpoints(kpoints_mp=args.kpoints_mp)
        task.set_run(mpi=args.mpi, server=server, jobname=args.jobname, nodes=args.nodes, ppn=args.ppn, queue=args.queue)
        task.set_llhpc(partition=args.partition, nodes=args.nodes, ntask=args.ntask, jobname=args.jobname, stdout=args.stdout, stderr=args.stderr)
        task.set_cdcloud(partition=args.partition, nodes=args.nodes, ntask=args.ntask, jobname=args.jobname, stdout=args.stdout, stderr=args.stderr)
        task.batch_a = args.batch_a     
        task.batch_c = args.batch_c            
        task.tetragonal(directory=args.directory, runopt=args.runopt, auto=args.auto, range_a=args.range_a, range_c=args.range_c)
    elif args.runtype == 5:
        # neb
        # we better set NSW manually in VTST neb calc. 
        # if not set, pymatflow.vasp.neb will set it to 100 automatically
        from pymatflow.vasp.neb import NebRun
        task = NebRun()
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
        from pymatflow.vasp.phonon import PhononRun
        task = PhononRun() 
        task.get_xyz(xyzfile)
        task.set_params(params=params, runtype="phonon")
        task.set_kpoints(kpoints_mp=args.kpoints_mp)
        task.supercell_n = args.supercell_n
        task.set_run(mpi=args.mpi, server=server, jobname=args.jobname, nodes=args.nodes, ppn=args.ppn, queue=args.queue)
        task.set_llhpc(partition=args.partition, nodes=args.nodes, ntask=args.ntask, jobname=args.jobname, stdout=args.stdout, stderr=args.stderr)
        task.phonon(directory=args.directory, runopt=args.runopt, auto=args.auto)
    elif args.runtype == 7:
        # phonopy
        from pymatflow.vasp.phonopy import PhonopyRun
        task = PhonopyRun()
        task.get_xyz(xyzfile)
        task.set_params(params=params, runtype="phonopy")
        task.set_kpoints(kpoints_mp=args.kpoints_mp)
        task.supercell_n = args.supercell_n
        task.set_run(mpi=args.mpi, server=server, jobname=args.jobname, nodes=args.nodes, ppn=args.ppn, queue=args.queue)
        task.set_llhpc(partition=args.partition, nodes=args.nodes, ntask=args.ntask, jobname=args.jobname, stdout=args.stdout, stderr=args.stderr)
        task.phonopy(directory=args.directory, runopt=args.runopt, auto=args.auto)
    elif args.runtype == 8:
        # sur pes
        from pymatflow.flow.surface_pes import VaspRun
        task = VaspRun()
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
        from pymatflow.vasp.opt import OptRun
        task = OptRun()
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
        from pymatflow.vasp.md import MdRun
        task = MdRun()
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
    elif args.runtype == 11:
        # vasp custom
        from pymatflow.vasp.custom import CustomRun
        task = CustomRun() 
        task.get_xyz(xyzfile)
        task.set_params(params=params, runtype="custom")
        task.set_kpoints(kpoints_mp=args.kpoints_mp)
        task.set_run(mpi=args.mpi, server=server, jobname=args.jobname, nodes=args.nodes, ppn=args.ppn, queue=args.queue)
        task.set_llhpc(partition=args.partition, nodes=args.nodes, ntask=args.ntask, jobname=args.jobname, stdout=args.stdout, stderr=args.stderr)
        task.set_cdcloud(partition=args.partition, nodes=args.nodes, ntask=args.ntask, jobname=args.jobname, stdout=args.stdout, stderr=args.stderr)
        task.custom(directory=args.directory, runopt=args.runopt, auto=args.auto)            
    else:
        pass                            