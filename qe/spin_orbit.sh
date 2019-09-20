#!/bin/bash

<<'COMMENT'
The calculation proceeds as follows:

1) make a self-consistent calculation for Pt (input=pt.scf.in,
   output=pt.scf.out). 

2) make a band structure calculation for Pt (input=pt.nscf.in,
   output=pt.nscf.out).
COMMENT

# environmental variables
PW_COMMAND=pw.x
PSEUDO_DIR=
TMP_DIR=


# self-consistent calculation
cat > pt.scf.in << EOF
Pt
Pt
 &control
    calculation = 'scf'
    restart_mode='from_scratch',
    prefix='Pt',
    tprnfor = .true.,
    tstress =.true.,
    pseudo_dir = '$PSEUDO_DIR/',
    outdir='$TMP_DIR/'
 /
 &system
    ibrav=  2, celldm(1) =7.42, nat=  1, ntyp= 1,
    lspinorb=.true.,
    noncolin=.true.,
    starting_magnetization=0.0,
    occupations='smearing',
    smearing = 'marzari-vanderbilt'
    degauss=0.02,
    ecutwfc =30.0,
    ecutrho =250.0,
 /
 &electrons
    mixing_beta = 0.7,
    conv_thr =  1.0d-8
 /
ATOMIC_SPECIES
Pt  0.0   Pt.rel-pz-n-rrkjus.UPF
ATOMIC_POSITIONS
Pt  0.0000000   0.00000000   0.0
K_POINTS AUTOMATIC
4 4 4 1 1 1
EOF
echo "  running the scf calculation for Pt with spin-orbit coupling...\c"
$PW_COMMAND < pt.scf.in > pt.scf.out
echo " done"

# a non self-consistent calculation
cat > pt.nscf.in << EOF
Pt
Pt
 &control
    calculation = 'nscf'
    restart_mode='from_scratch',
    prefix='Pt',
    tprnfor = .true.
    pseudo_dir = '$PSEUDO_DIR/',
    outdir='$TMP_DIR/'
 /
 &system
    ibrav=  2, celldm(1) =7.42, nat=  1, ntyp= 1,
    lspinorb=.true.,
    noncolin=.true.,
    starting_magnetization=0.0,
    occupations='smearing',
    smearing = 'marzari-vanderbilt'
    degauss=0.02,
    ecutwfc =30.0,
    ecutrho =250.0,
 /
 &electrons
    mixing_beta = 0.7,
    conv_thr =  1.0d-8
 /
ATOMIC_SPECIES
Pt  0.0   Pt.rel-pz-n-rrkjus.UPF 
ATOMIC_POSITIONS
Pt  0.0000000   0.00000000   0.0
K_POINTS
8
0.0 0.0 0.0 1.0
0.1 0.0 0.0 1.0
1.0 0.0 0.0 1.0
0.4 0.2 0.1 1.0
0.4 0.4 0.0 1.0
0.4 0.4 0.4 1.0
0.5 0.5 0.5 1.0
0.75 0.75 0.0 1.0
EOF
echo "  running the non-scf calculation for Pt with spin-orbit coupling...\c"
$PW_COMMAND < pt.nscf.in > pt.nscf.out
echo " done"
# clean TMP_DIR
echo "  cleaning $TMP_DIR...\c"
rm -rf $TMP_DIR/Pt*
echo " done"


echo
echo "done"
