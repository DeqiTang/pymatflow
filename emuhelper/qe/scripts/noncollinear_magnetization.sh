#!/bin/bash

<<'COMMENT'
1) make a self-consistent calculation for Fe (input=fe.scf.in,
   output=fe.scf.out). The number of computed bands is internally
   computed as equal to the number of electrons in the unit cell
   (16 in this case).

2) make a band structure calculation for Fe (input=fe.band.in,
   output=fe.band.out).
   The variable nbnd is explicitly set = 16.
   The list of k points given in input is the list of point where the
   bands are computed, the k-point weight is arbitrary and is not used.

3) make a self-consistent calculation for Fe with penalty functional
   where each component of the magnetization of the two atoms
   is constrained (input=fe.pen.in, output=fe.pen.out).
   Iron is a metal : the smearing technique is used for the 
   calculation of the Fermi energy (a value for the broadening
   degauss is provided).

4) make a self-consistent calculation for Fe with penalty functional
   where the angle between the direction of the magnetization of each atom
   and the z axis is constrained; mcons(1) = cosine of this angle.
   (input=fe.angl.in, output=fe.angl.out).

5) make a self-consistent calculation for Fe with penalty functional
   where each component of the total magnetization is constrained; 
   fixed_magnetization(ipol) = value of the magnetization.
   (input=fe.total.in, output=fe.total.out).

6) make a self-consistent calculation for Cu (input=cu.scf.in,
   output=cu.scf.out).
   Copper is also a metal. In this case the tetrahedron method is used
   for the calculation of the Fermi energy. K-points are automatically
   generated.

7) make a band structure calculation for Cu (input=cu.band.in,
   output=cu.band.out).
   The variable nbnd is explicitly set = 8.
   The list of k points given in input is the list of point where the
   bands are computed, the k-point weight is arbitrary and is not used.

8) make a self-consistent calculation for Cu (input=cu.cg.in,
   output=cu.cg.out) with cg diagonalization.

9) make a self-consistent calculation for Cu (input=cu.diis.in,
   output=cu.diis.out) with diis diagonalization.


10) make a self-consistent calculation for Ni (input=ni.scf.in,
   output=ni.scf.out).
   Nickel is a magnetic metal. A local-spin-density calculation is
   performed by specifying nspin=2 and an initial guess for the
   magnetization of each atomic species. This initial guess is used to
   build spin-up and spin-down starting charges from superposition of
   atomic charges.

11) make a band structure calculation for Ni (input=ni.band.in,
   output=ni.band.out).

12) make a scf calculation of molecular oxygen relaxing the atoms.
COMMENT



# environment variables
PW_COMMAND=pw.x
PSEUDO_DIR=
TMP_DIR=



# self-consistent calculation
cat > fe.scf.in << EOF
Fe
Iron
 &control
    calculation='scf'
    restart_mode='from_scratch',
    pseudo_dir = '$PSEUDO_DIR/',
    outdir='$TMP_DIR/',
    prefix='fe'
 /
 &system
    ibrav = 3, celldm(1) =5.217, nat= 1, ntyp= 1,
    ecutwfc = 25.0,ecutrho = 200.0,
    report=1,
    occupations='smearing', smearing='marzari-vanderbilt', degauss=0.05
    noncolin = .true.
    starting_magnetization(1) = 0.5
    angle1(1) = 90.0
    angle2(1) =  0.0
 /
 &electrons
    conv_thr = 1.0e-8
    mixing_beta = 0.2
 /
ATOMIC_SPECIES
 Fe 55.847 Fe.pz-nd-rrkjus.UPF
ATOMIC_POSITIONS
 Fe 0.0 0.0 0.0
K_POINTS
 11
   0.0625000  0.0625000  0.0625000   1.00
   0.0625000  0.0625000  0.1875000   3.00
   0.0625000  0.0625000  0.3125000   3.00
   0.0625000  0.0625000  0.4375000   3.00
   0.0625000  0.0625000  0.5625000   3.00
   0.0625000  0.0625000  0.6875000   3.00
   0.0625000  0.0625000  0.8125000   3.00
   0.0625000  0.0625000  0.9375000   3.00
   0.0625000  0.1875000  0.1875000   3.00
   0.0625000  0.1875000  0.3125000   6.00
   0.0625000  0.1875000  0.4375000   6.00
EOF
echo "  running the scf calculation for Fe...\c"
$PW_COMMAND < fe.scf.in > fe.scf.out
echo " done"

# band structure calculation
cat > fe.band.in << EOF
Fe
Iron
 &control
    calculation='bands'
    pseudo_dir = '$PSEUDO_DIR/',
    outdir='$TMP_DIR/',
    prefix='fe'
 /
 &system
    ibrav = 3, celldm(1) =5.217, nat= 1, ntyp= 1,
    ecutwfc = 25.0,ecutrho = 200.0,
    report=1, nbnd = 16
    noncolin = .true.
    starting_magnetization(1) = 0.5
    angle1(1) = 90.0
    angle2(1) =  0.0
 /
 &electrons
    conv_thr = 1.0e-8
    mixing_beta = 0.2
 /
ATOMIC_SPECIES
 Fe 55.847 Fe.pz-nd-rrkjus.UPF
ATOMIC_POSITIONS
 Fe 0.0 0.0 0.0
K_POINTS
 28
   0.0 0.0 0.0 1.0
   0.0 0.0 0.1 1.0
   0.0 0.0 0.2 1.0
   0.0 0.0 0.3 1.0
   0.0 0.0 0.4 1.0
   0.0 0.0 0.5 1.0
   0.0 0.0 0.6 1.0
   0.0 0.0 0.7 1.0
   0.0 0.0 0.8 1.0
   0.0 0.0 0.9 1.0
   0.0 0.0 1.0 1.0
   0.0 0.0 0.0 1.0
   0.0 0.1 0.1 1.0
   0.0 0.2 0.2 1.0
   0.0 0.3 0.3 1.0
   0.0 0.4 0.4 1.0
   0.0 0.5 0.5 1.0
   0.0 0.6 0.6 1.0
   0.0 0.7 0.7 1.0
   0.0 0.8 0.8 1.0
   0.0 0.9 0.9 1.0
   0.0 1.0 1.0 1.0
   0.0 0.0 0.0 1.0
   0.1 0.1 0.1 1.0
   0.2 0.2 0.2 1.0
   0.3 0.3 0.3 1.0
   0.4 0.4 0.4 1.0
   0.5 0.5 0.5 1.0
EOF
echo "  running band structure calculation for Fe...\c"
$PW_COMMAND < fe.band.in > fe.band.out
echo " done"

# clean TMP_DIR
echo "  cleaning $TMP_DIR...\c"
rm -rf $TMP_DIR/fe*
echo " done"

# self-consistent calculation with penalty functional
cat > fe.pen.in << EOF
 &control
    calculation='scf'
    restart_mode='from_scratch',
    pseudo_dir = '$PSEUDO_DIR/',
    outdir='$TMP_DIR/',
    prefix='fe'
 /
 &system
    ibrav = 3, celldm(1) =5.217, nat= 1, ntyp= 1,
    ecutwfc = 25.0,ecutrho = 200.0,
    report=1,
    occupations='smearing', smearing='marzari-vanderbilt', degauss=0.05
    noncolin = .true.
    starting_magnetization(1) = 0.5
    angle1(1) = 85.0
    angle2(1) =  0.0
    constrained_magnetization='atomic'
    lambda = 1
 /
 &electrons
    conv_thr = 1.0e-8
    mixing_beta = 0.2
 /
ATOMIC_SPECIES
 Fe 55.847 Fe.pz-nd-rrkjus.UPF
ATOMIC_POSITIONS
 Fe 0.0 0.0 0.0
K_POINTS
 11
   0.0625000  0.0625000  0.0625000   1.00
   0.0625000  0.0625000  0.1875000   3.00
   0.0625000  0.0625000  0.3125000   3.00
   0.0625000  0.0625000  0.4375000   3.00
   0.0625000  0.0625000  0.5625000   3.00
   0.0625000  0.0625000  0.6875000   3.00
   0.0625000  0.0625000  0.8125000   3.00
   0.0625000  0.0625000  0.9375000   3.00
   0.0625000  0.1875000  0.1875000   3.00
   0.0625000  0.1875000  0.3125000   6.00
   0.0625000  0.1875000  0.4375000   6.00
EOF
echo "  running scf calculation for Fe with penalty functional...\c"
$PW_COMMAND < fe.pen.in > fe.pen.out
echo " done"

# scf calculation with penalty functional (angle with z-axis constrained)
cat > fe.angl.in << EOF
 &control
    calculation='scf'
    restart_mode='from_scratch',
    pseudo_dir = '$PSEUDO_DIR/',
    outdir='$TMP_DIR/',
    prefix='fe'
 /
 &system
    ibrav = 3, celldm(1) =5.217, nat= 1, ntyp= 1,
    ecutwfc = 25.0,ecutrho = 200.0,
    report=1,
    occupations='smearing', smearing='marzari-vanderbilt', degauss=0.05
    noncolin = .true.
    starting_magnetization(1) = 0.5
    angle1(1) = 90.0
    angle2(1) =  0.0
    constrained_magnetization='atomic direction'
    lambda = 1
 /
 &electrons
    conv_thr = 1.0e-8
    mixing_beta = 0.2
 /
ATOMIC_SPECIES
 Fe 55.847 Fe.pz-nd-rrkjus.UPF
ATOMIC_POSITIONS
 Fe 0.0 0.0 0.0
K_POINTS
 11
   0.0625000  0.0625000  0.0625000   1.00
   0.0625000  0.0625000  0.1875000   3.00
   0.0625000  0.0625000  0.3125000   3.00
   0.0625000  0.0625000  0.4375000   3.00
   0.0625000  0.0625000  0.5625000   3.00
   0.0625000  0.0625000  0.6875000   3.00
   0.0625000  0.0625000  0.8125000   3.00
   0.0625000  0.0625000  0.9375000   3.00
   0.0625000  0.1875000  0.1875000   3.00
   0.0625000  0.1875000  0.3125000   6.00
   0.0625000  0.1875000  0.4375000   6.00
EOF
echo "  running the scf calculation for Fe with constrained angle...\c"
$PW_COMMAND < fe.angl.in > fe.angl.out
echo " done"

# clean TMP_DIR
echo "  cleaning $TMP_DIR...\c"
rm -rf $TMP_DIR/fe*
echo " done"

# scf calculation with penalty functional (total magnetization constrained)
cat > fe.total.in << EOF
 &control
    calculation='scf'
    restart_mode='from_scratch',
    pseudo_dir = '$PSEUDO_DIR/',
    outdir='$TMP_DIR/',
    prefix='fe'
 /
 &system
    ibrav = 3, celldm(1) =5.217, nat= 1, ntyp= 1,
    ecutwfc = 25.0,ecutrho = 200.0,
    report=1,
    occupations='smearing', smearing='marzari-vanderbilt', degauss=0.05
    noncolin = .true.
    starting_magnetization(1) = 0.5
    angle1(1) = 45.0
    angle2(1) = 53.0
    constrained_magnetization='total'
    fixed_magnetization(1)=0.3,
    fixed_magnetization(2)=0.4,
    fixed_magnetization(3)=0.5,
    lambda = 0.5
 /
 &electrons
    conv_thr = 1.0e-9
    mixing_beta = 0.3
 /
ATOMIC_SPECIES
 Fe 55.847 Fe.pz-nd-rrkjus.UPF
ATOMIC_POSITIONS
 Fe 0.0 0.0 0.0
K_POINTS AUTOMATIC
4 4 4 1 1 1
EOF
echo "  running the scf calculation for Fe with constrained magnetization...\c"
$PW_COMMAND < fe.total.in > fe.total.out
echo " done"

# clean TMP_DIR
echo "  cleaning $TMP_DIR...\c"
rm -rf $TMP_DIR/fe*
echo " done"

# self-consistent calculation
cat > cu.scf.in << EOF
 &control
    calculation='scf'
    restart_mode='from_scratch',
    pseudo_dir = '$PSEUDO_DIR/',
    outdir='$TMP_DIR/'
    prefix='cu'
 /
 &system
    ibrav = 2, celldm(1) =6.73, nat= 1, ntyp= 1,
    ecutwfc = 25.0, ecutrho = 300.0
    occupations='smearing', smearing='marzari-vanderbilt', degauss=0.02
    noncolin = .true.
    starting_magnetization(1) = 0.5
    angle1(1) = 90.0
    angle2(1) =  0.0
 /
 &electrons
    conv_thr = 1.0e-8
    mixing_beta = 0.7
 /
ATOMIC_SPECIES
 Cu 63.55 Cu.pz-d-rrkjus.UPF
ATOMIC_POSITIONS
 Cu 0.0 0.0 0.0
K_POINTS (automatic)
 8 8 8 0 0 0
EOF
echo "  running the scf calculation for Cu...\c"
$PW_COMMAND < cu.scf.in > cu.scf.out
echo " done"

# band structure calculation along delta, sigma and lambda lines
cat > cu.band.in << EOF
 &control
    calculation='bands'
    pseudo_dir = '$PSEUDO_DIR/',
    outdir='$TMP_DIR/',
    prefix='cu'
 /
 &system
    ibrav = 2, celldm(1) =6.73, nat= 1, ntyp= 1,
    ecutwfc = 25.0, ecutrho = 300.0, nbnd = 8
    noncolin = .true.
    starting_magnetization(1) = 0.5
    angle1(1) = 90.0
    angle2(1) =  0.0
 /
 &electrons
 /
ATOMIC_SPECIES
 Cu  63.55 Cu.pz-d-rrkjus.UPF
ATOMIC_POSITIONS
 Cu 0.0 0.0 0.0
K_POINTS
 28
   0.0 0.0 0.0 1.0
   0.0 0.0 0.1 1.0
   0.0 0.0 0.2 1.0
   0.0 0.0 0.3 1.0
   0.0 0.0 0.4 1.0
   0.0 0.0 0.5 1.0
   0.0 0.0 0.6 1.0
   0.0 0.0 0.7 1.0
   0.0 0.0 0.8 1.0
   0.0 0.0 0.9 1.0
   0.0 0.0 1.0 1.0
   0.0 0.0 0.0 1.0
   0.0 0.1 0.1 1.0
   0.0 0.2 0.2 1.0
   0.0 0.3 0.3 1.0
   0.0 0.4 0.4 1.0
   0.0 0.5 0.5 1.0
   0.0 0.6 0.6 1.0
   0.0 0.7 0.7 1.0
   0.0 0.8 0.8 1.0
   0.0 0.9 0.9 1.0
   0.0 1.0 1.0 1.0
   0.0 0.0 0.0 1.0
   0.1 0.1 0.1 1.0
   0.2 0.2 0.2 1.0
   0.3 0.3 0.3 1.0
   0.4 0.4 0.4 1.0
   0.5 0.5 0.5 1.0
EOF
echo "  running the band-structure calculation for Cu...\c"
$PW_COMMAND < cu.band.in > cu.band.out
echo " done"

# clean TMP_DIR
echo "  cleaning $TMP_DIR...\c"
rm -rf $TMP_DIR/cu*
echo " done"

# self-consistent calculation with cg diagonalization
cat > cu.cg.in << EOF
 &control
    calculation='scf'
    restart_mode='from_scratch',
    pseudo_dir = '$PSEUDO_DIR/',
    outdir='$TMP_DIR/'
    prefix='cu'
 /
 &system
    ibrav = 2, celldm(1) =6.73, nat= 1, ntyp= 1,
    ecutwfc = 25.0, ecutrho = 300.0
    occupations='smearing', smearing='marzari-vanderbilt', degauss=0.02
    noncolin = .true.
    starting_magnetization(1) = 0.5
    angle1(1) = 90.0
    angle2(1) =  0.0
 /
 &electrons
    conv_thr = 1.0e-8
    mixing_beta = 0.7
    diagonalization = 'cg'
 /
ATOMIC_SPECIES
 Cu 63.55 Cu.pz-d-rrkjus.UPF
ATOMIC_POSITIONS
 Cu 0.0 0.0 0.0
K_POINTS (automatic)
 8 8 8 0 0 0
EOF
echo "  running the scf calculation for Cu with cg diagonalization...\c"
$PW_COMMAND < cu.cg.in > cu.cg.out
echo " done"

# clean TMP_DIR
echo "  cleaning $TMP_DIR...\c"
rm -rf $TMP_DIR/cu*
echo " done"

# self-consistent calculation with diis diagonalization
#cat > cu.diis.in << EOF
# &control
#    calculation='scf'
#    restart_mode='from_scratch',
#    pseudo_dir = '$PSEUDO_DIR/',
#    outdir='$TMP_DIR/'
#    prefix='cu'
# /
# &system
#    ibrav = 2, celldm(1) =6.73, nat= 1, ntyp= 1,
#    ecutwfc = 25.0, ecutrho = 300.0
#    occupations='smearing', smearing='marzari-vanderbilt', degauss=0.02
#    noncolin = .true.
#    starting_magnetization(1) = 0.5
#    angle1(1) = 90.0
#    angle2(1) =  0.0
# /
# &electrons
#    conv_thr = 1.0e-8
#    mixing_beta = 0.7
#    diagonalization = 'diis'
# /
#ATOMIC_SPECIES
# Cu 63.55 Cu.pz-d-rrkjus.UPF
#ATOMIC_POSITIONS
# Cu 0.0 0.0 0.0
#K_POINTS (automatic)
# 8 8 8 0 0 0
#EOF
# echo "  running the scf calculation for Cu with diis diagonalization...\c"
# $PW_COMMAND < cu.diis.in > cu.diis.out
# echo " done"

# clean TMP_DIR
echo "  cleaning $TMP_DIR...\c"
rm -rf $TMP_DIR/cu*
echo " done"

# self-consistent calculation
cat > ni.scf.in << EOF
 &control
    calculation='scf'
    restart_mode='from_scratch',
    pseudo_dir = '$PSEUDO_DIR/',
    outdir='$TMP_DIR/'
    prefix='ni'
 /
 &system
    ibrav=2, celldm(1) =6.48, nat=1, ntyp=1,
    ecutwfc = 24.0, ecutrho = 288.0,
    occupations='smearing', smearing='marzari-vanderbilt', degauss=0.02
    noncolin = .true.
    starting_magnetization(1) = 0.1
    angle1(1) = 90.0
    angle2(1) =  0.0
 /
 &electrons
    conv_thr = 1.0e-8
    mixing_beta = 0.7
 /
ATOMIC_SPECIES
 Ni 58.69 Ni.pbe-nd-rrkjus.UPF
ATOMIC_POSITIONS
 Ni 0.0 0.0 0.0
K_POINTS
 60
   0.0625000  0.0625000  0.0625000   1.00
   0.0625000  0.0625000  0.1875000   3.00
   0.0625000  0.0625000  0.3125000   3.00
   0.0625000  0.0625000  0.4375000   3.00
   0.0625000  0.0625000  0.5625000   3.00
   0.0625000  0.0625000  0.6875000   3.00
   0.0625000  0.0625000  0.8125000   3.00
   0.0625000  0.0625000  0.9375000   3.00
   0.0625000  0.1875000  0.1875000   3.00
   0.0625000  0.1875000  0.3125000   6.00
   0.0625000  0.1875000  0.4375000   6.00
   0.0625000  0.1875000  0.5625000   6.00
   0.0625000  0.1875000  0.6875000   6.00
   0.0625000  0.1875000  0.8125000   6.00
   0.0625000  0.1875000  0.9375000   6.00
   0.0625000  0.3125000  0.3125000   3.00
   0.0625000  0.3125000  0.4375000   6.00
   0.0625000  0.3125000  0.5625000   6.00
   0.0625000  0.3125000  0.6875000   6.00
   0.0625000  0.3125000  0.8125000   6.00
   0.0625000  0.3125000  0.9375000   6.00
   0.0625000  0.4375000  0.4375000   3.00
   0.0625000  0.4375000  0.5625000   6.00
   0.0625000  0.4375000  0.6875000   6.00
   0.0625000  0.4375000  0.8125000   6.00
   0.0625000  0.4375000  0.9375000   6.00
   0.0625000  0.5625000  0.5625000   3.00
   0.0625000  0.5625000  0.6875000   6.00
   0.0625000  0.5625000  0.8125000   6.00
   0.0625000  0.6875000  0.6875000   3.00
   0.0625000  0.6875000  0.8125000   6.00
   0.0625000  0.8125000  0.8125000   3.00
   0.1875000  0.1875000  0.1875000   1.00
   0.1875000  0.1875000  0.3125000   3.00
   0.1875000  0.1875000  0.4375000   3.00
   0.1875000  0.1875000  0.5625000   3.00
   0.1875000  0.1875000  0.6875000   3.00
   0.1875000  0.1875000  0.8125000   3.00
   0.1875000  0.3125000  0.3125000   3.00
   0.1875000  0.3125000  0.4375000   6.00
   0.1875000  0.3125000  0.5625000   6.00
   0.1875000  0.3125000  0.6875000   6.00
   0.1875000  0.3125000  0.8125000   6.00
   0.1875000  0.4375000  0.4375000   3.00
   0.1875000  0.4375000  0.5625000   6.00
   0.1875000  0.4375000  0.6875000   6.00
   0.1875000  0.4375000  0.8125000   6.00
   0.1875000  0.5625000  0.5625000   3.00
   0.1875000  0.5625000  0.6875000   6.00
   0.1875000  0.6875000  0.6875000   3.00
   0.3125000  0.3125000  0.3125000   1.00
   0.3125000  0.3125000  0.4375000   3.00
   0.3125000  0.3125000  0.5625000   3.00
   0.3125000  0.3125000  0.6875000   3.00
   0.3125000  0.4375000  0.4375000   3.00
   0.3125000  0.4375000  0.5625000   6.00
   0.3125000  0.4375000  0.6875000   6.00
   0.3125000  0.5625000  0.5625000   3.00
   0.4375000  0.4375000  0.4375000   1.00
   0.4375000  0.4375000  0.5625000   3.00
EOF
echo "  running the scf calculation for Ni...\c"
$PW_COMMAND < ni.scf.in > ni.scf.out
echo " done"

# band structure calculation along delta, sigma and lambda lines
cat > ni.band.in << EOF
 &control
    calculation='bands'
    pseudo_dir = '$PSEUDO_DIR/',
    outdir='$TMP_DIR/'
    prefix='ni'
 /
 &system
    ibrav=2, celldm(1) =6.48, nat=1, ntyp=1,
    starting_magnetization(1)=0.7,
    ecutwfc = 24.0, ecutrho = 288.0, nbnd = 8
    noncolin = .true.
    starting_magnetization(1) = 0.5
    angle1(1) = 90.0
    angle2(1) =  0.0

 /
 &electrons
 /
ATOMIC_SPECIES
 Ni 58.69 Ni.pbe-nd-rrkjus.UPF
ATOMIC_POSITIONS
 Ni 0.0 0.0 0.0
K_POINTS
 28
   0.0 0.0 0.0 1.0
   0.0 0.0 0.1 1.0
   0.0 0.0 0.2 1.0
   0.0 0.0 0.3 1.0
   0.0 0.0 0.4 1.0
   0.0 0.0 0.5 1.0
   0.0 0.0 0.6 1.0
   0.0 0.0 0.7 1.0
   0.0 0.0 0.8 1.0
   0.0 0.0 0.9 1.0
   0.0 0.0 1.0 1.0
   0.0 0.0 0.0 1.0
   0.0 0.1 0.1 1.0
   0.0 0.2 0.2 1.0
   0.0 0.3 0.3 1.0
   0.0 0.4 0.4 1.0
   0.0 0.5 0.5 1.0
   0.0 0.6 0.6 1.0
   0.0 0.7 0.7 1.0
   0.0 0.8 0.8 1.0
   0.0 0.9 0.9 1.0
   0.0 1.0 1.0 1.0
   0.0 0.0 0.0 1.0
   0.1 0.1 0.1 1.0
   0.2 0.2 0.2 1.0
   0.3 0.3 0.3 1.0
   0.4 0.4 0.4 1.0
   0.5 0.5 0.5 1.0
EOF
echo "  running the band-structure calculation for Ni...\c"
$PW_COMMAND < ni.band.in > ni.band.out
echo " done"

# clean TMP_DIR
echo "  cleaning $TMP_DIR...\c"
rm -rf $TMP_DIR/ni*
echo " done"

# self-consistent calculation + relaxation of atoms
cat > o2.relax.in << EOF
 &control
    calculation='relax'
    restart_mode='from_scratch',!'restart', !
    pseudo_dir = '$PSEUDO_DIR/',
    outdir='$TMP_DIR/'
    prefix='o2'
 /
 &system
    ibrav = 1, celldm(1) =7.50, nat= 2, ntyp= 2,
    ecutwfc = 25.0,ecutrho = 200.0,
    report=1,
    occupations='smearing', smearing='marzari-vanderbilt', degauss=0.05
    noncolin = .true.
    starting_magnetization(1) = 0.5
    angle1(1) = 90.0
    angle2(1) =  0.0
    starting_magnetization(2) = 0.5
    angle1(2) = 90.0
    angle2(2) =  0.0
 /
 &electrons
    mixing_beta = 0.2
 /
 &ions
 /
ATOMIC_SPECIES
 O1 16.0 O.pbe-rrkjus.UPF
 O2 16.0 O.pbe-rrkjus.UPF
ATOMIC_POSITIONS
 O1 0.0 0.0 0.0
 O2 0.20 0.20 0.20
K_POINTS
 1
   0.0 0.0 0.0 1.00
EOF
echo "  running scf calculation with relax for oxygen molecule...\c"
$PW_COMMAND < o2.relax.in > o2.relax.out
echo " done"
# clean TMP_DIR
echo "  cleaning $TMP_DIR...\c"
rm -rf $TMP_DIR/o.*
echo " done"


echo
echo "done"
