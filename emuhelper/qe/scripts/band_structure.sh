#!/bin/bash

# 做能带计算以前需要做单点能计算: single point

# ===========
# File Names of Pseudopotential of Elements
Pse_B='B.pbe-n-kjpaw_psl.1.0.0.UPF'
Pse_N='N.pbe-n-kjpaw_psl.1.0.0.UPF'
# ===========
cat >BN.scf.in<<EOF
&CONTROL

calculation='bands'

title='BN'

prefix='BN'

restart_mode='from_scratch'

nstep=1000

outdir='./tmp'

pseudo_dir='./'

wf_collect=.true.

tstress=.true.

tprnfor=.true.

/

&SYSTEM

ibrav= 0

nat=2

ntyp=2

nbnd=

ecutwfc = 60.0,

ecutrho = 600.0,

input_DFT ='PBE',

occupations ='smearing',

degauss =1.0d-4,

smearing ='marzari-vanderbilt',

/

&ELECTRONS

electron_maxstep =1000,

conv_thr =1.0d-10,

mixing_mode = 'plain',

mixing_beta = 0.3d0,

scf_must_converge= .true.

/

ATOMIC_SPECIES

B 10.81 ${Pse_B}

N 14.01 ${Pse_N}

CELL_PARAMETERS (angstrom)

2.4700000286         0.0000000000         0.0000000000

-1.2350000143         2.1390827721         0.0000000000

0.0000000000         0.0000000000        10.0000000000

ATOMIC_POSITIONS (crystal)

B 0.666666687         0.333333343         0.500000000

N 0.333333313         0.666666627         0.500000000

K_POINTS automatic

12 12 1 0 0 0


EOF

cat >BN.relax.in<<EOF
&CONTROL

calculation='relax'

title='BN'

prefix='BN'

restart_mode='from_scratch'

nstep=1000

outdir='./tmp'

pseudo_dir='./'

wf_collect=.true.

tstress=.true.

tprnfor=.true.

/

&SYSTEM

ibrav= 0

nat=2

ntyp=2

ecutwfc = 60.0,

ecutrho = 600.0,

input_DFT ='PBE',

occupations = 'smearing',

degauss = 1.0d-4,

smearing = 'marzari-vanderbilt',

/

&ELECTRONS

electron_maxstep = 1000,

conv_thr = 1.0d-10,

mixing_mode = 'plain',

mixing_beta = 0.3d0,

scf_must_converge= .true.

/

&IONS

ion_dynamics='bfgs'

ion_positions = 'default'

/

ATOMIC_SPECIES

B 10.81 ${Pse_B}

N 14.01 ${Pse_N}

CELL_PARAMETERS (angstrom)

 2.511218514   0.000000000   0.000000000

-1.255609257   2.174779027   0.000000000

 0.000000000   0.000000000  10.000000000

ATOMIC_POSITIONS (crystal)

B       0.666666687   0.333333343   0.500000000

N       0.333333313   0.666666627   0.500000000

K_POINTS automatic

12 12 1 0 0 0
EOF


# pw.x < BN.vc-relax.in > BN.vc-relax.out
# pw.x < BN.relax.in > BN.relax.out

mpirun -np 3 pw.x < BN.vc-relax.in # > BN.vec-ralax.in
pw.x < BN.relax.in # > BN.relax.in


