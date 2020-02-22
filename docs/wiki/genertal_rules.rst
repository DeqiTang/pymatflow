General rules you should know while using Pymatflow
===================================================

High symmetry line used in Electronic Band and Phonon Band
----------------------------------------------------------
Research on solid state system usually revolves using the great tool of reciprocal
space. In Pymatflow, all the using of high symmetry point in reciprocal space goes
in a specific form.
The high symmetry k point path used in bands structure calculation or phonon band
structure calculation are always in format like this::

    [[kx, ky, kz, label, connect_indicator], ...] like [[0.0, 0.0, 0.0, 'GAMMA', 15], ...]

It is a list of high symmetry kpoint along with the connecting information.
Notes:
    * if ``connect_indicator`` in a kpoint is an integer, then it will connect to the following point through the number of kpoints defined by connect_indicator.
    * if connect_indicator in a kpoint is '|', then it will not connect to the following point.
    * [kx, ky, kz] are all in crystal format(fractional coordinates) in reciprocal space.

Besides, users can preparing the corresponding high symmmetry path in a file with
such a format::

    5
    0.0 0.0 0.0 #GAMMA 15
    x.x x.x x.x #XXX |
    x.x x.x x.x #XXX 10
    x.x x.x x.x #XXX 15
    x.x x.x x.x #XXX 20
When you provide such file to pymatflow(usually by ``--kpath-file`` or ``--qpath-file``),
Pymatflow will use it to bulid the corresponding data structure. Must be aware
that the ``x.x x.x x.x`` is in crystal coordinate(fractional coordinate)!

..
