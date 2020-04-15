The usage with --auto = 0
==========================

About the th*.py scripts
-------------------------

in order to use th*.py scripts and sz-*.py scripts and the function of --auto=2 in running scripts(like qe-relax.py), you have to prepare a server_yh.conf file and server_pbs.conf in ~/.pymatflow/, in format like this::

    [server]
    ip = x.x.x.x
    user = xxx
    serverdir = /xx/xx/xx # like /THFS/home/xxx/
    password = xxx

..
