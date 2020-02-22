import os

#from pymatflow.remote.ssh import ssh
#from pymatflow.remote.rsync import rsync

def server_handle(auto, server, directory, jobfilebase):
    """
    :param auto:
        0 do nothing
        1: copying files to server
        2: copying and executing
        3: pymatflow run in server, direct submit the job
        in order use auto=1, 2, you must make sure there is a working ~/.pymatflow/server_[pbs|yh].conf"
    :param server:
        can be 'pbs' or 'yh' or 'lsf_sz'
    :param jobfilebase:
        base name of submitting job script,
        like static-nscf and the coresponding job submit script
        would be static-nscf.pbs if server is pbs and sstatic-nscf.yh
        if server is yh
    """
    # server handle
    if auto == 0:
        pass
    elif auto == 1:
        from pymatflow.remote.rsync import rsync
        mover = rsync()
        if server == "pbs":
            mover.get_info(os.path.join(os.path.expanduser("~"), ".pymatflow/server_pbs.conf"))
        elif args.server == "yh":
            mover.get_info(os.path.join(os.path.expanduser("~"), ".pymatflow/server_yh.conf"))
        mover.copy_default(source=os.path.abspath(directory))
    elif auto == 2:
        from pymatflow.remote.ssh import ssh
        from pymatflow.remote.rsync import rsync
        mover = rsync()
        if server == "pbs":
            mover.get_info(os.path.join(os.path.expanduser("~"), ".pymatflow/server_pbs.conf"))
        elif server == "yh":
            mover.get_info(os.path.join(os.path.expanduser("~"), ".pymatflow/server_yh.conf"))
        mover.copy_default(source=os.path.abspath(directory))
        ctl = ssh()
        if server == "pbs":
            ctl.get_info(os.path.join(os.path.expanduser('~'), ".pymatflow/server_pbs.conf"))
            ctl.login()
            ctl.submit(workdir=directory, jobfile=jobfilebase+".pbs", server="pbs")
        elif server == "yh":
            ctl.get_info(os.path.join(os.path.expanduser('~'), ".pymatflow/server_yh.conf"))
            ctl.login()
            ctl.submit(workdir=directory, jobfile=jobfilebase+".sub", server="yh")
    elif auto == 3:
        os.chdir(directory)
        if server == "pbs":
            os.system("qsub %s" % jobfilebase+".pbs")
        elif server == "yh":
            os.system("yhbatch -p free %s" % jobfilebase+".sub")
        elif server == "lsf_sz":
            os.system("chmod 755 %s; bsub %s" % (jobfilebase+".lsf_sz", jobfilebase+".lsf_sz"))
        os.chdir("../")
