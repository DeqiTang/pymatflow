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
        can be 'pbs' or 'llhpc' or 'lsf_sz' or 'lsf_sustc
    :param jobfilebase:
        base name of submitting job script,
        like static-nscf and the coresponding job submit script
        would be static-nscf.pbs if server is pbs and sstatic-nscf.yh
        if server is llhpc
    """
    # server handle
    if auto == 0:
        pass
    elif auto == 1:
        from pymatflow.remote.rsync import Rsync
        mover = Rsync()
        if server == "pbs":
            mover.get_info(os.path.join(os.path.expanduser("~"), ".pymatflow/server_pbs.conf"))
        elif args.server == "llhpc":
            mover.get_info(os.path.join(os.path.expanduser("~"), ".pymatflow/server_llhpc.conf"))
        mover.copy_default(source=os.path.abspath(directory))
    elif auto == 2:
        from pymatflow.remote.ssh import Ssh
        from pymatflow.remote.rsync import Rsync
        mover = Rsync()
        if server == "pbs":
            mover.get_info(os.path.join(os.path.expanduser("~"), ".pymatflow/server_pbs.conf"))
        elif server == "llhpc":
            mover.get_info(os.path.join(os.path.expanduser("~"), ".pymatflow/server_llhpc.conf"))
        mover.copy_default(source=os.path.abspath(directory))
        ctl = Ssh()
        if server == "pbs":
            ctl.get_info(os.path.join(os.path.expanduser('~'), ".pymatflow/server_pbs.conf"))
            ctl.login()
            ctl.submit(workdir=directory, jobfile=jobfilebase+".pbs", server="pbs")
        elif server == "llhpc":
            ctl.get_info(os.path.join(os.path.expanduser('~'), ".pymatflow/server_llhpc.conf"))
            ctl.login()
            ctl.submit(workdir=directory, jobfile=jobfilebase+".slurm", server="llhpc")
    elif auto == 3:
        os.chdir(directory)
        if server == "pbs":
            os.system("qsub %s" % jobfilebase+".pbs")
        elif server == "llhpc":
            os.system("yhbatch %s" % jobfilebase+".slurm")
        elif server == "lsf_sz":
            os.system("chmod 755 %s; bsub %s" % (jobfilebase+".lsf_sz", jobfilebase+".lsf_sz"))
        elif server == "lsf_sustc":
            os.system("chmod 755 %s; bsub < %s" % (jobfilebase+".lsf_sustc", jobfilebase+".lsf_sustc"))            
        elif server == "cdcloud":
            os.system("sbatch < %s" % jobfilebase+".slurm_cd")
        os.chdir("../")
