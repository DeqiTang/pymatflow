"""
Molecular Dynamics calc
"""
import os
import re
import shutil


from pymatflow.remote.server import server_handle
from pymatflow.qe.pwscf import pwscf

class md_run(pwscf):
    """
    """
    def __init__(self):
        super().__init__()


    def md(self, directory="tmp-qe-md", inpname="md.in", output="md.out", runopt="gen", auto=0):
        """
        :param directory: a place for all the generated files
        """
        self.set_md()
        if runopt ==  "gen" or runopt == "genrun":
            if os.path.exists(directory):
                shutil.rmtree(directory)
            os.mkdir(directory)

            #os.system("cp *.UPF %s/" % directory)
            #os.system("cp %s %s/" % (self.arts.xyz.file, directory))

            # do not copy too many files at the same time or it will be slow
            # so we do not copy all UPF files in the directory but just copy
            # those used in the calculation.
            shutil.copyfile(self.arts.xyz.file, os.path.join(directory, os.path.basename(self.arts.xyz.file)))
            #all_upfs = [s for s in os.listdir() if s.split(".")[-1] == "UPF"]
            all_file = os.listdir()
            for element in self.arts.xyz.specie_labels:
                for item in all_file:
                    #if re.match("(%s)(.*)(upf)" % element, item, re.IGNORECASE):
                    if item.split(".")[0].lower() == element.lower() or item.split("_")[0].lower() == element.lower():
                        shutil.copyfile(item, os.path.join(directory, item))
                        break
            self.arts.pseudo.dir = os.path.abspath(directory)
            self.control.set_params({"pseudo_dir": os.path.abspath(directory)})
            #

            with open(os.path.join(directory, inpname), 'w') as fout:
                self.control.to_in(fout)
                self.system.to_in(fout)
                self.electrons.to_in(fout)
                self.ions.to_in(fout)
                self.arts.to_in(fout)
            # gen yhbatch script
            self.gen_llhpc(directory=directory, inpname=inpname, output=output, cmd="$PMF_PWX")
            # gen pbs script
            self.gen_pbs(directory=directory, inpname=inpname, output=output, cmd="$PMF_PWX", jobname=self.run_params["jobname"], nodes=self.run_params["nodes"], ppn=self.run_params["ppn"], queue=self.run_params["queue"])

        if runopt == "run" or runopt == "genrun":
            os.chdir(directory)
            os.system("%s $PMF_PWX < %s | tee %s" % (self.run_params["mpi"], inpname, output))
            os.chdir("../")
        server_handle(auto=auto, directory=directory, jobfilebase="md", server=self.run_params["server"])

    def vc_md(self, directory="tmp-qe-vc-md", inpname="vc-md.in", output="vc-md.out", runopt="gen", auto=0):
        """
        :param directory: a place for all the generated files
        """
        self.set_vc_md()
        if runopt ==  "gen" or runopt == "genrun":
            if os.path.exists(directory):
                shutil.rmtree(directory)
            os.mkdir(directory)

            #os.system("cp *.UPF %s/" % directory)
            #os.system("cp %s %s/" % (self.arts.xyz.file, directory))

            # do not copy too many files at the same time or it will be slow
            # so we do not copy all UPF files in the directory but just copy
            # those used in the calculation.
            shutil.copyfile(self.arts.xyz.file, os.path.join(directory, self.arts.xyz.file))
            all_upfs = [s for s in os.listdir() if s.split(".")[-1] == "UPF"]
            for element in self.arts.xyz.specie_labels:
                for upf in all_upfs:
                    #if upf.split(".")[0] == element:
                    if upf.split(".")[0].lower() == element.lower() or upf.split("_")[0].lower() == element.lower():
                        shutil.copyfile(upf, os.path.join(directory, upf))
                        break
            #
            self.arts.pseudo.dir = os.path.abspath(directory)
            self.control.pseudo_dir = os.path.abspath(directory)

            with open(os.path.join(directory, inpname), 'w') as fout:
                self.control.to_in(fout)
                self.system.to_in(fout)
                self.electrons.to_in(fout)
                self.ions.to_in(fout)
                self.cell.to_in(fout)
                self.arts.to_in(fout)

            # gen yhbatch script
            self.gen_yh(directory=directory, inpname=inpname, output=output, cmd="$PMF_PWX")
            # gen pbs script
            self.gen_pbs(directory=directory, inpname=inpname, output=output, cmd="$PMF_PWX", jobname=self.run_params["jobname"], nodes=self.run_params["nodes"], ppn=self.run_params["ppn"], queue=self.run_params["queue"])

        if runopt == "run" or runopt == "genrun":
            os.chdir(directory)
            os.system("%s $PMF_PWX < %s | tee %s" % (self.run_params["mpi"], inpname, output))
            os.chdir("../")
        server_handle(auto=auto, directory=directory, jobfilebase="vc-md", server=self.run_params["server"])

    def set_md(self):
        self.control.calculation('md')
        self.control.basic_setting("md")

        self.system.basic_setting(self.arts)
        self.electrons.basic_setting()
        self.ions.basic_setting('md')

    def set_vc_md(self):
        self.control.calculation('vc-md')
        self.control.basic_setting("vc-md")

        self.system.basic_setting(self.arts)
        self.electrons.basic_setting()
        self.ions.basic_setting('vc-md')
    #
