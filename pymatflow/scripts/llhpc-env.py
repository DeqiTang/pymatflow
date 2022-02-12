#!/bin/bash

import os


def main():
    with open(os.path.join(os.path.expanduser("~"), ".bashrc"), 'a') as fout:
        fout.write("# needed by pymatflow on TianHe II\n")
        fout.write("export LD_LIBRARY_PATH=/THFS/opt/intel/composer_xe_2013_sp1.3.174/mkl/lib/intel64:/THFS/opt/intel/composer_xe_2013_sp1.3.174/compiler/lib/intel64:$LD_LIBRARY_PATH")
        fout.write("export PMF_VASP_STD=/THFS/opt/vasp/5.4.4_neb/vasp.5.4.4/bin/vasp_std\n")
        fout.write("export PMF_VASP_GAM=/THFS/opt/vasp/5.4.4_neb/vasp.5.4.4/bin/vasp_gam\n")
        fout.write("export PMF_VASP_NCL=/THFS/opt/vasp/5.4.4_neb/vasp.5.4.4/bin/vasp_ncl\n")
        #fout.write("export PMF_ABINIT=")
        fout.write("export PMF_CP2K=/THFS/opt/cp2k/4.1/cp2k-4.1/exe/Linux-x86-64-intel-host/cp2k.popt\n")
        fout.write("export PMF_PWX=/THFS/opt/qe/qe-6.4/bin/pw.x\n")
        fout.write("export PMF_PROJWFCX=/THFS/opt/qe/qe-6.4/bin/projwfc.x\n")
        fout.write("export PMF_BANDSX=/THFS/opt/qe/qe-6.4/bin/bands.x\n")
        fout.write("export PMF_MOLECULARPDOSX=/THFS/opt/qe/qe-6.4/bin/molecularpdos.x\n")
        fout.write("export PMF_PPX=/THFS/opt/qe/qe-6.4/bin/pp.x\n")
        fout.write("export PMF_NEBX=/THFS/opt/qe/qe-6.4/bin/neb.x\n")
        fout.write("export PMF_PHX=/THFS/opt/qe/qe-6.4/bin/ph.x\n")
        fout.write("export PMF_Q2RX=/THFS/opt/qe/qe-6.4/bin/q2r.x\n")
        fout.write("export PMF_MATDYNX=/THFS/opt/qe/qe-6.4/bin/matdyn.x\n")
        fout.write("export PMF_PLOTBANDX=/THFS/opt/qe/qe-6.4/bin/plotband.x\n")
        fout.write("export PMF_DYNAMATX=/THFS/opt/qe/qe-6.4/bin/dynmat.x\n")
        fout.write("export PMF_SIESTA=/THFS/opt/siesta/siesta-4.1-b4/Obj/siesta\n")
    print("====================================================================\n")
    print("已经自动为你设置好了环境变量, 但是需要你执行一下命令来使其立即生效:\n")
    print("source ~/.bashrc\n")
    print("====================================================================\n")


if __name__ == "__main__":
    main()
