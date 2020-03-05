#!/bin/bash

import os


def main():
    with open(os.path.join(os.path.expanduser("~"), ".bashrc"), 'a') as fout:
        fout.write("# needed by pymatflow on TianHe II\n")
        fout.write("export LD_LIBRARY_PATH=/THFS/opt/intel/composer_xe_2013_sp1.3.174/mkl/lib/intel64:/THFS/opt/intel/composer_xe_2013_sp1.3.174/compiler/lib/intel64:$LD_LIBRARY_PATH")
        fout.write("export PMF_VASP_STD=/THFS/opt/vasp/5.4.4_neb/vasp.5.4.4/bin/vasp_std\n")
        fout.write("export PMF_VASP_GAM=/THFS/opt/vasp/5.4.4_neb/vasp.5.4.4/bin/vasp_gam\n")
        fout.write("export PMF_VASP_NCL=/THFS/opt/vasp/5.4.4_neb/vasp.5.4.4/bin/vasp_ncl\n")
    print("====================================================================\n")
    print("已经自动为你设置好了环境变量, 但是需要你执行一下命令来使其立即生效:\n")
    print("source ~/.bashrc\n")
    print("====================================================================\n")


if __name__ == "__main__":
    main()
