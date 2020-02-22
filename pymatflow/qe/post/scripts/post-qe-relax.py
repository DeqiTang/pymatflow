#!/usr/bin/env python
# _*_ coding: utf-8 _*_

import os
import argparse
from pymatflow.qe.post.opt import opt_post


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-d", "--directory", help="directory of relax running", type=str, default="tmp-qe-relax")
    parser.add_argument("-f", "--file", help="output of opt running", type=str, default="relax.out")
    parser.add_argument("-o", "--output", help="converted to xyz file", type=str, default="relaxed.xyz")
    parser.add_argument("--view-traj", type=str, default="yes",
            choices=["yes", "no"],
            help="whether to view the trajectory.")
    parser.add_argument("--engine", type=str, default="bash",
            choices=["bash", "python"],
            help="post processing engine to choose")

    args = parser.parse_args()
   
    if args.engine == "python":
        os.chdir(args.directory)
        task = opt_post(output=args.file, run_type='relax')
        task.export()
        if args.view_traj == "yes":
            task.view_trajectory()
        os.chdir("../")
    elif args.engine == "bash":
        os.chdir(args.directory)
        # touch postprocessing directory
        os.system("mkdir -p post-processing")
        with open("post-processing/post.sh", 'w') as fout:
            fout.write("#!/bin/bash\n")
            fout.write("\n")
            fout.write("# ==============================\n")
            fout.write("# define and get some variables\n")
            fout.write("opt_out=../%s\n" % args.file)
            fout.write("nat=`cat ${opt_out} | grep 'number of atoms/cel' | cut -d \'=\' -f 2`\n")
            fout.write("scf_cycles=`cat ${opt_out} | grep 'ATOMIC_POSITIONS (angstrom)' | wc -l`\n")
            fout.write("# if_optimized: 0 -> no, 1 -> yes\n")
            fout.write("if_optimized=`cat ${opt_out} | grep 'Begin final coordinates' | wc -l`\n")
            fout.write("final_struct_file=optimized.xyz\n")
            #fout.write("input_xyz=../%s\n")
            fout.write("\n")
            fout.write("\n")
            
            # get the trajectory
            fout.write("# ==============================\n")
            fout.write("# get trajectory\n")
            fout.write("\n") 
            fout.write("output_trajfile=trajectory.xyz\n")
            fout.write("for i in `seq -w 1 1 ${scf_cycles}`\n")
            fout.write("do\n")
            fout.write("  struct_begin=`cat ${opt_out} | grep -n 'ATOMIC_POSITIONS (angstrom)' | head -n +${i} | tail -n -1 | cut -d \':\' -f 1`\n")
            fout.write("  cat >> ${output_trajfile} <<EOF\n")
            fout.write("`echo \"${nat}\" | bc`\n")
            fout.write("scf scycles: ${i}, ion step: `echo \"${i} - 1\" | bc`\n")
            fout.write("EOF\n")
            fout.write("  cat ${opt_out} | head -n +`echo \"${struct_begin} + ${nat}\" | bc` | tail -n -`echo \"${nat}\" | bc` >> ${output_trajfile}\n")
            fout.write("done\n")
            fout.write("# end get the trajectory\n")
            fout.write("\n\n")
            # get the final structure
            fout.write("# get the final structure\n")
            fout.write("  begin_final=`cat ${opt_out} | grep -n \"Begin final\" | cut -d ':' -f1`\n")
            fout.write("  end_final=`cat ${opt_out} | grep -n \"End final\" | cut -d ':' -f1`\n")
            fout.write("  cat > ${final_struct_file}<<EOF\n")
            fout.write("`echo \"${nat}\" | bc`\n")
            fout.write("final structure\n")
            fout.write("EOF\n")
            fout.write("  cat ${opt_out} | head -n +`echo \"${end_final} - 1\" | bc` | tail -n -`echo \"${nat}\" | bc` >> ${final_struct_file}\n")
            #


        os.system("cd post-processing; bash post.sh; cd ../")
        os.chdir("../")
