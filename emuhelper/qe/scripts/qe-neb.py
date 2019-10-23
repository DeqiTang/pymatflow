#!/usr/bin/env python
# _*_ coding: utf-8 _*_

import argparse

from emuhelper.qe.neb import neb_run

"""
usage:
    qe-neb.py -f xxx.xyz -k '2 2 2 0 0 0' --ecutwfc 100
"""


control_params = {}
system_params = {}
electrons_params = {}
path_params = {}

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-d", "--directory", help="directory of the calculation", type=str, default="tmp-qe-neb")
    parser.add_argument("--runopt", help="gen, run, or genrun", type=str, default="genrun")
    parser.add_argument("--image1", help="the first image xyz file", type=str)
    parser.add_argument("--image2", help="the intermediate image xyz file", type=str)
    parser.add_argument("--image3", help="the last image xyz file", type=str)
    parser.add_argument("--ecutwfc", help="ecutwfc, default value: 100 Ry", type=int, default=100)
    parser.add_argument("-k", "--kpoints", help="set kpoints like '1 1 1 0 0 0'", type=str, default="1 1 1 0 0 0")
    parser.add_argument("--occupations", help="occupation type", type=str, default="smearing")
    parser.add_argument("--smearing", help="smearing type", type=str, default="gaussian")
    parser.add_argument("--degauss", help="value of the gaussian spreading (Ry) for brillouin-zone integration in metals.", type=float, default=0.001)
    # params for neb namelist &path
    parser.add_argument("--string-method", help="string_method", type=str, default="neb")
    parser.add_argument("--nstep-path", help="nstep_path", type=int, default=100)
    parser.add_argument("--opt-scheme", help="Specify the type of optimization scheme(sd, broyden, broyden2, quick-min, langevin)", type=str, default="broyden")
    parser.add_argument("--num-of-images", help="number of total images(including the initial and final image)", type=int, default=5)
    parser.add_argument("--k-max", help="Set them to use a Variable Elastic Constants scheme elastic constants are in the range [ k_min, k_max  ], this is useful to rise the resolution around the saddle point", type=float, default=0.3e0)
    parser.add_argument("--k-min", help="Set them to use a Variable Elastic Constants scheme elastic constants are in the range [ k_min, k_max  ], this is useful to rise the resolution around the saddle point", type=float, default=0.2e0)
    parser.add_argument("--ci-scheme", help="Specify the type of Climbing Image scheme(no-CI, auto, manual)", type=str, default="auto")
    parser.add_argument("--path_thr", help="path_thr", type=float, default=0.05)
    parser.add_argument("--ds", help="Optimisation step length ( Hartree atomic units )", type=float, default=1.e0)
    
    # ==========================================================
    # transfer parameters from the arg parser to opt_run setting
    # ==========================================================   
    args = parser.parse_args()
    directory = args.directory
    image1 = args.image1
    image2 = args.image2
    image3 = args.image3
    system_params["ecutwfc"] = args.ecutwfc
    system_params["occupations"] = args.occupations
    system_params["smearing"] = args.smearing
    system_params["degauss"] = args.degauss
    kpoints_mp = [int(args.kpoints.split()[i]) for i in range(6)]
    path_params["string_method"] = args.string_method
    path_params["nstep_path"] = args.nstep_path
    path_params["opt_scheme"] = args.opt_scheme
    path_params["num_of_images"] = args.num_of_images
    path_params["k_max"] = args.k_max
    path_params["k_min"] = args.k_min
    path_params["CI_scheme"] = args.ci_scheme
    path_params["path_thr"] = args.path_thr
    path_params["ds"] = args.ds

    task = neb_run(image1, image2, image3)
    task.neb(directory=directory, runopt=args.runopt, control=control_params, system=system_params, electrons=electrons_params, kpoints_mp=kpoints_mp, path=path_params)
