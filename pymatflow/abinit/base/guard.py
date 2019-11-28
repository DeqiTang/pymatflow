#!/usr/bin/env python
# _*_ coding: utf-8 _*_

import sys

class abinit_guard:
    """
    in the before, abinit_guard only cheeck the errors, but will not try to correct the
    error automatically.
    now, abinit_guard will automatically correct some errors for which there is only
    one solution, and if the number of solutions to some kind of erros is more that one
    it will not trying to correct automatically, but will stop the program.
    """
    def __init__(self, queen, electrons=None, ions=None, system=None, dfpt=None):
        """
        queen:
            type of abinit run, can be the following values:
            static, opt, md, neb, dfpt
            we will check the input argument according to value of queen.
            e.g. when queen == "opt", there are three arguments that can
            not be None, namely elctrons, ions, and system
        """
        self.queen = queen

        if queen == "static":
            if electrons == None or system == None:
                print("==============================================\n")
                print("                Warning !!!\n")
                print("==============================================\n")
                print("trouble while initialize abinit_guard:\n")
                print("if queen == static, you have to pass inelectrons\n")
                print("and system\n")
                sys.exit(1)
            self.electrons = electrons
            self.system = system

        if queen == "opt":
            if electrons == None or system == None or ions == None:
                print("==============================================\n")
                print("                Warning !!!\n")
                print("==============================================\n")
                print("trouble while initialize abinit_guard:\n")
                print("if queen == opt, you have to pass inelectrons\n")
                print("and system and ions\n")
                sys.exit(1)
            self.electrons = electrons
            self.system = system
            self.ions = ions

        if queen == "md":
            if electrons == None or system == None or ions == None:
                print("==============================================\n")
                print("                warning !!!\n")
                print("==============================================\n")
                print("trouble while initialize abinit_guard:\n")
                print("if queen == md, you have to pass inelectrons\n")
                print("and system and ions\n")
                sys.exit(1)
            self.electrons = electrons
            self.system = system
            self.ions = ions
        if queen == "neb":
            if electrons == None or system == None:
                print("==============================================\n")
                print("                warning !!!\n")
                print("==============================================\n")
                print("trouble while initialize abinit_guard:\n")
                print("if queen == neb, you have to pass inelectrons\n")
                print("and system\n")
                sys.exit(1)
            self.electrons = electrons
            self.system = system

    def check_all(self):
        """
        Note:
            we should not execute this function imediately after the construction
            of the object if we will do futher setting of parameters.
            it is recommended that you initialize the abinit_guard in __init__ of
            classes like static_run, and execute check_all when you finish setting
            in static_run, and decided to generate the input files.
        """
        if self.queen == "opt":
            self.check_optcell()

    def check_optcell(self):
        if "optcell" not in self.ions.params:
            return
        if self.ions.params["optcell"] == None or self.ions.params["optcell"] == 0:
            return
        # now optcell is non-zero, we have to make sure ecutsm is larger than 0
        # and the recommended value is 0.5 Ha.
        # see https://docs.abinit.org/variables/rlx/#ecutsm
        if "ecutsm" not in self.electrons.params or self.electrons.params["ecutsm"] == 0.0:
            print("==========================================================\n")
            print("                 Warning !!!\n")
            print("==========================================================\n")
            print("you set a non-zero optcell: %d\n" % self.ions.params["optcell"])
            print("so the value of ecutsm must be larger than 0 too!!!!\n")
            print("recommended value for ecutsm when optcell is non-zero is:\n")
            print("0.5 Ha\n")
            print("-----------------------------------------------------------\n")
            print("ecutsm is automatically set to 0.5 Ha\n")
            #sys.exit(1)
            self.electrons.params["ecutsm"] = 0.5
        #
