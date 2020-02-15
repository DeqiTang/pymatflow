import sys

class abinit_guard:
    """
    in the before, abinit_guard only cheeck the errors, but will not try to correct the
    error automatically.
    now, abinit_guard will automatically correct some errors for which there is only
    one solution, and if the number of solutions to some kind of erros is more that one
    it will not trying to correct automatically, but will stop the program.
    """
    def __init__(self):
        """
        """
        self.queen = None

    def set_queen(self, queen=None):
        self.queen = queen

    def check_all(self, electrons=None, ions=None, system=None, dfpt=None):
        """
        Note:
            we should not execute this function imediately after the construction
            of the object if we will do futher setting of parameters.
            it is recommended that you initialize the abinit_guard in __init__ of
            classes like static_run, and execute check_all when you finish setting
            in static_run, and decided to generate the input files.
        self.queen:
            type of abinit run, can be the following values:
            static, opt, md, neb, dfpt
            we will check the input argument according to value of self.queen.
            e.g. when queen == "opt", there are three arguments that can
            not be None, namely elctrons, ions, and system
        """
        if self.queen == None:
            print("=============================================================\n")
            print("                     WARNING!!!\n")
            print("-------------------------------------------------------------\n")
            print("abinit.base.guard.abinit_guard.check_all():\n")
            print("before execute check_all(), you should set the queen to some\n")
            print("value among -> opt, md, dfpt\n")
            sys.exit(1)
        #
        self.electrons = electrons
        self.ions = ions
        self.system = system
        self.dfpt = dfpt
        #
        if self.queen == "opt":
            self.check_optcell()
            self.check_ionmov()

        if self.queen == "md":
            self.check_optcell()
            self.check_ionmov()

        if self.queen == "dfpt":
            self.check_dfpt()

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
        if self.ions.params["optcell"] == 2:
            # then ionmov must be equal to one of the following: 2 3 13 15 22 25
            if self.ions.params["ionmov"] is not None and self.ions.params["ionmov"] not in [2, 3, 13, 15, 22, 25]:
                print("===========================================================\n")
                print("                        Warning !!!\n")
                print("===========================================================\n")
                print("the value of optcell is 2:\n")
                print("so the value of ionmov can only be in:\n")
                print("[2, 3, 13, 15, 22, 25]\n")
                print("however your ionmov is %d\n" % self.ions.params["ionmov"])
                print("-----------------------------------------------------------\n")
                sys.exit(1)

    def check_ionmov(self):
        if "ionmov" not in self.ions.params:
            return
        #
        if self.ions.params["ionmov"] == 13:
            # then nnos must be larger or equal to 1
            if "nnos" not in self.ions.params or self.ions.params["nnos"] == None or self.ions.params["nnos"] == 0:
                print("===========================================================\n")
                print("                        Warning !!!\n")
                print("===========================================================\n")
                print("the value of ionmov is 13:\n")
                print("so the value of nnos must be larger or equal to 1\n")
                print("so we automatically set nnos to 1 for you\n")
                print("-----------------------------------------------------------\n")
                #sys.exit(1)
                self.ions.params["nnos"] = 1
        #

    def check_dfpt(self):
        pass
