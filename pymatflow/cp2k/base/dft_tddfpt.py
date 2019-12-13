#!/usr/bin/env python
# _*_ coding: utf-8 _*_

# ====================
# ====================
class cp2k_dft_tddfpt:
    def __init__(self):
        self.params = {
                "CONVERGENCE": None,
                "DIAG_METHOD": None,
                "INVERT_S": None,
                "KERNEL": None,
                "LSD_SINGLETS": None,
                "MAX_KV": None,
                "NEV": None,
                "NLUMO": None,
                "NREORTHO": None,
                "OE_CORR": None,
                "PRECONDITIONER": None,
                "RESTARTS": None,
                "RES_ETYPE": None,
                }       
        self.status = False


