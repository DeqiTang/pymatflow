#!/usr/bin/env python
# _*_ coding: utf-8 _*_

# ==================================
# ==================================
class cp2k_dft_real_time_propagation:
    def __init__(self):
        self.params = {
                "ACCURACY_REFINEMENT": None,
                "APPLY_DELTA_PULSE": None,
                "ASPC_ORDER": None,
                "PERIODIC": None,
                "PROPAGATOR": None,
                "MAX_EXP": None,
                "MAX_ITER": None,
                "EPS_ITER": None,
                }
        self.status = False


