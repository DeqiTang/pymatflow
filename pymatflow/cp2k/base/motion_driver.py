#!/usr/bin/env python
# _*_ coding: utf-8 _*_

# ========================
# CP2K / MOTION / DRIVER
# ========================
class cp2k_motion_driver:
    def __init__(self):
        self.params = {

                }
        self.status = False

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 2:
                self.params[item.split("-")[-1]] = params[item]
            else:
                pass
