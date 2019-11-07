#!/usr/bin/env python
# _*_ coding: utf-8 _*_

class cp2k_motion_tmc:
    def __init__(self):
        pass

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 2:
                self.params[item.split("-")[-1]] = params[item]


