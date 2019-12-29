#!/usr/bin/env python
# _*_ coding: utf-8 _*_


def set_vdw(dft, usevdw=False):
    if usevdw == True:
        dft.xc.vdw_potential.status = True
    else:
        dft.xc.vdw_potential.status = False
