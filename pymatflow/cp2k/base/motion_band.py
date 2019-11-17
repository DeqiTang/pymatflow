#!/usr/bin/env python
# _*_ coding: utf-8 _*_

from pymatflow.base.xyz import base_xyz

# ====================
# CP2K / MOTION / BAND
# ====================
class cp2k_motion_band:
    def __init__(self):
        self.params = {
                "ALIGN_FRAMES": None,
                "BAND_TYPE": None, # CI-NEB, IT-NEB, SM
                "K_SPRING": None,
                "NPROC_REP": None,
                "NUMBER_OF_REPLICA": None,
                "POT_TYPE": None,
                "PROC_DIST_TYPE": None,
                "ROTATE_FRAMES": None,
                "USE_COLVARS": None,
                }
        self.images = []

        self.params["BAND_TYPE"] = "CI-NEB"
        self.params["NUMBER_OF_REPLICA"] = 10
        self.params["K_SPRING"] = 0.05
        self.params["ROTATE_FRAMES"] = "TRUE"
        self.params["ALIGN_FRAMES"] = "TRUE"

    def to_motion(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t&BAND\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t%s %s\n" % (item, str(self.params[item])))

        for image in self.images:
            fout.write("\t\t&REPLICA\n")
            fout.write("\t\t\tCOORD_FILE_NAME %s\n" % image.file)
            fout.write("\t\t&END REPLICA\n")
        fout.write("\t&END BAND\n")

    def get_images(self, images):
        for image in images:
            self.images.append(base_xyz(image))

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 2:
                self.params[item.split("-")[-1]] = params[item]


