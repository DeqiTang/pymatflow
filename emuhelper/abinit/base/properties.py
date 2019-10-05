#!/usr/bin/env python
# _*_ coding: utf-8 _*_


class properties:
    """
    """
    def __init__(self):
        self.params = {
                }
    def to_in(self, fout):
        # fout: a file stream for writing
        for item in self.params:
            if self.params[item] is not None:
                fout.write("%s %s\n" % (item ,str(self.params[item])))
        #
