# Usage:
# python -c "from pymatflow.test.test import *; unittest.main()"
#

from .unit.test_cpp_cp2k import *
#from .unit.test_cpp_cube_handle import *
#from .unit.test_fortran_cube_ciff_1d import *


if __name__ == "__main__":
    unittest.main(verbosity=1)