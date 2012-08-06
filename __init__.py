import sage.misc.preparser, os
from sage.all import *

_bord_prog_files = ["ChainCx.sage", "PMCAlg.sage", "TypeDStr.sage", "TypeDDStr.sage", "Plats.sage"]

for _bord_prog_file in _bord_prog_files:
    sage.misc.preparser.load(__path__[0] + os.sep + _bord_prog_file, globals(), attach=False)
