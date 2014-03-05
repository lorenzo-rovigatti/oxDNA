#!/usr/bin/env python

import base
import generators as gen
import numpy as np

s = gen.StrandGenerator()
str = s.generate(4, double=False, dir=np.array([0,1,0]), start_pos=np.array([0,0,1]))
syst = base.System([3, 3, 3])
syst.add_strand(str)
syst.print_lorenzo_output ("prova.conf", "prova.top")
