#!/usr/bin/env python3

import os

MAXIMUM_ENERGY = 3300 # cm^-1

#################################################### END OF INPUT

os.system("wget https://www.exomol.com/db/NH3/14N-1H3/CoYuTe/14N-1H3__CoYuTe.states.bz2")

Elow = 0
Ehigh = 100

while Ehigh <= MAXIMUM_ENERGY:
    os.system("wget -P coyute_trans/ https://www.exomol.com/db/NH3/14N-1H3/CoYuTe/14N-1H3__CoYuTe__" + "{:05d}".format(Elow)+ "-" + "{:05d}".format(Ehigh) + ".trans.bz2")
    Elow += 100
    Ehigh += 100