#!/usr/bin/env python3

import numpy
import math
import matplotlib
import matplotlib.pyplot as plt

class coll_fit:
    def __init__(self, m, c, T, dK, dP):
        self.m = m
        self.c = c
        self.T = T
        self.dK = dK
        if dP > 1:
            self.dP = 1
        else:
            self.dP = dP
    def get_fit(self, dE):
        return self.m * dE + self.c

class level:
    def __init__(self):
        self.E = 0.0
        self.g = -666.6
        self.J = -1
        self.K = -1
        self.p = -1
        self.qunt_nums = ""
        self.idx = -1

class coll_transition:
    def __init__(self, iup, ilow):
        self.idx_up = iup
        self.idx_low = ilow
        self.Cul_T = []
        self.Clu_T = []

def _get_lamda_levels_coll_rates(input_lamda_file=None):
    k = 1.3806488e-16       #Boltzman constant, erg/K
    h = 6.62607e-27         #Planck constant, erg*s
    c = 299792458.e2        #speed of light, sm/s
    h_c_k = h * c / k

    fin = open(input_lamda_file,'r')
    for _ in range(0,6):
        line = fin.readline()
    number_of_levels = int(line.split()[0])
    fin.readline()

    levels = []
    for _ in range(0, number_of_levels):
        line = fin.readline()
        line = line.split()
        temp_level = level()
        temp_level.idx = int(line[0])
        temp_level.E = float(line[1])
        temp_level.g = float(line[2])
        line = line[3].split("_")
        temp_level.J = int(line[0])
        temp_level.K = int(line[1])
        temp_level.p = int(line[2])
        levels.append(temp_level)
    
    fin.readline()
    line = fin.readline()
    number_of_rad_transitions = int(line.split()[0])
    fin.readline()
    for _ in range(0,number_of_rad_transitions):
        line = fin.readline()
            
    coll_rates = []

    for _ in range(0,6):
        line = fin.readline()
    number_of_coll_trans = int(line.split()[0])

    T = []
    line = fin.readline()
    line = fin.readline()
    number_of_coll_temps = int(line.split()[0])
    
    line = fin.readline()
    line = fin.readline().split()
    for j in range(0, number_of_coll_temps):
        T.append( float(line[j]) )
    line = fin.readline()

    for _ in range(0, number_of_coll_trans):
        line = fin.readline()
        line = line.split()
        lev_up = int(line[1]) - 1
        lev_low = int(line[2]) - 1
        temp_trans = coll_transition(lev_up, lev_low)
        for j in range(0, number_of_coll_temps):
                temp_trans.Cul_T.append( float(line[j+3]) )
                temp_trans.Clu_T.append( float(line[j+3])*levels[lev_up].g/levels[lev_low].g * math.exp( (levels[lev_low].E-levels[lev_up].E)*h_c_k/T[j] ) )
        temp_trans.Cul_T = numpy.array( temp_trans.Cul_T )
        temp_trans.Clu_T = numpy.array( temp_trans.Clu_T )
        coll_rates.append(temp_trans)
    fin.close()

    return levels, coll_rates, numpy.array(T)

def _make_model_for_coll_vs_T(levs, cols, T, dK, dP, iT):
    dE = []
    C_at_T = []
    if dP > 1:
        dP = 1
    for i in cols:
        dpar = abs(levs[i.idx_up].p - levs[i.idx_low].p)
        if dpar > 1:
            dpar = 1
        if abs(levs[i.idx_up].K - levs[i.idx_low].K) == dK and dpar == dP:
            dE.append( abs(levs[i.idx_up].E - levs[i.idx_low].E) )
            C_at_T.append( i.Clu_T[iT] )
    dE = numpy.array(dE)
    C_at_T = numpy.log10(C_at_T)
    A = numpy.vstack([dE, numpy.ones(len(dE))]).T
    m, c = numpy.linalg.lstsq(A, C_at_T)[0]
    return coll_fit(m, c, T[iT], dK, dP)

def get_coll_fitted_models(filename="o-nh3_LAMDA.dat"):
    levs, cols, T = _get_lamda_levels_coll_rates(filename)
    models = []
    if levs[0].K == 0: #ortho-NH3
        max_dK = 6
        for j in range(0, len(T)):
            models.append( _make_model_for_coll_vs_T(levs, cols, T, 0, 0, j) )
            models.append( _make_model_for_coll_vs_T(levs, cols, T, 3, 0, j) )
            models.append( _make_model_for_coll_vs_T(levs, cols, T, 6, 0, j) )
            models.append( _make_model_for_coll_vs_T(levs, cols, T, 0, 1, j) )
            models.append( _make_model_for_coll_vs_T(levs, cols, T, 3, 1, j) )
            models.append( _make_model_for_coll_vs_T(levs, cols, T, 6, 1, j) )
    else: #para-NH3
        max_dK = 4
        for j in range(0, len(T)):
            models.append( _make_model_for_coll_vs_T(levs, cols, T, 0, 0, j) )
            models.append( _make_model_for_coll_vs_T(levs, cols, T, 1, 0, j) )
            models.append( _make_model_for_coll_vs_T(levs, cols, T, 2, 0, j) )
            models.append( _make_model_for_coll_vs_T(levs, cols, T, 3, 0, j) )
            models.append( _make_model_for_coll_vs_T(levs, cols, T, 4, 0, j) )
            models.append( _make_model_for_coll_vs_T(levs, cols, T, 0, 1, j) )
            models.append( _make_model_for_coll_vs_T(levs, cols, T, 1, 1, j) )
            models.append( _make_model_for_coll_vs_T(levs, cols, T, 2, 1, j) )
            models.append( _make_model_for_coll_vs_T(levs, cols, T, 3, 1, j) )
            models.append( _make_model_for_coll_vs_T(levs, cols, T, 4, 1, j) )
    return levs, cols, T, models, max_dK

if __name__ == "__main__":
    levs, cols, T = _get_lamda_levels_coll_rates("p-nh3_LAMDA.dat")
    dE = []
    C_at_T = []
    dP = 1
    iT = 4
    dK = 4
    for i in cols:
        dpar = abs(levs[i.idx_up].p - levs[i.idx_low].p)
        if dpar > 1:
            dpar = 1
        if abs(levs[i.idx_up].K - levs[i.idx_low].K) == dK and dpar == dP:
            dE.append( abs(levs[i.idx_up].E - levs[i.idx_low].E) )
            C_at_T.append( i.Clu_T[iT] )
    cmodel = _make_model_for_coll_vs_T(levs, cols, T, dK, dP, iT)
    dE = numpy.array(dE)
    C_at_T = numpy.log10(C_at_T)
    print(cmodel.m, cmodel.c)
    plt.plot(dE, C_at_T, 'bo')
    plt.plot(dE, cmodel.m*dE + cmodel.c, 'r', label='Fitted line')
    plt.show()