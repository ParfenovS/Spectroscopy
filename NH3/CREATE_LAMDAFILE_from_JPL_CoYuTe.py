#!/usr/bin/env python3

import numpy
import math
import copy
import os
from dataclasses import dataclass
from Fit_C_vs_dE import get_coll_fitted_models

###########################################


# partition function at given T
Q_rs = 588.7816 # from JPL database
# temperature in [K] for which Q_rs is given
T = 300.0

# uncomment depending on which NH3 species is needed
molecule = "o-NH3"
#molecule = "p-NH3"

# maximum allowed values of quantum numbers
MAX_J = 15 # rotational quantum number J
# vibrational quantum numbers, see 14N-1H3__CoYuTe.readme
MAX_v1 = 0 # Symmetric stretch
MAX_v2 = 1 # Symmetric bend
MAX_v3 = 0 # Asymmetric stretch
MAX_v4 = 0 # Asymmetric bend

########################################### END OF INPUT

molecular_weight = "17."

input_file_jpl = "c017002.cat" # file with radiative transitions data from JPL database https://spec.jpl.nasa.gov/ftp/pub/catalog/c017002.cat
coyute_states_file = "14N-1H3__CoYuTe.states"
directory_with_coyute_transitions = "coyute_trans"

input_lamda_file = "p-nh3_LAMDA.dat" # file from LAMDA database https://home.strw.leidenuniv.nl/~moldata/NH3.html
output_file = "P-nh3.dat" # name of file with final results

if molecule == "o-NH3":
    input_lamda_file = "o-nh3_LAMDA.dat"
    output_file = "O-nh3.dat"

########################################### END OF ADDITIONAL SETTINGS

def check_symmetry(J, k):
    it_is_ortho = check_symmetry1(J, k)
    if molecule == "o-NH3": return it_is_ortho # ortho-NH3
    else: return not(it_is_ortho) # para-NH3

def check_symmetry1(J, k):
    it_is_ortho = False
    for n in range(0, J+2):
        if k == 3*n: # ortho-NH3
            it_is_ortho = True
            break
    return it_is_ortho

k = 1.3806488e-16       #Boltzman constant, erg/K
h = 6.62607e-27         #Planck constant, erg*s
c = 299792458.e2        #speed of light, sm/s
h_c_k = h * c / k

@dataclass
class level:
    E: float = -1.0
    g: int = -666
    J: int = -1
    K: int = -1
    p: int = -1
    v1: int = 0
    v2: int = 0
    v3: int = 0
    v4: int = 0
    gamma: int = -1
    idx: int = -1
    coyute_id: int = -1
    from_jpl: bool = False
    qunt_nums: str = ""
    def __eq__(self, other):
        equality = ( self.J == other.J and self.K == other.K and self.p == other.p and self.v1 == other.v1 and self.v2 == other.v2 and self.v3 == other.v3 and self.v4 == other.v4 )
        if self.g > 0 and other.g > 0:
            equality = equality and (self.g == other.g)
        '''if self.gamma > 0 and other.gamma > 0:
            equality = equality and (self.gamma == other.gamma)'''
        return equality

class transition:
    def __init__(self):
        self.freq = 0.0
        self.lgint = 0.0
        self.lev_up = level()
        self.lev_low = level()
        self.A = 0.0

def read_jpl_file(filename=None):
    transitions = []
    with open(filename,'r') as fin:
        for line in fin:
            trans = transition()
            trans.lev_low.from_jpl = True
            trans.lev_up.from_jpl = True
            trans.freq = float( line[0:13] )
            err = line[13:21]
            trans.lgint = float( line[21:29] )
            dr = line[29:31]
            trans.lev_low.E = float( line[31:41] )
            trans.lev_up.g = int( line[41:44] )
            tag = line[44:51]
            qnfmt = line[51:55]
            trans.lev_up.qunt_nums = line[55:66]
            trans.lev_up.K = int( (trans.lev_up.qunt_nums[2:4].replace("a", "1")).replace("b", "2").replace("c", "3").replace("-", " ") )
            trans.lev_up.qunt_nums = trans.lev_up.qunt_nums[0:2] + " " + (trans.lev_up.qunt_nums[2:4].replace("a", "1")).replace("b", "2").replace("c", "3").replace("-", " ") + " " + trans.lev_up.qunt_nums[4:6]
            trans.lev_up.p = int(trans.lev_up.qunt_nums[6:8])
            trans.lev_up.J = int(trans.lev_up.qunt_nums[0:2])
            trans.lev_low.qunt_nums = line[67:78]
            trans.lev_low.K = int( (trans.lev_low.qunt_nums[2:4].replace("a", "1")).replace("b", "2").replace("c", "3").replace("-", " ") )
            trans.lev_low.qunt_nums = trans.lev_low.qunt_nums[0:2] + " " + (trans.lev_low.qunt_nums[2:4].replace("a", "1")).replace("b", "2").replace("c", "3").replace("-", " ") + " " + trans.lev_low.qunt_nums[4:6]
            trans.lev_low.p = int(trans.lev_low.qunt_nums[6:8])
            trans.lev_low.J = int(trans.lev_low.qunt_nums[0:2])
            trans.A = 10.**trans.lgint * trans.freq**2. * Q_rs / trans.lev_up.g / ( math.exp(-trans.lev_low.E*h_c_k/T) * (1.-math.exp(-trans.freq*h*1.e6/k/T)) ) * 2.7964e-16
            trans.freq = trans.freq / 1000. #MHz -> GHz
            if trans.lev_up.J <= MAX_J and trans.lev_low.J <= MAX_J and check_symmetry(trans.lev_up.J, trans.lev_up.K) and check_symmetry(trans.lev_low.J, trans.lev_low.K):
                transitions.append( trans )
        for i in range(len(transitions)):
            for j in range(len(transitions)):
                if i != j and transitions[i].lev_up == transitions[j].lev_low:
                    transitions[i].lev_up.E = transitions[j].lev_low.E
                    break
            if transitions[i].lev_up.E == -1:
                transitions[i].lev_up.E = transitions[i].lev_low.E + transitions[i].freq * 1.e9 / c
    return transitions

# CoYuTe levels: awk '{if($2<200) print $2,$15,$16,$14}' 14N-1H3__CoYuTe.states | sort -n -k 1
def read_CoYuTe_states(filename):
    levels = []
    with open(filename, 'r') as fin:
        for line in fin:
            line = line.split()
            temp_lev = level()
            temp_lev.J = int(line[16-1])
            temp_lev.K = int(line[17-1])
            temp_lev.p = int(line[15-1])
            temp_lev.E = float(line[2-1])
            temp_lev.g = int(int(line[3-1]) / 3)
            temp_lev.coyute_id = int(line[1-1])
            v1 = int(line[9-1])
            temp_lev.v1 = v1
            v2 = int(line[10-1])
            temp_lev.v2 = v2
            v3 = int(line[11-1])
            temp_lev.v3 = v3
            v4 = int(line[12-1])
            temp_lev.v4 = v4
            #v5 = int(line[13-1])
            #v6 = int(line[14-1])
            gam = int(line[25-1])
            temp_lev.gamma = gam
            if temp_lev.J <= MAX_J and v1 <= MAX_v1 and v2 <= MAX_v2 and v3 <= MAX_v3 and v4 <= MAX_v4 and check_symmetry(temp_lev.J, temp_lev.K) and temp_lev.g > 0:
                levels.append(temp_lev)
    return levels


transitions = read_jpl_file(input_file_jpl)

levels = []
level_idx = 1
transitions[0].lev_low.idx = level_idx
levels.append( transitions[0].lev_low )
for trans in transitions:
    there_is_level = False
    for lev in levels:
        if (trans.lev_low == lev):
            there_is_level = True
            break
    if not(there_is_level):
        levels.append( trans.lev_low )
for trans in transitions:
    there_is_level = False
    for lev in levels:
        if (trans.lev_up == lev):
            there_is_level = True
            break
    if not(there_is_level):
        levels.append( trans.lev_up )

for lev in levels:
    for trans in transitions:
        if lev == trans.lev_up:
            lev.g = trans.lev_up.g
        if lev == trans.lev_low:
            lev.E = trans.lev_low.E


levels_coyute = read_CoYuTe_states(coyute_states_file)

levels_jpl = copy.deepcopy(levels)
for j in range(len(levels_coyute)):
    need_to_add_coyute_level = True
    for i in range(len(levels_jpl)):
        if levels_jpl[i] == levels_coyute[j]:
            need_to_add_coyute_level = False
            levels[i].coyute_id = levels_coyute[j].coyute_id
            levels[i].g = int(levels_coyute[j].g)
            break
    if need_to_add_coyute_level:
        levels.append(levels_coyute[j])

energies = []
for lev in levels:
    energies.append(lev.E)
E_sort_idx = numpy.argsort(energies)
ener_max = energies[len(E_sort_idx)-1]

levels0 = copy.deepcopy(levels)
levels = []
for i in range(len(E_sort_idx)):
    levels0[E_sort_idx[i]].idx = i + 1
    levels.append( levels0[E_sort_idx[i]] )

levels = numpy.array(levels)

for coyute_trans_file in os.listdir(directory_with_coyute_transitions):
    elow = int(coyute_trans_file[:-6].split("__")[2].split("-")[0])
    if (elow > ener_max): continue
    with open(directory_with_coyute_transitions + "/" + coyute_trans_file, 'r') as fin:
        for line in fin:
            line = line.split()
            id1 = int(line[0])
            id2 = int(line[1])
            A = float(line[2])
            first_lev_found = False
            second_lev_found = False
            trans = transition()
            for lev in levels:
                if lev.coyute_id == id1 and not first_lev_found:
                    trans.lev_up = copy.deepcopy(lev)
                    first_lev_found = True
                if lev.coyute_id == id2 and not second_lev_found:
                    trans.lev_low = copy.deepcopy(lev)
                    second_lev_found = True
                if first_lev_found and second_lev_found: break
            if first_lev_found and second_lev_found:
                if not trans.lev_low.from_jpl or not trans.lev_up.from_jpl:
                    trans.A = A
                    trans.freq = (trans.lev_up.E - trans.lev_low.E) * c * 1.e-9
                    transitions.append(trans)

freqs = []
for trans in transitions:
    freqs.append(trans.freq)
freq_sort_idx = numpy.argsort(freqs)
transitions0 = copy.deepcopy(transitions)
transitions = []
for i in range(len(freq_sort_idx)):
    transitions.append( transitions0[freq_sort_idx[i]] )

fout = open(output_file,'w')
fout.write('!MOLECULE\n'+molecule+'\n!MOLECULAR WEIGHT\n'+molecular_weight+'\n!NUMBER OF ENERGY LEVELS\n')
fout.write(str(len(levels)) + '\n')
fout.write('! LEVEL + ENERGY(CM-1) + WEIGHT + QUANTUM NOS.  J,K,inv-sym,v1,v2,v3,v4\n')
for i in range(0,len(levels)):
    fout.write(str(levels[i].idx)+'  '+str(levels[i].E)+'  '+str(levels[i].g)+' '+str(levels[i].J)+' '+str(levels[i].K)+' '+str(levels[i].p)+' '+str(levels[i].v1)+' '+str(levels[i].v2)+' '+str(levels[i].v3)+' '+str(levels[i].v4)+'\n')

fout.write('!NUMBER OF RADIATIVE TRANSITIONS\n')
fout.write( str(len(transitions)) )
fout.write('\n!TRANS + UP + LOW + EINSTEINA(s^-1) + FREQ(GHz) + E_low(K)\n')
i1 = 1
for trans in transitions:
    up_lev_found = False
    low_lev_found = False
    id_up = -1
    id_low = -1
    for lev in levels:
        if lev == trans.lev_up and not up_lev_found:
            up_lev_found = True
            id_up = lev.idx
        if lev == trans.lev_low and not low_lev_found:
            low_lev_found = True
            id_low = lev.idx
    if up_lev_found and low_lev_found:
        fout.write( str(i1)+"  "+str(id_up)+"  "+str(id_low)+"  " )
        fout.write( str(trans.A)+"  "+str(trans.freq)+" "+str(trans.lev_low.E*h*c/k)+'\n' )
        if id_up <= id_low: print("problematic rad. trans: " + str(i1) + " " + str(id_up) + " " + str(id_low))
        i1 = i1 + 1

lev_lamda, coll_lamda, T, coll_models, max_dK = get_coll_fitted_models(input_lamda_file)

fout.write("!NUMBER OF COLL PARTNERS\n1\n!COLLISIONS BETWEEN\n2                 "+molecule+" coll with p-H2\n")
fout.write("!NUMBER OF COLL TRANS\n")
fout.write(str( int(len(levels)*(len(levels)-1)/2)) )
#fout.write("3")
fout.write("\n!NUMBER OF COLL TEMPS\n")
fout.write( str(len(T)) + "\n")
fout.write("!COLL TEMPS\n")
for i in range(0, len(T)):
    fout.write( " " + str(T[i]) )
fout.write("\n !TRANS + UP + LOW + COLLRATES(cm^3 s^-1)\n")

trans_idx = 1
for i in range(2, len(lev_lamda)):
    for j in range(1, i):
        fout.write(str(trans_idx) + " " + str(i) + " " + str(j))
        for it in range(len(T)):
            fout.write( ' {:.1E}'.format(coll_lamda[trans_idx-1].Cul_T[it]) )
        fout.write("\n")
        trans_idx = trans_idx + 1
        if trans_idx > len(coll_lamda):
            break
    if trans_idx > len(coll_lamda):
            break
i1 = i + 1

for i in range(i1, len(levels)+1):
    for j in range(1, i):
        fout.write(str(trans_idx) + " " + str(i) + " " + str(j))
        dK = abs( levels[i-1].K - levels[j-1].K )
        dE = levels[i-1].E - levels[j-1].E
        dP = abs(levels[i-1].p - levels[j-1].p)
        inter_vibrational_trans = False
        if abs(levels[i-1].v2 - levels[j-1].v2) > 0 or abs(levels[i-1].v4 - levels[j-1].v4) > 0: inter_vibrational_trans = True
        if dP > 1:
            dP = 1
        gg = levels[j-1].g / levels[i-1].g
        for it in range(len(T)):
            rate = 0.0
            if dK <= max_dK:
                for im in coll_models:
                    if im.dK == dK and im.T == T[it] and im.dP == dP:
                        rate = math.pow(10., im.get_fit(dE)) * gg * math.exp(dE * h_c_k / T[it])
                        if inter_vibrational_trans: rate *= 1.e-3
                        break
            fout.write( ' {:.1E}'.format(rate) )
        fout.write("\n")
        trans_idx = trans_idx + 1

fout.close()
