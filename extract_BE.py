from aiida.orm import DataFactory
from aiida.orm.querybuilder import QueryBuilder
import sys
import numpy as np
import pandas as pd

def BE(calc_id,shift): # extra term for the energy shift between nscf and scf
    import xml.etree.ElementTree as ET
    import glob

    c = load_node(calc_id)
    occu = c.out.output_band.get_array('occupations')
    bands = c.out.output_band.get_array('bands') #--unit in ev.

    CBM = np.inf
    VBM = - np.inf
    for i in range(occu.shape[0]):
        for j in range(len(occu[i])):
            if float(occu[i,j+1]) <= 1e-4: break
        CBMt = float(bands[i,j+1])
        VBMt = float(bands[i,j])
        if CBMt < CBM: CBM = CBMt
        if VBMt > VBM: VBM = VBMt
    #band_gap = (CBM - VBM)
    VBM = VBM - shift
    CBM = CBM - shift
    return VBM,CBM

def site_term_atm(calc_id):
    c = load_node(calc_id)
    desc = c.label
    term_atm = desc.split()[8] # extract termination from label
    
    formula = c.inp.structure.get_kind_names() # extract elemets from label, and remove oxygen.
    
    formula.remove('O')

    #print(formula)

    atm_num = sum(c.inp.structure.get_composition().values())

    if atm_num ==26: # A termination contains 26 atoms
        A_site = term_atm
        formula.remove(A_site)
        B_site = formula[0]

    elif atm_num == 25:# B termination is 25 atoms
        B_site = term_atm
        formula.remove(B_site)
        A_site = formula[0]

    else:
        print(str(calc_id)+" doesn't have the right structure")
    return A_site,B_site,term_atm

def shift_fermi(calc_id):
    
    c = load_node(calc_id)
    nscf_fermi = c.res.fermi_energy # nscf fermi energy

    par_c = c.inp.structure.inp.output_structure #parent relaxation structure fermi energy
    par_fermi = par_c.res.fermi_energy

    shift_energy = nscf_fermi - par_fermi

    return shift_energy

def vac_BE():

    qb = QueryBuilder()
    
    qb.append(Group, tag = "group",project = "id",filters={"id":{"in":[263, 265, 266, 267, 268, 278, ]}})
    qb.append(JobCalculation, member_of="group", tag = "calculation",filters={"state":{"==":"FINISHED"}},project=["id"])
    calc_list = qb.dict()
    
    print "Total slab structures in vacuum %s . . ."%len(calc_list)

    with open("vacuum_slab_BE.txt","w") as f:

        for bulk_calc in calc_list:
            shift_energy = shift_fermi(bulk_calc['calculation']['id'])
            VBM, CBM = BE(bulk_calc['calculation']['id'],shift_energy)
            A_site, B_site, term_site = site_term_atm(bulk_calc['calculation']['id'])
            f.write(str(A_site) + "    " + str(B_site) + "    " + str(term_site) + "    " + str(VBM) + "    " + str(CBM) + "    " + str(bulk_calc['calculation']['id'])+"\n")

def water_BE():

    qb = QueryBuilder()
    
    qb.append(Group, tag = "group",project = "id",filters={"id":{"in":[269, 271, 272, 273, 274, 275, 277, 279, 280, ]}})
    qb.append(JobCalculation, member_of="group", tag = "calculation",filters={"state":{"==":"FINISHED"}},project=["id"])
    calc_list = qb.dict()
    
    print "Total slab structures in water %s . . ."%len(calc_list)

    with open("water_slab_BE.txt","w") as f:

        for bulk_calc in calc_list:
            shift_energy = shift_fermi(bulk_calc['calculation']['id'])
            VBM, CBM = BE(bulk_calc['calculation']['id'],shift_energy)
            A_site, B_site, term_site = site_term_atm(bulk_calc['calculation']['id'])
            f.write(str(A_site) + "    " + str(B_site) + "    " + str(term_site) + "    " + str(VBM) + "    " + str(CBM) + "    " + str(bulk_calc['calculation']['id'])+"\n")


def main():

    vac_BE()
    water_BE()

if __name__=='__main__':
    main()
