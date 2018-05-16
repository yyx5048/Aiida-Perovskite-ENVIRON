#=========================================================================#
#--Workflow for restart ENVIRON water conditions -------------------------#
#=========================================================================#
#-- import Aiida module 
from aiida import load_dbenv, scheduler
from aiida.orm import DataFactory
from aiida.orm import Code
from aiida.orm.querybuilder import QueryBuilder
from aiida.work.run import run
from aiida.work.workfunction import workfunction as wf
from aiida.orm.data.base import Str, Float, List, Int
import numpy as np
import json

StructureData = DataFactory('structure')

ParameterData = DataFactory('parameter')

Ang2Bohr = 1.88973 #--angstrom to bohr constant

def load_failed_environ_calc(group_node):
    """
    load output structure of previous relaxed structures
    """
    qb = QueryBuilder()
    qb.append(Group,tag = 'group',filters={'id':{'in':[group_node]}})
    qb.append(JobCalculation, member_of = 'group',project = ['id'], tag = 'calculation',filters={
        "or":[
            #{'state':{'==':'SUBMISSIONFAILED'}},
            {'state':{'==':'FAILED'}}]})

    id_list = qb.dict()
    res = [x['calculation']['id'] for x in id_list]
    return res 

#def load_parent_calc(child_calc):
#    """
#    trace back to the parent calculation
#    """
#    parent_calc = child_calc.inp.structure.inp.output_structure
#    parent_struc = parent_calc.out.output_structure # finished calculation structure is used in others
#    return parent_calc,parent_struc

#def push_in_structure_list():
#    """
#    slab structure with pushed in oxygen atoms
#    are stored in this list 
#    """
#    qb = QueryBuilder()
#    qb.append(Group,tag = 'group',filters={'id':{'in':[209,214,215,216,218,224]}})
#    qb.append(JobCalculation, member_of = 'group',tag = 'calculation',filters={
#        "or":[
#            {'state':{'==':'FINISHED'}}]})
#    qb.append(StructureData,output_of = 'calculation',tag = 'shrink_structure',project= ['id'])

#    id_list = qb.dict()
#    res = [x['shrink_structure']['id'] for x in id_list]
#    return res


def slab_structure_list():
    """
    load the slab structure list and label in a list of dictionaries
    for structure matching
    """
    group_node = 199 # original slab structure generation without push-in oxygens
    qb = QueryBuilder()
    qb.append(Group,tag = 'group',filters={'id':{'in':[group_node]}})
    qb.append(StructureData, member_of = 'group', tag = 'slab_structure',filters={'label':{'ilike':'%slab%'}}, project = ['id','label'])
    res = [x['slab_structure'] for x in qb.dict()]#slab structure id and their label
    return res

@wf
def environ_dielectric_regions(slab_struct):
    z_len = [sit['position'][2]for sit in slab_struct.get_attrs()['sites']]#--grep z direction coordinates

    z_mid = (max(z_len)+min(z_len))/2*Ang2Bohr
    z_spread = (max(z_len)-min(z_len))/2*Ang2Bohr #--a lttle bit more initial cavity
    z_xy_mid1 = slab_struct.get_ase().cell[0][0]/2*Ang2Bohr #-- mid point of x direction
    z_xy_mid2 = slab_struct.get_ase().cell[1][1]/2*Ang2Bohr #-- mid point of y direction
    diel_region = List()
    diel_region.append([1, 1, z_xy_mid1,z_xy_mid2,z_mid,z_spread,0.5,2,3]) #--dielectric regions card for ENVIRON
    diel_region.store()

    return diel_region

def main():
    #===load code#===#
    codename = 'QE-SNB-6.1_20_CPUs_update@aci_ENVIRON_update'

    code = Code.get_from_string(codename)

    #pushin_list = push_in_structure_list() # exclude the structure with pushed in oxygen atoms
    
    par_s_info = slab_structure_list() #load parent slab structure info (id and label)

    par_s_id = [x['id'] for x in par_s_info]

    par_s_label = [x['label'] for x in par_s_info]

    group_node = [213,221,225] #-- group number that contains faield calculations   

    KpointsData = DataFactory('array.kpoints')
    kpoints = KpointsData()
    kpoints.set_kpoints_mesh([4,4,1])

    for group_id in group_node:

        failed_environ_list = load_failed_environ_calc(group_id)
        
        for cid in failed_environ_list: #-- load calculation id (cid)

            c1 = load_node(cid)

            lab = c1.label
            
            #================#
            #slab label match#
            #================#
            # better query: using label searching the structure

            identifier_label = " ".join(lab.split()[:5])+" "

            identi_idx = par_s_label.index(identifier_label)
            
            s = load_node(par_s_id[identi_idx])

            #par_c, par_s = load_parent_calc(c1) #load parent calculation and structure

            old_parameters = c1.inp.parameters
            parameters_dict = old_parameters.get_dict()
            old_settings = c1.inp.settings
            settings_dict = old_settings.get_dict()

            parameters_dict['CONTROL'].update({'restart_mode':'from_scratch'})
            #parameters_dict['ELECTRONS'].update({'mixing_beta':0.2})
            parameters_dict['IONS'].update({'ion_dynamics':'damp'})


            settings_dict['ENVIRON']['environ_type'] = 'water' #--change environ type to water
            settings_dict['ENVIRON']['environ_thr'] = 5
            settings_dict['ENVIRON']['verbose'] = 0 #--set to lower verbosity
            settings_dict['ENVIRON'].update({'environ_restart':False}) #--restart from previous environ calculation
            settings_dict['ENVIRON'].update({'env_dielectric_regions':1}) #--diel region card
            diel_region = environ_dielectric_regions(s)
            settings_dict['ENVIRON'].update({'DIELECTRIC_REGIONS':diel_region[0]}) #--need to update environ card

            settings_dict.update({'cmdline':['-ndiag','1']})


            c2 = code.new_calc()
            c2.label = str(lab)
            c2.set_resources({"num_machines": 1})
            c2.set_max_wallclock_seconds(24*60*60)
            c2.use_structure(s)
            c2.use_pseudos_from_family('SSSP_eff_aiida')
            c2.use_kpoints(kpoints)
            c2.use_parameters(ParameterData(dict=parameters_dict))
            c2.use_settings(ParameterData(dict=settings_dict))
            c2.set_custom_scheduler_commands('#PBS -A ixd4_d_g_sc_default \n#PBS -l pmem=8gb')
            c2.store_all()
        
            c2.submit()
    return

if __name__ == '__main__':
    main()
