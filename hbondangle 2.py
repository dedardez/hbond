

#! /usr/bin/env python

__author__ = "Suman Sirimulla"
__authoraffliation__ = "Department of Chemistry & Biochemsitry, Northern Arizona University"
__authoremail___= "suman.sirimulla@nau.edu"

__co_author__='Erik M. Chavez'
__co_authoraffliation__ = "Undergraduate Researcher,Chemistry, Northern Arizona University"
__co_authoremail___ = "emc253@nau.edu"


import re
import math
import numpy as np
import itertools
import sys
import os
from Tkinter import Tk
from Tkinter import *
import ttk
from tkFileDialog import askdirectory

Tk().withdraw()


userpdb=askdirectory(title='Select PDB containing folder')
pdbfiles=os.listdir(userpdb)
useroutpath=askdirectory(title='Where would you like the ouput to be saved?')

useroutname=raw_input('Enter desired name for output file: ')+'.txt'

#the desired name and filepath of output file
rawoutfile=open(useroutpath+'/'+useroutname,'w')
rawoutfile.close()
rawoutfile=open(useroutpath+'/'+useroutname,'a')
cycnum=0


#outward loop
#repeats searh for each file in pdb folder
print("Checking pdb files...")
header='{:4s} {:6s} {:3s} {:3s} {:2s} {:6s} {:4s} {:3s} {:3s} {:1s} {:6s} {:3s} {:3s} {:4s} {:2s} {:6s} {:6s} {:6s}'.format('PDB ','Het#   ',
                                                    'Het','Res',
                                                     'Ty','Het#  ',
                                                     'Het ','Res',
                                                     'Ty ','|','Pro#   ',
                                                     'Atm','Res',
                                                     'R# ','Ty',
                                                     'C-X--D','  Dis ','C-X--'+chr(227))
rawoutfile.write(header+ '\n')



for fil in pdbfiles:

    InFileName = userpdb +'/'+ fil
    #OutFileName = "results_01_23.txt"

    InFile = open(InFileName, 'r')
    #OutFile=open(OutFileName,'a')
    dict_hethalogen = { }
    dict_hethal={}
    dict_hetcarbon= { }
    dict_conect = { }
    dict_atom = { }

    i = 0
    for line in InFile:
            record_type=line[0:6]
            atom_type=line[76:78] 
            if record_type == 'HETATM' and (atom_type == 'BR' or atom_type == 'CL' or atom_type == ' F' or atom_type == ' I' ):
                    dict_hethalogen[i] = {'record_type':line[0:6], 'serial_number':line[6:11],'fullname':line[12:16], 'altloc':line[16:17], 'resname':line[17:20], 'chainid':line[21:22], 'resseq':int(line[22:26]), 'icode':line[26:27], 'x':float(line[30:38]), 'y':float(line[38:46]), 'z':float(line[46:54]), 'atom_type':line[76:78]}
                    #print dict_hethalogen[i]
                    #print line
                    i = i + 1
    f = i
    #print (dict_hethalogen)
    InFile.close()

    InFile = open(InFileName, 'r')

    #Benzene Dictionaries
    modres=[]
    dict_tyr = { }
    dict_phe = { }
    dict_trp = { }
    dict_cen_benz={ }

    ty=0
    ph=0
    tr=0
    for line in InFile:
            record_type=line[0:6]
            res_name=line[17:20]
            side_chain=line[12:16]
            if record_type == 'ATOM  ' and res_name == 'TYR' and (side_chain==' CG 'or side_chain==' CD1' or side_chain==' CD2' or side_chain==' CE1' or side_chain==' CE2' or side_chain==' CZ '):
                    dict_tyr[ty] = {'record_type':line[0:6], 'serial_number':line[6:11],'fullname':line[12:16], 'altloc':line[16:17], 'resname':line[17:20], 
                    'chainid':line[21:22], 'resseq':int(line[22:26]), 'icode':line[26:27], 'x':float(line[30:38]), 'y':float(line[38:46]), 'z':float(line[46:54]), 'atom_type':line[76:78]}
                    #print(dict_tyr[ty]['resname'], dict_tyr[ty]['fullname'])
                    ty += 1
            elif record_type == 'ATOM  ' and res_name == 'PHE' and (side_chain==' CG 'or side_chain==' CD1' or side_chain==' CD2' or side_chain==' CE1' or side_chain==' CE2' or side_chain==' CZ '):
                    dict_phe[ph] = {'record_type':line[0:6], 'serial_number':line[6:11],'fullname':line[12:16], 'altloc':line[16:17], 'resname':line[17:20], 
                    'chainid':line[21:22], 'resseq':int(line[22:26]), 'icode':line[26:27], 'x':float(line[30:38]), 'y':float(line[38:46]), 'z':float(line[46:54]), 'atom_type':line[76:78]}
                    #print(dict_phe[ph]['resname'], dict_phe[ph]['fullname'])
                    ph += 1
            elif record_type == 'ATOM  ' and res_name == 'TRP' and (side_chain==' CD2'or side_chain==' CE2' or side_chain==' CE3' or side_chain==' CZ2' or side_chain==' CZ3' or side_chain==' CH2'):
                    dict_trp[tr] = {'record_type':line[0:6], 'serial_number':line[6:11],'fullname':line[12:16], 'altloc':line[16:17], 'resname':line[17:20], 
                    'chainid':line[21:22], 'resseq':int(line[22:26]), 'icode':line[26:27], 'x':float(line[30:38]), 'y':float(line[38:46]), 'z':float(line[46:54]), 'atom_type':line[76:78]}
                    #print(dict_trp[tr]['resname'], dict_trp[tr]['fullname'])
                    tr += 1
                

    InFile.close()

    #Finding the center and normal vector of benzene rings
    def benz_centerx(x1,x2,x3,x4,x5,x6):
            return (x1+x2+x3+x4+x5+x6)/6
    def benz_centery(y1,y2,y3,y4,y5,y6):
            return (y1+y2+y3+y4+y5+y6)/6
    def benz_centerz(z1,z2,z3,z4,z5,z6):
            return (z1+z2+z3+z4+z5+z6)/6


    def normal_v(x1,x2,x3,y1,y2,y3,z1,z2,z3):
            v1=[x2-x1,y2-y1,z2-z1]
            v2=[x3-x1,y3-y1,z3-z1]
            return np.cross(v1,v2)
        
    #dictionary with only one index. used for loops
    dict_cycle={}

    ty+=1
    benzall=0
    #print("")
    #print("Tyrosine rings")
    #print("")
    for carb in range(0,ty//6):
            dict_cycle[0]={'cenx': benz_centerx(dict_tyr[6*carb]['x'],dict_tyr[(6*carb)+1]['x'],dict_tyr[(6*carb)+2]['x'],dict_tyr[(6*carb)+3]['x'],dict_tyr[(6*carb)+4]['x'],dict_tyr[(6*carb)+5]['x']),
                           'ceny': benz_centery(dict_tyr[6*carb]['y'],dict_tyr[(6*carb)+1]['y'],dict_tyr[(6*carb)+2]['y'],dict_tyr[(6*carb)+3]['y'],dict_tyr[(6*carb)+4]['y'],dict_tyr[(6*carb)+5]['y']),
                           'cenz': benz_centerz(dict_tyr[6*carb]['z'],dict_tyr[(6*carb)+1]['z'],dict_tyr[(6*carb)+2]['z'],dict_tyr[(6*carb)+3]['z'],dict_tyr[(6*carb)+4]['z'],dict_tyr[(6*carb)+5]['z'])
                           }
            dict_cen_benz[benzall]={'cenx': dict_cycle[0]['cenx'],
                                    'ceny': dict_cycle[0]['ceny'] ,
                                    'cenz': dict_cycle[0]['cenz'],
                                    'norm_v': normal_v(dict_cycle[0]['cenx'],dict_tyr[6*carb]['x'],dict_tyr[(6*carb)+1]['x'],dict_cycle[0]['ceny'],dict_tyr[carb*6]['y'],dict_tyr[(carb*6)+1]['y'],
                                                   dict_cycle[0]['cenz'],dict_tyr[carb*6]['z'],dict_tyr[(carb*6)+1]['z']),
                                    'resseq':dict_tyr[6*carb]['resseq'],
                                    'resname':dict_tyr[6*carb]['resname'],
                                    'sn_range': [dict_tyr[6*carb]['serial_number'],dict_tyr[(6*carb)+5]['serial_number']]}


            benzall+=1
    ph+=1
    for carb in range(0,ph//6):
            dict_cycle[0]={'cenx': benz_centerx(dict_phe[6*carb]['x'],dict_phe[(6*carb)+1]['x'],dict_phe[(6*carb)+2]['x'],dict_phe[(6*carb)+3]['x'],dict_phe[(6*carb)+4]['x'],dict_phe[(6*carb)+5]['x']),
                           'ceny': benz_centery(dict_phe[6*carb]['y'],dict_phe[(6*carb)+1]['y'],dict_phe[(6*carb)+2]['y'],dict_phe[(6*carb)+3]['y'],dict_phe[(6*carb)+4]['y'],dict_phe[(6*carb)+5]['y']),
                           'cenz': benz_centerz(dict_phe[6*carb]['z'],dict_phe[(6*carb)+1]['z'],dict_phe[(6*carb)+2]['z'],dict_phe[(6*carb)+3]['z'],dict_phe[(6*carb)+4]['z'],dict_phe[(6*carb)+5]['z'])
                           }
            dict_cen_benz[benzall]={'cenx': dict_cycle[0]['cenx'],
                                    'ceny': dict_cycle[0]['ceny'] ,
                                    'cenz': dict_cycle[0]['cenz'],
                                    'norm_v': normal_v(dict_cycle[0]['cenx'],dict_phe[6*carb]['x'],dict_phe[(6*carb)+1]['x'],dict_cycle[0]['ceny'],dict_phe[carb*6]['y'],dict_phe[(carb*6)+1]['y'],
                                                   dict_cycle[0]['cenz'],dict_phe[carb*6]['z'],dict_phe[(carb*6)+1]['z']),
                                    'resseq':dict_phe[6*carb]['resseq'],
                                    'resname':dict_phe[6*carb]['resname'],
                                    'sn_range': [dict_phe[6*carb]['serial_number'],dict_phe[(6*carb)+5]['serial_number']]}


            benzall+=1

    tr+=1

    for carb in range(0,tr//6):
            dict_cycle[0]={'cenx': benz_centerx(dict_trp[6*carb]['x'],dict_trp[(6*carb)+1]['x'],dict_trp[(6*carb)+2]['x'],dict_trp[(6*carb)+3]['x'],dict_trp[(6*carb)+4]['x'],dict_trp[(6*carb)+5]['x']),
                           'ceny': benz_centery(dict_trp[6*carb]['y'],dict_trp[(6*carb)+1]['y'],dict_trp[(6*carb)+2]['y'],dict_trp[(6*carb)+3]['y'],dict_trp[(6*carb)+4]['y'],dict_trp[(6*carb)+5]['y']),
                           'cenz': benz_centerz(dict_trp[6*carb]['z'],dict_trp[(6*carb)+1]['z'],dict_trp[(6*carb)+2]['z'],dict_trp[(6*carb)+3]['z'],dict_trp[(6*carb)+4]['z'],dict_trp[(6*carb)+5]['z'])
                           }
            dict_cen_benz[benzall]={'cenx': dict_cycle[0]['cenx'],
                                    'ceny': dict_cycle[0]['ceny'] ,
                                    'cenz': dict_cycle[0]['cenz'],
                                    'norm_v': normal_v(dict_cycle[0]['cenx'],dict_trp[6*carb]['x'],dict_trp[(6*carb)+1]['x'],dict_cycle[0]['ceny'],dict_trp[carb*6]['y'],dict_trp[(carb*6)+1]['y'],
                                                   dict_cycle[0]['cenz'],dict_trp[carb*6]['z'],dict_trp[(carb*6)+1]['z']),
                                    'resseq':dict_trp[6*carb]['resseq'],
                                    'resname':dict_trp[6*carb]['resname'],
                                    'sn_range': [dict_trp[6*carb]['serial_number'],dict_trp[(6*carb)+5]['serial_number']]}

        
            benzall+=1
        
        
    InFile = open(InFileName, 'r')
    n=0
    for line in InFile:
        record_type=line[0:6]
        serial_number=line[6:11]
        if record_type == 'CONECT':
            for q in dict_hethalogen:
                if serial_number==dict_hethalogen[q]['serial_number']:
                    dict_conect[n] = {'serial_number1':line[6:11], 'serial_number2':line[11:16]}    
                    n = n + 1  
                    
    InFile.close()
    i = 0
    n= 0
    InFile = open(InFileName, 'r')
    for line in InFile:
        record_type=line[0:6]
        serial_number=line[6:11]
        atom_type=line[76:78]
        res_name=line[17:20]
        found=False
        if record_type == 'HETATM' and atom_type == ' C' and (res_name not in modres):
            for con in dict_conect:
                if serial_number==dict_conect[con]['serial_number2']:
                    for hal in range(len(dict_hethalogen)):
                        if dict_hethalogen[hal]['serial_number']==dict_conect[con]['serial_number1']:
                            dict_hethal[n]=dict_hethalogen[hal]
                                
                            dict_hetcarbon[n] = {'record_type':line[0:6], 'serial_number':line[6:11],'fullname':line[12:16], 'altloc':line[16:17], 
                                                 'resname':line[17:20], 'chainid':line[21:22], 'resseq':int(line[22:26]), 'icode':line[26:27],
                                                 'x':float(line[30:38]), 'y':float(line[38:46]), 'z':float(line[46:54]), 'atom_type':line[76:78],
                                                 'H1':None,'H2':None,'H3':None}  
                                             

                            n = n + 1
                            break
    w = n 
    InFile.close()

    i = 0
    d = 0
    s = 0
    for j in range(0,w):
            InFile = open(InFileName, 'r')
            for line in InFile:
                    record_type=line[0:6]
                    if record_type == 'ATOM  ':
                            #print line
                            try:
                                distance = math.sqrt(((dict_hethal[j]['x'] - float(line[30:38]))**2) + ((dict_hethal[j]['y'] - float(line[38:46]))**2) + ((dict_hethal[j]['z'] - float(line[46:54]))**2)) # calculate length between halogen and atoms in protein
                            except(ZeroDivisionError,ValueError):
                                continue
                            #print distance
                            if distance < 4.0: # if length is less than 4.0 
                                    dict_atom[s] = {'record_type':line[0:6], 'serial_number':line[6:11],'fullname':line[12:16], 'altloc':line[16:17], 'resname':line[17:20], 'chainid':line[21:22], 'resseq':int(line[22:26]), 'icode':line[26:27], 'x':float(line[30:38]), 'y':float(line[38:46]), 'z':float(line[46:54]), 'atom_type':line[76:78]}
                                            #print distance
                                    #print dict_atom[s]
                      
                                    v1 = [(dict_hethal[j]['x']-dict_hetcarbon[j]['x']),(dict_hethal[j]['y']-dict_hetcarbon[j]['y']),(dict_hethal[j]['z']-dict_hetcarbon[j]['z'])] # vector v1 (hethalogen - hetcarbon)
                                    v2 = [(dict_hethal[j]['x']-dict_atom[s]['x']),(dict_hethal[j]['y']-dict_atom[s]['y']),(dict_hethal[j]['z']-dict_atom[s]['z'])] # vector v2 (hethalogen- atom)

                                    def dotproduct(v1, v2):
                                            return sum((a*b) for a, b in zip(v1, v2))
                                    def length(v):
                                            return math.sqrt(dotproduct(v, v))
                                
                                    length_v1 = length(v1) #calculate length of vector 1
                                    length_v2 = length(v2) #calculate length of vector 2
                                
                                    def angle(v1, v2): 
                                            return math.acos(dotproduct(v1, v2) / (length_v1 * length_v2))  
                                    try:
                                        Angle = angle(v1, v2) *180/3.141592
                                    except(ZeroDivisionError,ValueError):
                                        continue
                  
                                    rawout='{:4s} {:6d} {:3s} {:3s} {:2s} {:6d} {:4s} {:3s} {:3s} {:1s} {:6d} {:3s} {:3s} {:4d} {:2s} {:6.2f} {:6.2f}'.format(InFileName[-8:-4],int(dict_hetcarbon[j]['serial_number']),
                                                     dict_hetcarbon[j]['fullname'],dict_hetcarbon[j]['resname'],
                                                     dict_hetcarbon[j]['atom_type'],int(dict_hethal[j]['serial_number']),
                                                     dict_hethal[j]['fullname'],dict_hethal[j]['resname'],
                                                     dict_hethal[j]['atom_type'],'|',int(dict_atom[s]['serial_number']),
                                                     dict_atom[s]['fullname'],dict_atom[s]['resname'],
                                                     int(dict_atom[s]['resseq']),dict_atom[s]['atom_type'],
                                                     round(Angle,2),round(distance,2))
                                    rawoutfile.write(rawout+'\n')
                                
                    
    InFile.close()
    
    #Finding benzine rings within 6a
    rings=True
    if benzall==0:
        rings=False

    def benz_dis(x1,x2,y1,y2,z1,z2):
        return math.sqrt((x2-x1)**2 +(y2-y1)**2 +(z2-z1)**2)

    def plane_angle(x1,x2,y1,y2,z1,z2,normp):
        hv=[x2-x1,y2-y1,z2-z1]
        return round(math.degrees(math.asin((np.dot(hv,normp))/(benz_dis(hv[0],0,hv[1],0,hv[2],0)*benz_dis(normp[0],0,normp[1],0,normp[2],0)))),2)
    def bh_bond_angle(hv,bv):
        return round(math.degrees(math.acos((np.dot(hv, bv) / ((benz_dis(hv[0],0,hv[1],0,hv[2],0) * benz_dis(bv[0],0,bv[1],0,bv[2],0)))))),2)


    for ent in range(0,w):
        for b in range(0,benzall):
            try:
                bdis=benz_dis(dict_cen_benz[b]['cenx'],dict_hethal[ent]['x'],
                              dict_cen_benz[b]['ceny'],dict_hethal[ent]['y'],
                              dict_cen_benz[b]['cenz'],dict_hethal[ent]['z'])
                              
                bpangle=plane_angle(dict_cen_benz[b]['cenx'],dict_hethal[ent]['x'],
                                    dict_cen_benz[b]['ceny'],dict_hethal[ent]['y'],
                                    dict_cen_benz[b]['cenz'],dict_hethal[ent]['z'],
                                    dict_cen_benz[b]['norm_v'])
                baangle=bh_bond_angle([dict_hetcarbon[ent]['x']-dict_hethal[ent]['x'],
                                      dict_hetcarbon[ent]['y']-dict_hethal[ent]['y'],
                                      dict_hetcarbon[ent]['z']-dict_hethal[ent]['z']],
                                     [dict_cen_benz[b]['cenx']-dict_hethal[ent]['x'],
                                      dict_cen_benz[b]['ceny']-dict_hethal[ent]['y'],
                                      dict_cen_benz[b]['cenz']-dict_hethal[ent]['z']])
            except(ZeroDivisionError,ValueError):
                continue
            if rings and bdis<= 6:
                rawout='{:4s} {:6d} {:3s} {:3s} {:2s} {:6d} {:4s} {:3s} {:3s} {:1s} {:6s} {:3s} {:3s} {:4d} {:2s} {:6.2f} {:6.2f} {:6.2f}'.format(InFileName[-8:-4],int(dict_hetcarbon[j]['serial_number']),
                                                     dict_hetcarbon[ent]['fullname'],dict_hetcarbon[ent]['resname'],
                                                     dict_hetcarbon[ent]['atom_type'],int(dict_hethal[ent]['serial_number']),
                                                     dict_hethal[ent]['fullname'],dict_hethal[ent]['resname'],
                                                     dict_hethal[ent]['atom_type'],'|','      ',
                                                     '  R ',dict_cen_benz[b]['resname'],
                                                     int(dict_cen_benz[b]['resseq']),' ',
                                                     round(baangle,2),round(bdis,2),round(bpangle,2))                                
                rawoutfile.write(rawout+'\n')
                    

rawoutfile.close() 
print('Halogen interaction Analysis Complete')