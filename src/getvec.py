#!/usr/bin/env python3
import re
import sys
import os
import routines
import xmvb_orb as xmlib
import numpy as np  
import sys

#Core function which does the job
def main(input_file):   
#    print(len(sys.argv)," arguments",sys.argv[0] )
    print("getvec.py - SH 2024                --")
    return input_file


if __name__ == "__main__": # permet d'utiliser comme une librairie qu'on importe
    if len(sys.argv) >= 2: 
        zeta=1
        numatoms=[]
        coeffs=[]
        indices=[]
        input_file = sys.argv[1]
        input_file_name, input_file_ext = os.path.splitext(input_file)
        output_file_name=input_file_name+".orbb"
        main(input_file)
        input_file_name, input_file_ext = routines.file_extension(input_file)
#        print(input_file_name, input_file_ext)
        pos,line=routines.detect_keyword(input_file, "VEC", 0)
        coeffs,nvect = routines.read_vec(input_file,coeffs,pos+1)
        vect=routines.make_table(coeffs)#
#        routines.write_orbs("screen",vect,0,len(vect))
        print()
        fin=int(input("highest MO:"))
        print("Select MOs 1-",fin, ' among ', len(vect)," MOs from ",input_file,"to ", output_file_name," in .orb format")
        routines.write_orbs(output_file_name,vect,0,fin) 
#        routines.write_orbs("screen",vect,0,fin)
#------
        #indices=routines.compte_AO(vect) # compte_AO est fait dans make_orb
#        print('##############VECT ##','\n',"{:8.5f}".format(float(vect[0][1])))
#        print('ifin',coeffs)
#        routines.write_vec("screen",coeffs,1, nvect)
#        routines.write_vec(output_file_name, coeffs,1, nvect)
#        indices=routines.compte_AO(vect) # compte_AO est fait dans make_orb
        #print(vec
        #routines.write_orb("screen",vect,indices,0,len(vect))
#        print('len(indices)',len(indices),end=' ' )
        ao_orb,coeffs_orb=routines.make_orb(vect,indices)
#        print('len(indices_orb)',len(indices_orb),'indices_orb[46]',indices_orb[46])
#        print(len(coeffs_orb)," orbitals")#,'coeffs_orb[46]', "{:5.3f}".format(coeffs_orb[46][1]))
#        print(indices)
#        print('---ao     _orb)',ao_orb)
#        print('---coeffs_orb)',len(coeffs_orb),coeffs_orb[0])
#        routines.write_vec("screen",coeffs_orb,1,len(coeffs_orb) )
#        routines.write_orb("screen",vect,ao_orb,0,len(vect))
        #routines.write_orb("screen",vect,ao_orb)
#        print(ao_orb)
        #copy to ccoeffs only the coeff to print aligned to indices
#        ccoeffs=[]
#        for i in range(len(indices)):
#             print('-----------')
#             print('bcl indices(',i,')',indices[i],len(indices[i]))    
#             for j in range(len(indices[i])):
                #ccoeffs.append(coeffs[indices[i][i]])
#                print(len(indices[i]),'indices(',i,',',j,' )=',indices[i][j],end=' ,')
            #   coeffs[i]=ccoeffs
            #   print('coeff  (',i,',',j,')',coeffs[i][j],end='')
            #    print(' coeffs(',i,',',j,')=',coeffs[i][1][indices[i][j]],end='')
            #    ccoeffs.append(coeffs[i][1][indices[i][j]])
#             print(i,'-----------')
#             #print('ccoeffs',ccoeffs)    
        if len(sys.argv) == 3:
            xmvb_input_file=sys.argv[2]
#            xmvb_orb_file=sys.argv[3]
            #print('xmvb_input_file',xmvb_input_file,'xmvb_orb_file',xmvb_orb_file)
#            print('+++++++++++++++++++++++++++++++-',ao_orb)
            vb_orb_coeffs,vb_orb_indices=routines.read_orb(output_file_name)
#            print()
#            print('--------------------------------',ao_orb)
#            print('vb_orb_coeffs',vb_orb_coeffs,'indices',vb_orb_indices)
#            #routines.write_orb("screen",vb_orb_coeffs,vb_orb_indices)
            pos, line=routines.detect_keyword(xmvb_input_file, "bfi", 0)
            if pos == -1:
                print('bfi not found, use lower case bfi')
                print('########.... .STOP.  ######## bfi')
                sys.exit()
            print('bfi found in ',xmvb_input_file,'at line',pos)
            bfi_nom,bfi_noa,list_om,list_oa=routines.read_bfi(xmvb_input_file,pos)
            print(bfi_nom,'MOs to freeze',': ',list_om)
            print(bfi_noa,'OAs to keep  ',': ',list_oa)
            for p in list_oa:  # test inutile en principe
                try:
                     t=list_oa.index(p)
                     #print(t, end=', ')
                except ValueError:
                     print('error',p,' not il the AOs list',list_oa,list_om)
            print()
            new_coeffs = []
            new_MOs = []
            print('>>',len(vb_orb_coeffs),' vb vectors and list_om=',list_om)
            for i in range(len(vb_orb_coeffs)):
                 if i+1 not in list_om: # i est forcement + petit que fin car relu 
                      new_coeffs.append(vb_orb_coeffs[i])
                      #print('\n##',i+1, vb_orb_coeffs[i])
                      print(i+1, end=', ')
                      #print('\n__',i+1,vb_orb_indices[i] )
                      for k in vb_orb_indices[i]:
                                #print(list_oa.index(k)+1,end=' ')
                                continue
                      #print('--',vb_orb_indices[list_oa[i]])   
                 else : # i est dans les MO's a garder
                      new_MOs.append(vb_orb_coeffs[i])
                      #print('\n##',i+1, vb_orb_coeffs[i])
                      #print(i+1, end=', ')
#
#            print('new_coeffs',new_coeffs[2])
            fin=len(new_coeffs)
            finmos=len(new_MOs)
            xmvb_input_file_name, xmvb_input_file_ext = os.path.splitext(xmvb_input_file)
            xmvb_output_file=xmvb_input_file_name+'.orbb'
            gamess_output_file=xmvb_input_file_name+'.vecc'
            routines.write_orbs(xmvb_output_file,new_coeffs,0,fin) 
            routines.write_vec(gamess_output_file,new_MOs,0,finmos-1) 
                    # new_row = []
                    # print('new_row',new_row,list_om[1], 'i=',i, len(vb_orb_coeffs[list_om[i]]))
                    # new_row.append(coeffs[list_om[i-1]])
                    # new_coeffs.append(new_row)
                    # print('new_coeffs',new_coeffs[list_om[i]])
#            routines.write_orb("screen",new_coeffs,indices)
     #        coeffs = new_coeffs
    #        for i in range(len(indices)):
    #             for j in range(len(indices[i])):
    #                  print('i',i,'j',j,indices[i][j],end=' ')
    #                  print('list_oa.index(indices[i][j])',list_oa.index(indices[i][j]+1), end=" .  ")   
    #                  new_indices[i][j] = list_oa.index(indices[i][j])
    #        print()
        #routines.write_orb("screen",coeffs,indices)routines.write_orb("screen",coeffs,indices)
    #    routines.write_orb("screen",coeffs,new_indices)
        sys.exit() 
