#!/usr/bin/env python3
import re
import sys
import os
import routines 
#import xmvb_orb as xmlib
import numpy as np  
import sys

#Core function which does the job
def main(input_file):   
#    print(len(sys.argv)," arguments",sys.argv[0] )
    print("+--------getvec.py - SH 2024 ---------------------")
    print('| getvec.py file_w_VEC [file_w_bfi]')
    return input_file


if __name__ == "__main__": # permet d'utiliser comme une librairie qu'on importe
    NORB_RHF=10000
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
        if input_file_ext == ".dat":
            pos,line=routines.detect_keyword(input_file, "$DATA", 0)
            pos2,line2=routines.detect_keyword(input_file, "$END", pos)
# met dans lines les lignes entre DATA et END
            with open(input_file, 'r') as file:
                in_data=False
                N_protons=[]
                for line_num, line in enumerate(file, 1):  
                    if '$END' in line:
                        in_data=False
                        break
                    if '$DATA' in line.upper():
                        in_data=True
              #          print('|  ',line.upper(), end=' ')
                        continue    
                    if in_data and len(line.split()) ==5: # 5 mots sur la ligne, c'est un atome et le second champ est le numero atomique
              #          print('*',line.split()[1],end='')
                        N_protons.append(float(line.split()[1]))
              #      print('|  ',line.upper(),end='_')
              #      print('|  ',len(line.split()),line.split())
            NORB_RHF =sum(N_protons)/2
            print('|  there are ',NORB_RHF*2 ,' electrons in the neutral molecule hence best choice is ',NORB_RHF,' orbitals')
#        print(input_file_name, input_file_ext)
        pos,line=routines.detect_keyword(input_file, "VEC", 0)
        coeffs,nvect = routines.read_vec(input_file,coeffs,pos+1)
        vect=routines.make_table(coeffs)#
        NORB_RHF=min(len(vect),NORB_RHF)
#        routines.write_orbs("screen",vect,0,len(vect))
        str1="| ==>   highest MO: (hit return for all MOs)  : ["+str(NORB_RHF)+"] max= "+str(len(vect))+" :"
        print("|  The following number MUST be consistent with the bfi (if provided) ")
        try:
            fin=int(input(str1))    
        except:
            print("|             I will take all ",NORB_RHF,' orbitals')
            fin = int(NORB_RHF)
        print("|     Select MOs 1-",fin, ' among ', len(vect)," MOs from ",input_file,"to ", output_file_name," in .orb format")
        print("                 -_________-  ") 
        print('|  all electron bfi:')
        routines.make_bfi(vect)
        print()
#        routines.write_orbs("screen",vect,0,fin)
        routines.write_orbs(output_file_name,vect,0,fin) 
        ao_orb,coeffs_orb=routines.make_orb(vect,indices)
        # Analyse des orbitales, proposition de bfi
        tab=[]
        pi_orb=[]
        sigma_orb=[]
        ratio_sigma_pi=3.0
        for i in range(len(ao_orb)):
            tab.append(len(ao_orb[i]))   
        mini=min(tab)
        maxi=max(tab)
        systeme_pi=int(maxi/ratio_sigma_pi)
        print('|  max  and min number of AO\'s :',maxi,'/',mini,'=',maxi/mini)
        print('|  systeme_pi for MOs with NOAs <',systeme_pi,'aos (ratio sigma/pi=',ratio_sigma_pi)
        if mini < systeme_pi: #  arbitrary ratio sigma/pi
            list_pi_aos=[]
            for i in range(fin):
                    if tab[i] <= systeme_pi:
                        pi_orb.append(i+1)
#                        print(i+1,end=' , ')
                    if tab[i] > systeme_pi:
                        sigma_orb.append(i+1)
            print('|  sigma orbitals:',sigma_orb)
            print('|  pi orbitals:',pi_orb)
#            print('|  oa_orb                                   :',ao_orb)
            for i in range(len(pi_orb)):
                index=i
#                print('pi_orb[',index,']=',pi_orb[index],':' ,ao_orb[pi_orb[index]-1])
                for j in range(len(ao_orb[pi_orb[index]-1])):
#                    print(ao_orb[pi_orb[index]-1][j]+1,end=' ')
                    list_pi_aos.append(ao_orb[pi_orb[index]-1][j]+1)
               # print(ao_orb[pi_orb[i]],end=' ')
#                'keep only unique values in list_pi_aos'    
#            list_pi_aos = list(dict.fromkeys(list_pi_aos))
            longue_pi_orb=tab.index(mini)
            print()
            print('|  temptative bfi:')
            print(' $bfi')
            print('  ',len(sigma_orb),' ',len(ao_orb[longue_pi_orb]     ))
            print(routines.makeSTR(sigma_orb,0))
#            print ((list_pi_aos))
            #print(routines.makeSTR(list_pi_aos))
            longue_pi_orb=tab.index(mini)
            #print(longue_pi_orb+1)
            print(routines.makeSTR(ao_orb[longue_pi_orb],1))
            print(' $end')
            print('| ------------    .')



        if len(sys.argv) == 3:
            xmvb_input_file=sys.argv[2]
            if xmvb_input_file=='orb':
                if os.path.isfile(xmvb_input_file):
                    print('|  orb is a reserved word, but a file named orb does exist')
                ao_orb, coeff_orb = routines.make_orb(vect,indices)
                output='orb.orb'
                routines.make_dollarorb_file(ao_orb,fin, output)
                print('|  $orb section written in ', output)
                print('|  for an all electron calculation and ',fin,' orbitals')
                sys.exit()  
            if not os.path.isfile(xmvb_input_file):
                    print('|  file ',xmvb_input_file,' not found')
                    sys.exit()
#            xmvb_orb_file=sys.argv[3]
            vb_orb_coeffs,vb_orb_indices=routines.read_orb(output_file_name)
#            routines.write_orb("screen",vb_orb_coeffs,vb_orb_indices)
            pos, line=routines.detect_keyword(xmvb_input_file, "bfi", 0)
            if pos == -1:
                print('|  bfi not found, use lower case bfi')
                print('   ########.... .STOP.  ######## bfi')
                print('+-     ------------------------------------------------------')
                sys.exit()
            print('|  bfi found in ',xmvb_input_file)
            bfi_nom,bfi_noa,list_om,list_oa=routines.read_bfi(xmvb_input_file,pos)
            print('| ', bfi_nom,' MOs to freeze',': ',list_om)
            print('| ',bfi_noa,' OAs to keep  ',': ',list_oa)
            for p in list_oa:            # test inutile en principe
                try:
                     t=list_oa.index(p)
                     #print(t, end=', ')
                except ValueError:
                     print('|   error',p,' not in the AOs list',list_oa,list_om)
            print('|  ')
            new_coeffs = []
            new_MOs = []
            print('| > among ',len(vb_orb_coeffs),' orbitals, ')# and list_om=',list_om)
            print('| > here is the VB orb list       : ', end=' ')
            compteorbvb=0
            for i in range(len(vb_orb_coeffs)):
                 if i+1 not in list_om: # i est forcement + petit que fin car relu 
                      new_coeffs.append(vb_orb_coeffs[i])
                      compteorbvb+=1
                      #print('\n##',i+1, vb_orb_coeffs[i])
                      print(i+1, end=' ')
                      #print('\n__',i+1,vb_orb_indices[i] )
                      for k in vb_orb_indices[i]:
                                #print(list_oa.index(k)+1,end=' ')
                                continue
                      #print('--',vb_orb_indices[list_oa[i]])   
                 else : # i est dans les MO's a garder
                      new_MOs.append(vect[i])
                      #print('\n##',i+1, vb_orb_coeffs[i])
                      #print(i+1, end=', ')
            # now add the rest of the MOs to the VEC
            for i in range(len(vb_orb_coeffs)):
                 if i+1 not in list_om : # i est forcement + petit que fin car relu 
                    print(' complement',i+1, end=' ')
                    new_MOs.append(vect[i])
           
            print('| ')
#            print('new_coeffs',new_coeffs[2])
            fin=len(new_coeffs)
            finmos=len(new_MOs)
            xmvb_input_file_name, xmvb_input_file_ext = os.path.splitext(xmvb_input_file)
            xmvb_output_file=xmvb_input_file_name+'.orbb'
            gamess_output_file=xmvb_input_file_name+'.vecc'
            print('#|-   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   ------')
            print('#| writes the',finmos,'VECs to', gamess_output_file)
            routines.wwrite_vec(gamess_output_file,new_MOs,0,finmos) 
            print('#| and those ',compteorbvb,'VB orbs to ', xmvb_output_file)
            routines.write_orbs(xmvb_output_file,new_coeffs,0,fin) 
            print("#+-  -----------------------------------------------------------------------")
            ao_orb, coeff_orb = routines.make_orb(new_coeffs,indices)
            routines.make_dollarorb(ao_orb,fin)
            print("#+- - ----------------------------------------------------------------------")
#            print("|  ")
#            print("|    ____  _  __                _____ ______ _________      ________ _____ ")
#            print("|   / __ \| |/ /               / ____|  ____|__   __\ \    / /  ____/ ____|")
#            print("|  | |  | | : /               | |  __| |__     | |   \ \  / /| |__ | |     ")
#            print("|  | |  | |  <                | | |_ |  __|    | |    \ \/ / |  __|| |     ")
#            print("|  | |__| | : \               | |__| | |____   | |     \  /  | |___| |____ ")
#            print("|   \____/|_|\_\               \_____|______|  |_|      \/   |______\_____|")
#            print("|                                                                          ")
#            print("|                                                                          ")
#            print("+- -    -------------------------------------------------------------------")
 #           print('finmos:',finmos,'new_MOs',new_MOs[1])
            #new_MOs_array=routines.to_array(new_MOs)
#            routines.wwrite_vec('screen',new_MOs,1,finmos) 
    else:
            print('+-  -----------------------------------------------------------------------')
            print('| Usage: getvec.py file_w_VEC [file_w_bfi]')
            print('|  -> file_w_VEC (.dat or .inp) $VEC is read and put in file_w_VEC.orbb ')
            print('|                                               ')
            print('|  if [file_bfi]== "orb"  then make the $orb section ')
            print('|                                               ')
            print('| with [file_w_bfi] the VEC is splitted between ')
            print('|  -> sigma ($VEC written in file_w_bfi.vecc) to use un gamess')
            print('|  -> pi (orb written in file_w_bfi.orbb, to use with xmvb') 
            print('+-----------------------------------------------------------------------')
            print('+ needs routines.py -----------------------------------------------------')
            print('+-----------------------------------------------------------------------')
            sys.exit() 

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
