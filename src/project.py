#!/usr/bin/env python3
# write a code that diagonalise a matrix with an appropraite numpy routine
# take the MCSCF from 
## $MCSCF CISTEP=GUGA $END
## $DRT    GROUP=CS  FORS=.TRUE. NMCC=13 NDOC=2 NVAL=2 $END
## $GUGDIA NSTATE=2 $END
## $GUGDM2 WSTATE(1)=1.0,0.0 $END

import os
import sys
import re
import routines
import numpy as np
def CS_weight(i,Ci,Sij):
    wi=Ci[i]*np.dot(Sij,Ci)[i]
#    print('CS_weight(',i,')=',wi)
    return wi

def Get_CIVECT(CAS_file_name, state):
    offset=4
    if state == -1:
        offset=2   
        print("xmo file")
        pos_CASvect,line=routines.detect_keyword(CAS_file_name,"******  COEFFICIENTS OF STRUCTURES", 0)
    else:
        pos_CASvect,line=routines.detect_keyword(CAS_file_name, "LAGRANGIAN CONVERGED", 0)
        search_state=" STATE #    "+str(state)
        pos_CASvect,line=routines.detect_keyword(CAS_file_name, search_state, pos_CASvect)
    if (pos_CASvect == -1):
        print("Get_CIVECT  no converged Ci coeffs in  ",CAS_file_name)
        print("for state :",state)
        quit(   )
    pos_CASvectfin = routines.detect_blank(CAS_file_name,pos_CASvect+offset)
    CI_SIZE=pos_CASvectfin-1-(pos_CASvect+offset-1)
#    print('pos_CASvect',pos_CASvect+offset, CI_SIZE, pos_CASvectfin-1)
    # on decale de 3 lignes pour .log et de 2 lignes pour xmo 
    CI_conf, CI_vect=Read_CIVECT(CAS_file_name, pos_CASvect+offset-1,pos_CASvectfin-1,offset)
#    print('GET_CIVEC CI_vect:',CI_vect)
    return CI_conf,CI_vect

def Read_CIVECT(file, pos,fin,offset):
    #offset=4 pour .log et 2 pour .xmo
    # indique la position a lire
    with open(file, 'r') as f:
        lines = f.readlines()[pos:fin]
        n=fin-pos
#        print('n',pos,n,fin)
        CI_vect = np.zeros(n)
        CI_conf = []
        for i in range(n):
# Det,       CI_vect[i] = float(re.split('\|+| +|\n', lines[i])[7])
# GUGA
            CI_vect[i] = float(re.split(' +',lines[i])[2])
            if offset == 4: # log
                str=re.split(' +|\n',lines[i])[3]
                CI_conf.append(str)
            else: # xmo
                #CI_conf[i] = re.split(' +|\n',lines[i])[3:]
                str=''
                k=4
                while True:
                    if re.split(' +|\n',lines[i])[k]!= '':
                        str+=' '+re.split(' +|\n',lines[i])[k]   
#                        print(str,end=' ')
                    else:
                        break
                    k+=1
                CI_conf.append(str)
                continue
    return CI_conf, CI_vect

print("+--------project.py - SH 2024 ---------------------")
print('| project.py file_CI ')
if len(sys.argv) <= 1:
    print('| file_CI can be either a .log or a .xmo file')
    input_file = "proj_xm.xmo"
else: 
    input_file = sys.argv[1]
CAS_file_name = input_file
if not os.path.exists(CAS_file_name):
    log_files = [file for file in os.listdir() if file.endswith(".log")]
    xmo_files = [file for file in os.listdir() if file.endswith(".xmo")]
    print("Available .log .xmo files:")
    for file in xmo_files:
        print(file, end=' ')
    for file in log_files:
        print(file, end=' ')
    print()
    input_file = input("Enter the file name: ")
input_file_name, input_file_ext = os.path.splitext(input_file)
#    CAS_file_name = input_file_name + ".log"
#input_file_name, input_file_ext = os.path.splitext(input_file)
###        pos=pos+2
###        posfin=(routines.detect_blank(input_file,pos))
###        print('pos & posfin',pos,posfin, posfin-(pos)), 
###        with open(input_file, 'r') as f:
###            lines = f.readlines()[pos-1:posfin]
###            n=posfin-pos
###            print('n',pos,n,posfin,lines)
###            CI_vect = np.zeros(n)
###            CI_conf = np.zeros(n)
###            for i in range(n):
###                print(i,'||',re.split(' +|\n',lines[i])[2])
###                CI_vect[i]=float(re.split(' +',lines[i])[2])
###                CI_conf[i] = re.split(' +|\n',lines[i])[4]
###            print('CI_vect:',CI_vect)

# Le calcul CAS est dans nom.log
if input_file_ext == ".log":
    state = int(input("Enter the state number: "))
else:
    state = -1  # xmo file only 1 state
# Define your matrix
CI_conf,CI_vect=Get_CIVECT(CAS_file_name, state)
print('CI_vect:',CI_vect,CI_conf)
CAS_file_namevb = input_file_name+"_vb"+input_file_ext
CI_confvb,CI_vectvb=Get_CIVECT(CAS_file_namevb, state)
print('CI_vect VB:',CI_vectvb,CI_confvb)
print(len(CI_vect),len(CI_vectvb))
k=0
ttab=[]
tab=[]
while True:
    if k >= len(CI_confvb):
        break
    else:
        l=0
        ttab=re.split(' +|:|\n',CI_confvb[k])
        print(len(ttab),'ttab',ttab,end=' ')
        while True:
            #print('\n ttab[',k,']=',end=' ')
            if l >= len(ttab):
                break
            else:
                if ttab[l] != '':
                    #print(ttab[l],end=' ')
                    tab.append(int(ttab[l]))
                l+=1
        k+=1
print('-----------------------')
print(tab,'max=',max(tab))
# offset the MCSCF conf by max(tab) to write the conf
OFFSET=max(tab)
k=0
ttab=[]
tab=[]
inidi =[]
new_indices=[]
while True:
    if k >= len(CI_conf):
        break
    else:
        l=0
        indices=re.split(' +|\n',CI_conf[k])  
        print(len(indices),'indices',indices,end=' ')
        llist=''
        while True:
            if l >= len(indices):
                break
            else:
                if indices[l] != '':
                    if ':' in  indices[l] :
                        indic=re.split(':',indices[l])  
                        debut=int(indic[0])+OFFSET
                        fin=int(indic[1])+OFFSET
                        print('debut',debut,'fin',fin)
                        inidi.append(str(debut)+':' +str(fin)  )
                        llist+=str(debut)+':' +str(fin)+' '
                    else :
                        inidi.append(str(int(indices[l])+OFFSET))
                        llist+=str(int(indices[l])+OFFSET)+' '
            print(inidi)
#                    tab.append(int(indices[l]))
            l+=1
        new_indices.append(llist)
#        CI_conf[k]=(int(CI_conf[k])+max(tab))
        k+=1
    print('new_indices',new_indices)
quit()
val=int(re.split(' +|:|\n',CI_confvb[k]))
print(re.split(' +|\n',CI_confvb[len(CI_vectvb)-1])[1])
###state_number=1
###search_state=" STATE #    "+str(state_number)
###pos_CASvect=routines.detect_keyword(CAS_file_name, "LAGRANGIAN CONVERGED", 0)
###pos_CASvect=routines.detect_keyword(CAS_file_name, search_state, pos_CASvect)
###pos_CASvectfin = routines.detect_blank(CAS_file_name,pos_CASvect + 3)
###
###CI_SIZE=pos_CASvectfin-1-(pos_CASvect+3)
###print('pos_CASvect',pos_CASvect, CI_SIZE)
#### on decale de 2 ligne
###CI_vect=Read_CIVECT(CAS_file_name, pos_CASvect+3,pos_CASvectfin-1)
# by convention the CAS comes last
#                 ******  OVERLAP OF VB STRUCTURES  ******  
#
#
#                  1            2            3            4            5
#       1       1.000000    -0.572037    -0.571954    -0.980288     0.155406
#       2      -0.572037     1.000000     0.673136     0.502923     0.000015
#       3      -0.571954     0.673136     1.000000     0.502812     0.000016
#       4      -0.980288     0.502923     0.502812     1.000000     0.000000
#       5       0.155406     0.000015     0.000016     0.000000     1.000000
#
# Remove the line "stop"
# Rest of the code...
NFZC= routines.Read_INTEGER(CAS_file_name, "NFZC=",4)  
NDOC= routines.Read_INTEGER(CAS_file_name, "NDOC=",4)  
NVAL= routines.Read_INTEGER(CAS_file_name, "NVAL=",4)  
NMCC= routines.Read_INTEGER(CAS_file_name, "NMCC=",4)  
numel= routines.Read_INTEGER(CAS_file_name, "NUMBER OF ELECTRONS                          =",5)  
print('-#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#__#_#_')
print('NFZC',NFZC,'NDOC',NDOC,'NVAL',NVAL,end=', ')
print('NMCC',NMCC,end=', ')
print('numel',numel)    
quit()
Svbvb = np.array([[1., .177, 0.07299], [0.177, 1., 0.0], [0.07299, 0.00, 1.]])
SOMvb=np.array([[.8768], [0.441327], [0.108482]])
normed=np.array([[.0000], [0.0000], [1.0]])
sol=np.linalg.solve(Svbvb,SOMvb)
#sigma=np.array([[1.0], [1.0]])
#ss=np.dot(sigma.T,sigma)
#print('<sigma|sigma>',ss)
#write the matrix
Svbvbsol=np.dot(Svbvb,sol)  
norm_sol=np.dot(sol.T,Svbvbsol)  
print('vecteur sol', sol)
print("norme sol",norm_sol)
print('-----')
sol = sol / np.sqrt(norm_sol)
norm_sol = np.dot(sol.T, np.dot(Svbvb, sol))
print('vecteur sol', sol)
print("norme sol", norm_sol)
print('-----')
wcs = np.zeros(len(sol))  # Initialize the variable "wcs" with zeros
print('wcs',     wcs)
print('sol et Svbvb',sol, Svbvb)
for i in range(0, len(sol)):
    wcs[i]=CS_weight(i, sol, Svbvb)
    print('wcs_', i, '=', round(wcs[i], 3))
#   w[i]=w[i]+Svbvb[i][j]*sol[j]
#print('w_',i,'=',w[i])
quit() 
numel= routines.Read_INTEGER(CAS_file_name, "NUMBER OF ELECTRONS")  
# Rest of the code...
Svbvb = np.array([[1., .177, 0.07299], [0.177, 1., 0.0], [0.07299, 0.00, 1.]])
SOMvb=np.array([[.8768], [0.441327], [0.108482]])
normed=np.array([[.0000], [0.0000], [1.0]])
sol=np.linalg.solve(Svbvb,SOMvb)
#sigma=np.array([[1.0], [1.0]])
#ss=np.dot(sigma.T,sigma)
#print('<sigma|sigma>',ss)
#write the matrix
Svbvbsol=np.dot(Svbvb,sol)  
norm_sol=np.dot(sol.T,Svbvbsol)  
print('vecteur sol', sol)
print("norme sol",norm_sol)
print('-----')
sol = sol / np.sqrt(norm_sol)
norm_sol = np.dot(sol.T, np.dot(Svbvb, sol))
print('vecteur sol', sol)
print("norme sol", norm_sol)
print('-----')
wcs = np.zeros(len(sol))  # Initialize the variable "wcs" with zeros
print('wcs',     wcs)
print('sol et Svbvb',sol, Svbvb)
for i in range(0, len(sol)):
    wcs[i]=CS_weight(i, sol, Svbvb)
    print('wcs_', i, '=', round(wcs[i], 3))
#   w[i]=w[i]+Svbvb[i][j]*sol[j]
#print('w_',i,'=',w[i])
