#!/usr/bin/env python3
# write a code that diagonalise a matrix with an appropraite numpy routine
# take the MCSCF from 
## $MCSCF CISTEP=GUGA $END
## $DRT    GROUP=CS  FORS=.TRUE. NMCC=13 NDOC=2 NVAL=2 $END
## $GUGDIA NSTATE=2 $END
## $GUGDM2 WSTATE(1)=1.0,0.0 $END

import os
import cclib
import sys
import re
import routines
import numpy as np


def print_lin_matrix(title,matrix):
    print('- ',title,end=' ')
    print('- mtrx printed in lines ------------')
    for i in range(len(matrix)):
        print(title,'[',f"{i:4d}",']',end="\t")
        for j in range(len(matrix[i])):
            print(f"{matrix[i][j]:7.4f}", end="\t")
            if (j + 1) % 5 == 0:
                print()
        print()

def print_matrix(title,matrix):
    print('- ',title,end=' ')
    print('- mtrx printed in lines ------------')
    for i in range(len(matrix)):
        print(title,'[',f"{i:4d}",']',end="\t")
        for j in range(len(matrix[i])):
            print(f"{matrix[i][j]:7.4f}", end="\t")
            if (j + 1) % 5 == 0:
                print()
        print()


def CS_weight(i,Ci,Sij):
    wi=Ci[i]*np.dot(Sij,Ci)[i]
#    print('CS_weight(',i,')=',wi)
    return wi
def Offset_conf(conf,offset):
    # offset the conf by the largest MO in VB conf
    # cares of the : in the confs
    # usefull to detect the max/min of the orbitals.
    k=0
    l=0
    inidi =[]
    new_indices=[]
    while True:
        if k >= len(conf):
            break
        else:
            l=0
        indices=re.split(' +|\n',conf[k])  
        #print(len(indices),'indices',indices,end=' ')
        llist=''
        while True:
            if l >= len(indices):
                break
            else:
                if indices[l] != '':
                    if ':' in  indices[l] :
                        indic=re.split(':',indices[l])  
                        debut=int(indic[0])+offset
                        fin=int(indic[1])+offset
                       # print('debut',debut,'fin',fin)
                        inidi.append(str(debut)+':' +str(fin)  )
                        llist+=str(debut)+':' +str(fin)+' '
                    else :
                        inidi.append(str(int(indices[l])+offset))
                        llist+=str(int(indices[l])+offset)+' '
            #print(inidi)
            l+=1
        new_indices.append(llist)
        k+=1
    #print('new_indices',new_indices)
    return new_indices

def Get_CIVECT(CAS_file_name, state):
    offset=4
    if state == -1:  #  xmo file
        offset=2   
        pos_CASvect,line=routines.detect_keyword(CAS_file_name,"******  COEFFICIENTS OF STRUCTURES", 0)
    else:            # log file
        pos_CASvect,line=routines.detect_keyword(CAS_file_name, "LAGRANGIAN CONVERGED", 0)
        search_state=" STATE #    "+str(state)
        pos_CASvect,line=routines.detect_keyword(CAS_file_name, search_state, pos_CASvect)
    if (pos_CASvect == -1):
        print("Get_CIVECT  no converged Ci coeffs in  ",CAS_file_name)
        print("for state :",state)
        quit(   )
    pos_CASvectfin = routines.detect_blank(CAS_file_name,pos_CASvect+offset)
    CI_SIZE=pos_CASvectfin-1-(pos_CASvect+offset-1)
    # on decale de 3 lignes pour .log et de 2 lignes pour xmo 
    CI_conf, CI_vect=Read_CIVECT(CAS_file_name, pos_CASvect+offset-1,pos_CASvectfin-1,offset)
    return CI_conf,CI_vect

def est_reel(chaine):
    return '.' in chaine

 # check the lenght & size
def check_size(file_name,pos,fin):
    with open(file_name, 'r') as file:
        k=0 # line counter
        noa=0
        nom=0
        NOA=-1  
        NOM=0
        #print('--')
        while k < pos: # skip the first pos lines
            file.readline() 
            k+=1
        newblock=False
        while k < fin: # will stop 
            line = file.readline() 
            if len(line.split()) >=1:
                noa+=1
                newblock=True
                nom=len(line.split())-1
                #print('',len(line.split()),end=' ')
            elif newblock and len(line.split()) == 0:
                if NOA==-1:
                    NOA=noa-1
                    noa=0
#                    print('NOA=',NOA,end=' ')
                NOM=NOM+nom
#                print('size=',NOA,'x',NOM,end=' ')
                newblock=False
            k+=1
#        print('NOM=',NOM,end=' ')
    file.close()
    return NOA,NOM

def read5cols(file_name,length,size,pos,fin):# read a file with 3 blank lines, +1 to skip, then blocks of columns of length lines
                                         #  reading starts at pos. end reading when size vectors are read    
                                         # pos  A                 ******  OVERLAP OF VB STRUCTURES  ******
                                         # .    B
                                         #      C
                                         #      D           1            2            3            4            5
                                         #       1       1.000000     0.283693     0.222389    -0.370804    -0.341846
                                         #       2       0.283693     1.000000     0.283693    -0.149626    -0.132904
                                         #       3       0.222389     0.283693     1.000000    -0.370806    -0.122930
                                         #       4      -0.370804    -0.149626    -0.370806     1.000000     0.262162
    if size==-1:
        # should be able to read undetermined size and length
        # and the end of the section 
        size=length
        #continue
  #  nblock= (fin-pos -5)/length
  #  print('nblock= environ',nblock, (fin-pos -5-3*(nblock-1))/length, (fin-pos -5-3*(nblock-1))//length)
  #  print('size= environ',  (fin-pos -5-3*(nblock-1))/nblock)
    Stot=np.zeros((length,size))
    ioa=0
    nblock=0
    line='blanck'
    space=0
    firstom=0
    lenlue=0
    with open(file_name, 'r') as file:
        k=0 # line counter
        print('--')
        while k < pos: # skip the first pos lines
            file.readline() 
            k+=1
        while k < fin: # will stop 
            line = file.readline() 
            k+=1
            if len(line.split()) == 0: #line is empty=space between block
                if space==0: 
    #                print('new block starts, nblock',nblock,end=' ')
                    firstom=firstom+lenlue
                    ioa=0
                line = file.readline() 
                k+=1
                space+=1
            else: 
                if space !=0: # new block
                    space=0
                    nblock+=1
                    line = file.readline() 
                    k+=1
                lenlue=len(line.split() )-1
                for i, val in enumerate(line.split()[1:]):
                    iom=i+firstom
    #                print('S(',iom,',',ioa,')','=',val,end=' ')
                    Stot[iom][ioa] = float(val)
                ioa+=1
        print("| ------",iom+1,' columns have been read')
        file.close()
        return Stot

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
    f.close()
    return CI_conf, CI_vect
def ls_dir(beginning,end,k):
    files = [file for file in os.listdir() if file.endswith(end) and file .startswith(beginning)]
    for file in files:
        print(file, end=' ')
        k+=1
        if k%5 == 0:
            print() 
    return  k
def collect_confs(CI_conf):
    # collect the confs from the CI_conf in a single table of integers
    # cares of the : in the confs
    # usefull to detect the max/min of the orbitals.
    k=0
    tab=[]
    ttab=[]
    while True:
        if k >= len(CI_conf):
            break
        else:
            l=0
            ttab=re.split(' +|:|\n',CI_conf[k])
      #      print(len(ttab),'ttab',ttab,end='M ')
            while True:
                if l >= len(ttab):
                    break
                else:
                    if ttab[l] != '':
                        tab.append(int(ttab[l]))
                    l+=1
            k+=1
    return tab

# =========------------------------------------------------------
# == main =------------------------------------------------------
# =========------------------------------------------------------
print("+--------project.py - SH 2024 ---------------------")
print('| project.py file_CI ',len(sys.argv), 'arguments')
print("+--------                     ---------------------")
if len(sys.argv) <= 1:
    print('| file_CI can be either a .log or a .xmo file')
    print('| usually file_xm.xmo + a file_xm_vb to make the _vbp')
    
input_file = sys.argv[1]
if not os.path.exists(input_file):
    print('ls ',input_file,'*.log or *.xmo')
    k=ls_dir(input_file,'.log',0)
    k=ls_dir(input_file,'.xmo',k)
    input_file = input("Enter the file name: ")
input_file_name, input_file_ext = os.path.splitext(input_file)
# Le calcul CAS est dans nom.log ?
if input_file_ext == ".log":
    reponse=input("Enter the state number: [1] ")
    try:
        state = int(reponse)    
    except:
        state = 1

else:
    state = -1  # xmo file only 1 state
VB_inp_file = input_file_name+'_vb.inp'
OVERL_inp_file = input_file_name+'_vbp.inp'
must_write_OVERL=False  
if not os.path.exists(OVERL_inp_file): # if that one exists the xmo should also be done
    must_write_OVERL=True
    if os.path.exists(VB_inp_file):
       command='cp '+VB_inp_file+' '+OVERL_inp_file 
       print('|  ',command,' the inp file for the OVERLAP calculation is mandatory')
       os.system(command)
    else:
        print('Error : no file ',VB_inp_file,' to copy')
OVERL_output_file = input_file_name+'_vbp.xmo'
OVERL_input_file = input_file_name+'_vbp.xmii' #for .xmi, to avoid erasing
OVERL_file_orb=input_file_name+'_vbp.orbb'     #for .orb, to avoid erasing
    
print('|  read files :\n| ',input_file,end=':')
CI_conf,CI_vect=Get_CIVECT(input_file, state)
lenCI=len(CI_vect)
print('',lenCI,end=' CI vect, ')
if len(sys.argv) >= 3:
    VB_file  = sys.argv[2]
else:
    if len(sys.argv)==2  :
        print('You need to build the _xm.* files')
        tableau=collect_confs(CI_conf)
#        print('apres collect_conf',tableau) 
        toprint=routines.make_conf_from_gamess(CI_conf)
        print('Make a ',input_file_name,'_xm.xmi file     with ')
        print("$ctrl \n vbftyp=det WFNTYP=struc iscf=3 itmax=1000 bprep nstate=0 iprint=3 guess=read \n nmul=1 nstr=",lenCI,"\n  $end")
        print("$str ")
        for ii in range(len(toprint)):
            print(toprint[ii],' ; ', CI_conf[ii],ii+1,'... ',CI_vect[ii])    
        print("$end \n")
        if os.path.exists(input_file_name+'.dat'):
           print('now type :  getvec.py ',input_file_name+'.dat \n \n')
    VB_file  = input_file_name+"_vb.xmo"
print(VB_file,end=':')
try:
    NVBCONF= routines.Read_INTEGER(VB_file, " Number of Structures:",12)  
except:   
    print('no valid file ',VB_file)
    NVBCONF=len(CI_conf)
    quit()
VB_conf,VB_vect=Get_CIVECT(VB_file, state)
if (len(VB_vect) != NVBCONF):
    print('Error : should read ',NVBCONF,' in ', VB_file,' but read ', len(VB_file))
else:
    print('',len(VB_vect),end='VB vect, ')
file_orb_VB=input_file_name+'_vb.orb'
print(file_orb_VB,end=':')
VB_orb_coeffs,VB_orb_aos=routines.read_orb(file_orb_VB) 
print(len(VB_orb_coeffs),end=' VB orbs, ')
file_orb_CI=input_file_name+'.orb'
print(file_orb_CI,end=':')
CI_orb_coeffs,CI_orb_aos=routines.read_orb(file_orb_CI)
print(len(CI_orb_coeffs),end=' CI orbs')
print()
if must_write_OVERL:
# offset the MCSCF conf by the largest MO in VB conf
    OFFSET=max(collect_confs(VB_conf))  
    print('|  Largest VB orb number, ',OFFSET,', is used as offset for the ',min(collect_confs(CI_conf)),'-', max(collect_confs(CI_conf)) ,' CI orbitals      ')
    print('|  hence CI orb are now numbered from ',OFFSET+min(collect_confs(CI_conf)),' to ',OFFSET+max(collect_confs(CI_conf)),' and written in ',OVERL_file_orb )
    dec_CI_conf=Offset_conf(CI_conf,OFFSET)
    print('-------------------------------------------------')
    with open(VB_file , "r") as file:
        for line in file:
            MULT=1
            if routines.is_in(line, "MULTIPLICITY"):
               MULT= routines.Read_INT(line, "MULTIPLICITY")
               print('| Multiplicity is ',MULT)
    print('|  filling the ', OVERL_input_file,' xmi input file ' )
    with open(OVERL_input_file,'a') as file:
        line='made by project.py \n $ctrl   ; ============='+str(NVBCONF)+' + '+str(lenCI)+'=============.======+\n'
        line=line+'  vbftyp=det WFNTYP=struc iscf=3 itmax=10' 
        line=line+'   iprint=-1, nmul='+str(MULT)+' nstr='+str(NVBCONF+lenCI)+' guess=read \n $end \n'
        file.write(line)    
        line= ' $struc  ; ============= \n '
        file.write(line)   
        #ttt=collect_confs(VB_conf)
        for i in range(len(VB_conf)):
            line= '  '+ VB_conf[i]+'    ; VB '+str(i+1) + str( VB_vect[i])+'\n'  
            #print(*['%4.0f' % int(val) for val in VB_conf[i].split()],end=' ')
            #print(line)
            file.write(line)   
        for i in range(lenCI):
            line= '  '+ dec_CI_conf[i]+'    ; CI '+str(i+1) + str( CI_vect[i])+'\n'  
            #print(*['%4.0f' % int(val) for val in dec_CI_conf[i].split()],end=' ')
            #print('  ',dec_CI_conf[i],'       ; CI ',i+1 , CI_vect[i])
            file.write(line)   
        line= ' $end  ; ============= \n '
        file.write(line)   
        pos, line=routines.detect_keyword(VB_file, "$bfi", 0)
        bfi_nom,bfi_noa,list_om,list_oa=routines.read_bfi(VB_file,pos)
        line= ' $bfi  ; ============= \n '
        line=line+' '+str(bfi_nom)+' '+str(bfi_noa)+ '\n   '+routines.makeSTR(list_om)+'\n   '+routines.makeSTR(list_oa)+'\n $end \n'
        file.write(line)   
    OVERL_coeffs=[]
    OVERL_aos=[]
    for k in range(len(VB_orb_coeffs)):
   #     print('VBcoeffs', k, len(VB_orb_coeffs[k]),end=' ')
        OVERL_coeffs.append(VB_orb_coeffs[k])
        OVERL_aos.append(VB_orb_aos[k])
    for k in range(len(CI_orb_coeffs)):
   #     print('ICcoeffs', k, len(VB_orb_coeffs[k]),end=' ')
        OVERL_coeffs.append(CI_orb_coeffs[k])
        OVERL_aos.append(CI_orb_aos[k])
    routines.write_orb(OVERL_file_orb,OVERL_coeffs,OVERL_aos,0,len(OVERL_coeffs))
    print('  ',OVERL_file_orb,' written')   
    
    routines.write_DOLLARORB(OVERL_input_file,OVERL_aos,0,len(OVERL_aos))
    print('| $orb section, with ',len(OVERL_aos),' orbitals is to get in ORBB    ')
    print('| $end ')

OVERL_output_file = input_file_name+'_vbp.xmo'
if not os.path.exists(OVERL_output_file):
    print('- submit the calculation:',OVERL_output_file, ' is required to continue')
    quit()

OVERLP_size= routines.Read_INTEGER(OVERL_output_file, " Number of Structures:",12)  
if OVERLP_size != NVBCONF+lenCI           :
    print('Error in OVERLP_size',OVERLP_size,'instead of ',NVBCONF+lenCI           )
    quit()
else:
    print('OVERLP_size',OVERLP_size,' is OK with VB+CI ')
posit,line=routines.detect_keyword(OVERL_output_file,"******  OVERLAP OF VB STRUCTURES ",0)
posfin,line=routines.detect_keyword(OVERL_output_file,"******  HAMILTONIAN ",posit)
# read  The  Matrix at pos (3 lines between each blocks)
# 5 intergers first (5 cols) then each line as integer (line) then 5 reals 
#
n,m=check_size(OVERL_output_file,posit+1,posfin)
Stot=[]
print('| In ', OVERL_output_file, ' the overlap matrix is (',n,'x',m,')')
Stot=read5cols(OVERL_output_file,OVERLP_size,OVERLP_size,posit,posfin) # read a file with 3 blank lines, +1 to skip, then blocks of columns of length lines
#print_matrix('Stot',Stot)
#print('CI_vect[1]',CI_vect[1])
Svbvb=np.zeros((NVBCONF,NVBCONF))
part_SOM=np.zeros((NVBCONF,len(CI_vect)))
SOMvb=np.zeros(NVBCONF)
print('| ' )
print('| ------------------------------------------------------------')
print('|  Get the overlaps between each of the CI\'s CSF of ',CI_vect,'with each of the ',NVBCONF,'VB conf')
for ivb in range(len(VB_conf)):
    for jvb in range(len(VB_vect)):
        Svbvb[ivb][jvb]=Stot[ivb][jvb]
    for jmo in range(len(CI_vect)):
        part_SOM[ivb][jmo]=Stot[ivb][jmo+NVBCONF] # part_SOM is the part that concerns the overlap between each MO configurations of the CI and each VB configurations.
#print_lin_matrix('S_CI/VB^T',part_SOM.T)
print('| ' )
#print('| ------------------------------------------------------------')
#print('| The vector of the overlap between ',NVBCONF,' VB CSF\' (structure) ')
#print_lin_matrix('S_VB/VB',Svbvb)
print('| ' )
#print('| ------------------------------------------------------------')
#print('| The CI vector of this state ',CI_vect,' has the following overlaps with each VB\'s structures')
SOMvb=np.dot(part_SOM,CI_vect)
#SOMvb=SOMvb.T
#print('SOMvb',SOMvb)    
sol=np.linalg.solve(Svbvb,SOMvb)
Svbvbsol=np.dot(Svbvb,sol)  
norm_sol=np.dot(sol.T,Svbvbsol)  
print('| projected  : norme before=',f"{norm_sol:4.3f}")
sol = sol / np.sqrt(norm_sol)
norm_sol = np.dot(sol.T, np.dot(Svbvb, sol))
#print("|  norm sol after............",f"{norm_sol:4.3f}")
#print('vecteur sol', sol)
print('-----')
w=[]
for i in range(len(sol)):
    w.append(CS_weight(i, sol, Svbvb))
#    print('CS_weight(',i,')=',CS_weight(i, sol, Svbvb))

trust=np.dot(SOMvb,sol)
trust2=np.dot(SOMvb,VB_vect)
print('|   Solution of the projection (C[i] are ab initio, while K[i] are projected) : trust/C[i]//K[i]',f"{trust*100:7.2f}/{trust2*100:7.2f}",'%')
print('|  VB[i]     C[i]   k[i]     w[i]        <VB i| ψ_CI>   ')
for i in range(len(sol)):
    print('|  ',f"{i+1:3d}",'',f"{VB_vect[i]:7.3f}",f"{sol[i]:7.3f}",f"{w[i]*100:7.2f}",'%     ',f"{SOMvb[i]:7.3f}",'  ')

#print(f"{matrix[i][j]:7.4f}", end="\t")
yes=input('Do you want to print the CI/VB overlaps ? (y/n)')
if yes != 'n':
    print('|  VB[i]     C[i]   ',end=' ')   
    for i in range(len(CI_vect)):
        print('   CI',f"{i:2d}",end='  ')
    print()
    for i in range(len(sol)):
        print('|  ',f"{i+1:3d}",'  ',f"{VB_vect[i]:7.3f}",end=' ')
        for j in range(len(CI_vect)):
            print(' ',f"{part_SOM[i][j]:7.4f}",end=' ') 
        print()

quit()
print(f"{matrix[i][j]:7.4f}", end="\t")
   # teste ' a la main ' sur le resultat HULIS de l'acrolein:
   # a faire gaffe que HuLiS range l'OM puis VB alors que ici c'est l'inverse
   #          	ψtot	   ψ1   	ψ2	   ψ3
   #         ψt	1.00000	0.87683	0.44133	0.10848
   #         ψ1	0.87683	1.00000	0.17701	0.07299
   #         ψ2	0.44133	0.17701	1.00000	0.00000
   #         ψ3	0.10848	0.07299	0.00000	1.00000
   #         
   #         
   #         2- Eigenvectors (Ci) and weights (Wi) of the mesomery
   #         Method			HLP
   #         Trust factor			92.50319%
   #         ψi	      Ei	      Ci	 WiHLP
   #         ψ1	4α + 5.30137β	0.88741	84.11657%
   #         ψ2	4α + 3.94000β	0.32002	15.26774%
   #         ψ3	4α + 2.00000β	0.05250	0.61569%
   #         

print('-----')
print('a la main@@')
Svbvb = np.array([[1., .177, 0.07299], [0.177, 1., 0.0], [0.07299, 0.00, 1.]])
SOMvb=np.array([[.8768], [0.441327], [0.108482]])
normed=np.array([[.0000], [0.0000], [1.0]])
sol=np.linalg.solve(Svbvb,SOMvb)
Svbvbsol=np.dot(Svbvb,sol)  
norm_sol=np.dot(sol.T,Svbvbsol)  
print('vecteur sol', sol)
print("norme sol",norm_sol)
print('-----')
quit()
for j in range(len(VB_conf)):
    print_matrix('SOMVb',SOMvb)
print_matrix('SOMvb',SOMvb)
quit()
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
