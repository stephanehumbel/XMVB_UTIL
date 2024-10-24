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
    i=0
    tampon=[]
    if offset>=1:   
        str_ofset="1:"+str(offset)+'   '
        #print('str_ofset',str_ofset )
    else:
       str_ofset =" "
    while i < len(new_indices):
        tampon.append(str_ofset+new_indices[i])
        i+=1
    #print('new_indices',new_indices)
    return tampon

def Get_CIVECT(CAS_file_name, state):
    offset=4
    if state == -1:  #  xmo file
        offset=2   
        pos_CASvect,line=routines.detect_keyword(CAS_file_name,"******  COEFFICIENTS OF STRUCTURES", 0)
    else:            # log file
        pos_CASvect,line=routines.detect_keyword(CAS_file_name, "LAGRANGIAN CONVERGED", 0)
        search_state=""+str(state)+'  ENERGY'
        pos_CASvect,line=routines.detect_keyword(CAS_file_name, search_state, pos_CASvect)
    if (pos_CASvect == -1):
        print("Get_CIVECT  no converged Ci coeffs in  ",CAS_file_name)
        print("for state :",state)
        quit(   )
    pos_CASvectfin = routines.detect_blank(CAS_file_name,pos_CASvect+offset)
    CI_SIZE=pos_CASvectfin-1-(pos_CASvect+offset-1)
    # on decale de 3 lignes pour .log et de 2 lignes pour xmo 
    CAS_conf, CAS_vect=Read_CIVECT(CAS_file_name, pos_CASvect+offset-1,pos_CASvectfin-1,offset)
    return CAS_conf,CAS_vect

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
        CAS_vect = np.zeros(n)
        CAS_conf = []
        for i in range(n):
# Det,       CAS_vect[i] = float(re.split('\|+| +|\n', lines[i])[7])
# GUGA
            CAS_vect[i] = float(re.split(' +',lines[i])[2])
            if offset == 4: # log
                str=re.split(' +|\n',lines[i])[3]
                CAS_conf.append(str)
            else: # xmo
                #CAS_conf[i] = re.split(' +|\n',lines[i])[3:]
                str=''
                k=4
                while True:
                    if re.split(' +|\n',lines[i])[k]!= '':
                        str+=' '+re.split(' +|\n',lines[i])[k]   
#                        print(str,end=' ')
                    else:
                        break
                    k+=1
                CAS_conf.append(str)
                continue
    f.close()
    return CAS_conf, CAS_vect
def ls_dir(beginning,end,k):
    files = [file for file in os.listdir() if file.endswith(end) and file .startswith(beginning)]
    for file in files:
        print(file, end=' ')
        k+=1
        if k%5 == 0:
            print() 
    return  k
def collect_confs(CAS_conf):
    # collect the confs from the CAS_conf in a single table of integers
    # cares of the : in the confs
    # usefull to detect the max/min of the orbitals.
    k=0
    tab=[]
    ttab=[]
    while True:
        if k >= len(CAS_conf):
            break
        else:
            l=0
            ttab=re.split(' +|:|\n',CAS_conf[k])
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
print('| CAS_LEW.py file_CAS file_VB ',len(sys.argv), 'arguments')
print("+- must have .orb and .xmo of both              ---")
print("+--------                     ---------------------")
if len(sys.argv) == 1:
    print('| file_CI can be either a .log or a .xmo file')
    print('| usually file_xm.xmo + a file_xm_vb to make the _vbp')
    print('| at final stage cas.xmo + vb.xmo  and will search $cas+_overl.xmo ')
    quit()
if len(sys.argv) == 2:
#For an xmo file')
    CAS_file = sys.argv[1]
    CAS_file_name, CAS_file_ext = os.path.splitext(CAS_file)
#For a  GAMESS log file')
    if CAS_file_ext == '.log' and os.path.exists(CAS_file):
        print (CAS_file_name+'.log exists must be like $MCSCF CISTEP=GUGA MAXIT=200 QUAD=.F. $END $GUGDIA NSTATE=2 $END $GUGDM2 WSTATE(1)=1.0,0.0 $END')
        CAS_conf,CAS_vect=Get_CIVECT(CAS_file_name+".log",1)
        toprint=routines.make_conf_from_gamess(CAS_conf)
        num,line=routines.detect_keyword(CAS_file_name+".log", "SPIN MULTIPLICITY", 0)    
        MULT=routines.Read_INT(line,"MULTIPLICITY")
        num,line=routines.detect_keyword(CAS_file_name+".log", "NUMBER OF CORE MOL", 0)    
        offset=routines.Read_INT(line,"ORBITALS")
        num,line=routines.detect_keyword(CAS_file_name+".log", "  NMCC", 0)    
        nmcc=routines.Read_INT(line,"NMCC")
        num,line=routines.detect_keyword(CAS_file_name+".log", "  NDOC", 0)    
        ndoc=routines.Read_INT(line,"NDOC")
        num,line=routines.detect_keyword(CAS_file_name+".log", "  NALP", 0)    
        NALP=routines.Read_INT(line,"NALP")
        nae=ndoc*2+NALP
        num,line=routines.detect_keyword(CAS_file_name+".log", "              NVAL", 0)    
        nval=routines.Read_INT(line,"NVAL")
        nae=ndoc*2+NALP
        nao=ndoc+NALP+nval
#        num,line=routines.detect_keyword(CAS_file_name+".log", "TOTAL NUMBER OF ATOMS", 0)    
#        natoms=routines.Read_INT(line,"ATOMS")
        print('number of core Orbs=',offset,CAS_conf,CAS_vect)
        print('other numbers      =',nmcc,ndoc,nval)
        #get geom
        symbol,zat, x,y,z,natoms=routines.read_geom(CAS_file_name+".log")
        basis_set = input('give the basis set (6-31G):')
        if basis_set == '':
            basis_set='6-31G'
        if offset >= 0:
            CAS_conf=Offset_conf(toprint,offset)
        print('')
        print('from ',CAS_file_name+'.log ; ',basis_set)
        print('$ctrl \n str=full nao='+str(nao)+ '  nae='+str(nae)+'   nmul='+str(MULT),' orbtyp=hao frgtyp=sao int=libcint ; nstr='+str(len(CAS_conf)))
        print(' basis='+basis_set, ' iscf=5 iprint=3 guess=read itmax=0 \n $end')
        print('$frag \n ', natoms,'\n spxyzdxxyyzzxyxzyzz 1-'+str(natoms)+'\n $end') 
        norb=nmcc+nao
        print('$orb \n 1*'+str(norb))
        for i in range(norb):
            print(' 1')
        print('$end')
        print('$stru  ;',len(CAS_conf),' confs')
        routines.write_conf("screen",CAS_conf,CAS_vect)
        print('$end') 
        print('$geo   ; ',natoms,' atoms') 
        for i in range(natoms):
            print(f"   {symbol[i]:5}{x[i]:18.9f}{y[i]:18.9f}{z[i]:18.9f}")
        print('$end') 
        if os.path.exists(CAS_file_name+'.dat'):
            input_file=CAS_file_name+'.dat'
            log_file=CAS_file_name+'.log'
            # GET the MO's from dat
            coeffs=[]
            pos,line=routines.detect_keyword(input_file, "VEC", 0)
#            print('|  read files :\n| ',input_file,end=':')
#            print (pos+1)
            coeffs,nvect = routines.read_vec(input_file,coeffs,pos+1)
            vect=routines.make_table(coeffs)#
            norb=nmcc+nval+ndoc+NALP
            type_OA=routines.read_basis(log_file)
            reord_OA=[]
            #print(type_OA,reord_OA)
            reord_OA=routines.reorder_OA(type_OA)
            #print(reord_OA)
            new_vect=[]
            for j in range(len(vect)):
              new_orb=[]
              for i in range(len(reord_OA)):
                new_orb.append(vect[j][reord_OA[i]])
              new_vect.append(new_orb)
#            routines.write_orbs("screen",new_vect,norb-2,norb) 
            print('$gus') 
            routines.write_orbs("screen",new_vect,0,norb) 
            print('$end') 

            print('')
        else:
            print('  ', CAS_file_name+'.dat not found ')
        print('                 ---===---===:::========')
        quit()
    elif CAS_file_ext == '.xmo' and os.path.exists(CAS_file):
        print('|  read files :\n| ',CAS_file_name+".xmo",end=':')
        CAS_conf,CAS_vect=Get_CIVECT(CAS_file_name+".xmo", -1)
        lenCI=len(CAS_vect)
        print('',lenCI,end=' CI vect.')
        print('')    
        print('                 ---===---===:::========')    
        routines.write_conf("screen",CAS_conf,CAS_vect)
        routines.write_conf(".conf"  ,CAS_conf,CAS_vect)
        print('                 ---===---===:::========')    
    else:   
        print('### >>> ',CAS_file,'  not found <<<<< ##')
    quit()
#k=1    
if len(sys.argv) >= 3: 
    VB_file = sys.argv[1]
    CAS_file = sys.argv[2]
    VB_file_name, VB_file_ext = os.path.splitext(VB_file)
    CAS_file_name, CAS_file_ext = os.path.splitext(CAS_file)
    if not os.path.exists(VB_file_name+'.xmo'):
        ext_error='### >>> ',VB_file_name+'.xmo not found <<<<< ## '
        quit(ext_error) 
    if not os.path.exists(CAS_file_name+'.xmo'):
        ext_error='### >>> ',CAS_file_name+'.xmo not found <<<<< ## '
        quit(ext_error)
    state = -1  # xmo file only 1 state
    if len(sys.argv) ==4:
        must_write_OVERL=False  
        OVERL_file_xmo=sys.argv[3]
        print('|OVERLAP files :\n| ',OVERL_file_xmo,end='')
    if len(sys.argv) ==3:
        must_write_OVERL=True
        OVERL_file_orb=CAS_file_name+'_overl.orb'
        OVERL_file_xmi=CAS_file_name+'_overl.xmi'
        OVERL_file_xmo=CAS_file_name+'_overl.xmo'
# comes here only if more than 2 arguments
# reads the xmo files
print('|  read files :\n| ',CAS_file,end=':')
CAS_conf,CAS_vect=Get_CIVECT(CAS_file, state)
lenCI=len(CAS_vect)
print('',lenCI,end=' CI vect, ')
print(VB_file,end=':')
try:
    NVBCONF= routines.Read_INTEGER(VB_file, " Number of Structures:",12)  
except:   
    print('no valid file ',VB_file)
    NVBCONF=len(CAS_conf)
    quit()
VB_conf,VB_vect=Get_CIVECT(VB_file, state)
if (len(VB_vect) != NVBCONF):
    print('Error : should read ',NVBCONF,' in ', VB_file,' but read ', len(VB_file))
else:
    print('',len(VB_vect),end='VB vect, ')
file_orb_CAS=CAS_file_name+'.orb'
print(file_orb_CAS,end=':')
CAS_orb_coeffs,CAS_orb_aos=routines.read_orb(file_orb_CAS) 
print(len(CAS_orb_coeffs),end=' CAS orbs, ')
file_orb_VB=VB_file_name+'.orb'
print(file_orb_VB,end=':')
VB_orb_coeffs,VB_orb_aos=routines.read_orb(file_orb_VB)
print(len(VB_orb_coeffs),end=' VB orbs')
print()
##
##
if must_write_OVERL:
# offset the MCSCF conf by the largest MO in VB conf
    OFFSET=max(collect_confs(VB_conf))  
    print('|  Largest VB orb number, ',OFFSET,', is used as offset for the ',min(collect_confs(CAS_conf)),'-', max(collect_confs(CAS_conf)) ,' CI orbitals      ')
    print('|  hence CI orb are now numbered from ',OFFSET+min(collect_confs(CAS_conf)),' to ',OFFSET+max(collect_confs(CAS_conf)),' and written in ',OVERL_file_orb )
    dec_CI_conf=Offset_conf(CAS_conf,OFFSET)
    print('-------------------------------------------------')
    with open(VB_file , "r") as file:
        for line in file:
            MULT=1
            if routines.is_in(line, "MULTIPLICITY"):
               MULT= routines.Read_INT(line, "MULTIPLICITY")
               print('| Multiplicity is ',MULT)
    print('|  filling the ', OVERL_file_xmi,' xmi input file ' )
    with open(OVERL_file_xmi,'w') as file:
        line='made by CAS_LEW.py \n $ctrl   ; ============='+str(NVBCONF)+' + '+str(lenCI)+'=============.======+\n'
        line=line+'  vbftyp=det WFNTYP=struc iscf=1 itmax=0 int=libcint basis=6-31G ' 
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
            line= '  '+ dec_CI_conf[i]+'    ; CI '+str(i+1) + str( CAS_vect[i])+'\n'  
            #print(*['%4.0f' % int(val) for val in dec_CI_conf[i].split()],end=' ')
            #print('  ',dec_CI_conf[i],'       ; CI ',i+1 , CAS_vect[i])
            file.write(line)   
        line= ' $end  ; ============= \n '
        file.write(line)   
     ##   pos, line=routines.detect_keyword(VB_file, "$bfi", 0)
     ##   bfi_nom,bfi_noa,list_om,list_oa=routines.read_bfi(VB_file,pos)
     ##   line= ' $bfi  ; ============= \n '
     ##   line=line+' '+str(bfi_nom)+' '+str(bfi_noa)+ '\n   '+routines.makeSTR(list_om)+'\n   '+routines.makeSTR(list_oa)+'\n $end \n'
     ##   file.write(line)   
    OVERL_coeffs=[]
    OVERL_aos=[]
    for k in range(len(VB_orb_coeffs)):
   #     print('VBcoeffs', k, len(VB_orb_coeffs[k]),end=' ')
        OVERL_coeffs.append(VB_orb_coeffs[k])
        OVERL_aos.append(VB_orb_aos[k])
    for k in range(len(VB_orb_coeffs)):
   #     print('ICcoeffs', k, len(VB_orb_coeffs[k]),end=' ')
        OVERL_coeffs.append(CAS_orb_coeffs[k])
        OVERL_aos.append(CAS_orb_aos[k])
    routines.write_orb(OVERL_file_orb,OVERL_coeffs,OVERL_aos,0,len(OVERL_coeffs))
    print('  ',OVERL_file_orb,' written')   
    
    routines.write_DOLLARORB(OVERL_file_xmi,OVERL_aos,0,len(OVERL_aos))
    print('| $orb section, with ',len(OVERL_aos),' orbitals is to get in ORBB    ')
    print('| $end ')
    print('- ------------------------------------------------------------------')
    print('- submit the calculation:',OVERL_file_xmo, ' is required to continue')
    print('- ------------------------------------------------------------------')
    quit()  
if not os.path.exists(OVERL_file_xmo):
    print('- submit the calculation:',OVERL_file_xmo, ' is required to continue')
    quit()

OVERLP_size= routines.Read_INTEGER(OVERL_file_xmo, " Number of Structures:",12)  
if OVERLP_size != NVBCONF+lenCI           :
    print('Error in OVERLP_size',OVERLP_size,'instead of ',NVBCONF+lenCI           )
    quit()
else:
    print('|  OVERLP_size',OVERLP_size,' is OK with VB+CI ')
posit,line=routines.detect_keyword(OVERL_file_xmo,"******  OVERLAP OF VB STRUCTURES ",0)
posfin,line=routines.detect_keyword(OVERL_file_xmo,"******  HAMILTONIAN ",posit)
# read  The  Matrix at pos (3 lines between each blocks)
# 5 intergers first (5 cols) then each line as integer (line) then 5 reals 
#
n,m=check_size(OVERL_file_xmo,posit+1,posfin)
Stot=[]
print('| In ', OVERL_file_xmo, ' the overlap matrix is (',n,'x',m,')')
Stot=read5cols(OVERL_file_xmo,OVERLP_size,OVERLP_size,posit,posfin) # read a file with 3 blank lines, +1 to skip, then blocks of columns of length lines
#print_matrix('Stot',Stot)
#print('CAS_vect[1]',CAS_vect[1])
Smomo=np.zeros((lenCI,lenCI))
Svbvb=np.zeros((NVBCONF,NVBCONF))
part_SOM=np.zeros((NVBCONF,len(CAS_vect)))
SOMvb=np.zeros(NVBCONF)
print('| ' )
print('| ------------------------------------------------------------')
print('|  Get the overlaps between each of the CI\'s CSF of ',lenCI,'CAS CSFs with each of the ',NVBCONF,'VB conf')
for imo in range(lenCI):
    for jmo in range(lenCI):
        Smomo[imo][jmo]=Stot[imo+NVBCONF][jmo+NVBCONF] # Smomo is the part that concerns the overlap between each MO configurations of the CI.

for ivb in range(len(VB_conf)):
    for jvb in range(len(VB_vect)):
        Svbvb[ivb][jvb]=Stot[ivb][jvb]
    for jmo in range(len(CAS_vect)):
        part_SOM[ivb][jmo]=Stot[ivb][jmo+NVBCONF] # part_SOM is the part that concerns the overlap between each MO configurations of the CI and each VB configurations.

#print_lin_matrix('S_CI/CI',Smomo)
#print_lin_matrix('S_CI/VB^T',part_SOM.T)
#print('| The vector of the overlap between ',NVBCONF,' VB CSF\' (structure) ')
#print_lin_matrix('S_VB/VB',Svbvb)
print('| ' )
#print('| ------------------------------------------------------------')
#print('| The CI vector of this state ',CAS_vect,' has the following overlaps with each VB\'s structures')
##
# just check the norm of the CAS vector ------------------
norm_CAS_vect=np.dot(CAS_vect.T,np.dot(Smomo,CAS_vect))
if abs(norm_CAS_vect-1.0) > 0.01: 
    print('>>>>>>> .  Error in the norm of the CAS wf: ',f"{norm_CAS_vect:4.3f} <<<<<<<<<")
    quit()
print('| Norm of the CAS wf is OK: ',f"{norm_CAS_vect:4.3f}")
print('| solving the projection of the CAS wf on the VB wf')
# define & solv the HLP equations       ------------------
SOMvb=np.dot(part_SOM,CAS_vect)
sol=np.linalg.solve(Svbvb,SOMvb)
# align VB_vect on the CAS_vect         ------------------
if np.dot(VB_vect.T,SOMvb) < 0:
    print('| turn the vectors so they better aligns on the CAS:'," VB_vect=", VB_vect, 'CAS in the VB space=', SOMvb)
    VB_vect=-VB_vect
if np.dot(VB_vect.T,sol) < 0:
    sol=-sol
    Svbvbsol=np.dot(Svbvb,sol)  
    print('|  the solution is inverted to have >0  overlaps')
Svbvbsol=np.dot(Svbvb,sol)  
norm_sol=np.dot(sol.T,Svbvbsol)  
print('| projected wf_K obtained. Norme wf_K=',f"{norm_sol:4.3f}")
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
    for i in range(len(CAS_vect)):
        print('   CI',f"{i:2d}",end='  ')
    print()
    for i in range(len(sol)):
        print('|  ',f"{i+1:3d}",'  ',f"{VB_vect[i]:7.3f}",end=' ')
        for j in range(len(CAS_vect)):
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
