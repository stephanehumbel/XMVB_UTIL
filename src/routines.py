import sys
import os
import copy 
import numpy as np 
import re       
                        

class Orb:
    def __init__(self, zeta : int=0, coeffs : list[list[float]]=[], indices : list[list[int]]=[], numatoms : list[list[int]]=[]):
        self.zeta = zeta
        self.coeffs = coeffs
        self.indices = indices
        self.numatoms = numatoms

def file_extension(file_name):
    fname, extension = os.path.splitext(file_name)
    return fname, extension

def is_in(line, keyword):
    if keyword in line:
        return True
    return False

def Read_INT(line, keyword):
    words = re.split(' +|=+|\n',line)
    strin=' '.join(words) # twice to remove all the empty strings 
    words = re.split(' +|=+|\n',strin)
    index = words.index(keyword)
    integer = int(words[index + 1])
    return integer

def detect_keyword(file_path, keyword, start_line):
    '''
    Description : get the line number of the *last* occurence of keyword in file_path.
    Args:
        parameter (file_path): name of the file (possibly the path?).
        parameter (keyword): string to find
        parameter (start_line): the line the search starts from

    Returns:
        type: returns the line number (lin_num) of the *last* occurence of keyword or -1 if keyword is not found
    '''
#    print('detect_keyword',file_path,keyword,start_line)
    with open(file_path, 'r') as file:
        ret_num = -1                # not found
        ret_line=''                      # not found
        for line_num, line in enumerate(file, 1):  
            if line_num >= start_line and keyword in line:
                ret_line = line
                ret_num = line_num
    return ret_num, ret_line 

def Read_INTEGER(file_path, STRING, size):
    '''
    Description : get the integer value of the *last* occurence of keyword in file_path.
    Args:
        parameter (file_path): name of the file (possibly the path?).
        parameter (STRING): string to find. It must contain the = character, if so
        parameter (size): the size of the integer to typically 4 or 5 characters

    Returns:
        type: returns the integer value of the *last* occurence of keyword or -1 if keyword is not found
    '''
    numline, line_to_read = detect_keyword(file_path, STRING, 1)
    #print('Read_INTEGER',numline, line_to_read) 
    if numline > 0:
        col1=line_to_read.find(STRING)
        num= line_to_read[col1+len(STRING):col1+len(STRING)+size]
        #print(STRING,':',col1,':::',num)
        return int(num)



def detect_blank(file_path, start_line):
    '''
    Description: Get the position of the next blank line in the file.
    Args:
        parameter (file_path): Path to the file.
        parameter (start_line): The line number to start the search from.

    Returns:
        type: Returns the line number of the next blank line or -1 if not found.
    '''
    with open(file_path, 'r') as file:
        for line_num, line in enumerate(file, start=1):
            if line_num >= start_line and line.strip() == '':
                return line_num
    return -1  # not found

def read_orb(file_name):
    '''
    Description :  Reads coef and ao in  an .orb file 
    Args:
        (file_name):  the file .orb

    Returns:
        the tables all_coeffs and all_aos
    '''
    all_coeffs = []
    all_aos = []
    coef=[]
    ao=[]
    iom=0
    with open(file_name, 'r') as file:
      noa = []
      for  line in file:
          values = []
          if ("#" in line):
             if iom > 0:
                 all_coeffs.append(coef)
                 all_aos.append(ao)
                 #print(iom,'zz',len(all_coeffs),all_aos,all_coeffs)
                 iom+=1
                 coef=[]
                 ao=[]
             else:
                 iom+=1

          if not("#" in line):
             values=line.split()
             if iom==0:
                 noa.append(values)
             else:
                 #print("||",values[0], len(values))
                 toread=len(values)//2
                 for i in range(toread):
                 #   print(i,end='')
                    coef.append(float(values[2*i]))  # add the last to vectors
                    ao.append(int(values[2*i+1]))

      # last orb must be updated

    all_coeffs.append(coef)
    all_aos.append(ao)
    #print('rr',len(all_coeffs),all_aos,all_coeffs)
    return all_coeffs, all_aos

def make_dollarorb_file(ao_orb,fin,filename): # writes the $orb of vect
    if filename == 'screen':
        print( '$orb')
        for numorb in range(fin):
            print(len(ao_orb[numorb]),end=' ')
        for numorb in range(fin):
            print()
            for j in range(len(ao_orb[numorb])):
                print(' ',ao_orb[numorb][j]+1,end='')     
        print()
        print('$end')
    else:   
        with open(filename, 'w') as f:
            f.write('$orb')
            f.write("\n")
            for numorb in range(fin):
                f.write(f"{len(ao_orb[numorb]):4d}")
            for numorb in range(fin):
                f.write("\n")
                for j in range(len(ao_orb[numorb])):
                    f.write(f"{ao_orb[numorb][j]+1:4d}")
            f.write("\n")
            f.write("$end")
        f.close()   


def make_dollarorb(ao_orb,fin): # writes the $orb of vect
    print( '$orb')
    for numorb in range(fin):
        print(len(ao_orb[numorb]),end=' ')
    for numorb in range(fin):
        print()
        for j in range(len(ao_orb[numorb])):
            print(' ',ao_orb[numorb][j]+1,end='')     
    print()
    print('$end')

def make_dolorb(vect,fin): # writes the $orb of vect
    print( '$orb')
    for numorb in range(fin):
        make_oneorb(vect,numorb)    
        print()
    print( '$end')


def make_oneorb(vect,numorb): # writes the ao's of an orbital for the $orb
    for i in range(len(vect)):
        if vect[numorb][i] != 0:
            print(f"{i+1} ",end='') 

def make_bfi(vect): # make an all electron bfi
    print(' $bfi')
    print('  0    '  ,len(vect)   )
    print(' ')
    #imprime la liste des entiers entre 1 et len(ao_orb)
    tab_entier=range(1,len(vect)+1)
    for i in tab_entier:
        print (i,end=' ')
    print(' ')
    print(' $end')
#   
 
def compte_AO(vect):
    '''
    Description : compte et retourne la liste des AO non nuls dans un vecteur MO
    Args:
        parameter (vect): le vecteur MO (array de NAO)
        parameter (nao): nombre de AO non nuls dans le vecteur
        parameter (tab_ao): liste des AO NAO
    
    Returns:
        nao, tab_ao
    ''' 
    #tableau est le tableau des OM, sur les OA, mais tout demarre à tableau[0][0]
    tableau=vect
#    tableau=[[val for val in pair[1]] for pair in vect]
    # Initialisation du compteur de AO non nuls et de la liste des AO
    tab_ao = []
#    print('compte_AO',len(tableau),vect[3],len(tableau[0]))
    # Parcours du vecteur MO  tous les vecteurs du tableau on la meme taille?
    for i in range(len(tableau)):
        tab_ao.append([])
        for j in range(len(tableau[i])):
            nao = 0
            #print('bcl',i,len(tableau[i]),end='')     
            if tableau[i][j] != 0.:
                # Ajouter l'indice du AO a la liste des AO non nuls
                tab_ao[i].append(j)
#                tab_ao[i][nao]=j
                # Incrémenter le compteur de AO non nuls
                nao += 1
    #print('l',nao,tab_ao,len(tab_ao[0]),len(tab_ao[1]))
      
    # Retourner le nombre total de AO non nuls et la liste des AO non nuls
    print('fin compte_AO:',end='' )
    for i in range(len(tab_ao)):
        print(len(tab_ao[i]),end=' ')   
    print()
    return tab_ao  
        
# READ_geom
def read_geom(file_name):
    '''
    read geometry from a gamess .log file
    '''
    a0=0.529177 # Bohr radius
    num,line=detect_keyword(file_name, "TOTAL NUMBER OF ATOMS", 0) 
    natoms=Read_INT(line,"ATOMS")
    num,line=detect_keyword(file_name, "ATOM      ATOMIC                      COORDINATES", 0) 
    with open(file_name, 'r') as file:
        symbol=[]
        zat=[]
        x=[]
        y=[]
        z=[]
        line_num=0
        for line in file:
#            print(line[0:2])
            values=[]
            line_num+=1
            if line_num >= num+2:
                if line_num >= num+natoms+2:
                    break   
#                print(line_num, line,num)
                values = re.split(' +|\n',line)
#                print('#--',values, len(values),values[2])
                symbol.append(values[1])
                zat.append(values[2])
                x.append(float(values[3])*a0)   
                y.append(float(values[4])*a0)   
                z.append(float(values[5])*a0)   
#    print('read_geom',x,y,z,natoms)
    return symbol,zat,x,y,z,natoms
# READ_VEC
def read_vec(file_path,vectors,start_line):
    '''
    Description : reads MO's from either a .inp or a .dat file).
    Args:
        parameter (file_path): name of the file (possibly the path?).
        parameter (vectors): the MO's as two tuples:vectors[i][0]=N, vectors[i][1]=Coef np.array
        parameter (start_line): the line the reading starts from: must be the line after $VEC

    Returns:
        type: returns the  MO's as vectors, (until $END is reached)
    '''
#    def read15(line,values):
#        values=[]
#        toread = (len(line) - 1 - 5) // 15 # skip 5 (+1) digits and get n (number of float nF15.8)
#        for i in range(toread):            # so partially filled lines are read
#            start_index = 5 + i * 15       # skip 5 digits and the already read floats
#            end_index = start_index + 15   # field as nF15.8
#            values.append(float(line[start_index:end_index]))
    
    prev_vector_number=1
    item=1
   # print('|  read_vec from:  ',file_path,' OM ', item, end=' ')
    with open(file_path, 'r') as file:
        values = []
        for line_num, line in enumerate(file, start=1):
            if line_num < start_line:
                continue  # skip the first lines

            if line.strip() == '$END':
                vectors.append((item, np.array(values)))  # add the last to vectors
#                print('>last',item,line[0:2],len(values),values[0:12])
                break  # stop reading after $END

            # get from format (I2,I3,5F15.8)
            vector_number = int(line[0:2])
#            print('values',vector_number,')',prev_vector_number,values)
            if vector_number == prev_vector_number:
               toread = (len(line) - 1 - 5) // 15 # skip 5 (+1) digits and get n (number of float nF15.8)
               for i in range(toread):            # so partially filled lines are read
                   start_index = 5 + i * 15       # skip 5 digits and the already read floats
                   end_index = start_index + 15   # field as nF15.8
                   values.append(float(line[start_index:end_index]))
#               print('######', vector_number, 'stored', values)
            else:
               vectors.append((item, np.array(values)))  # add to vectors
               item=item+1
               prev_vector_number=vector_number
               values=[]
               toread = (len(line) - 1 - 5) // 15 # skip 5 (+1) digits and get n (number of float nF15.8)
               #print(item,end=' ')   
               for i in range(toread):            # so partially filled lines are read
                   start_index = 5 + i * 15       # skip 5 digits and the already read floats
                   end_index = start_index + 15   # field as nF15.8
                   values.append(float(line[start_index:end_index]))
#        print('read_vec',vectors)
#        print('-',item,end=' ')   
#        print()
    return vectors,vector_number
# 
def to_array(liste):
    ''' convert a list of list to a numpy array'''
    vectors=[]  
    item=0
    print('to_array',end=': ')
    for values in liste:
        print(item,values)
        v=[]
        for i in range(len(values)):
          v.append(float(values[i]))    
        vectors.append(np.array(v))
        item+=1
    return vectors

def make_table(coeffs):
    """ returns a table of MO's from  this bizarre tuple thing tableau[i]=coeffs[i][1] from read_vec:
    that I should re write
    """
#    print(" routines.make_table for data conversion (should be rewrite someday)")
#    print('make_table', end=': ')
    tableau = []
    for i in range(len(coeffs)): 
        orbital = []
        value=[]
#        print(coeffs[i][1])
        #print(i,end=' ')
        for value in coeffs[i][1]:
            #print("{:5.3e}".format(value), end=' ')
            orbital.append(value)
            if abs(value)>1e-6:
                continue #print("{:5.3f}".format(value),i, end=' ')
        tableau.append(orbital)
#        print("tableau d orbital",tableau[i])
    return tableau

def make_orb(vect, indices):
    ''' Returns a reduced vector of MOs with only the non-zero AOs and the list of AO'''
    ao_orb = []
    coeffs_orb = []
    orb_coeffs = []
    orb_ao = []
#    print('make_orb',len(vect),len(indices))
    for i in range(len(vect)):
        #print()
        #print('coeffs(',i,',[1])',end='')
        #value=[]
        noa = 0  
        ioa = 0  
        orb_coeffs = []
        orb_ao = []
#        print()
#        print('orbitale n°',i,len(vect[i]),": ",end='')
        for val in vect[i]:
            #value=float(val)
            #print("{:5.3f}".format(value), end=' ')
            #continue
            if val != 0:
#                print(ioa,'',end=' ')
                orb_coeffs.append(val)
                orb_ao.append(ioa)
                noa += 1
            ioa += 1
            #for j in range(len(vect[i])):
            #    if vect[i][j] != 0:
            #        print('A',j,'  ',vect[i][j],end='')
            #        orb_coeffs.append(vect[i][j])
            #    else:
            #        continue# print('Z',j,'|',end='')
        coeffs_orb.append(orb_coeffs)
        ao_orb.append(orb_ao)
#        print(noa,end='Z')

    #print('make_orb',len(coeffs_orb),len(coeffs_orb[0]))#,coeffs_orb) 
    return ao_orb,coeffs_orb



# PRINT_VEC
def wwrite_vec(file_path, vectors, deb, fin):
    '''
    Description : write MO's as $VEC to a file.
    Args:
        parameter (vectors): the MO's vectors(array(NAO's))
        parameter (deb, fin): beginning and end of the printing
        parameter (file_path): path to the output file ; "screen" for screen output

    Returns:
        None
    '''
    if deb < 0 or fin > len(vectors) + 1:
        print('Error: Invalid limits in write_vec', deb, fin, len(vectors))
        sys.exit(1)
    numvec=0
    if file_path == 'screen':
        print('$VEC ',end='')
        for values in vectors:
            line_number = 1
            if deb <= numvec <= fin:
                for i in range(0, len(values)):
                    if i % 5 == 0:
                        print(f"\n{numvec % 100:2d}{line_number:3d}",end='')
                        line_number += 1
                    print(f"{values[i]:15.8E}",end='')
            numvec+=1
        print('\n$END')
    else:
        with open(file_path, 'w') as output_file:
             output_file.write(' $VEC')
             for values in vectors:
                line_number = 1
                if deb <= numvec <= fin:
                     for i in range(0, len(values)):
                         if i % 5 == 0:
                             output_file.write('\n')
                             output_file.write(f"{numvec+1 % 100:2d}{line_number:3d}")
                             line_number += 1
                         output_file.write(f"{values[i]:15.8E}")

                numvec+=1
             output_file.write('\n $END')
             output_file.write('\n')



def write_vec(file_path, vectors, deb, fin):
    '''
    Description : write MO's as $VEC to a file.
    Args:
        parameter (vectors): the MO's vectors(i,array(NAO's))
        parameter (deb, fin): beginning and end of the printing
        parameter (file_path): path to the output file ; "screen" for screen output

    Returns:
        None
    '''
    if deb < 1 or fin > len(vectors) + 1:
        print('Error: Invalid limits in write_vec', deb, fin, len(vectors))
        sys.exit(1)
    if file_path == 'screen':
        print('$VEC ',end='')
        for vector_number, values in enumerate(vectors):
            line_number = 1
            indice = vectors[vector_number][0]
            if deb <= indice <= fin:
                for i in range(0, len(values[1])):
                    if i % 5 == 0:
                        print(f"\n{indice % 100:2d}{line_number:3d}",end='')
                        line_number += 1
                    print(f"{values[1][i]:15.8E}",end='')
        print('\n$END')
    else:
        with open(file_path, 'w') as output_file:
             output_file.write(' $VEC')
             for vector_number, values in enumerate(vectors):
                 line_number = 1
                 indice = vectors[vector_number][0]
                 if deb <= indice <= fin:
                     for i in range(0, len(values[1])):
                         if i % 5 == 0:
                             output_file.write('\n')
                             output_file.write(f"{indice % 100:2d}{line_number:3d}")
                             line_number += 1
                         output_file.write(f"{values[1][i]:15.8E}")

             output_file.write('\n $END')
             output_file.write('\n')


def write_orbs(filename, phis, deb, fin):
    """ write the table phis table of MO's to a file or screen from deb to fin"""
    temponao=[]
#    print('|   write_orbs',end='._._.')
#    print('write_orbs',filename,len(phis),deb,phis[deb])
    if filename == 'screen':
        for i in range(deb , fin):
            compte=0
            for j in range(0  , len(phis[i])):
                if phis[i][j] != 0:
                    compte+=1
            print(f"{compte:4d}",end='')
            temponao.append(compte)
        for i in range(deb , fin):
            print() 
            compte=0
            print("# _._  Orbital  :   ",i+1," write_orbs---",temponao[i], "//" ,end='')
            for j in range(0  , len(phis[i])):
                if phis[i][j] != 0:
                    if (compte) % 4 == 0:
                            print()
                    print(f"{float(phis[i][j]):13.10f}{j+1:4d}  ",end='')
                    compte+=1
                #print(f"{float(phis[i][j]):13.10f}{j+1:4d}  ",end='')
            #print()
        print() 
    else:        
        with open(filename, 'w') as f: 
            for i in range(deb , fin):
                compte=0
                for j in range(0  , len(phis[i])):
                    if phis[i][j] != 0:
                        compte+=1
                f.write(f"{compte:4d}")
                temponao.append(compte)
                #print(f"{compte:4d}",end='')
            for i in range(deb , fin):
                f.write("\n")
                compte=0
                f.write(f"# __  Orbital  :   {i+1:4d} write_orbs----NAO={temponao[i]:4d}")
                for j in range(0  , len(phis[i])):
                    if phis[i][j] != 0:
                        if (compte) % 4 == 0:
                                f.write("\n")
                        #print(f"{float(phis[i][j]):13.10f}{j+1:4d}  ",end='')
                        f.write(f"{float(phis[i][j]):13.10f}{j+1:4d} ")
                        compte+=1
                    #print(f"{float(phis[i][j]):13.10f}{j+1:4d}  ",end='')
            #for i in range(len(indices)):
#    print("-end write_orbs-----") 
def write_orb(filename, coeffs, indices, deb, fin):
#    print('write_orb',filename,len(coeffs),len(indices))
    if filename == 'screen':
        print('write_orb:',end=''  )
        for i in range(deb , fin):
            print(f"{len(indices[i]):4d}",end='')
        print()
        for i in range(deb , fin):
            print(f"# ORBITAL {i+1:4d}  NAO = {len(indices[i]):4d}")
            count = 0
            for j in range(len(indices[i])):
                print(f"{float(coeffs[i][j]):13.10f}{(indices[i][j]):4d}  ",end='')
                #print(f"{float(coeffs[i][indices[i][j]]):13.10f}{(indices[i][j])+1:4d}  ",end='')
                count += 1
                if (j+1) % 4 == 0 and j != len(indices[i])-1:
                    print()
            print()
    else:        
        with open(filename, 'w') as f: 
            for i in range ( deb , fin):
                f.write(f"{len(indices[i]):4d}")
            f.write("\n")
            for i in range (deb , fin):
                f.write(f"# ORBITAL {i+1:4d}  NAO = {len(indices[i]):4d}      routines.write_orb({filename})\n")
                count = 0   
                for j in range(len(indices [i])):
#                    print(' coeffs(',i,',',j,')=',coeffs[i][1][indices[i][j]],end='')
                    #f.write(f"{coeffs[i][1][indices[i][j]]:13.10f}{(indices[i][j])+1:4d}  ")
                    f.write(f"{coeffs[i][j]:13.10f}{(indices[i][j]):6d}  ")
                    count += 1
                    if (j+1) % 4 == 0 and j != len(coeffs[i])-1:
                        f.write("\n")
                f.write("\n")
        f.close()

def write_conf(filename, CONF,COEF):
    ssize=len(CONF)
#    print('write_conf',filename,len(CONF),ssize)
    if filename == 'screen':
       for ii in range(ssize):
            print(' ', CONF[ii],';',ii+1,'... ',COEF[ii])    
#       print(" \n")
    else:
       print(ssize, "confs",end='')
       with open(filename, 'w') as f: 
            for ii in range(ssize):  
#                print(' ', CONF[ii],';',ii+1,'... ',COEF[ii])    
                f.write(f"{CONF[ii]}  ; {ii+1} ... {COEF[ii]:13.10f}")
                f.write("\n")
       print(" written on ",filename)    


def write_DOLLARORB(filename, AOS, deb, fin):
    if filename == 'screen':
        print('$orb')
        for i in range(len(AOS)):
            #print(*['%4.0f' % int(val) for val in VB_conf[i].split()],end=' ')
            print('',len(AOS[i]),end='')
        for i in range(len(AOS)):
            print()
            for j in range(len(AOS[i])):
                print(' ',AOS[i][j],end='')
        print()
        print('$end')
    else:        
        with open(filename, 'a') as f: 
            f.write("\n")
            f.write('$orb')
            f.write("\n")
            for i in range ( deb , fin):
                f.write(f"{len(AOS[i]):4d}")
            for i in range(len(AOS)):
                f.write("\n")
                for j in range(len(AOS[i])):
                    f.write(f"{AOS[i][j]:4d}")
            f.write("\n")
            f.write("$end")
        f.close()



def makeSTR(list_of_int,offset):
        ''' convert a list of integers to a string, replacing consecutive integers with a range '''
        ''' offset is the value to add to the integers in the list (corrextion de zero)'''
        list_of_int.sort()
        str_list = []
        i = 0
        while i < len(list_of_int):
            start = list_of_int[i]+offset
            end = start
            while i + 1 < len(list_of_int) and list_of_int[i + 1] == end + 1:
                end = list_of_int[i + 1]+offset
                i += 1
            if start == end:
                str_list.append(str(start))
            else:
              #  print('start',start,end)
                str_list.append(f"{start}-{end}")
            i += 1
        return ' '.join(str_list)
    

def read_bfi(file_path, start_line):
    bfi_nom = 0
    bfi_noa = 0
    list_om = []
    list_oa = []
    with open(file_path, 'r') as file:
        lines = file.readlines()[start_line:]
        bfi_nom, bfi_noa = map(int, lines[0].split())
        counter = 0
        index = 1
        while  'end' not in lines[index]:
            line = lines[index]
            ranges = re.split(r'\s+', line.strip())
            for r in ranges:
                if '-' in r:
                    start, end = map(int, r.split('-'))
                    for num in range(start, end+1):
                        if counter < bfi_nom:
                            list_om.append(num)
                        elif counter < bfi_nom + bfi_noa:
                            list_oa.append(num)
                        else:
                            break
                        counter += 1
                else: 
                    if r=='':   # empty line
                       break
                    if counter < bfi_nom:
                        list_om.append(int(r))
                    elif counter < bfi_nom + bfi_noa and index < len(lines):
                        list_oa.append(int(r))
                    else:
                        break
                    counter += 1
            index += 1
    file.close()
    return bfi_nom, bfi_noa, list_om, list_oa




def make_conf_from_gamess(CI_conf):
    ''' 
    takes a table of strings like 22201210000 and convert to a table of string for xmvb conf
    like 1:3 6 6 5 7
    '''
    tab=[]
#    print('make_conf_from_gamess')
    for i in range(len(CI_conf)):
#        print(CI_conf[i],end='')
        k=0
        TheXMconf = ''
        Theconf = CI_conf[i]
#        print(Theconf,' . ')
        while True: 
            if k >= len(Theconf):
                break
            else:
                if Theconf[k] == '0':
                    TheXMconf=TheXMconf
                if Theconf[k] == '1':
                    TheXMconf=TheXMconf+' '+str(k+1)    
                if Theconf[k] == '2':
                    TheXMconf=TheXMconf+' '+str(k+1)+' '+str(k+1)    
            k+=1
#            print()
#        print('TheXMconf',TheXMconf,'CI_conf',CI_conf[i])   
        tab.append(TheXMconf)
 #   print('in the routines',tab)
    return tab


def Make_conf_from_gamess(CI_conf):
    ''' 
    takes a table of strings like 22201210000 and convert to a table of string for xmvb conf
    like 1:3 6 6 5 7
    '''
    tab=[]
    print('make_conf_from_gamess')
    for i in range(len(CI_conf)):
        print(CI_conf[i],end='')
        k=0
        TheXMconf = ''
        alasuite=False
        Theconf = CI_conf[i]
        print(Theconf,' . ')
        while True: 
            if k >= len(Theconf)-1:
                break
            else:
                print('.k=',k,Theconf[k],'|',TheXMconf,'',end='')
                if Theconf[k] == '0':
                    print('vide :',str(k+1),'X',end='')
                    print()
                    alasuite=False
                    continue
                if Theconf[k] == '1':
                    TheXMconf=TheXMconf+' '+str(k+1)
                    alasuite=False
                    print ()
                    continue
                if Theconf[k] == '2':
                    if not alasuite:
                        TheXMconf=TheXMconf+' '+str(k+1)
                        try:
                            if Theconf[k+1] == '2':
                                alasuite=True
                            else:
                                print('fin de serie de 2',str(k+1),')))',end='')
                                print ()
                                alasuite=False
                                TheXMconf=TheXMconf+':'+str(k+1)
                        except:
                            TheXMconf+=' '+str(k+1)
                            break
                    else:
                        try:
                            if Theconf[k+1] == '2':
                                break
                            else:
                                TheXMconf=TheXMconf+':'+str(k+2)
                                alasuite=False
                        except:
                                TheXMconf=TheXMconf+':'+str(k+1)
                                alasuite=False
                                break
            k+=1
        print()
        print('TheXMconf',TheXMconf,'CI_conf',CI_conf[i])   
        tab.append(TheXMconf)
    return tab

