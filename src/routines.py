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
        for line_num, line in enumerate(file, 1):  
            if line_num >= start_line and keyword in line:
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
    tableau=[[val for val in pair[1]] for pair in vect]
    # Initialisation du compteur de AO non nuls et de la liste des AO
    tab_ao = []
    ##print('compte_AO',len(tableau),vect[3],len(tableau[0]))
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
    return tab_ao  
        

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
#    print('read_vec',file_path,start_line)
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
               for i in range(toread):            # so partially filled lines are read
                   start_index = 5 + i * 15       # skip 5 digits and the already read floats
                   end_index = start_index + 15   # field as nF15.8
                   values.append(float(line[start_index:end_index]))
#        print('read_vec',vectors)
    return vectors,vector_number
# 
def make_table(coeffs):
    tableau = []
    for i in range(len(coeffs)): 
        print("##----##",i)
        for value in coeffs[i][1]:
            if abs(value)>1e-6:
                print("{:5.3f}".format(value),i, end=' ')
            tableau.append(value)
    return tableau

def make_orb(coeffs, indices):
    ''' Returns a reduced vector of MOs with only the non-zero AOs'''
    coeffs_orb = coeffs
    coeffs_orb = []
    orb_coeffs = []
    print('make_orb',len(coeffs),len(indices))#,coeffs_orb[0][1])
    for i in range(len(coeffs)):
        #print()
        #print('coeffs(',i,',[1])',end='')
        for value in coeffs[i][1]:
            continue#print("{:5.3f}".format(value), end=' ')
        for j in range(len(coeffs[i][1])):
            if coeffs[i][1][j] != 0:
                orb_coeffs.append(coeffs[i][1][j])
            #    print('A',j,'  ',end='')#,orb_coeffs,end='')
            else:
                continue# print('Z',j,'|',end='')
        coeffs_orb.append(orb_coeffs)
    print('make_orb',len(coeffs_orb),len(coeffs_orb[0]),coeffs_orb) 
    return coeffs_orb



# PRINT_VEC
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


def write_orb(filename, coeffs, indices):
    #new_orb_values = new_orb_data.coeffs
    #new_indices = new_orb_data.indices
    #print (coeffs)
#    print('write_orb',filename,len(coeffs),len(indices))
    if filename == 'screen':
        for i in range(len(indices)):
            print(f"{len(indices[i]):4d}",end='')
        print()
        for i in range(len(indices)):
            print(f"# ORBITAL {i+1:4d}  NAO = {len(indices[i]):4d}")
            count = 0
            for j in range(len(indices[i])):
                print(f"{coeffs[i][1][indices[i][j]]:13.10f}{(indices[i][j])+1:4d}  ",end='')
                count += 1
                if (j+1) % 4 == 0 and j != len(indices[i])-1:
                    print()
            print()
    else:        
        with open(filename, 'w') as f: 
            for i in range(len(indices)):
                f.write(f"{len(indices[i]):4d}")
            f.write("\n")
            for i in range(len(indices)):
                f.write(f"# ORBITAL {i+1:4d}  NAO = {len(indices[i]):4d}\n")
                count = 0   
                for j in range(len(indices [i])):
#                    print(' coeffs(',i,',',j,')=',coeffs[i][1][indices[i][j]],end='')
                    f.write(f"{coeffs[i][1][indices[i][j]]:13.10f}{(indices[i][j])+1:4d}  ")
                    count += 1
                    if (j+1) % 4 == 0 and j != len(coeffs[i][1])-1:
                        f.write("\n")
                f.write("\n")

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
    return bfi_nom, bfi_noa, list_om, list_oa