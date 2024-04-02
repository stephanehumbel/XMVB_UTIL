import sys
import os
import copy 
import numpy as np 
import re       
                        

class Orb:
    def __init__(self, zeta : int=0, coeffs : list[list[float]]=[], indices : list[list[int]]=[], numatoms : list(list([int]))=[]):
        self.zeta = zeta
        self.coeffs = coeffs
        self.indices = indices
        self.numatoms = numatoms

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
    # Initialisation du compteur de AO non nuls et de la liste des AO
    nao = 0
    tab_ao = []
        
    # Parcours du vecteur MO
    for i in range(len(vect)):
        # Si le coefficient du AO est non nul
        if vect[i] != 0.:
            # Ajouter l'indice du AO à la liste des AO non nuls
            tab_ao.append(i)
            # Incrémenter le compteur de AO non nuls
            nao += 1
      
    # Retourner le nombre total de AO non nuls et la liste des AO non nuls
    return nao, tab_ao  
        

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
    with open(file_path, 'r') as file:
        values = []
        for line_num, line in enumerate(file, start=1):
            if line_num < start_line:
                continue  # skip the first lines

            if line.strip() == '$END':
                vectors.append((item, np.array(values)))  # add the last to vectors
                print('>last',item,line[0:2],len(values),values[0:12])
                break  # stop reading after $END

            # get from format (I2,I3,5F15.8)
            vector_number = int(line[0:2])
            print('values',vector_number,')',prev_vector_number,values)
            if vector_number == prev_vector_number:
               toread = (len(line) - 1 - 5) // 15 # skip 5 (+1) digits and get n (number of float nF15.8)
               for i in range(toread):            # so partially filled lines are read
                   start_index = 5 + i * 15       # skip 5 digits and the already read floats
                   end_index = start_index + 15   # field as nF15.8
                   values.append(float(line[start_index:end_index]))
               print('######',vector_number,'stored',values ))')
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
    return vectors



# PRINT_VEC
def write_vec(vectors, deb, fin, file_path):
    '''
    Description : write MO's as $VEC to a file.
    Args:
        parameter (vectors): the MO's vectors(i,array(NAO's))
        parameter (deb, fin): beginning and end of the printing
        parameter (file_path): path to the output file

    Returns:
        None
    '''
    if deb < 1 or fin > len(vectors) + 1:
        print('Error: Invalid limits in write_vec', deb, fin, len(vectors))
        sys.exit(1)

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


