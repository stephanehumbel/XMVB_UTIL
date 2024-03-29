import sys
import os
import copy 
import numpy as np 
import re       
                        

class Orb:
    def __init__(self, zeta : int=0, coeffs : list[list[str]]=[], indices : list[list[int]]=[], numatoms : list(list([int]))=[]):
        self.zeta = zeta
        self.coeffs = coeffs
        self.indices = indices
        self.numatoms = numatoms


def readorb(file_name):
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
                    coef.append(str(values[2*i]))  # add the last to vectors
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
        
def writorb(vectors,outfile):
    '''          
    Description : write MO's as .orb to a file.
    Args:        
        parameter (vectors): the MO's vectors(i,array(of coef for each AO's))
        parameter (outfile): path to the output file
    
    Returns:
        None

    Uses: compte_AO
    '''
    nao=0
    tab_ao=[]
    with open(outfile, 'w') as output_file:
    # first line : the nao_count 
         with vect in vectors:
              nao,tab_ao=compte_AO(vect)
              output_file.write(f"{nao:4d}")

