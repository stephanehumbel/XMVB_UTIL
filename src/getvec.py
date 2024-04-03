#!/usr/bin/env python3
import re
import sys
import os
import routines
import xmvb_orb as xmlib
import numpy as np  

class Orb:
    def __init__(self, zeta : int=0, coeffs : list[list[float]]=[], indices : list[list[int]]=[], numatoms : list(list([int]))=[]):
        self.zeta = zeta
        self.coeffs = coeffs
        self.indices = indices
        self.numatoms = numatoms

#Core function which does the job
def main(input_file):   
    return input_file


if __name__ == "__main__": # permet d'utiliser comme une librairie qu'on importe

   if len(sys.argv) == 2:
       zeta=1
       numatoms=[]
       coeffs=[]
       indices=[]
       input_file = sys.argv[1]
       input_file_name, input_file_ext = os.path.splitext(input_file)
       output_file_name=input_file_name+"_2"+input_file_ext
       main(input_file)
       input_file_name, input_file_ext = routines.file_extension(input_file)
#       print(input_file_name, input_file_ext)
       pos=routines.detect_keyword(input_file, "VEC", 0)
       coeffs,nvect = routines.read_vec(input_file,coeffs,pos+1)
#       print('ifin',coeffs)
       routines.write_vec(coeffs,1, nvect, "screen")
       routines.write_vec(coeffs,1, nvect, output_file_name)
       indices=routines.compte_AO(coeffs) 
       print('len(indices)',len(indices),indices)
       #copy to ccoeffs only the coeff to print aligned to indices
       ccoeffs=[]
       for i in range(len(indices)):
            print('-----------')
            print('bcl indices(',i,')',indices[i],len(indices[i]))    
            for j in range(len(indices[i])):
               #ccoeffs.append(coeffs[indices[i][i]])
               print(len(indices[i]),'indices(',i,',',j,' )=',indices[i][j],end=' ,')
           #   coeffs[i]=ccoeffs
           #   print('coeff  (',i,',',j,')',coeffs[i][j],end='')
               print(' coeffs(',i,',',j,')=',coeffs[i][1][indices[i][j]],end='')
           #    ccoeffs.append(coeffs[i][1][indices[i][j]])
            print('')
            print('-----------')
            #print('ccoeffs',ccoeffs)    
       routines.write_orb(output_file_name,coeffs,indices)
       sys.exit() 