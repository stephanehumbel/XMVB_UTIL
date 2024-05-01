#!/usr/bin/env python3
# write a code that diagonalise a matrix with an appropraite numpy routine
import os
import sys
import re
import routines
import numpy as np
def CS_weight(i,Ci,Sij):
    wi=Ci[i]*np.dot(Sij,Ci)[i]
#    print('CS_weight(',i,')=',wi)
    return wi
def Get_CIVECT(file, state ):
    search_state=" STATE #    "+str(state)
    pos_CASvect=routines.detect_keyword(CAS_file_name, "LAGRANGIAN CONVERGED", 0)
    pos_CASvect=routines.detect_keyword(CAS_file_name, search_state, pos_CASvect)
    pos_CASvectfin = routines.detect_blank(CAS_file_name,pos_CASvect + 3)
    CI_SIZE=pos_CASvectfin-1-(pos_CASvect+3)
    print('pos_CASvect',pos_CASvect, CI_SIZE)
    # on decale de 2 ligne
    CI_vect=Read_CIVECT(CAS_file_name, pos_CASvect+3,pos_CASvectfin-1)
    return CI_vect
    

def Read_CIVECT(file, pos,fin):
    with open(file, 'r') as f:
        lines = f.readlines()[pos:fin]
        n=fin-pos
        print('n',pos,n,fin)
        CI_vect = np.zeros(n)
        for i in range(n):
# Det,       CI_vect[i] = float(re.split('\|+| +|\n', lines[i])[7])
# GUGA
            CI_vect[i] = float(re.split(' +',lines[i])[2])
        print('CI_vect',CI_vect)
    

input_file = sys.argv[1]
state = sys.argv[2]
input_file_name, input_file_ext = os.path.splitext(input_file)
# Le calcul CAS est dans nom.out
CAS_file_name = input_file_name + ".log"
if not os.path.exists(CAS_file_name):
    log_files = [file for file in os.listdir() if file.endswith(".log")]
    print("Available .log files:")
    for file in log_files:
        print(file, end=' ')
    print()
    input_file = input("Enter the file name: ")
    input_file_name, input_file_ext = os.path.splitext(input_file)
    CAS_file_name = input_file_name + ".log"
# Define your matrix
CI_vect=Get_CIVECT(CAS_file_name, state)
print('CI_vect',CI_vect)
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

Svbvb = np.array([[1., .177, 0.07299], [0.177, 1., 0.0], [0.07299, 0.00, 1.]])
SOMvb = np.array([[.8768], [0.441327], [0.108482]])
normed = np.array([[.0000], [0.0000], [1.0]])
sol = np.linalg.solve(Svbvb, SOMvb)
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
