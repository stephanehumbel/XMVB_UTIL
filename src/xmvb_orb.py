#!/usr/bin/python3
import re
import sys
import os

class Orb:
    def __init__(self, zeta : int=0, coeffs : list[list[float]]=[], indices : list[list[int]]=[], numatoms : list(list([int]))=[]):
        self.zeta = zeta
        self.coeffs = coeffs
        self.indices = indices
        self.numatoms = numatoms

def readorb(file_name):
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

#def read_inp(file_name):


def make_numatoms(aos, zeta):
    numatoms = []
    zeta = int(zeta)
    for i in range(0, len(aos), 1):
        nums = []
        if zeta == 1:
            for j in range (0, len(aos[i]), 1):
                nums.append(aos[i][j])
            numatoms.append(nums)
        else:
            for j in range (0, len(aos[i]), zeta):
                num = aos[i][j+zeta-1]//zeta
                nums.append(num)
            numatoms.append(nums)
    return numatoms


def change_zeta_coeffs(orb_data, new_zeta):
    new_zeta = int (new_zeta)
    old_zeta = int(orb_data.zeta)
    old_coeffs = orb_data.coeffs
    #print (old_coeffs)
    new_coeffs = []
    print(f"Old Zeta : {old_zeta}",f"New Zeta : {new_zeta}")
    if new_zeta == old_zeta:
        print("The orbital data is already in the desired format.")
        return orb_data
    if new_zeta > old_zeta:
        for i in range (0, len(old_coeffs), 1):
            compteur = 0
            stock = []
            for j in range (0, len (old_coeffs[i]), 1):
                if compteur == old_zeta:
                    while compteur < new_zeta:
                        stock.append(str(0.0000000000))
                        compteur += 1
                    stock.append(old_coeffs[i][j])
                    compteur = 1
                else:
                    stock.append(old_coeffs[i][j])
                    compteur += 1
            stock.append(str(0.0000000000))
            new_coeffs.append(stock)
    if new_zeta < old_zeta:
        for i in range (0, len(old_coeffs), 1):
            boucle = 0
            compteur = 0
            while boucle < len(old_coeffs[i]):
                while compteur != new_zeta:
                    stock.append(old_coeffs[i][boucle])
                    compteur += 1
                    boucle += 1
                compteur = 0
                boucle += old_zeta - new_zeta
            new_coeffs.append(stock)
    return new_coeffs

def change_zeta_aos(orb_data, new_zeta):
    new_zeta = int (new_zeta)
    numatoms = orb_data.numatoms
    new_aos = []
    for i in range (0, len(numatoms), 1):
        stock = []
        for j in range (0, len(numatoms[i]), 1):
            for k in range (1, new_zeta+1, 1):
                stock.append(new_zeta*(numatoms[i][j]-1)+k)
            #print (stock)
        new_aos.append(stock)
    return new_aos

    if new_zeta > old_zeta:
        for j in range(len(atomsnum)):
            for k in range(1, new_zeta+1, 1):
                new_indices.append(new_zeta*(atomsnum[j]-1)+k)
        #print (new_indices)
        for l, item in enumerate(old_coeffs):
            new_coeffs.append(item)
            if (l+1) % old_zeta == 0:
                for k in range(1, new_zeta-old_zeta+1, 1):
                    new_coeffs.append(str(0.0000000000))
        new_orb_data = Orb(new_zeta, numorb, new_nao_count, new_coeffs, new_indices, atomsnum)
        return new_orb_data
    if new_zeta < old_zeta:
        duke = 1
        for j in range(len(atomsnum)):
            for k in range(1, new_zeta+1, 1):
                new_indices.append(new_zeta*(atomsnum[j]-1)+k)
        while duke < len(old_coeffs):
            if duke % (old_zeta) == 0:
                duke += (old_zeta-new_zeta)
                #print ('bloup')
                continue
            else:
                new_coeffs.append(old_coeffs[duke-1])
                #print ('blip')
                duke += 1
        #print (new_coeffs)
        new_orb_data = Orb(new_zeta, numorb, new_nao_count, new_coeffs, new_indices, atomsnum)
    #else:
        #print("The orbital data is already in TZ format.")
        return new_orb_data

def symm_numatoms(file_path,orb_data):
    numatoms = orb_data.numatoms
    new_numatoms = []
    with open(file_path) as f:
        for line in f:
            numbers = line.split()
            for i in range(0, len(numbers), 2):
                num1 = numbers[i]
                num2 = numbers[i+1]
                for j in range(0, len(numatoms), 1):
                    new_nums = []
                    for k in range(0, len(numatoms[j]), 1):
                        if numatoms[j][k] == int(num1):
                            new_nums.append(int(num2))
                        elif numatoms[j][k] == int(num2):
                            new_nums.append(int(num1))
                    new_numatoms.append(new_numatoms)
    return new_numatoms

def write_new_indices(numatoms,zeta):
    new_indices = []
    for i in range (0, len(numatoms), 1):
        for j in range (0, len(numatoms[i]), 1):
            indi = []
            for k in range(1, zeta+1, 1):
                indi.append(zeta*(numatoms[i][j]-1)+k)  
            new_indices.append(indi)
    return new_indices  

def write_orb_file(filename, new_orb_data):
    new_zeta = int(new_orb_data.zeta)
    new_orb_values = new_orb_data.coeffs
    new_indices = new_orb_data.indices
    with open(filename, 'w') as f: 
        for i in range(len(new_indices)):
            f.write(f"   {len(new_orb_values[i])}   ")
        f.write("\n")
        for i in range(len(new_indices)):
            f.write(f"# ORBITAL          {i}  NAO =      {len(new_orb_values[i])}\n")
            for j in range(len(new_orb_values[i])):
                f.write(f"   {new_orb_values[i][j]}   {str(new_indices[i][j])}   ")
                if (j + 1) % new_zeta == 0:
                    f.write("\n")


#function which allows the user to input the name of the file to be read
def get_file_extension(file_name):
    fname, extension = os.path.splitext(file_name)
    return fname, extension

#Core function which does the job
def main(input_file, zeta, new_zeta, sym):   
    input_file_name, input_file_ext = get_file_extension(input_file)
    #print(input_file_name, input_file_ext)
    if input_file_ext == ".orb":
        #Create a new blank output file with the same name as the input file added with a suffix containing new_zeta
        output_file = input_file_name + f"_zeta{new_zeta}{sym}.orb"
        coeffs, indices = readorb(input_file)
        #print (coeffs)
        numatoms = make_numatoms(indices, zeta)
        #print (numatoms)
        orb_data = Orb(zeta, coeffs, indices, numatoms)
        if sym is not None:
            new_numatoms = symm_numatoms(sym, orb_data)
            new_indices = write_new_indices(new_numatoms, zeta)
            new_orb_data = Orb(zeta, coeffs, new_indices, new_numatoms)
        if new_zeta is not None or new_zeta != zeta:
            new_coeffs = change_zeta_coeffs(orb_data, new_zeta)
            #print (new_coeffs)
            new_aos = change_zeta_aos(orb_data, new_zeta)
            #print (new_aos)
            new_orb_data = Orb(new_zeta, new_coeffs, new_aos, numatoms)

        write_orb_file(output_file, new_orb_data)
    else:
        print("The input file must be a .orb file.")
    return output_file


if len(sys.argv) == 4:
    file = sys.argv[1]
    if sys.argv[2] .isdigit():
        zeta = sys.argv[2]
        new_zeta = sys.argv[3]
        symmetry = None
    else:
        zeta = sys.argv[2]
        symmetry = sys.argv[2]
        new_zeta = None
    output_file = main(file, zeta, new_zeta, symmetry)
    print (f"New file {output_file} created.")
elif len(sys.argv) == 5:
    file = sys.argv[1]
    zeta = sys.argv[2]
    new_zeta = sys.argv[3]
    symmetry = sys.argv[4]
    output_file = main(file, zeta, new_zeta, symmetry)
    print (f"New file {output_file} created.")
else:
    print("No argument provided. Please provide an argument.")
    print("Usage: python3 xmvbtz.py <input_file> <zeta> <new_zeta> <symmetry>")
    print("Example: python3 xmvbtz.py Tetra_1_1.orb 3 4 sym.txt")
    print("Example: python3 xmvbtz.py Tetra_1_1.orb 3 4")
    print("sym.txt is a file containing the symmetry of the molecule. It is optional.")
    print("If sym.txt is not provided, the program will only change the zeta value.")
    print("Example of sym.txt for permutation: 1  2  3  4 ")
    print("The first number is the old atom number. The second number is the new atom number.")