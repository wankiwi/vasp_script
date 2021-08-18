#!/usr/bin/env python
# coding: utf-8
__author__ = 'wankw (wankaiweii@gmail.com)' 

from ase.neb import NEB
from ase.io import read
from ase.io import write
import numpy as np
np.set_printoptions(suppress=True)
import os

def check_move_far(initial,final):
    is_pos = initial.get_positions()
    fs_pos = final.get_positions()
    cell_vec = initial.get_cell()
    is_scaled_pos = is_pos @ np.linalg.inv(cell_vec)#get_scaled_positions输出没有负号
    fs_scaled_pos = fs_pos @ np.linalg.inv(cell_vec)
    dist = np.linalg.norm(list(map(np.linalg.norm, (is_pos - fs_pos))))
    print (f'Distance between is&fs is {dist}!')

    n = 0
    far_list = []
    for i in abs(is_scaled_pos - fs_scaled_pos) >  0.5:
        if i.any():
            print (f'Atom no.{n+1} moved more than half of cell length! vector_diff: {abs(is_scaled_pos-fs_scaled_pos)[n]}')
            far_list.append(n)
        n += 1
        
    if far_list != []:
        temp = 'atoms' if len(far_list)>1 else 'atom'
        print (f'{len(far_list)} {temp} moved more than half of cell length, fix it automatically?(y/n)')
        answer = input('')
    #is_scaled_pos_new = is_scaled_pos.copy()
    if answer == 'y':
        for i in far_list:
            a = abs(is_scaled_pos[i] - fs_scaled_pos[i]) >  0.5
            b = is_scaled_pos[i] > [0.5,0.5,0.5]
            c = np.array([-1,-1,-1]) ** b * a
            is_scaled_pos[i] = is_scaled_pos[i] + c
            
    return is_scaled_pos

def linear_interpolation(ini,fin):
    images = [ini.copy() for i in range(n_image+1)] + [fin]
    #linearly_interpolates
    neb = NEB(images)
    neb.interpolate()
    #generate POSCAR and movie.xyz
    movie = []
    for i in range(len(images)):
        if not os.path.isdir(f'{i:02}'):
            os.mkdir(f'{i:02}')
        images[i].write(f'{i:02}/POSCAR',vasp5=True,direct=True)
        images[i].write(f'{i:02}/POSCAR.xyz')
        with open(f'{i:02}/POSCAR.xyz') as xyz_file:
            line = xyz_file.readline()
            while line:
                movie.append(line)
                line = xyz_file.readline()
                
    with open('movie_line.xyz','w') as movie_xyz:
        movie_xyz.writelines(movie)

#Reference: S. Smidstrup, A. Pedersen, K. Stokbro and H. Jonsson, Improved initial guess for minimum energy path calculations, J. Chem. Phys. 140, 214106 (2014).
def idpp_interpolation(ini,fin):
    images = [ini.copy() for i in range(n_image+1)] + [fin]
    #idpp_interpolates
    neb = NEB(images)
    neb.interpolate('idpp')
    #generate POSCAR and movie.xyz
    movie = []
    for i in range(len(images)):
        if not os.path.isdir(f'{i:02}'):
            os.mkdir(f'{i:02}')
        images[i].write(f'{i:02}/POSCAR',vasp5=True,direct=True)
        images[i].write(f'{i:02}/POSCAR.xyz')
        with open(f'{i:02}/POSCAR.xyz') as xyz_file:
            line = xyz_file.readline()
            while line:
                movie.append(line)
                line = xyz_file.readline()
                
    with open('movie_idpp.xyz','w') as movie_xyz:
        movie_xyz.writelines(movie)
        
def get_version():
    return '1.1 (2021.8.18, wankaiweii@gmail.com)'

if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description='Takes initial and final CARfiles, generate the initial guess '
                                                 'images between them by linear interpolation or image dependent '
                                                 'pair potential (idpp) interpolation. The initial guess files are '
                                                 'written to the directories 00 to NI+1, where NI is the number of specified images.')
    
    parser.add_argument('-v', '--version', action='version', version=get_version(),help='Display version')
    parser.add_argument('-i','--initial_state_carfile', type=str, action='store', default=r'is/CONTCAR',
                        help='The storage location of initial state CARfile. [Optional] [default="is/CONTCAR"]')   
    parser.add_argument('-f','--final_state_carfile', type=str, action='store', default=r'fs/CONTCAR',
                        help='The storage location of final state CARfile. [Optional] [default="fs/CONTCAR"]')
    parser.add_argument('-m','--interpolation_method', type=str, action='store',choices={'idpp','line'}, default='idpp',
                        help='The method of interpolation. [Optional] [default="idpp"]')
    parser.add_argument('-n','--number_of_images', type=int, action='store', default=5,
                        help='The number of interpolation. [Optional] [default=5]')
    
    args = parser.parse_args()
    
    initial = read(rf'{args.initial_state_carfile}')
    final = read(rf'{args.final_state_carfile}')
    n_image = args.number_of_images
    interpolation_method = args.interpolation_method
    
    initial.set_scaled_positions(check_move_far(initial,final))
    
    if interpolation_method == 'line':
        linear_interpolation(initial,final)
    elif interpolation_method == 'idpp':
        idpp_interpolation(initial,final)
        
    print (f'Generate {n_image} images between {args.initial_state_carfile} & {args.final_state_carfile} by {interpolation_method} interpolation!')

