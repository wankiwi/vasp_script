#!/usr/bin/env python
# coding: utf-8
__author__ = 'wankw (wankaiweii@gmail.com)' 

from ase.neb import NEB
from ase.io import read
from ase.io import write
import numpy as np
from numpy.linalg import norm
import os

def check_move_far(initial,final):
    is_pos = initial.get_scaled_positions(wrap=False)
    fs_pos = final.get_scaled_positions(wrap=False)
    dist = norm(list(map(norm, (is_pos - fs_pos)@initial.get_cell()[:])))
    print (f'Distance between is&fs is {dist}!')

    n = 0
    far_id_list = np.argwhere(abs(is_pos-fs_pos)>0.5)

    for i in far_id_list:
        print (f'Movement of atom {i[0]+1} is too far! Check it! is: {is_pos[i[0]]} --> fs: {fs_pos[i[0]]}')
    return far_id_list

def wrap_atoms_by_id(far_id_list, init_atoms, fin_atoms):
    pos_is = init_atoms.get_scaled_positions(wrap=False)
    pos_fs = fin_atoms.get_scaled_positions(wrap=False)
    
    for i in far_id_list:
        if pos_is[i[0]][i[1]] > pos_fs[i[0]][i[1]]:
            pos_fs[i[0]][i[1]] = pos_fs[i[0]][i[1]]+1
        elif pos_is[i[0]][i[1]] < pos_fs[i[0]][i[1]]:
            pos_is[i[0]][i[1]] = pos_is[i[0]][i[1]]+1

    init_wrap_atoms = init_atoms.copy()
    init_wrap_atoms.set_scaled_positions(pos_is)

    fin_wrap_atoms = fin_atoms.copy()
    fin_wrap_atoms.set_scaled_positions(pos_fs)
    return init_wrap_atoms, fin_wrap_atoms

def linear_interpolation():
    images = [initial.copy() for i in range(n_image+1)] + [final]
    #linearly_interpolates
    neb = NEB(images)
    neb.interpolate()
    #generate POSCAR and movie.xyz
    movie = []
    for i in range(len(images)):
        if not os.path.isdir(f'{i:02}'):
            os.mkdir(f'{i:02}')
        images[i].write(f'{i:02}/POSCAR',vasp5=True,direct=True)
        images[i].write(f'{i:02}/POSCAR.xyz', format='xyz')
        with open(f'{i:02}/POSCAR.xyz') as xyz_file:
            line = xyz_file.readline()
            while line:
                movie.append(line)
                line = xyz_file.readline()
                
    with open('movie_line.xyz','w') as movie_xyz:
        movie_xyz.writelines(movie)

#Reference: S. Smidstrup, A. Pedersen, K. Stokbro and H. Jonsson, Improved initial guess for minimum energy path calculations, J. Chem. Phys. 140, 214106 (2014).
def idpp_interpolation():
    images = [initial.copy() for i in range(n_image+1)] + [final]
    #idpp_interpolates
    neb = NEB(images)
    neb.interpolate('idpp')
    #generate POSCAR and movie.xyz
    movie = []
    for i in range(len(images)):
        if not os.path.isdir(f'{i:02}'):
            os.mkdir(f'{i:02}')
        images[i].write(f'{i:02}/POSCAR',vasp5=True,direct=True)
        images[i].write(f'{i:02}/POSCAR.xyz', format='xyz')
        with open(f'{i:02}/POSCAR.xyz') as xyz_file:
            line = xyz_file.readline()
            while line:
                movie.append(line)
                line = xyz_file.readline()
                
    with open('movie_idpp.xyz','w') as movie_xyz:
        movie_xyz.writelines(movie)
        
def get_version():
    return '1.0 (2020.12.8, wankaiweii@gmail.com)'

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
    parser.add_argument('-w','--wrap_tag', type=bool, action='store', default=True,
                        help='Wrap positions or not. [Optional] [default=True]')
    
    
    args = parser.parse_args()
    
    initial_ini = read(rf'{args.initial_state_carfile}')
    final_ini = read(rf'{args.final_state_carfile}')
    n_image = args.number_of_images
    interpolation_method = args.interpolation_method
    wrap_tag = args.wrap_tag
   
    far_id_list = check_move_far(initial_ini, final_ini)
    if far_id_list.any() and wrap_tag:
        initial, final = wrap_atoms_by_id(far_id_list, initial_ini, final_ini)
        print ('\nAfter wrap!')
        far_id_list = check_move_far(initial, final)
    else:
        initial = initial_ini
        final = final_ini
 
    if interpolation_method == 'line':
        linear_interpolation()
    elif interpolation_method == 'idpp':
        idpp_interpolation()
        
    print (f'Generate {n_image} images between {args.initial_state_carfile} & {args.final_state_carfile} by {interpolation_method} interpolation!')
