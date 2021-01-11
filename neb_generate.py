#!/usr/bin/env python
# coding: utf-8
__author__ = 'wankw (wankaiweii@gamil.com)' 

from ase.neb import NEB
from ase.io import read
from ase.io import write
import os

def check_move_far(initial,final):
    is_pos = initial.get_positions()
    fs_pos = final.get_positions()
    cell_vec = initial.get_cell_lengths_and_angles()[0:3]
    n = 0
    for i in abs(is_pos - fs_pos) >  0.5 * cell_vec:
        if i.any():
            print (f'Movement of atom {n+1} is too far! Check it! is: {is_pos[n]} --> fs: {fs_pos[n]}')
        n += 1

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
        images[i].write(f'{i:02}/POSCAR.xyz')
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
        images[i].write(f'{i:02}/POSCAR.xyz')
        with open(f'{i:02}/POSCAR.xyz') as xyz_file:
            line = xyz_file.readline()
            while line:
                movie.append(line)
                line = xyz_file.readline()
                
    with open('movie_idpp.xyz','w') as movie_xyz:
        movie_xyz.writelines(movie)
        
def get_version():
    return '1.0 (2020.12.8, wankaiweii@gamil.com)'

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
    
    check_move_far(initial,final)

    if interpolation_method == 'line':
        linear_interpolation()
    elif interpolation_method == 'idpp':
        idpp_interpolation()
        
    print (f'Generate {n_image} images between {args.initial_state_carfile} & {args.final_state_carfile} by {interpolation_method} interpolation!')
