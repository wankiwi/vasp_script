#!/usr/bin/env python
# coding: utf-8

from ase.neb import NEB
from ase.io import read
from ase.io import write
import os

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
    return '1.1 (2020.12.14, wankaiweii@gamil.com)'

if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description='Takes initial and final CARfiles, and linear or idpp method interpolation' 
                                     'the specified number of images between them. The interpolated files are written to the'
                                     'directories 00 to NI+1, where NI is the number of specified images.')
    
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
    
    if interpolation_method == 'line':
        linear_interpolation()
    elif interpolation_method == 'idpp':
        idpp_interpolation()
        
    print (f'Generate {n_image} images between {args.initial_state_carfile} & {args.final_state_carfile} by {interpolation_method} interpolation!')

