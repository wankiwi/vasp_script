#!/usr/bin/env python
# coding: utf-8

__author__ = 'wankw (wankaiweii@gmail.com)' 

import re
import numpy as np
from plotille import *
from matplotlib import pyplot as plt
import matplotlib.ticker as ticker
import argparse

## plot convergence log in terminate
def term_plot(x_data,y_data,x_label,y_label):
    #formate: float -> int
    def _num_formatter(val, chars, delta, left=True):
        align = '<' if left else ''
        return '{:{}{}d}'.format(int(val), align, chars)
    #fig.register_label_formatter(float, _num_formatter)
    #fig.register_label_formatter(int, _num_formatter)
    
    fig = Figure()
    fig.width = 50
    fig.height = 25
    fig.set_x_limits(min(x_data), max(x_data))
    fig.set_y_limits(min(y_data), max(y_data))
    fig.x_label = x_label
    fig.y_label = y_label
    fig.color_mode = 'byte'
    #lc: line color
    fig.plot(x_data,y_data,lc=155)
    print(fig.show())


def mp_plot(x_data,y_data,x_label,y_label):        
    plt.xlabel(x_label) 
    plt.ylabel(y_label)
    plt.plot(x_data,y_data,'k--o',mfc='r')
    plt.gca().yaxis.set_major_formatter(ticker.FormatStrFormatter('%.5f'))
    plt.tight_layout()
    plt.savefig(f'conv-{y_label.split()[0]}.png',dpi=300)


def get_version():
    return '1.5 (2021.2.10, wankaiweii@gmail.com)'

parser = argparse.ArgumentParser(description='Plot the convergence curve in VASP calculation.')

parser.add_argument('-v', '--version', action='version', version=get_version(),help='display version')
parser.add_argument("-y", "--y_variable", default='f', choices=['f', 'e'],
                    help='the variable (force/energy) you want to plot [default=force]')
parser.add_argument("-m", "--plot_method", default='term', choices=['term', 'mp'],
                    help='plot method (in terminal/by matplotlib) you want to use [default=term]')
parser.add_argument("-n", "--last_n", default=0, type=int,
                    help='the number of last steps you want to plot [default=all]')
parser.add_argument("-l", "--log_mode", action='store_true', default=False,
                    help='generate the check_conv.log file [default=Flase]')
args = parser.parse_args()

y_variable = args.y_variable
plot_method = args.plot_method
last_n = args.last_n
log_mode = args.log_mode


#atom information list & selective dynamics list 
with open('POSCAR') as file:
    poscar = file.readlines()

a_vector = np.array(poscar[2].split(),dtype=float)
b_vector = np.array(poscar[3].split(),dtype=float)
c_vector = np.array(poscar[4].split(),dtype=float)

tmp_list1 = poscar[5].split()
tmp_list2 = list(map(int, poscar[6].split()))
atominfo_list = ''
for i in range(len(tmp_list1)):
    atominfo_list+=(f'{tmp_list1[i]} '*tmp_list2[i])
atominfo_list = atominfo_list.split()

selective_list = []
selec_tag = re.compile(r'[TF]')
if poscar[7].split()[0][0] == ('s' or 'S'):
    for i in poscar[9:9+len(atominfo_list)]:
        selective_list.append(selec_tag.findall(i))
else:
    for i in range(sum(tmp_list2)):
        selective_list.append(['T', 'T', 'T'])
selective_list_array = np.array(selective_list)
selective_list_array[selective_list_array=='T'] = 1
selective_list_array[selective_list_array=='F'] = 0
selective_list_array = np.array(selective_list_array, dtype=int)


with open ('OUTCAR') as file:
    outcar = file.read()
    
#EDIFFG
EDIFFG=re.findall('EDIFFG = (.*)   stopping-criterion for IOM', outcar)[0]
EDIFFG=float(EDIFFG)

#The number of atoms fixed in x y z direction
n_fix = len([i for i in selective_list_array if all(i == np.array([0,0,0]))])
n = np.shape(selective_list_array)[0] - n_fix

#force & position list
force_tag = re.compile(f'TOTAL-FORCE \(eV/Angst\)\n {"-"*83}\n(.*?)\n -', re.DOTALL)
force_array = [i.split() for i in force_tag.findall(outcar)]
force_array = np.array(force_array,dtype=float)
dim1,dim2 = force_array.shape
force_array = force_array[:, :, np.newaxis]
force_array = force_array.reshape((dim1, int(dim2/6), 6))

total_force_list = []
max_force_list = []
max_force_atom_seq_list = []
ave_force_list = []
for i in force_array:
    total_f_ionic_step = [(j[0]**2 + j[1]**2 + j[2]**2)**0.5 for j in (i[:,3:6] *  selective_list_array)]
    total_force_list.append(total_f_ionic_step)
    
    max_f = max(total_f_ionic_step)
    max_force_list.append(max_f)
    
    max_force_atom_seq = total_f_ionic_step.index(max_f) + 1
    max_force_atom_seq_list.append(max_force_atom_seq)
    
    ave_f = sum(total_f_ionic_step)/n
    ave_force_list.append(ave_f)
    
    
total_force_array = np.array(total_force_list, dtype=float)

#energy list
energy_tag = re.compile(f'total drift.*?energy\(sigma->0\) =(.*?)\n', re.DOTALL)
energy_array = np.array(energy_tag.findall(outcar), dtype=float)

if y_variable == 'f':
    y_data = max_force_list[-last_n:]
    x_data = list(range(1,len(max_force_list)+1))[-last_n:]
    x_label = 'Step'
    y_label = 'max_F (eV/A)'
elif y_variable == 'e':
    y_data = energy_array[-last_n:]
    x_data = list(range(1,len(energy_array)+1))[-last_n:]
    x_label = 'Step'
    y_label = 'Energy (eV)'
   

print(f"Last {len(x_data)} steps were plotted!")
if plot_method == 'term':
    term_plot(x_data,y_data,x_label,y_label)
elif plot_method == 'mp':
    mp_plot(x_data,y_data,x_label,y_label)
    print (f'conv-{y_label.split()[0]}.png generated by matplotlib!')    
print(f'{"="*90}')
if max_force_list[-1] >= abs(EDIFFG):
    print(f'[set EDIFFG:{EDIFFG:.2e} (eV/A), not converged!]')
else:
    print(f'[set EDIFFG:{EDIFFG:.2e} (eV/A), converged!]')    
print(f'After {len(max_force_list)} ionic steps, max force converged to  {max_force_list[-1]:.6f} at atom {max_force_atom_seq_list[-1]}, energy(sigma->0): {energy_array[-1]}.')
print(f'{n_fix} of {len(atominfo_list)} atoms were fixed, average force: {ave_force_list[-1]:.6f}.')


(step_num, atom_num) = np.shape(total_force_array)
with open('check_conv.log','w+') as log:
    for i in range(step_num):
        log.write(f'Ionic step: {i+1:>4}\n')
        log.write('-------atom-------||---------Position x y z----------------||-------------------Force x y z total------------\n')
        for j in range(atom_num):
            log.write(f'{j+1:>4}{atominfo_list[j]:>4}  ')
            [log.write(f'{k:3}') for k in selective_list[j]]
            [log.write(f'{k: 13.7f}') for k in force_array[i][j]]
            log.write(f'{total_force_array[i][j]: 12.7f}')
            log.write('\n')
        log.write('-------------------------------------------------------------------------------------------------------------\n')
        log.write(f'After {len(max_force_list)} ionic steps, max force converged to  {max_force_list[i]:.6f} at atom {max_force_atom_seq_list[i]}, energy(sigma->0): {energy_array[i]}.\n')
        log.write(f'{n_fix} of {len(atominfo_list)} atoms were fixed, average force: {ave_force_list[i]:.6f}.\n\n')
