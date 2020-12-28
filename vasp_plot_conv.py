#!/usr/bin/env python
# coding: utf-8

__author__ = 'wankw (wankaiweii@gmail.com)' 

import xml.etree.cElementTree as et
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
    return '1.2 (2020.12.28, wankaiweii@gmail.com)'

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

last_n = args.last_n

#parse vasprun.xml
tree = et.parse('vasprun.xml')
root = tree.getroot()

#incar
EDIFFG = [i.text for i in root[[i.tag for i in root].index('incar')] 
            if i.attrib['name'] == 'EDIFFG'][0]       
EDIFFG = float(EDIFFG)

#selective dynamics list
selective_list = []
selective_tag_array = []
for j in [i for i in root.findall('structure/varray') if i.attrib['name'] == 'selective'][0]:
    selective_list.append(j.text)
    selective_tag_array.append(j.text.replace('T','1').replace('F','0').split()) 
selective_tag_array = np.array(selective_tag_array,dtype=int)

#The number of atoms fixed in x y z direction
n_fix = len([i for i in selective_tag_array if all(i == np.array([0,0,0]))])
n = np.shape(selective_tag_array)[0] - n_fix

#The number of atoms fixed in x y z direction
n_fix = len([i for i in selective_tag_array if all(i == np.array([0,0,0]))])
n = np.shape(selective_tag_array)[0] - n_fix

#force list
force_list = []
selective_force_list = []
max_force_list = []
max_force_atom_list = []
ave_force_list = []
for j in [i for i in root.findall('calculation/varray') if i.attrib['name'] == 'forces']:    
    force_array = [k.text.split() for k in j]
    force_array = np.array(force_array,dtype=float)
    selective_force_array = selective_tag_array * force_array      
    selective_force = [(i[0]**2 + i[1]**2 + i[2]**2)**0.5  for i in selective_force_array]

    max_f = max(selective_force)
    max_force_list.append(max_f)
    
    max_f_atom = selective_force.index(max_f) + 1
    max_force_atom_list.append(max_f_atom)
    
    ave_f = sum(selective_force)/n
    ave_force_list.append(ave_f)
    
    force_list.append(force_array)
    selective_force_list.append(selective_force)
    
force_array = np.array(force_list)
selective_force_array = np.array(selective_force_list)
selective_force_array = selective_force_array[:, :, np.newaxis]

force_log_array = np.c_[force_array, selective_force_array]

#energy list
energy_list = []
for j in [i for i in root.findall('calculation')]:
    energy_list.append(float([k.text for k in j.findall('scstep/energy/i')][-1]))


if args.y_variable == 'f':
    y_data = max_force_list[-last_n:]
    x_data = list(range(1,len(max_force_list)+1))[-last_n:]
    x_label = 'Step'
    y_label = 'max_F (eV/A)'
elif args.y_variable == 'e':
    y_data = energy_list[-last_n:]
    x_data = list(range(1,len(energy_list)+1))[-last_n:]
    x_label = 'Step'
    y_label = 'Energy (eV)'
   

print(f"Last {len(x_data)} steps were plotted!")
if args.plot_method == 'term':
    term_plot(x_data,y_data,x_label,y_label)
elif args.plot_method == 'mp':
    mp_plot(x_data,y_data,x_label,y_label)

print(f'{"="*90}')
if max_force_list[-1] >= abs(EDIFFG):
    print(f'[set EDIFFG:{EDIFFG:.2e} (eV/A), not converged!]')
else:
    print(f'[set EDIFFG:{EDIFFG:.2e} (eV/A), converged!]')    
print(f'After {len(max_force_list)} ionic steps, max force converged to  {max_force_list[-1]:.6f} at atom {max_force_atom_list[-1]}, energy(sigma->0): {energy_list[-1]}.')
print(f'{n_fix} of {np.shape(selective_tag_array)[0]} atoms were fixed, average force: {ave_force_list[-1]:.6f}.')


if args.log_mode:
    #force&energy log for every ionic step
    atominfo_list = [j[0].text for j in [i for i in root.findall('atominfo/array/set')][0]]
    
    pos_list = []
    for j in [i for i in root.findall('calculation/structure/varray')]:
        pos_list.append([k.text.split() for k in j])
        
    pos_array = np.array(pos_list, dtype=float)
    
    (step_num, atom_num, temp) = np.shape(pos_array)
    with open('check_conv.log','w+') as log:
        for i in range(step_num):
            log.write(f'Ionic step: {i+1:>4}\n')
            log.write('-------atom-------||---------Position x y z--------------||-------------------Force x y z total-----------\n')
            for j in range(atom_num):
                log.write(f'{j+1:>4}{atominfo_list[j]:>4}{selective_list[j]}  ')
                [log.write(f'{k: 12.7f}') for k in pos_array[i][j]]
                log.write('  ')
                [log.write(f'{k: 12.7f}') for k in force_log_array[i][j]]
                log.write('\n')
            log.write('----------------------------------------------------------------------------------------------------------\n')
            log.write(f'After {len(max_force_list)} ionic steps, max force converged to  {max_force_list[-1]:.6f} at atom {max_force_atom_list[-1]}, energy(sigma->0): {energy_list[-1]}.\n')
            log.write(f'{n_fix} of {np.shape(selective_tag_array)[0]} atoms were fixed, average force: {ave_force_list[-1]:.6f}.\n\n')
    print ('\ncheck_conv.log file generated successfully!')