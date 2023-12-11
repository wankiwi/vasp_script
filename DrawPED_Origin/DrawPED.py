#!/usr/bin/python3
# -*- coding: UTF-8 -*-
"""
@Description :   
@Author      :   KiwiWan 
@Email       :   wankaiweii@gmail.com
@Version     :   v1.0
@Time        :   2023/12/11 13:17:49
"""

import originpro as op
import numpy as np
import os
import sys
import matplotlib.pyplot as plt

bar_width = 2
spacing_between_bar = 4
nbins = 100*spacing_between_bar

class EnergyTagPair(object):
    class EnergyTag(object):
        def __init__(self, energy, tag):
            self.energy = energy
            self.tag = tag

        def __repr__(self):
            return f"{self.energy}:'{self.tag}'"

    def __init__(self, energy_list, tag_list):
        if len(energy_list) != len(tag_list):
            raise ValueError("Energy list and tag list must be of the same length.")
        self.pairs = [self.EnergyTag(energy, tag) for energy, tag in zip(energy_list, tag_list)]

    def __getitem__(self, index):
        return self.pairs[index]

    def __repr__(self):
        return '[' + ', '.join(repr(pair) for pair in self.pairs) + ']'

class PEDGenerator(object):
    def __init__(self, file_name):
        print ('Loading data...')
        self.file_name = file_name
        self.bar_width = bar_width
        self.spacing_between_bar = spacing_between_bar
        self.color_dict = {'r': '#FF0000', 'y': '#FFFF00', 'b': '#0000FF',
                           'g': '#00FF00', 'k': '#000000'}
        self.axis_name = {'x_axis':'Reaction coordinate',
                          'y_axis': 'Free Energy (eV)'}
        self.connect_type = 'cubic'
        
        # Extract data from file
        self.data = self.load_data(file_name)
        self.legend, self.color, self.energy_tag_list = self.parse_data()

        self.book = None
        self.gl_l = None
        self.gl_s = None

    def load_data(self, file_name):
        return np.loadtxt(file_name, delimiter=',', dtype=str)

    def set_x_axis_name(self, name_list):
        self.axis_name['x_axis'] = name_list[0]
        self.axis_name['y_axis'] = name_list[1]

    def set_connect_type(self, connect_type):
        self.connect_type = connect_type
    
    def parse_data(self):
        legend = []
        color = []
        energy_tag_list = []

        for i in enumerate(self.data):
            if i[0] % 2 == 0:
                if i[1][0] in self.color_dict:
                    color.append(self.color_dict[i[1][0]])
                else:
                    color.append(i[1][0])
                energy_list = i[1][1:].astype(float)
            else:
                legend.append(i[1][0])
                tag_list = i[1][1:]
                energy_tag_list.append(EnergyTagPair(energy_list, tag_list))
        return legend, color, energy_tag_list
    
    @staticmethod
    def cubic_fit(rcoord, energy, nbins=nbins):
        """ Build a cubic polynomial that satisfies
            f(x0) = y0
            f(x1) = y1
            f'(x0) = 0
            f'(x1) = 0.

        Args:
            rcoord (_type_): [x0, x1]
            energy (_type_): [y0, y1]
            nbins (int, optional): The number of insert points between x0 to x1. Defaults to 100.

        Returns:
            xfit, yfit: the points of fitted curve.
        """
        
        x0, y0 = rcoord[0], energy[0]
        x1, y1 = rcoord[1], energy[1]
        A = np.array([
                [x0**3,   x0**2, x0, 1],
                [x1**3,   x1**2, x1, 1],
                [3*x0**2,  2*x0,  1, 0],
                [3*x1**2,  2*x1,  1, 0]])
        b = np.array([y0, y1, 0, 0])
        coeff = np.linalg.solve(A, b)
        xfit = np.linspace(x0, x1, nbins)
        yfit = np.polyval(coeff, xfit)
        return xfit, yfit
    
    def gen_line_scatter_data(self, energy_tag):
        """_summary_

        Args:
            connect_type (str, optional): 'cubic' or 'line'. Defaults to 'cubic'.

        Returns:
            line: data for line and scatter plot.
        """
        energy_list = [i.energy for i in energy_tag]
        tag_list = [i.tag for i in energy_tag]
        scatter_data = np.zeros((len(energy_list), 2))
            
        if self.connect_type == 'line':
            for i in enumerate(energy_list):
                scatter_data[i[0]] = [i[0]*(self.spacing_between_bar+bar_width)+
                                    1/2*bar_width, i[1]]
            
            line_data = np.zeros((len(energy_list)*2, 2))
            for i in enumerate(energy_list):
                line_data[i[0]*2] = [i[0]*(self.spacing_between_bar+bar_width), i[1]]
                line_data[i[0]*2+1] = [(i[0])*(self.spacing_between_bar+bar_width)+
                                       bar_width, i[1]]
        
        elif self.connect_type == 'cubic':
            for i in enumerate(energy_list):
                scatter_data[i[0]] = [i[0]*(self.spacing_between_bar), i[1]]
            
            _connect_list = ['line'] * (len(tag_list)-1)
            for i in range(1, len(tag_list)-1):
                if 'ts' in tag_list[i].lower():
                    _connect_list[i-1] = 'cubic'
                    _connect_list[i] = 'cubic'
            if 'ts' in tag_list[-1].lower():
                raise ValueError("The last point can't be a transition state.")
            
            line_data = []
            for i in enumerate(_connect_list):
                _rcood = [scatter_data[i[0],0], scatter_data[i[0]+1,0]]
                _energy = [scatter_data[i[0],1], scatter_data[i[0]+1,1]]
                if i[1] == 'cubic':
                    xfit, yfit = self.cubic_fit(_rcood, _energy)
                    line_data.extend([i for i in zip(xfit, yfit)])
                elif i[1] == 'line':
                    line_data.extend([i for i in zip(_rcood, _energy)])
            line_data = np.array(line_data).reshape(-1,2)
        else:
            raise ValueError("connet_type must be 'cubic' or 'line'.")
        
        return line_data, scatter_data

    def save_data(self, file_name='PED.data'):
        print (f'Saving {file_name}...')
        data_file = open(file_name, 'w')
        for i in zip(self.legend, self.color, self.energy_tag_list):
            data_file.write(f'#{i[0]} {i[1]}\n')            
            axis_name = np.array(list([self.axis_name['x_axis'], self.axis_name['y_axis']]))
            axis_name = axis_name.reshape(1,2)
            tag_name = np.array([i.tag for i in i[2]])
            tag_name = np.hstack((['Tag name'], tag_name))
            tag_name = tag_name.reshape(-1,1)

            line_data, scatter_data = self.gen_line_scatter_data(i[2])
            line_data = np.vstack((axis_name, line_data))
            scatter_data = np.vstack((axis_name, scatter_data))
            scatter_data = np.hstack((scatter_data, tag_name))
            
            # Fill the line and scatter data to the same length
            max_rows = max(line_data.shape[0], scatter_data.shape[0])
            filled_line_data = np.full((max_rows, line_data.shape[1]), '', dtype=object)
            filled_scatter_data = np.full((max_rows, scatter_data.shape[1]), '', dtype=object)
            filled_line_data[:line_data.shape[0], :line_data.shape[1]] = line_data
            filled_scatter_data[:scatter_data.shape[0], :scatter_data.shape[1]] = scatter_data
            
            line_scatter_data = np.hstack((filled_line_data, filled_scatter_data))
            np.savetxt(data_file, line_scatter_data, delimiter=',', fmt='%s')
            data_file.write(f'{"#"*50}\n')
        data_file.close

    def plot_with_mpl(self, plot_kwargs={}, scatter_kwargs={}, file_name='PED.png'):
        print (f'Saving {file_name}...')
        plt.figure(figsize=(10, 6))
        plot_kwargs_default = {'linewidth': 2}
        scatter_kwargs_default = {'s': 600, 'marker': '_', 'c': 'k'}
        plot_kwargs = {**plot_kwargs_default, **plot_kwargs}
        scatter_kwargs = {**scatter_kwargs_default, **scatter_kwargs}
        
        for i in zip(self.legend, self.color, self.energy_tag_list):
            line_data, scatter_data = self.gen_line_scatter_data(i[2])            
            plt.plot(line_data[:, 0], line_data[:, 1], i[1], label=i[0], **plot_kwargs)          
            plt.scatter(scatter_data[:, 0], scatter_data[:, 1], **scatter_kwargs)
            
        plt.title('Potential Energy Diagram')
        plt.xlabel(self.axis_name['x_axis'])
        plt.ylabel(self.axis_name['y_axis'])
        plt.legend()

        plt.savefig(f'{file_name}', dpi=600)
  
    @staticmethod
    def origin_shutdown_exception_hook(exctype, value, traceback):
        op.exit()
        sys.__excepthook__(exctype, value, traceback)
              
    def plot_with_origin(self, file_name='PED.opju'):
        print (f'Saving {file_name}...')
        path = os.getcwd()
        if os.path.exists(rf'{path}\{file_name}'):
            os.remove(rf'{path}\{file_name}')
        
        if op and op.oext:
            sys.excepthook = self.origin_shutdown_exception_hook
        op.set_show(False)
        
        self.book = op.new_book(type='w', lname='PED data')
        sheet_list = []
        _n = 0
        for i in zip(self.legend, self.energy_tag_list):
            line_data, scatter_data = self.gen_line_scatter_data(i[1])
            if _n == 0:
                wks = self.book[0]
                wks.name = f'{i[0]}'
            elif _n > 0:
                wks = self.book.add_sheet(f'{i[0]}', active=False)
            wks.from_list(0, line_data[:, 0], lname=self.axis_name['x_axis'], axis='X')
            wks.from_list(1, line_data[:, 1], lname=self.axis_name['y_axis'], axis='Y')
            wks.from_list(2, scatter_data[:, 0], lname=self.axis_name['x_axis'], axis='X')
            wks.from_list(3, scatter_data[:, 1], lname=self.axis_name['y_axis'], axis='Y')
            wks.from_list(4, [i.tag for i in i[1]], lname='Tag name', axis='Y')
            sheet_list.append(wks)
            _n += 1
        
        self.graph = op.new_graph(lname='PED', template=op.path('e') + 'doubley.otp')
        self.gl_l = self.graph[0]
        self.gl_s = self.graph[1]
        
        for i in zip(self.color, sheet_list):
            p1 = self.gl_l.add_plot(i[1], colx='A', coly='B', type='line')
            p1.color = i[0]
            p2 = self.gl_s.add_plot(i[1], colx='C', coly='D', type='scatter')
            p2.color = i[0]
        self.gl_l.rescale()
        self.gl_s.rescale()
        
        legend = self.gl_l.label('Legend')
        legend_text = ''
        for i in range(len(self.legend)):
            legend_text += f'\l({i+1}) %({i+1},@WS)\n'
        legend_text = legend_text[:-2]
        legend.text = legend_text
        legend.set_int('left', 0)
        legend.set_int('showframe',0)
        
        op.save(rf'{path}\{file_name}')
        op.exit()

if __name__ == '__main__':
    print ('Initializing...')
    ped_generator = PEDGenerator('test.csv')
    ped_generator.save_data('PED.data')
    # ped_generator.plot_with_mpl(plot_kwargs={'linestyle': ':'}, scatter_kwargs={'marker': '_'}, file_name='PED.png')
    ped_generator.plot_with_origin('PED.opju')
    print ('DrawPED Finished!')
