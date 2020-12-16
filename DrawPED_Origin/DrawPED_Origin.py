#!/usr/bin/env python
# coding: utf-8

import originpro as op
import pandas as pd
import os


def line_scateer_plot(i, path_line_pd, path_scatter_pd, color):
    #generate worksheet
    wks = book[i]
    wks.name = f'path{i+1}'
    #superscript: \+(text), subcript: \-(text)
    wks.from_list(0, path_line_pd.iloc[:,0], lname='Reaction coordinate', axis = 'X')
    wks.from_list(1, path_line_pd.iloc[:,1], lname='Energy/eV', axis = 'Y')
    wks.from_list(2, path_scatter_pd.iloc[:,0], lname='Reaction coordinate', axis = 'X')
    wks.from_list(3, path_scatter_pd.iloc[:,1], lname='Energy/eV', axis = 'Y')
    wks.from_list(4, path_scatter_pd.iloc[:,2], lname='Name', axis = 'Y')      
    
    #plot potential energy diagram
    p1 = gl_1.add_plot(wks, colx='A', coly='B', type='l')
    p1.color = color   
    gl_1.rescale()
    p2 = gl_2.add_plot(wks, colx='C', coly='D', type='s')
    p2.color = color
    gl_2.rescale()


def generate_line_scateer_pd(i, data_pd):
    color = data_pd.iloc[(i+1)*2-2, 0]
    energy_pd = data_pd.iloc[(i+1)*2-2, 1:]
    name_pd = data_pd.iloc[(i+1)*2-1, 1:]
    step_num = int(len(energy_pd))
    path_line_pd = pd.DataFrame(columns=['react_cord', 'energy'])
    path_scatter_pd = pd.DataFrame(columns=['react_cord', 'energy','name'])
    for j in range(step_num):
        path_line_pd.loc[2*j] = [j*(bar_width+spacing_between_bar), energy_pd.iloc[j]]
        path_line_pd.loc[2*j+1] = [j*(bar_width+spacing_between_bar)+bar_width, energy_pd.iloc[j]]
        path_scatter_pd.loc[j] = [j*(bar_width+spacing_between_bar)+(bar_width)/2, energy_pd.iloc[j], name_pd.iloc[j]]
    return path_line_pd, path_scatter_pd, color


if __name__ == '__main__':
    bar_width = 2
    spacing_between_bar = 4
    data_pd = pd.read_excel('test.xlsx',header=None)
    
    #support hex color code, https://www.sojson.com/rgb.html
    color_dict = {'r':'#FF0000', 'y':'#FFFF00' , 'b':'#0000FF', 'g':'#00FF00', 'k':'#000000' }    
    
    # Very useful, especially during development, when you are
    # liable to have a few uncaught exceptions.
    # Ensures that the Origin instance gets shut down properly.
    import sys
    def origin_shutdown_exception_hook(exctype, value, traceback):
        '''Ensures Origin gets shut down if an uncaught exception'''
        op.exit()
        sys.__excepthook__(exctype, value, traceback)
    if op and op.oext:
        sys.excepthook = origin_shutdown_exception_hook
    
    
    # Set Origin instance visibility.
    op.set_show(False)
    
    path_num = int(data_pd.shape[0]/2)
    
    graph = op.new_graph(lname='PED', template=op.path('e') + 'doubley.otp')
    #graph layer 
    gl_1 = graph[0]
    gl_2 = graph[1]
    
    book = op.new_book(type='w', lname='PED data')
    
    for i in range(path_num):
        if i > 0:
            book.add_sheet(f'path{i+1}_scatter', active = False)
        path_line_pd, path_scatter_pd, color = generate_line_scateer_pd(i, data_pd)
        if color in color_dict.keys(): 
            color = color_dict[color]
        line_scateer_plot(i, path_line_pd, path_scatter_pd, color)
    
    # Save the opju to your UFF.
    path = os.getcwd()
    op.save(rf'{path}\PED.opju')
    op.exit()
    
    print ('PED.opju generated successfully!')
