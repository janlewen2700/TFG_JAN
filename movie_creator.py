#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#"""
#Created on Sun Dec 24 12:45:47 2023

#@author: jan
#"""

import matplotlib as mpl
from matplotlib import pyplot
import numpy as np
import imageio.v2 as imageio
import time

time_step = 100001;
time=0;
data = []

# make a color map of fixed colors
cmap = mpl.colors.ListedColormap(['yellow','blue','white'])
bounds=[0.5,1.5,2.5,3.5]
norm = mpl.colors.BoundaryNorm(bounds, cmap.N)

img_color = [];

# make values from -5 to 5, for this example
while time<time_step:
    name = "./output/time" + str(int(time)) + ".txt";
    file = open(name, "r")
    comp = []
    for row in file:
        comp.append([int(x) for x in row.split()])
    # tell imshow about color map so that only set colors are used
    pyplot.pcolor(comp,cmap = cmap,norm=norm);
    # make a color bar
    #pyplot.colorbar(cmap=cmap,norm=norm,boundaries=bounds,ticks=[1,2,3,4]);
    name = "./movie/time" + str(int(time));
    pyplot.savefig( name + '.png');
    img_color.append(imageio.imread(name+'.png'));
    pyplot.close();
    time = time+1000;

print('\n only pitures, no movie output\n');
imageio.mimsave('./movie/time_evolution.gif', img_color, format='GIF', duration = 100);



