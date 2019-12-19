#!/usr/bin/env python
# coding: utf-8


import matplotlib as plt
import os, sys
import numpy as np
import pandas as pd
import csv
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D


#load spatial data
path = "/Users/Marie/Desktop/"
dirs = os.listdir(path)
for file in dirs:
    if file.endswith(".csv"):
        raw_data = pd.read_csv(path + file, delimiter = '\t')
      


#name the column for the genes
raw_data.rename( columns={'Unnamed: 0':'genes'}, inplace=True )


#reshape data
raw_data = raw_data.transpose()
raw_data.columns = raw_data.iloc[0]
raw_data = raw_data.drop("genes")



#make one vector with all genenames + coordinates
genenames = raw_data.columns



#make list with coordinates
coords = []
for i in raw_data.index: 
    coords.append(i)


#seperate coordinates at x and store new coordinates in new arrays

for i in range(len(coords)):
    array = np.zeros(shape=(len(coords),1))
    zcoords = [float(i.split('x')[0]) for i in coords]

for i in range(len(coords)):
    array = np.zeros(shape=(len(coords),1))
    xcoords = [float(i.split('x')[1]) for i in coords]

for i in range(len(coords)):
    array = np.zeros(shape=(len(coords),1))
    ycoords = [float(i.split('x')[2]) for i in coords]
    


#make dataframe with coordinates
coords = pd.DataFrame(list(zip(xcoords, ycoords, zcoords)), 
               columns =['x', 'y', 'z'])



#plot multiple slices in 3d
fig = plt.figure()
ax = Axes3D(fig)

col = raw_data["DPM1"]

x = xcoords
y = ycoords
z = zcoords

p = ax.scatter(x, y, z, c = col)
fig.colorbar(p)


#subset raw_data and coords depending on developmental stage
data1 = raw_data.iloc[0:238]
data2 = raw_data.iloc[238:1753]
data3 = raw_data.iloc[1753:3111]

coords1 = coords.iloc[0:238]
coords2 = coords.iloc[238:1753]
coords3 = coords.iloc[1753:3111]

#function to rotate the layers

def rotate_quadrant(x):
    if x > np.pi:
        return 2 * np.pi - x
    elif x < -np.pi:
        return 2 * np.pi + x
    else:
        return x
    
def invert_x(tu):
    if tu[0] > np.pi / 2 or tu[0] < -np.pi / 2:
        return -tu[1]
    else:
        return tu[1]
    
def invert_y(tu):
    if tu[0] < 0:
        return -tu[1]
    else:
        return tu[1]


def rotate_layer(xcoords, ycoords, degrees, translateX=0, translateY=0):
    radians = degrees / 180 * np.pi
    meanx = np.mean(xcoords)
    meany = np.mean(ycoords)
    x = xcoords - meanx
    y = ycoords - meany
    r = np.sqrt(np.power(x, 2) + np.power(y, 2))
    t = np.arctan2(y,x)-radians
    x = r/np.sqrt(np.power(np.tan(t), 2) + 1)
    y = r/np.sqrt(1/np.power(np.tan(t), 2) + 1)
    t = np.array(list(map(rotate_quadrant, t)))
    x = np.array(list(map(invert_x, zip(t, x))))
    y = np.array(list(map(invert_y, zip(t, y))))
    x += meanx + translateX 
    y += meany + translateY
    return x, y



#plot multiple slices in 2D with colour gradient
#rearrange graphs for first developmental stage
x_1, y_1 = coords["x"][0:55], coords["y"][0:55]
x_2, y_2 = rotate_layer(coords["x"][55:114], coords["y"][55:114], 0, 8, 2)
x_3, y_3 = rotate_layer(coords["x"][114:183], coords["y"][114:183], 0, 4, 3)
x_4, y_4 = rotate_layer(coords["x"][183:237], coords["y"][183:237], 0, 6, 4)

i = list(x_1) + list(x_2) + list(x_3) + list(x_4) 
j = list(y_1) + list(y_2) + list(y_3) + list(y_4)
col = ['grey'] * len(coords["x"][0:55]) + ['red'] * len(coords["x"][55:114]) + ['blue'] * len(coords["x"][114:183]) + ['yellow'] * len(coords["x"][183:237])
#col = data1["DPM1"]
plt.scatter(i, j, s = 10, c = col)


#plot multiple slices in 2D with colour gradient
#rearrange graphs for second developmental stage
x_5, y_5 = rotate_layer(coords["x"][238:343], coords["y"][238:343], 80)
x_6, y_6 = rotate_layer(coords["x"][343:447], coords["y"][343:447], 80, -2, 1.5)
x_7, y_7 = rotate_layer(coords["x"][447:606], coords["y"][447:606], 80, 2, 6)
x_8, y_8 = rotate_layer(coords["x"][606:792], coords["y"][606:792], 80, -0.5, 4)
x_9, y_9 = rotate_layer(coords["x"][792:1006], coords["y"][792:1006], 78, -6, 11)
x_10, y_10 = rotate_layer(coords["x"][1006:1215], coords["y"][1006:1215], 80, 2.5, 10)
x_11, y_11 = rotate_layer(coords["x"][1215:1393], coords["y"][1215:1393], 78, -7, 12)
x_12, y_12 = rotate_layer(coords["x"][1393:1575], coords["y"][1393:1575], 78, -7, 13)
x_13, y_13 = rotate_layer(coords["x"][1575:1752], coords["y"][1575:1752], 70, -3, 13)

i = list(x_5) + list(x_6) + list(x_7) + list(x_8) + list(x_9) + list(x_10) + list(x_11) + list(x_12) + list(x_13)
j = list(y_5) + list(y_6) + list(y_7) + list(y_8) + list(y_9) + list(y_10) + list(y_11) + list(y_12) + list(y_13)
col = ['grey'] * len(coords["x"][238:343]) + ['grey'] * len(coords["x"][343:447]) + ['grey'] * len(coords["x"][447:606]) + ['grey'] * len(coords["x"][606:792]) + ['grey'] * len(coords["x"][792:1006]) + ['grey'] * len(coords["x"][1006:1215]) + ['grey'] * len(coords["x"][1215:1393]) + ['grey'] * len(coords["x"][1393:1575]) + ['red'] * len(coords["x"][1575:1752])
#col = data2["DPM1"]
plt.scatter(i, j, s = 10, c = col)



#plot multiple slices in 2D with colour gradient
#rearrange graphs for second developmental stage
x_14, y_14 = rotate_layer(coords["x"][1754:1992], coords["y"][1754:1992], 60)
x_15, y_15 = rotate_layer(coords["x"][1992:2234], coords["y"][1992:2234], 58, 1, -1.5)
x_16, y_16 = rotate_layer(coords["x"][2234:2477], coords["y"][2234:2477], 63, 3, -3)
x_17, y_17 = rotate_layer(coords["x"][2477:2703], coords["y"][2477:2703], 60, 2, -3.5)
x_18, y_18 = rotate_layer(coords["x"][2703:2901], coords["y"][2703:2901], 60, 6, 1)
x_19, y_19 = rotate_layer(coords["x"][2901:3111], coords["y"][2901:3111], 65, 6, -2)

i = list(x_14) + list(x_15) + list(x_16) + list(x_17) + list(x_18) + list(x_19) 
j = list(y_14) + list(y_15) + list(y_16) + list(y_17) + list(y_18) + list(y_19) 
col = ['grey'] * len(coords["x"][1754:1992]) + ['grey'] * len(coords["x"][1992:2234]) + ['grey'] * len(coords["x"][2234:2477]) + ['grey'] * len(coords["x"][2477:2703]) + ['grey'] * len(coords["x"][2703:2901]) + ['red'] * len(coords["x"][2901:3111]) 
#col = data2["DPM1"]
plt.scatter(i, j, s = 10, c = col)


#make new lists with coordinates 
new_x1 = list(x_1) + list(x_2) + list(x_3) + list(x_4) 
new_y1 = list(y_1) + list(y_2) + list(y_3) + list(y_4) 
new_z1 = [1] * len(coords["x"][0:55]) + [2] * len(coords["x"][55:114]) + [3] * len(coords["x"][114:183]) + [4] * len(coords["x"][183:237])

new_x2 = list(x_5) + list(x_6) + list(x_7) + list(x_8) + list(x_9) + list(x_10) + list(x_11) + list(x_12) + list(x_13) 
new_y2 = list(y_5) + list(y_6) + list(y_7) + list(y_8) + list(y_9) + list(y_10) + list(y_11) + list(y_12) + list(y_13) 
new_z2 = [1] * len(coords["x"][238:343]) + [2] * len(coords["x"][343:447]) + [3] * len(coords["x"][447:606]) + [4] * len(coords["x"][606:792]) + [5] * len(coords["x"][792:1006]) + [6] * len(coords["x"][1006:1215]) + [7] * len(coords["x"][1215:1393]) + [8] * len(coords["x"][1393:1575]) + [9] * len(coords["x"][1575:1752])

new_x3 = list(x_14) + list(x_15) + list(x_16) + list(x_17) + list(x_18) + list(x_19) 
new_y3 = list(y_14) + list(y_15) + list(y_16) + list(y_17) + list(y_18) + list(y_19) 
new_z3 = [1] * len(coords["x"][1754:1992]) + [2] * len(coords["x"][1992:2234]) + [3] * len(coords["x"][2234:2477]) + [4] * len(coords["x"][2477:2703]) + [5] * len(coords["x"][2703:2901]) + [6] * len(coords["x"][2901:3111]) 


#create sqlite database for genedata in VR
import sqlite3
import os
os.remove('databaseHeart.sqlite')
conn = sqlite3.connect('databaseHeart.sqlite')
c = conn.cursor()
c.execute("CREATE TABLE datavalues ('gene_id' REAL, 'cell_id' REAL, 'value' REAL);")
c.execute("CREATE TABLE genes ('id' INTEGER NOT NULL UNIQUE,'gname' varchar(20) COLLATE NOCASE);")
c.execute("CREATE TABLE cells ('id' INTEGER NOT NULL UNIQUE,'cname' varchar(20) COLLATE NOCASE);")
c.execute("CREATE UNIQUE INDEX gnameIDX ON genes (gname);")
c.execute("CREATE UNIQUE INDEX cnameIDX ON cells (cname);")
c.execute("CREATE INDEX gene_id_data ON datavalues ('gene_id');")
c.execute("CREATE INDEX cell_id_data ON datavalues ('cell_id');")
conn.commit()



#fill database with genenames
for gname in zip(genenames, range(len(genenames))):
    c.execute("insert into genes values ("+ str(gname[1]) + ",'" + gname[0] + "');") # insert into genes values (0, 'Gm25252');
conn.commit()


#fill database with cellnames
for cname in range(3111):
    c.execute("insert into cells values ("+ str(cname) + ",'" + str(cname) + "');")
conn.commit() 


#fill database with expression values
for cellindex in range(3111):
    row = raw_data.iloc[cellindex]
    for gname in zip(genenames, range(len(genenames))):
        if gname[0] in row:
            expression = row[gname[0]]
            if expression > 0:
                c.execute("insert into datavalues values ("+ str(gname[1]) + "," + str(cellindex) + "," + str(expression) +");")
           
            
conn.commit()


CellID1 = raw_data.index[range(0,237)]
CellID2 = raw_data.index[range(238,1752)]
CellID3 = raw_data.index[range(1753,3110)]

data1 = {'cell_id': CellID1, 'dim_1': new_x1, 'dim_2': new_y1, 'dim_3': new_z1}
data2 = {'cell_id': CellID2, 'dim_1': new_x2, 'dim_2': new_y2, 'dim_3': new_z2}
data3 = {'cell_id': CellID3, 'dim_1': new_x3, 'dim_2': new_y3, 'dim_3': new_z3}

Coordinates1 = pd.DataFrame(data1)
Coordinates2 = pd.DataFrame(data2)
Coordinates3 = pd.DataFrame(data3)
Coordinates1.to_csv(r'/Users/Marie/Desktop/export_coordinates1.csv', index = None, header=True, sep = " ")
Coordinates2.to_csv(r'/Users/Marie/Desktop/export_coordinates2.csv', index = None, header=True, sep = " ")
Coordinates3.to_csv(r'/Users/Marie/Desktop/export_coordinates3.csv', index = None, header=True, sep = " ")




