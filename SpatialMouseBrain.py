#!/usr/bin/env python
# coding: utf-8



import matplotlib as plt
import os, sys
import numpy as np
import pandas as pd
import csv
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D



#display spatial data
path = "/Users/Marie/Desktop/MOBDaten/"
dirs = os.listdir(path)
for file in dirs:
    if file.endswith(".csv"):
        print(file)


#load spatial data into one list
raw_data = []

for file in dirs:
    if file.endswith(".csv"):
        data = pd.read_csv(path + file, delimiter = '\t')
        raw_data.append(data)
        print(raw_data)


#name the column for the coordinates
length = len(raw_data)

for i in range(length):
    raw_data[i].rename( columns={'Unnamed: 0':'coordinates'}, inplace=True )



#make one vector with all genenames + coordinates
genenames = []

for i in range(length):
    names = raw_data[i].columns.tolist()
    genenames.extend(names)

genenames = list(set(genenames))



#seperate coordinates at x and store new coordinates in new arrays
xcoords_list = []

for i in range(length):
    array = np.zeros(shape=(len(raw_data[i]),1))
    array = [float(i.split('x')[0]) for i in raw_data[i]["coordinates"]]
    xcoords_list.append(array)

ycoords_list = []
    
for i in range(length):
    array = np.zeros(shape=(len(raw_data[i]),1))
    array = [float(i.split('x')[1]) for i in raw_data[i]["coordinates"]]
    ycoords_list.append(array)




#plot one slice in 2D

colors = ['red' if value > 0 else 'grey' for value in raw_data[10]["Penk"]]
plt.scatter(xcoords_list[10], ycoords_list[10], c = colors)



#plot one slice in 2D with colour gradient
plt.scatter(xcoords_list[10], ycoords_list[10], c=raw_data[10]["Penk"])
plt.colorbar()



#plot one slice in 3d
fig = plt.figure()
ax = Axes3D(fig)

p = ax.scatter(xcoords_list[10], ycoords_list[10], c=raw_data[10]["Penk"])
fig.colorbar(p)



#plot multiple slices in 3d
fig = plt.figure()
ax = Axes3D(fig)

s = raw_data[10]["Penk"]
t = raw_data[7]["Penk"]
r = raw_data[3]["Penk"]
col = s.append(t, ignore_index = True).append(r, ignore_index = True)

x = xcoords_list[10] + xcoords_list[7] + xcoords_list[3]
y = ycoords_list[10] + ycoords_list[7] + ycoords_list[3]
z = [1] * len(xcoords_list[10]) + [2] * len(xcoords_list[7]) + [3] * len(xcoords_list[3])

p = ax.scatter(x, y, z, c = col)
fig.colorbar(p)


#function for rotation and rescaling the layers

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
#rearrange graphs
x_1, y_1 = xcoords_list[10], ycoords_list[10]
x_2, y_2 = rotate_layer(xcoords_list[7], ycoords_list[7], 30, 2.5, -3.75)
x_3, y_3 = rotate_layer(xcoords_list[3], ycoords_list[3], 85)
x_3, y_3 = rotate_layer(x_3, y_3, 32, 3, -3)
x_4, y_4 = rotate_layer(xcoords_list[4], ycoords_list[4], 20, 0.75, -2.75)
x_5, y_5 = rotate_layer(xcoords_list[5], ycoords_list[5], 15, -1, 2)
x_6, y_6 = rotate_layer(xcoords_list[11], ycoords_list[11], 85)
x_6, y_6 = rotate_layer(x_6, y_6, 45, 1, -2)
x_7, y_7 = rotate_layer(xcoords_list[0], ycoords_list[0], 85)
x_7, y_7 = rotate_layer(x_7, y_7, 60, 1, -0.5)
x_8, y_8 = rotate_layer(xcoords_list[8], ycoords_list[8], 85)
x_8, y_8 = rotate_layer(x_8, y_8, 48, 1, 0)
x_9, y_9 = rotate_layer(xcoords_list[2], ycoords_list[2], 85)
x_9, y_9 = rotate_layer(x_9, y_9, 45, 2.5, -2)
x_10, y_10 = rotate_layer(xcoords_list[1], ycoords_list[1], 85)
x_10, y_10 = rotate_layer(x_10, y_10, 85)
x_10, y_10 = rotate_layer(x_10, y_10, 85)
x_10, y_10 = rotate_layer(x_10, y_10, 65, -5, -3.5)
x_11, y_11 = rotate_layer(xcoords_list[9], ycoords_list[9], 85)
x_11, y_11 = rotate_layer(x_11, y_11, 85)
x_11, y_11 = rotate_layer(x_11, y_11, 85)
x_11, y_11 = rotate_layer(x_11, y_11, 68, -1.5, -0.5)
x_12, y_12 = rotate_layer(xcoords_list[6], ycoords_list[6], 85)
x_12, y_12 = rotate_layer(x_12, y_12, 85)
x_12, y_12 = rotate_layer(x_12, y_12, 85)
x_12, y_12 = rotate_layer(x_12, y_12, 64, -2.75, -5.75)


i = xcoords_list[10] + list(x_2) + list(x_3) + list(x_4) + list(x_5) + list(x_6) + list(x_7) + list(x_8) + list(x_9) + list(x_10) + list(x_11) + list(x_12)
j = ycoords_list[10] + list(y_2) + list(y_3) + list(y_4) + list(y_5) + list(y_6) + list(y_7) + list(y_8) + list(y_9) + list(y_10) + list(y_11) + list(y_12)
#col = ['grey'] * len(xcoords_list[10]) + ['grey'] * len(xcoords_list[7]) + ['grey'] * len(xcoords_list[3]) + ['grey'] * len(xcoords_list[4]) + ['grey'] * len(xcoords_list[5]) + ['grey'] * len(xcoords_list[11]) +['grey'] * len(xcoords_list[0]) + ['grey'] * len(xcoords_list[8]) + ['grey'] * len(xcoords_list[2]) + ['grey'] * len(xcoords_list[1]) + ['grey'] * len(xcoords_list[9]) + ['red'] * len(xcoords_list[6])
s = raw_data[10]["Penk"]
col = s.append(raw_data[7]["Penk"], ignore_index = True).append(raw_data[3]["Penk"], ignore_index = True).append(raw_data[4]["Penk"], ignore_index = True).append(raw_data[5]["Penk"], ignore_index = True).append(raw_data[11]["Penk"], ignore_index = True).append(raw_data[0]["Penk"], ignore_index = True).append(raw_data[8]["Penk"], ignore_index = True).append(raw_data[2]["Penk"], ignore_index = True).append(raw_data[1]["Penk"], ignore_index = True).append(raw_data[9]["Penk"], ignore_index = True).append(raw_data[6]["Penk"], ignore_index = True)
plt.scatter(i, j, s = 10, c = col)



#plot multiple slices in 3d
fig = plt.figure()
ax = Axes3D(fig)

s = raw_data[10]["Penk"]
col = s.append(raw_data[7]["Penk"], ignore_index = True).append(raw_data[3]["Penk"], ignore_index = True).append(raw_data[4]["Penk"], ignore_index = True).append(raw_data[5]["Penk"], ignore_index = True).append(raw_data[11]["Penk"], ignore_index = True).append(raw_data[0]["Penk"], ignore_index = True).append(raw_data[8]["Penk"], ignore_index = True).append(raw_data[2]["Penk"], ignore_index = True).append(raw_data[1]["Penk"], ignore_index = True).append(raw_data[9]["Penk"], ignore_index = True).append(raw_data[6]["Penk"], ignore_index = True)

i = xcoords_list[10] + list(x_2) + list(x_3) + list(x_4) + list(x_5) + list(x_6) + list(x_7) + list(x_8) + list(x_9) + list(x_10) + list(x_11) + list(x_12)
j = ycoords_list[10] + list(y_2) + list(y_3) + list(y_4) + list(y_5) + list(y_6) + list(y_7) + list(y_8) + list(y_9) + list(y_10) + list(y_11) + list(y_12)
z = [1] * len(xcoords_list[10]) + [2] * len(xcoords_list[7]) + [3] * len(xcoords_list[3]) + [4] * len(xcoords_list[4]) + [5] * len(xcoords_list[5]) + [6] * len(xcoords_list[11]) + [7] * len(xcoords_list[0]) + [8] * len(xcoords_list[8]) + [9] * len(xcoords_list[2]) + [10] * len(xcoords_list[1]) + [11] * len(xcoords_list[9]) + [12] * len(xcoords_list[6])

p = ax.scatter(i, j, z, c = col)
fig.colorbar(p)



#make new lists with coordinates 
new_x = [x_1, x_2, x_3, x_4, x_5, x_6, x_7, x_8, x_9, x_10, x_11, x_12]
new_y = [y_1, y_2, y_3, y_4, y_5, y_6, y_7, y_8, y_9, y_10, y_11, y_12]

#rearrange dataframes
new_data = [pd.DataFrame(raw_data[10]), pd.DataFrame(raw_data[7]), pd.DataFrame(raw_data[3]), pd.DataFrame(raw_data[4]), pd.DataFrame(raw_data[5]), pd.DataFrame(raw_data[11]), pd.DataFrame(raw_data[0]), pd.DataFrame(raw_data[8]), pd.DataFrame(raw_data[2]), pd.DataFrame(raw_data[1]), pd.DataFrame(raw_data[9]),pd.DataFrame(raw_data[6])]

genenames.remove("coordinates")



#create sqlite database for genedata in VR
import sqlite3
import os
os.remove('database.sqlite')
conn = sqlite3.connect('database.sqlite')
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
for cname in range(3116):
    c.execute("insert into cells values ("+ str(cname) + ",'" + str(cname) + "');")
conn.commit() 


#calculate total cellnumber
totlength = 0
for i in range(length):
    totlength += len(new_x[i])
print(totlength)    


#make a dictionary
cellindices = dict()
index = 0
for i in range(length):
    for j in range(len(new_data[i])):
        cellindices[index] = (i,j)
        index += 1
        


#fill database with expression values
for cellindex in range(3116):
    d = cellindices[cellindex]
    #print(type(correctDataFrame))
    correctDataFrame = new_data[d[0]]
    row = correctDataFrame.iloc[d[1]]
    for gname in zip(genenames, range(len(genenames))):
        if gname[0] in row:
            expression = row[gname[0]]
            if expression > 0:
                c.execute("insert into datavalues values ("+ str(gname[1]) + "," + str(cellindex) + "," + str(expression) +");")
           
            
conn.commit()



CellID = range(3116)
complete_x = x_1 + list(x_2) + list(x_3) + list(x_4) + list(x_5) + list(x_6) + list(x_7) + list(x_8) + list(x_9) + list(x_10) + list(x_11) + list(x_12)
complete_y = y_1 + list(y_2) + list(y_3) + list(y_4) + list(y_5) + list(y_6) + list(y_7) + list(y_8) + list(y_9) + list(y_10) + list(y_11) + list(y_12)
complete_z = [1] * len(xcoords_list[10]) + [2] * len(xcoords_list[7]) + [3] * len(xcoords_list[3]) + [4] * len(xcoords_list[4]) + [5] * len(xcoords_list[5]) + [6] * len(xcoords_list[11]) + [7] * len(xcoords_list[0]) + [8] * len(xcoords_list[8]) + [9] * len(xcoords_list[2]) + [10] * len(xcoords_list[1]) + [11] * len(xcoords_list[9]) + [12] * len(xcoords_list[6])
data = {'cell_id': CellID, 'dim_1': complete_x, 'dim_2': complete_y, 'dim_3': complete_z}
#data = [CellID, complete_x, complete_y, complete_z]
Coordinates = pd.DataFrame(data)
Coordinates.to_csv(r'/Users/Marie/Desktop/export_coordinates.csv', index = None, header=True, sep = " ")


