# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""
import sys
import numpy
from shapely.geometry import LineString, MultiLineString

from matplotlib import pyplot as plt
from shapely.ops import polygonize, polygonize_full
from shapely.geometry import Polygon, MultiPolygon, mapping
from shapely.ops import linemerge, cascaded_union
from descartes import PolygonPatch

# References
# http://stackoverflow.com/questions/34475431/plot-unions-of-polygons-in-matplotlib
# http://toblerity.org/shapely/manual.html#MultiLineString
# http://stackoverflow.com/questions/21824157/how-to-extract-interior-polygon-coordinates-using-shapely

'''
1	0	1	4	1 
2	0	3	4	3 
3 	0	1	2	4
4	2	0	0	3
5	2	4	4	1
6 	2	0	4	3

1	0	0	1	0
2	1	0	2	0 
3 	0	1	1	1
4	1	1	2	1
5	0	0	0	1
6 	1	0	1	1
7 	2	0	2	1
'''

fin = open('geometry.txt', 'r')

# ==================================================
# Read points
# ==================================================
n_point = int(fin.readline())

points = []
for i in range(n_point):
    str = fin.readline()
    point = str.split("\t")
    points.append([float(point[1]), float(point[2])])

# Print points
print 'Points'
for i in range(n_point):
    print i, ' th points : ', points[i]
print '\n'

# ==================================================
# Read lines
# ==================================================
n_line = int(fin.readline())

lines = []
for i in range(n_line):
    str = fin.readline()
    line = str.split("\t")
    lines.append([int(line[1]), int(line[2])])

# Print lines
print 'Lines'
for i in range(n_line):
    print i, ' th lines : ', lines[i]
print '\n'

# ==================================================
# Make lines with points
# ==================================================
linepoints = []
for i in range(n_line):
    poi1, poi2 = lines[i]
    linepoint = [((points[poi1-1]),(points[poi2-1]))]
    linepoints.extend(linepoint)

# Print line list with points
print 'Lines with points'
for i in range(n_line):
    print i, ' th lines with points : ', linepoints[i]
print '\n'

'''
# Read geometry file
a = numpy.loadtxt('geometry.txt', usecols = [1,2,3,4])

# Make line list 
lines = []
for i in range(len(a)):
    line = [((a[i,0],a[i,1]),(a[i,2],a[i,3]))]
    lines.extend(line)
    print i, ' th line : ', line
'''

# ==================================================
# Make mutilinestring from line list
# ==================================================
multilines = MultiLineString(linepoints)

# Polygonize
result, dangles, cuts, invalids = polygonize_full(multilines)

result = MultiPolygon(result)
polygon = cascaded_union(result)

# Make mutilinestring from line list
multilines = MultiLineString(linepoints)

# Polygonize
result, dangles, cuts, invalids = polygonize_full(multilines)

result = MultiPolygon(result)
polygon = cascaded_union(result)

##################################################
multilines = polygon.boundary.union(result.boundary)

# Polygonize
result, dangles, cuts, invalids = polygonize_full(multilines)

polygon = MultiPolygon(result)
##################################################

# Print polygon
print 'Polygons - ', len(polygon)
for i in range(len(polygon)):
    print i, ' th point : ', polygon[i]
print '\n'

# ==================================================
# Extract points
# ==================================================
points = []
for i in range(len(polygon)):
    point = list(polygon[i].exterior.coords)
    for j in range(len(point)):
        x = point[j]
        if x not in points:
            points.extend([(x[0], x[1])])

# Print points
print 'Points - ', len(points)
for i in range(len(points)):
    print i, ' th point : ', points[i]
print '\n'

# ==================================================
# Face connectivity
# ==================================================
conns = []
for i in range(len(polygon)):
    conn  = polygon[i].exterior.coords
    count = len(polygon[i].exterior.coords)
    face  = []
    for j in range(count-1):
        face.append(points.index(conn[j]))

    # Make connectivity
    conns.append(face)

# Print connectivity
print 'Connectivity - ', len(conns)
for i in range(len(conns)):
    print i, ' th connectivity : ', conns[i]
print '\n'

# ==================================================
# ==================================================
# ==================================================
sys.exit()
# ==================================================

# Print points and lines
print '#########################################'

points = []
for i in range(len(polygon)):
    #print i, polygon[i]
    x, y = polygon[i].exterior.coords.xy
    poly = polygon[i]
    interior = list(polygon[i].interiors)
    exterior = list(poly.exterior.coords)
    print i, ' th polygon : ', exterior
    aaa = list(exterior[0])
    #print aaa[0], aaa[1]

#m = mapping(polygon)
#print m['coordinates']
sys.exit()

# Plot polygon graph
# http://stackoverflow.com/questions/26935701/ploting-filled-polygons-in-python

ring_patch = PolygonPatch(result)

xrange = [-5, 10]
yrange = [-5, 10]

fig = plt.figure(1, figsize=(5,5), dpi=90)
ax  = fig.add_subplot(111)

ax.add_patch(ring_patch)
ax.set_title('Filled Polygon')
ax.set_xlim(*xrange)
ax.set_ylim(*yrange)
ax.set_aspect(1)