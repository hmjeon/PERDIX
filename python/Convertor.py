# -*- coding: utf-8 -*-
"""
Convertor - From points and lines to faces

"""

from shapely.geometry import MultiLineString
from shapely.geometry import MultiPolygon
from shapely.ops import polygonize_full
from shapely.ops import cascaded_union

# Open file stream
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
print 'Points: ', n_point
for i in range(n_point):
    print i+1, ' th points : ', points[i]
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
print 'Lines:',  n_line
for i in range(n_line):
    print i+1, ' th lines : ', lines[i]
print '\n'

fin.close()

# ==================================================
# Make lines with points
# ==================================================
linepoints = []
for i in range(n_line):
    poi1, poi2 = lines[i]
    linepoint = [((points[poi1-1]),(points[poi2-1]))]
    linepoints.extend(linepoint)

# Print line list with points
print 'Lines with points: ', n_line
for i in range(n_line):
    print i+1, ' th lines with points : ', linepoints[i]
print '\n'

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
##################################################

polygon = MultiPolygon(result)

# Print polygon
print 'Polygons: ', len(polygon)
for i in range(len(polygon)):
    print i+1, ' th', polygon[i]
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
print 'Points: ', len(points)
for i in range(len(points)):
    print i+1, ' th point : ', points[i]
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
print 'Face connectivity: ', len(conns)
for i in range(len(conns)):
    print i+1, ' th connectivity : ', len(conns[i]), ' - ', conns[i]
print '\n'

# ==================================================
# Write file
# ==================================================
fout = open('output.txt', 'w')

fout.write('%d\t' % len(points))
fout.write('%d\n' % len(conns))
for i in range(len(points)):
    fout.write('%s \t %s \t' % points[i])
    fout.write('0.0\n')

for i in range(len(conns)):
    fout.write('%d \t' % len(conns[i]))
    for j in range(len(conns[i])):
        entity = conns[i][j] + 1
        fout.write('%d \t' %entity)
    
    fout.write('\n')
fout.close()