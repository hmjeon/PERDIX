"""
!
! =============================================================================
!
! Module - Shapely
! Last Updated : 04/28/2018, by Hyungmin Jun (hyungminjun@outlook.com)
!
! =============================================================================
!
! This is part of PERDIX-2L, which allows scientists to build and solve
! the sequence design of complex DNAnanostructures.
! Copyright 2018 Hyungmin Jun. All rights reserved.
!
! License - GPL version 3
! PERDIX-2L is free software: you can redistribute it and/or modify it under
! the terms of the GNU General Public License as published by the Free Software
! Foundation, either version 3 of the License, or any later version.
! PERDIX-2L is distributed in the hope that it will be useful, but WITHOUT
! ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
! FOR A PARTICULAR PURPOSE. See the GNU General Public License
! for more details.
! You should have received a copy of the GNU General Public License along with
! this program. If not, see <http://www.gnu.org/licenses/>.
!
! -----------------------------------------------------------------------------
!
"""

import sys
from shapely.geometry import MultiLineString
from shapely.geometry import MultiPolygon
from shapely.ops import polygonize_full
from shapely.ops import cascaded_union

# Open file stream
if len(sys.argv) is 1:
    fin = open('test.igs', 'r')
    filename = 'test'
    filetype = 'igs'
if len(sys.argv) is 2:
    fin = open(sys.argv[1], 'r')
    filename, filetype = sys.argv[1].split('.')

print '\nFilename: ', filename, '\nFiletype: ', filetype, '\n'

# ==================================================
# Read GEO data
# ==================================================
if filetype == 'geo':
    str = fin.readline()
    str = str.split("\t")

    n_point = int(str[0])
    n_line  = int(str[1])
    n_face  = int(str[2])

    # ==================================================
    # Read points
    # ==================================================
    points = []
    for i in range(n_point):
        str = fin.readline()
        point = str.split()
        points.append([float(point[1]), -float(point[2])])

    # Print points
    print 'Points: ', n_point
    for i in range(n_point):
        print i+1, ' th points : ', points[i]
    print '\n'

    # ==================================================
    # Read lines
    # ==================================================
    lines = []
    for i in range(n_line):
        str = fin.readline()
        line = str.split()
        lines.append([int(line[1]), int(line[2])])
        
    # Print lines
    print 'Lines:',  n_line
    for i in range(n_line):
        print i+1, ' th lines : ', lines[i]
    print '\n'
    
    # ==================================================
    # Make lines with points
    # ==================================================
    linepoints = []
    for i in range(n_line):
        poi1, poi2 = lines[i]
        linepoint = [((points[poi1-1]),(points[poi2-1]))]
        linepoints.extend(linepoint)

# ==================================================
# Read IGES file
# ==================================================
if filetype == 'igs' or filetype == 'iges':

    linepoints = []
    while True:
        str = fin.readline()
        if str == '':
            break
        str = str.split(',')

        # Add lines
        points = []
        if str[0] == '110':

            # no count first and last items
            index = len(str) - 2
            for i in (range(1,len(str)-1)):
                if i == 5:
                    split0 = str[5].split(';') 
                    points.append(float(split0[0]))
                else:
                    points.append(float(str[i]))
            else:
                if index is 3:
                    # 2 items are below
                    str = fin.readline()
                    str = str.split(',')
                    split0 = str[2].split(';')
                    points.append(float(str[0]))
                    points.append(float(str[1]))
                    points.append(float(split0[0]))
                if index is 4:
                    # 1 item is below
                    str = fin.readline()
                    str = str.split(',')
                    split0 = str[1].split(';')
                    points.append(float(str[0]))
                    points.append(float(split0[0]))

            # Set zero
            if abs(points[0]) < 0.0000001:
                points[0] = 0.0
            if abs(points[1]) < 0.0000001:
                points[1] = 0.0
            if abs(points[2]) < 0.0000001:
                points[2] = 0.0
            if abs(points[3]) < 0.0000001:
                points[3] = 0.0

            linepoint = [((points[0], points[1]), (points[3], points[4]))]
            linepoints.extend(linepoint)

    n_line = len(linepoints)

fin.close()

# Print line list with points
print 'Lines with points: ', n_line
for i in range(n_line):
    print i+1, ' th lines with points : ', linepoints[i]
print '\n'

# ==================================================
# Make mutilinestring from line list
# ==================================================
multilines = MultiLineString(linepoints)

x = multilines.intersection(multilines)

# Polygonize
result, dangles, cuts, invalids = polygonize_full(x)

result = MultiPolygon(result)
polygon = cascaded_union(result)

# Make mutilinestring from line list
#multilines = MultiLineString(linepoints)

# Polygonize
#result, dangles, cuts, invalids = polygonize_full(multilines)

#result = MultiPolygon(result)
#polygon = cascaded_union(result)

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
# Open file stream
fout = open(filename+'_shapely.geo', 'w')

fout.write('%d\t' % len(points))
fout.write('0\t')
fout.write('%d\n' % len(conns))
for i in range(len(points)):
    fout.write('%d \t' % (i+1))
    fout.write('%s \t %s \t\n' % points[i])

for i in range(len(conns)):
    fout.write('%d \t' % (i+1))
    fout.write('%d \t' % len(conns[i]))
    for j in range(len(conns[i])-1, -1, -1):
        entity = conns[i][j] + 1
        fout.write('%d \t' %entity)
    
    fout.write('\n')
fout.close()