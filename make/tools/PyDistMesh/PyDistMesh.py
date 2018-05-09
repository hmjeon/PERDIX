"""
!
! =============================================================================
!
! Module - PyDistMesh
! Last Updated : 05/04/2018, by Hyungmin Jun (hyungminjun@outlook.com)
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

import numpy as np
import scipy.spatial as spspatial

# Local imports
import distmesh
import distmesh.mlcompat as ml
import distmesh.utils as dmutils

__all__ = ['distmesh2d']

# Open file stream
if len(sys.argv) is 3:
    fin = open(sys.argv[1], 'r')
    filename, filetype = sys.argv[1].split('.')

# sys.argv[0] - Python runnung file name
# sys.argv[1] - First inputs
# sys.argv[2] - Second inputs, mesh spacing parameter
#print " * File name              -> ", filename
#print " * File type              -> ", filetype
#print " * Mesh spacing parameter -> ", sys.argv[2]

# Open file stream
str  = fin.readline()
str  = str.split("\t")
n_point = int(str[0])
n_line  = int(str[1])
n_face  = int(str[2])

#print " * The number of points   -> ", n_point
#print " * The number of lines    -> ", n_line
#print " * The number of faces    -> ", n_face

point_id = []
point_x  = []
point_y  = []
for i in range(n_point):
	str = fin.readline()
	str = str.split()
	point_id.append(int(str[0]))
	point_x.append(float(str[1]))
	point_y.append(float(str[2]))

face = []
str  = fin.readline()
str  = str.split()

for i in range(n_point):
	face.append((point_x[int(str[i+2])-1], point_y[int(str[i+2])-1]))
fin.close()

face.append((point_x[int(str[2])-1], point_y[int(str[2])-1]))

pv   = face
posx = []
posy = []

for i in range(len(pv)):
	posx.append(pv[i][0])
	posy.append(pv[i][1])

minx = min(posx)
miny = min(posy)
maxx = max(posx)
maxy = max(posy)
#print " * Min x and y            -> ", minx, ",  ", miny
#print " * Max x and y            -> ", maxx, ",  ", maxy

max = max([abs(minx), abs(miny), abs(maxx), abs(maxy)])

for i in range(len(pv)):
	pv[i] = (pv[i][0] / max, pv[i][1] / max)

fd = lambda p: distmesh.dpoly(p, pv)
p, t = distmesh.distmesh2d(fd, distmesh.huniform, float(sys.argv[2]), (minx,miny,maxx,maxy), pv)

#for i in range(len(p)):
#	print p[i,0], p[i,1]

# Open file stream
fout = open(filename+'_distmesh.geo', 'w')
fout.write('%d\t' % len(p))
fout.write('0\t')
fout.write('%d\n' % len(t))

for i in range(len(p)):
    fout.write('%5d %14.5f %14.5f \n' % (i, p[i,0], p[i,1]))

for i in range(len(t)):
    fout.write('%5d %5d %5d %5d %5d \n' % (i, 3, t[i,0]+1, t[i,1]+1, t[i,2]+1))

fout.close()
quit()