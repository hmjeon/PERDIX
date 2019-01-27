%
% =============================================================================
%
% Module - DistMesh
% Last Updated : 05/04/2018, by Hyungmin Jun (hyungminjun@outlook.com)
%
% =============================================================================
%
% This is part of PERDIX-2L, which allows scientists to build and solve
% the sequence design of complex DNAnanostructures.
% Copyright 2018 Hyungmin Jun. All rights reserved.
%
% License - GPL version 3
% PERDIX-2L is free software: you can redistribute it and/or modify it under
% the terms of the GNU General Public License as published by the Free Software
% Foundation, either version 3 of the License, or any later version.
% PERDIX-2L is distributed in the hope that it will be useful, but WITHOUT
% ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
% FOR A PARTICULAR PURPOSE. See the GNU General Public License
% for more details.
% You should have received a copy of the GNU General Public License along with
% this program. If not, see <http://www.gnu.org/licenses/>.
%
% -----------------------------------------------------------------------------
%

function [] = meshing(filename, n_mesh)

% clear all;
close all;

[filepath,name,ext] = fileparts(filename);

fid    = fopen(filename, 'r');
tline  = fgetl(fid);
buffer = sscanf(tline, '%d %d %d%');

n_point = buffer(1);
n_line  = buffer(2);
n_face  = buffer(3);
    
for i = 1: n_point
    tline  = fgetl(fid);
    buffer = sscanf(tline, '%d %f %f');
    
    point(i).id  = buffer(1);
    point(i).x(1)= buffer(2);
    point(i).x(2)= buffer(3);
end

for i = 1: n_face
    tline  = fgetl(fid);
end

face = str2num(tline);

% Set line connectivity
for i = 1: n_point
    line(i,1) = point(face(i+2)).x(1);
    line(i,2) = point(face(i+2)).x(2);
end

line(n_point+1,1) = point(face(3)).x(1);
line(n_point+1,2) = point(face(3)).x(2);

% Scale
A    = [abs(min(line(:,1))); abs(min(line(:,2))); abs(max(line(:,1))); abs(max(line(:,2)))];
A    = max(A);
line = line / A;

% Clean up
fclose(fid);

% [ p, t ] = distmesh_2d ( fd, fh, h, box, iteration_max, fixed );
% where:
% fd       : the name of a distance function defining the region;
% fh       : the name of a mesh density function;
% h        : the nominal mesh spacing;
% box      : defining a box that contains the region;
% iter_max : the maximum number of iterations;
% fixed    : a list of points which must be included in the mesh, or '[]', if no fixed points are given.
% p(output): a list of node coordinates;
% t(output): a list of node indices forming triangles;

% For deploytool
if(isnumeric(n_mesh) == 0)
    n_mesh = str2num(n_mesh)
end

[p,t]=distmesh2d(@dpoly, @huniform, n_mesh, [min(line(:,1)),min(line(:,2)); max(line(:,1)),max(line(:,2))], line, line);

% writing file list
fid = fopen(strcat(name,'_distmesh.geo'), 'w');

fprintf(fid,'%d %d %d\n', size(p,1), 0, size(t,1));
for i=1:size(p,1)
    fprintf(fid,'%d %f %f\n', i, p(i,1), p(i,2));
end

for i=1:size(t,1)
    fprintf(fid,'%d %d %d %d %d\n', i, 3, t(i,1), t(i,2), t(i,3));
end

% Clean up
fclose(fid);
exit;
end