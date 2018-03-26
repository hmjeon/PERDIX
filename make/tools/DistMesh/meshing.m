%
% ---------------------------------------------------------------------------------------
%
%                           Wrapper moduel for DistMesh
%
%                                                                    Updated : 2017/03/27
%
% Comments: This module is to generate the triangular meshes using DistMesh
%
% Script written by Hyungmin Jun (hyungminjun@outlook.com)
% Copyright Hyungmin Jun, 2018. All rights reserved.
%
% ---------------------------------------------------------------------------------------
%

function [] = meshing(filename, n_mesh)

% clear all;
close all;
addpath src

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
A    = [abs(max(line(1,:))); abs(max(line(2,:))); abs(min(line(1,:))); abs(min(line(1,:)))];
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

[p,t]=distmesh2d(@dpoly, @huniform, n_mesh, [min(line(:,1)),min(line(:,2)); max(line(:,1)),max(line(:,2))], line, line);

% writing file list
fid = fopen(strcat('input\',name,'_distmesh.geo'), 'w');

fprintf(fid,'%d %d %d\n', size(p,1), 0, size(t,1));
for i=1:size(p,1)
    fprintf(fid,'%d %f %f\n', i, p(i,1), p(i,2));
end

for i=1:size(t,1)
    fprintf(fid,'%d %d %d %d %d\n', i, 3, t(i,1), t(i,2), t(i,3));
end

% Clean up
fclose(fid);

end