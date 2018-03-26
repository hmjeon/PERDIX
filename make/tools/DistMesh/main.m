% DistMesh - A simple Mesh Generator
% http://persson.berkeley.edu/distmesh/
% http://people.sc.fsu.edu/~jburkardt/m_src/distmesh/distmesh.html

clear all;
close all;
addpath src

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

%pv=[0.0 0.0; 2.0 0.0; 2.0 1.0; 1.0 1.0; 1.0 2.0; 0.0 2.0; 0.0 0.0];
%[p,t]=distmesh2d(@dpoly, @huniform, 0.4, [-1,-1; 2,2], pv, pv);

%pv=[3.2078 -15.3108; 5.4321 -15.8051; 5.7177 -18.0599; 2.3053 -18.342; 3.2078 -15.3108];
%[p,t]=distmesh2d(@dpoly, @huniform, 0.5, [0,-20; 20,0], pv, pv);

%pv=[0.00 -1.06; 12.38 -0.00; 17.04 -14.66; 30.26 -11.99; 24.59 -30.50; 2.99 -28.67; 0.00 -1.06];
%[p,t]=distmesh2d(@dpoly, @huniform, 7.0, [0,-31; 31,0], pv, pv);

pv=[1.0 0.0; 3.0 0.0; 3.0 1.0; 4.0 1.0; 4.0 3.0; 3.0 3.0; 3.0 4.0; 1.0 4.0; 1.0 3.0; 0.0 3.0; 0.0 1.0; 1.0 1.0; 1.0 0.0];
[p,t]=distmesh2d(@dpoly, @huniform, 0.5*3, [0,0; 3,3], pv, pv);

filename = 'output.geo'

% writing file list
fid = fopen(filename, 'w');


fprintf(fid,'%d %d %d\n', size(p,1), 0, size(t,1));
for i=1:size(p,1)
    fprintf(fid,'%d %f %f\n', i, p(i,1), p(i,2));
end

for i=1:size(t,1)
    fprintf(fid,'%d %d %d %d %d\n', i, 3, t(i,1), t(i,2), t(i,3));
end

% Clean up
fclose(fid);