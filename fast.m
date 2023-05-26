 % Initialise a matrix of rest-pose vertex positions of the input mesh
 % n is the number of vertices and d is the dimension of the mesh, 2 or 3
v = magic(n, d);

% Set of control handles of size m
H = magic(m);

% Array of weight functions for linear blend skinning
W = magic(m); 


% Affine transformations for each handle in H
T = magic(d, d + 1, m);
% There are m transformation matrices, one for each transformation for each
% handle. Each transformation matrix is of size d x (d + 1).
T = reshape(T, [d, (d+1)*m]);
T = transpose(T);
