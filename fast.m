 % Initialise a matrix of rest-pose vertex positions of the input mesh
 % n is the number of vertices and d is the dimension of the mesh, 2 or 3
v = zeros(n, d);

% Set of control handles of size m
H = zeros(m);

% Array of weight functions for linear blend skinning
W = zeros(m); 

% Constructing the Matrix M for calculation of new vertex positions through
% LBS

vp = [v ones(n, 1)];
M = [];
for j = 1:m
    mult = rowfun(@W[j], v);
    M = [M vp.*mult];
end

% Affine transformations for each handle in H
T = zeros(d, d + 1, m);
% There are m transformation matrices, one for each transformation for each
% handle. Each transformation matrix is of size d x (d + 1).
T = reshape(T, [d, (d+1)*m]);
T = transpose(T);

V_prime = M*T;
