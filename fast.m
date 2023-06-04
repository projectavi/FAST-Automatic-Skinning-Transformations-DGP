function FAST = fast(V, F, W, C, T)
   n = size(V, 1);
   m = size(C, 1);
   d = size(V, 2);

   M = zeros(n, (d+1)*m);
   
   for i = 1:n
       for j = 1:m
           M(i, (j-1)*(d+1)+(1:d)) = W(i, j) * V(i, :);
           M(i, j*(d+1)) = W(i, j);
       end
   end

   T = reshape(T, [d, (d+1)*m]);
   T = transpose(T);

   L = cotangent(V, F);

   r = size(W, 2);

   K = calculateK(d, r, n, C, L, W);

   % Local-Global Optimization

   % Initialising the transformed vertices as the rest-positions
   T = V;

   % Set an epsilon tolerance value
   epsilon = 1e-6;
   delta = epsilon + 1;

   while delta > epsilon
       T_prev = T;

       R = zeros(r, d, d);
        
       % Local Step
       for i = 1:r
           % Weighted sum under the transformed vertices
           Q_i = transpose(T) * diag(W(:, i));

           % SVD
           [U_svd, sigma_svd, V_svd] = svd(Q_i);

           R(i,:,:) = U_svd * transpose(V_svd);
       end

       % Global Step
       
       
   end

function K = calculateK(d, r, n, C, L, W)
K = zeroes(d*r, n);

for i = 1:r
    C_i = C(i, :, :);
    
    K_temp = L * reshape(C_i, d, d) .* W(:, i);
    
    K((i-1)*(d+1):i*d, :) = K_temp;
end

%% TEMP

%% Initialise a matrix of rest-pose vertex positions of the input mesh
% % n is the number of vertices and d is the dimension of the mesh, 2 or 3
%v = zeros(n, d);

%% Set of control handles of size m
%H = zeros(m);

%% Array of weight functions for linear blend skinning
%W = zeros(m); 

%% Constructing the Matrix M for calculation of new vertex positions through
%% LBS

%vp = [v ones(n, 1)];
%M = [];
%for j = 1:m
%    mult = rowfun(@W[j], v);
%    M = [M vp.*mult];
%end

%% Affine transformations for each handle in H
%T = zeros(d, d + 1, m);
%% There are m transformation matrices, one for each transformation for each
%% handle. Each transformation matrix is of size d x (d + 1).
%T = reshape(T, [d, (d+1)*m]);
%T = transpose(T);

%V_prime = M*T;