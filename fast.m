% Load model
[V, F] = readOBJ("Tests\cow.obj");

%t = tsurf(F, V);

[P, H] = farthest_points(V, 10);

W = biharmonic_bounded(V, F, H, eye(size(H,1)));

%vis_weights(V, F, W);

% Rotation matrix for 45 degrees about the x-axis
R = [1, 0, 0; 
     0, cosd(45), -sind(45); 
     0, sind(45), cosd(45)];

% Translation vector
t = transpose([10, 20, 30]);

% Constructing the affine transformation matrix
T = [R, t; 
     0, 0, 0, 1];

m = 10;

C = zeros(16, m);

for j = 1:m
    T_j = T;
    c_j = reshape(T_j, 16, 1);

    C(:, j) = c_j;
end

V_prime = FAST(V, F, W, C);

tsurf(F, V_prime);

function vis_weights(V, F, W)
    t = tsurf(F, V, 'CData', W(:, 1));

    for i = 2:size(W, 2)
        t.CData = W(:, i);
    end
end

function V_prime = FAST(V, F, W, C)
   n = size(V, 1);
   m = size(C, 2);
   d = size(V, 2);

   M = zeros(n, (d+1)*m);
   
   for i = 1:n
       for j = 1:m
           M(i, (j-1)*(d+1)+(1:d)) = W(i, j) * V(i, :);
           M(i, j*(d+1)) = W(i, j);
       end
   end

   L = cotangent(V, F);

   r = size(W, 2);

   K = zeros(d*r, n);

   %for i = 1:r
       %C_i = C(i, :, :);

       %K_temp = L * reshape(C_i, [d, d]) .* W(:, i);

       %K((i-1)*(d+1):i*d, :) = K_temp;
   %end

   % Local-Global Optimization

   % Initialising the transformed vertices as the rest-positions
   T = V;

   % Set an epsilon tolerance value
   epsilon = 1e-6;
   delta = epsilon + 1;

   while delta > epsilon
       T_prev = T;

       R = zeros(d, d, m);

       rotations = {};

        % Loop over each bone
        for j = 1:m
            % Compute Qj
            Q_j = zeros(3, 3);
            for i = 1:n
                diff_V = V(i, :)' - reshape(C(1:12, j), 3, 4) * [1; 0; 0; 1];
                W_ij = W(i, j);
                start_idx = (j-1)*d + 1;
                end_idx = j*d;
                T_ik = transpose(T_prev(start_idx:end_idx));
                Q_j = Q_j + W_ij * (diff_V * transpose(T_ik) - transpose(diff_V) * T_ik);
            end
            
            % Perform SVD
            [U, ~, V_svd] = svd(Q_j);
            
            % Ensure we have a proper rotation
            D = eye(3);
            D(3, 3) = sign(det(U * V_svd'));
            
            % Compute the rotation matrix R_j
            R_j = U * D * V_svd';
            
            % Assign R_j to the j-th rotation matrix in R
            R(:, :, j) = R_j;

            rotations{end + 1} = R_j;
        end

       % Global Step
        A = zeros(n*d, m*(d+1));

        for i = 1:n
            for j = 1:m
                idx = (j-1)*(d+1) + (1:d);
                A((i-1)*d + (1:d), idx) = W(i, j) * R(:, :, j); % Using rotation matrix R for bone j
                A((i-1)*d + (1:d), j*(d+1)) = W(i, j) * eye(d);
            end
        end

        Q = A' * (L' * L) * A;

        
        b = zeros((d + 1) * m, 1); % The vector b in the linear system
        
        % Compute the right-hand side b
        for j = 1:m
            sum_j = zeros(d+1, 1);
            for i = 1:n
                e = V(i, :)' - C(:, j);
                e = [e; 1];
                sum_j = sum_j + W(i, j) * R(:, :, j) * e;
            end
            b((j-1)*(d+1) + 1 : j*(d+1)) = sum_j;
        end
        
        % Solve the linear system
        T = Q \ (M' * L * b);

        % Compute the change to check for convergence
        delta = norm(T - T_prev, 'fro');
   end  

   V_prime = M * T;
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