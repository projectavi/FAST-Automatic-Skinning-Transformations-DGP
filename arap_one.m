%% ARAP

clear;

% Load model
[V_rest, F] = readOBJ("Tests\cow.obj");

n = size(V_rest, 1);
d = size(V_rest, 2);
m = size(F, 1);

% Construct the N(i), the set of vertices connected to vertex i

A = adjacency_matrix(F);

t = tsurf(F, V_rest);
axis equal;
hold on;
xlim([-6, 6]);
ylim([-6, 6]);
frame = getframe(gcf);
v = VideoWriter('old.mp4', 'MPEG-4');
open(v);
writeVideo(v, frame);

V = V_rest;
theta = 0;
for i = 1:41
    b = [42, 100];
    bc = zeros(size(b, 2), d);
    bc(1, :) = [-sin(theta)*1.5, cos(theta)*1.5, 0];
    V = my_arap(V_rest, V_rest, F, b, bc, 100);
    theta = theta + 0.05*pi;
    t.Vertices = V;
    frame = getframe(gcf);
    writeVideo(v, frame);
end

close(v);

% After the deformation, if it was rigid then for each cell i there must be
% a rotation matrix R_i such that
% Ni = N(i, A);
% for j = 1:size(Ni)
%   V_prime(i, :) - V_prime(Ni(j), :) = R_i * (V(i, :) - V(Ni(j), :));
% end
% 
% When the deformation is not rigid, the best approximating rotation matrix
% can still be found by finding R_i that minimizes an energy where w is a
% matrix of weights
% energy = 0;
% for j = 1:size(Ni)
%   energy = energy + w(i, Ni(j)) * square(abs((V_prime(i, :) - V_prime(Ni(j), :)) - R_i*(V(i, :) - V(Ni(j), :))))
% end

% From now on, e(i, j)  = V(i, :) - V(j, :), similarly for e_prime with
% V_prime

% Since terms not including R_i are constant in the minimisation, they are
% not required, thus the the optimisation is reduced to
% energy = 0
% for jj = 1:size(Ni)
%   j = Ni(jj);
%   energy = energy + w(i, j) * e(i, j) * transpose(e_prime(i, j));  
% argmax_Ri (Trace(R_i * energy)) = argmax_Ri (Trace(R_i * S_i))

% Find the covariance matrix which will be represented by S_i
% S_i = P_i * D_i * transpose(Pprime_i)
% D_i is a diagonal matrix containing the weights w(i, j)
% Ni = N(i, A);
% P_i = zeros(3, size(Ni))
% Pprime_i = zeros(3, size(Ni))
% for j = 1:size(Ni)
%   P_i(Ni(j), :) = e(i, Ni(j));
%   Pprime_i(Ni(j), :) = e_prime(i, Ni(j));
% end
% 
% The rotation matrix that maximises the Trace(R_i * S_i) is symmetric
% positive semi-definite. Thus it can be obtained from the SVD of S_i
% [U_i, Sig_i, V_iT] = svd(S_i)
% R_i = transpose(V_iT) * transpose(U_i)

% The total energy measures the total rigidity of the deformation which is
% obtained by summing over the deviations from rigidity per cell.
% E = 0
% for i = 1:n
%   cell_energy = 0;
%   for j = 1:size(Ni)
%       cell_energy = cell_energy + w(i, Ni(j)) * square(abs((V_prime(i, :) - V_prime(Ni(j), :)) - R_i*(V(i, :) - V(Ni(j), :))))
%   end
%   cell_energy = cell_energy * w(i)
%   E = E + cell_energy
% end
% In this case, w(i) is a cell energy, fixed.
% The only variables in this definition are the deformed vertex positions
% because even all the R_i are well-defined functions of the deformed
% vertex positions through SVD.

% To compute weights for edges so that they are as mesh-independent as
% possible, cotangent weights are used.
% w(i, j) = 0.5 * (cot(a_ij) + cot(b_ij))
% a_ij and b_ij are the angles opposite of the mesh edge (i, j). For a
% boundary edge only one such angle exists.
% w(i) = 1 because energy is an integrated quantity and thus cell energy is
% proportional to cell area.
% The cotangent weights are likely in the cotangent Laplacian.

% The idea is to solve for the deformed vertex positions that minimise the
% energy of the deformed shape under user defined modelling constraints.
% Thus, the rigid transformations R_i need to be solved for along with the
% deformed vertex positions. Thus, the energy is viewed as a function of
% these two and we seek the minimum energy under variation of both
% quantities.

% To solve for the next local minimum energy state with a start of initial
% position vectors and rotations we alternate optimisation. Thus, for a
% fixed set of rigid transformations, find deformed vertex positions that
% minimise the energy. Then, fix those, and find the set of rigid
% transformations that minimise that energy. Continue alternating until a
% local minimum energy is reached.

% To find rigid deformations to minimise the energy given a set of deformed
% vertex positions, we seek the deformations which minimise per cell
% energy. This is done using the SVD of the covariance matrix as discussed
% above.

% To compute the optimal deformed vertex positions from the given
% rotations, produce the gradient of the energy wrt the deformed vertex
% positions. Then compute partial derivatives wrt each vertex. 
% The only terms for which the derivative does not vanish are those
% involving neighbouring vertices.
% Thus the partial derivative of the energy wrt to a specific vertex's
% deformed positions is equal to (for a fixed i)
% grad = 0
% Ni = N(i, A)
% for jj = 1:size(Ni)
%   j = Ni(jj)
%   grad = grad + (4 * w(i, j) * ((V_prime(i, :) - V_prime(j, :)) - 0.5 * (R_i + R_j)(V(i, :) - V(j, :)))
% end

% The next step is to set all of the partial derivatives equal to 0 to
% solve for a minimum. This results in a sparse linear system of equations.
% sum over all j in Ni (w(i, j) * (V_prime(i, :) - V_prime(j, :))) = sum
% over all j in Ni (0.5 * w(i, j) * (R_i + R_j) * (V(i, :) - V(j, :)))
% The left side of this system is the discrete Laplace-Beltrami operator
% applied to V_prime. Thus the system can be written as L * V_prime = b
% where L is the aforementioned operator and b is an n-dimensional vector
% whose ith row/entry is the right hand side of the linear system.
% The modelling constraints can be embedded into this system using an
% additional constraint of V_prime(k, :) = c(k) for all k in a set of
% indices of the constrained vertices (H), ie the handles.

% This is incorportated into L * V_prime = b by erasing the corresponding
% rows and columns from L and updating b with the constrained values from
% c.

% The sparse Cholesky factorization with fill-reducing reordering is an
% efficient solver since L is symmetric and positive-definite.

% Overall, first the w(i, j) are precomputed and the system matrix is
% precomputed (L).
% Given an initial guess of V_prime (original V with the constrained
% positions), the local rotations are estimated.
% Then, fixing those local rotations, new positions V_prime are obtained by
% solving L * V_prime = b by doing V_prime = inv(L) * b. This is repeated
% until the energy is at an acceptable level. 

%{
% Scrambled original position
V = V_rest;
%V = rand(n, d);

%t1 = tsurf(F, V);
%hold on;
t = tsurf(F, V);
hold on;
axis equal;
frame = getframe(gcf);
v = VideoWriter('araph4.avi');
open(v);
writeVideo(v, frame);
writeVideo(v, frame);
writeVideo(v, frame);
writeVideo(v, frame);
writeVideo(v, frame);

L = cotmatrix(V_rest, F) + 1e-9*speye(n);
L = -L; % Make the matrix positive definite

% Determine the transformation of the 'handles' for the constraint.

%H = [2612, 49]; % Handle indices
H = [2612, 49];

h = size(H, 2); % number of handles, for this purpose we will use every 290th vertex

C = zeros(h, d);

% Initial Guess
V_prime = V;

% Set the target shape to extract handle values from
target_shape = V_rest;

o = -1;

for i = 1:h
    mod = [0, 0, 0];
    original_pos = target_shape(H(i), :);
    modded = original_pos + mod;
    C(i, :) = modded;
    scatter3(modded(1), modded(2), modded(3), 'filled');
    V_prime(H(i), :) = original_pos;
    o = 1;
end



max_iter = 5;

% Initialise array of local rotation matrices
R = zeros(n, d, d);

offsetx = 5;

%% Computing ARAP Energy to Check
energy = 0;
for i = 1:n
    Ni = N(i, A);
    for j = Ni
        cot_weight = L(i, j);
        
        edge_diff_rest = V_rest(i, :) - V_rest(j, :);

        rotated_diff_rest = squeeze(R(i, :, :)) * edge_diff_rest';

        edge_diff_deformed = V_prime(i, :) - V_prime(j, :);

        %squared_diff = cot_weight * norm(rotated_diff_rest - edge_diff_deformed')^2;
        squared_diff = norm(rotated_diff_rest - edge_diff_deformed')^2;

        energy = energy + squared_diff;
    end
end

fprintf('The ARAP Energy is: %d\n', energy);

% Compute the covariance matrix S
CSM = covariance_scatter_matrix(V_prime, F, 'Energy', 'spokes');

% Optimisation Loops
for iter = 1:max_iter
    
    %% Fixed Deformed Positions, Finding Rotations
    
    S = CSM * repmat(V_prime, d, 1);
    % dim by dim by n list of covariance matrices
    S = permute(reshape(S,[n d d]),[2 3 1]);

    % Iterate through vertices i
    for i = 1:n
        Si = S(:, :, i);

        %Testing new function
        %Si = covariance(V_rest, V_prime, F, L, i);


        % Compute the SVD of Si
        [Ui, Sigi, ViT] = svd(Si);

        Ri = ViT * Ui';
        if det(Ri) <= 0
            Ui(:, 1) = Ui(:, 1) * -1;
            Ri = ViT * Ui';
        end

        if det(Ri - eye(d)) < 1e-12
            Ri = eye(d);
        end

        %[UR, sigR, VR] = svd(Ri);

        %Ri = UR * VR';

        %Ri = Ui * ViT';
        %if det(Ri) <= 0
        %    Ui(:, 1) = Ui(:, 1) * -1;
        %    Ri = Ui * ViT';
        %end
        
        %Testing new function
        %Ri = ViT * diag([ones(1, d-1) det(ViT*Ui')]) * Ui';

        % For testing purposes
        %Ri = eye(d);

        R(i, :, :) = Ri;
    end

    %S = reshape(S, [d d n]);

    %R = fit_rotations(S);

    %% Fixed Rotations, Finding Deformed Positions
    
    % Construction of the b vector
    b = zeros(n, d);

    % b has to be n x d so the construction is a multiplication.
    % R is n x d x d
    % CSM is nd x nd

    for i = 1:n
        RHS = 0;
        Ni = N(i, A);
        for j = Ni
            cot_weight = 2*L(i, j); % Compensating for how minimisation works

            edge_diff_rest = V_rest(i, :) - V_rest(j, :);

            avg_rotated_diff_rest = (squeeze(R(i, :, :)) + squeeze(R(j, :, :)));

            rotated_diff_rest = avg_rotated_diff_rest * edge_diff_rest';

            RHS = RHS + cot_weight * rotated_diff_rest';
        end
        b(i, :, :) = RHS;
    end
    
    %R = reshape(R, [d d n]);
    %[b, K] = arap_rhs(V_rest, F, R);
    
    %Rcol = reshape(permute(R,[3 1 2]),size(K,2),1);
    %Bcol = K * Rcol;
    %b = reshape(Bcol,[n d]);

    %b = -L * V_rest;

    V_prime = min_quad_with_fixed(2*L, b, H, C); % Additional Compensation

    % Compute the covariance matrix S
    CSM = covariance_scatter_matrix(V_prime, F, 'Energy', 'spokes');

    R = reshape(R, [n d d]);

    %% Computing ARAP Energy to Check
    energy = 0;
    for i = 1:n
        Ni = N(i, A);
        for j = Ni
            cot_weight = L(i, j);
            
            edge_diff_rest = V_rest(i, :) - V_rest(j, :);

            rotated_diff_rest = squeeze(R(i, :, :)) * edge_diff_rest';

            if (rotated_diff_rest ~= squeeze(eye(d)) * edge_diff_rest')
                [i, j];
            end

            edge_diff_deformed = V_prime(i, :) - V_prime(j, :);

            %squared_diff = cot_weight * norm(rotated_diff_rest - edge_diff_deformed')^2;
            squared_diff = norm(rotated_diff_rest - edge_diff_deformed')^2;

            energy = energy + squared_diff;
        end
    end

    fprintf('The ARAP Energy is: %d\n', energy);

    offsetx = 10 * iter;

    t.Vertices = V_prime;
    axis equal;
    %drawnow;
    frame = getframe(gcf);
    writeVideo(v, frame);
end

offsetx = 12;

%tsurf(F, V_prime + [offsetx, 0, 0]);

close(v);

% t.Vertices = V_prime;
%}

function V_prime = my_arap(V_in, V_rest, F_in, H, bc, iterations)
    V = V_in;
    F = F_in;

    n = size(V_in, 1);
    d = size(V_in, 2);

    A = adjacency_matrix(F);
    
    h = size(H, 2); % number of handles, for this purpose we will use every 290th vertex

    C = zeros(h, d);
    
    % Initial Guess
    V_prime = V;

    for i = 1:h
        handle_position = V(H(i), :) + bc(i, :);
        C(i, :) = handle_position;
        V_prime(H(i), :) = handle_position;
        scatter3(handle_position(1), handle_position(2), handle_position(3), 'filled');
    end

    max_iter = iterations;
    
    % Initialise array of local rotation matrices
    R = zeros(n, d, d);
    
    L = cotmatrix(V_rest, F) + 1e-9*speye(n);
    L = -L; % Make the matrix positive definite
    
    %% Computing ARAP Energy to Check
    %{
    energy = 0;
    for i = 1:n
        Ni = N(i, A);
        for j = Ni
            cot_weight = L(i, j);
            
            edge_diff_rest = V_rest(i, :) - V_rest(j, :);
    
            rotated_diff_rest = squeeze(R(i, :, :)) * edge_diff_rest';
    
            edge_diff_deformed = V_prime(i, :) - V_prime(j, :);
    
            %squared_diff = cot_weight * norm(rotated_diff_rest - edge_diff_deformed')^2;
            squared_diff = norm(rotated_diff_rest - edge_diff_deformed')^2;
    
            energy = energy + squared_diff;
        end
    end
    
    fprintf('The ARAP Energy is: %d\n', energy);
    %}
    
    % Compute the covariance matrix S
    CSM = covariance_scatter_matrix(V_rest, F, 'Energy', 'spokes');
    
    % Optimisation Loops
    for iter = 1:max_iter
        
        %% Fixed Deformed Positions, Finding Rotations
        
        S = CSM * repmat(V_prime, d, 1);
        % dim by dim by n list of covariance matrices
        S = permute(reshape(S,[n d d]),[2 3 1]);
    
        % Iterate through vertices i
        
        for i = 1:n
            Si = S(:, :, i);
    
            % Compute the SVD of Si
            [Ui, Sigi, ViT] = svd(Si);
    
            Ri = ViT * Ui';
            if det(Ri) <= 0
                Ui(:, 1) = Ui(:, 1) * -1;
                Ri = ViT * Ui';
            end
    
            %if det(Ri - eye(d)) < 1e-12
            %    Ri = eye(d);
            %end
    
            R(i, :, :) = Ri';
        end

        V_prime(H, :) = V_rest(H,:) + bc;

        %{
        SS = zeros(size(CSM,1), d);
        SS(:,1:d) = CSM*repmat(V_prime,d,1);
        SS = permute(reshape(SS,[size(CSM,1)/d d d]),[2 3 1]);

        R = fit_rotations(SS, 'SinglePrecision' ,false);

        R = reshape(R, [n d d]);
        

        V_prime(H, :) = V_rest(H,:) + bc;
        %}
    
        %% Fixed Rotations, Finding Deformed Positions
        
        % Construction of the b vector
        b = zeros(n, d);
    
        for i = 1:n
            RHS = 0;
            Ni = N(i, A);
            for j = Ni
                cot_weight = 2*L(i, j); % Compensating for how minimisation works
                edge_diff_rest = V_rest(i, :) - V_rest(j, :);
                avg_rotated_diff_rest = (squeeze(R(i, :, :)) + squeeze(R(j, :, :)));
                rotated_diff_rest = avg_rotated_diff_rest * edge_diff_rest';
                RHS = RHS + cot_weight * rotated_diff_rest';
            end
            b(i, :, :) = RHS;
        end
    
        V_prime = min_quad_with_fixed(2*L, b, H, C); % Additional Compensation
    
        % Compute the covariance matrix S
        %CSM = covariance_scatter_matrix(V_prime, F, 'Energy', 'spokes');

        %fprintf('The ARAP Energy is: %d\n', arap_energy(V_rest, F, V_prime));
    
        R = reshape(R, [n d d]);
    
        %% Computing ARAP Energy to Check
        %{
        energy = 0;
        for i = 1:n
            Ni = N(i, A);
            for j = Ni           
                edge_diff_rest = V_rest(i, :) - V_rest(j, :);
                rotated_diff_rest = squeeze(R(i, :, :)) * edge_diff_rest';
                edge_diff_deformed = V_prime(i, :) - V_prime(j, :);
                squared_diff = norm(rotated_diff_rest - edge_diff_deformed')^2;
                energy = energy + squared_diff;
            end
        end
    
        fprintf('The ARAP Energy is: %d\n', energy);
        %}
    end
end

% Function to get triangles incident upon vertex i
function Ci = cell(i, F)
    Ci = F(any(F == i, 2), :);
end

% N(i) function for one-ring nearest neighbours
function Ni = N(i, A)
    Ni = find(A(i, :));
end