%% ARAP

% Load model
[V, F] = readOBJ("Tests\cow.obj");

n = size(V, 1);
d = size(V, 2);
m = size(F, 1);

% Construct the N(i), the set of vertices connected to vertex i

A = adjacency_matrix(F);

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
    
L = cotmatrix(V, F);

% Determine the transformation of the 'handles' for the constraint.
h = 10; % number of handles, for this purpose we will use every 290th vertex

H = []; % Handle indices
for i = 1:h
    H = [H, 290 * i];
end

C = zeros(h, d);

% Initial Guess
V_prime = V;

for i = 1:h
    original_pos = V(H(i), :);
    original_pos(3) = original_pos(3) + 10.5;
    C(i, :) = original_pos;
    V_prime(H(i), :) = original_pos;
end

max_iter = 25;

% Compute the covariance matrix S
CSM = covariance_scatter_matrix(V, F, 'Energy', 'spokes');
S = CSM * repmat(V, d, 1);
% dim by dim by n list of covariance matrices
S = permute(reshape(S,[n d d]),[2 3 1]);

% Initialise array of local rotation matrices
R = zeros(n, d, d);

% Optimisation Loop
for iter = 1:max_iter
    % Iterate through vertices i
    for i = 1:n
        Si = S(:, :, i);

        % Compute the SVD of Si
        [Ui, Sigi, ViT] = svd(Si);

        Ri = transpose(ViT) * transpose(Ui);

        R(i, :, :) = Ri;
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