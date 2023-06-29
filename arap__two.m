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
v = VideoWriter('my_rotation_dynamic.mp4', 'MPEG-4');
open(v);
writeVideo(v, frame);

V = V_rest;
theta = 0;
for i = 1:41
    b = [42, 100];
    bc = zeros(size(b, 2), d);
    bc(1, :) = [-sin(theta)*1.5, cos(theta)*1.5, 0];
    C = V_rest(b, :) + bc;
    %V = arap(V_rest, F, b, C, 'Energy', 'spokes');
    V = arap_two(V_rest, F, b, C, 100);
    theta = theta + 0.05*pi;
    t.Vertices = V;
    frame = getframe(gcf);
    writeVideo(v, frame);
end

close(v);

function V_prime = arap_two(V, F, H, C, max_iterations)
    % Get number of vertices
    n = size(V, 1);
    d = size(V, 2);

    V_prime = [];

    energy = [];

    Aeq = [];
    Beq = [];

    if isempty(C)
        C = sparse(0, d);
    end

    if isempty(V_prime)
        V_prime = V;
        V_prime = V_prime(:, 1:d);
    end

    % This is for dynamic
    M = massmatrix(V, F);

    L = cotmatrix(V, F);

    CSM = covariance_scatter_matrix(V, F, 'Energy', 'spokes');

    m = size(F, 1);
    mesh_size = size(F, 2);

    [KX, KY, KZ] = covariance_block(V, F, L, n);
    Z = sparse(size(V, 1), n);

    K = [KX Z Z KY Z Z KZ Z Z ;
         Z KX Z Z KY Z Z KZ Z];
    K = [K; Z Z KX Z Z KY Z Z KZ];

    % Precompute the multiplier
    %[~, K] = arap_rhs(V, F, [], 'Energy', 'spokes');

    V_prime_prev = V;
    energy = inf;
    
    for iteration = 1:max_iterations
        change = max(abs(V_prime(:) - V_prime_prev(:)));

        V_prime_prev = V_prime;
    
        V_prime(H, :) = C;
    
        % Compute the covariance matrices
        S = zeros(n*d, d);
        S(:, 1:d) = CSM * repmat(V_prime, d, 1);
    
        S = permute(reshape(S, [n d d]), [2 3 1]);
        
        R = zeros([d d n]);
        SS = zeros([d d n]);

        for i = 1:n
            [US, SigS, VS] = svd(S(:, :, i));

            Ri = VS * US';

            SS(:, :, i) = SigS;

            if(det(Ri) < 0) 
                US(:, end) = -US(:, end);
                Ri = VS * US';
            end

            R(:, :, i) = Ri;
        end
    
        V_prime(H, :) = C;
    
        preF = [];
    
        Rcol = reshape(permute(R, [3 1 2]), size(K, 2), 1);
        Bcol = K * Rcol;
        B = reshape(Bcol, [size(Bcol, 1)/d d]);
    
        LHS = -0.5*L + 0.5*M;
        RHS = -1*B + M*(-V);

        %zQ = -0.5*dyn_alpha*data.L+DQ+alpha_tik*speye(size(data.L));
        %        zL = -dyn_alpha*B+Dl;
    
        [V_prime, preF] = min_quad_with_fixed(LHS, RHS, H, C, Aeq, Beq, preF);
    
        energy_prev = energy;
        energy = trace(V_prime' * (LHS) * V_prime + V_prime' * (RHS)) + trace(V' * (0.5 * L * V));
        if energy > energy_prev
            fprintf('arap: energy (%g) increasing (over %g) (iter %d)\n', ...
              energy,energy_prev,iteration);
        end
    end
    V_prime = V_prime(:, 1:d);
end

function [KX, KY, KZ] = covariance_block(V, F, L, n)
    E = edges(F);

    % Build sparse adjacency matrix with position differences in each
    % dimension instead of binary existentials

    AX = sparse(E(:, 1), E(:,2), V(E(:, 1),1) - V(E(:, 2), 1), n, n);
    AY = sparse(E(:, 1), E(:,2), V(E(:, 1),2) - V(E(:, 2), 2), n, n);
    AZ = sparse(E(:, 1), E(:,2), V(E(:, 1),3) - V(E(:, 2), 3), n, n);

    AX = AX-AX';
    AY = AY-AY';
    AZ = AZ-AZ';

    KX = L.*AX;
    KX = 0.5*(KX + diag(sum(KX, 2)));

    KY = L.*AY;
    KY = 0.5*(KY + diag(sum(KY, 2)));

    KZ = L.*AZ;
    KZ = 0.5*(KZ + diag(sum(KZ, 2)));
end





    

