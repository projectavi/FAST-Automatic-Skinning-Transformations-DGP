[V, F] = readOBJ("Tests\cow.obj");
%tsurf(F, V)
n = size(V, 1);
W = ones(n, 2);
t_x = 2;
t_y = 2;
t_z = 2;
T1 = [1, 0, 0, t_x;
     0, 1, 0, t_y;
     0, 0, 1, t_z];
T2 = [1, 0, 0, t_y*5;
     0, 1, 0, t_z*9;
     0, 0, 1, t_x*3];
T = [T1, T2];

V_prime = linear_blend_skinning(V, F, W, T);

%tsurf(F, V_prime)

t = tsurf(F, V_prime);
axis([0 1000 0 1000 0 1000]);
v = VideoWriter('translation2.avi');
open(v);

for i=1:360
    t_x = 2 * i;
    t_y = 2 * i;
    t_z = 2 * i;
    T(:, 4) = [t_x; t_y; t_z];
    V_prime = linear_blend_skinning(V, F, W, T);
    t.Vertices = V_prime;
    drawnow;
    frame = getframe(gcf);
    writeVideo(v, frame);
    %pause(0.0005);
end

close(v);


function V_prime = linear_blend_skinning(V, F, W, T)
    n = size(V, 1);
    d = size(V, 2);

    m = 2; % for now

    M = zeros(n, (d+1)*m);
    
    for i = 1:n
       for j = 1:m
           M(i, (j-1)*(d+1)+(1:d)) = W(i, j) * V(i, :);
           M(i, j*(d+1)) = W(i, j);
       end
    end

    % There are m transformation matrices, one for each transformation for each
    % handle. Each transformation matrix is of size d x (d + 1).
    T = reshape(T, [d, (d+1)*m]);
    T = transpose(T);

    V_prime = M * T;
end