m = n - size(P2PVtxIds,1); % n - m is those flexible vertices

[x0,y0,z0,w0] = deal(rot(1),rot(2),rot(3),rot(4));
R = [1-2*y0*y0-2*z0*z0  2*x0*y0 - 2*w0*z0   2*x0*z0 + 2*w0*y0
    2*x0*y0 + 2*w0*z0   1-2*x0*x0-2*z0*z0   2*y0*z0 - 2*w0*x0
    2*x0*z0 - 2*w0*y0   2*y0*z0 + 2*w0*x0   1-2*x0*x0-2*y0*y0]; % rotate matrix
[X, Y, Weight] = find(L);
X = [X X+n X+2*n];
Y = [Y Y+n Y+2*n];
L_mat = sparse(X, Y, repmat(Weight, 3, 1));
delta = full(L)*x; %x is actual coordinate of mesh 
delta(P2PVtxIds,:) = delta(P2PVtxIds,:)*R'; % rotation is here to apply

for i=1:n
    % get neighbors of i
    Neighbors=[NeighborVex(i, f)', i];
    n_i = length(Neighbors);
    V_i = x(Neighbors, :)';
    
    %the matrix to solve T
    A = zeros(3 * n_i, 7);
    for r=1:length(Neighbors)
        A(r,:)       =  [V_i(1,r) 0 V_i(3,r) (-1)*V_i(2,r) 1 0 0];
        A(r+n_i,:)   =  [V_i(2,r) (-1)*V_i(3,r) 0 V_i(1,r) 0 1 0];
        A(r+2*n_i,:) =  [V_i(3,r) V_i(2,r) (-1)*V_i(1,r) 0 0 0 1];
    end
    
    delta_i = delta(i,:)';
    delta_ix = delta_i(1);
    delta_iy = delta_i(2);
    delta_iz = delta_i(3);
    
    Ainv = pinv(A);
    s = Ainv(1,:);
    h1 = Ainv(2,:);
    h2 = Ainv(3,:);
    h3 = Ainv(4,:);

    T_delta = [s * delta_ix         - h3 * delta_iy + h2 * delta_iz
               h3 * delta_ix        + s * delta_iy  - h1 * delta_iz
               (-1) * h2 * delta_ix + h1 * delta_iy + s * delta_iz];
    
    % update weights
    r_xyz = [Neighbors Neighbors+n Neighbors+2*n];
    L_mat(i,r_xyz) =      (-1)*L_mat(i,r_xyz) + T_delta(1,:);
    L_mat(i+n,r_xyz) =    (-1)*L_mat(i+n,r_xyz) + T_delta(2,:);
    L_mat(i+2*n,r_xyz) =  (-1)*L_mat(i+2*n,r_xyz) + T_delta(3,:);
end
L_mat=L_mat'*L_mat;

y = LaplacianEditing(L_mat, x, Bound, P2PVtxIds', p2pDsts, 1); % get answers


%% Laplacian editing
% the implementation is based on the demo given by paper writer
% input L the laplacian matrix
% x the coordinate
% static_index the vertices needs to be fixed, AKA boundaries
% handle_index the vertices needs to be moved
% new_pos the new coordinate of moved vertices
% w the required weight
% output V the new position of whole mesh
function V = LaplacianEditing(L, x, static_index, handle_index, new_pos, w)

n = length(x); % vertice num
x(handle_index,:) = new_pos; % apply new pos
V = reshape(x, 3*n, 1);

constrain_points = [static_index handle_index]; % constrains need to be fufilled
vec = sparse(1, 3*n);
vec([constrain_points, constrain_points+n, constrain_points+2*n]) = w;
D = diag(vec);
B = D' * D;
rhs = B * V;
A = L + B;

V_new_pos = A\rhs;
V = [V_new_pos(1:n) V_new_pos((n+1):(2*n)) V_new_pos((2*n+1):(3*n))];
end
%% find neighbor
function vj = NeighborVex(vi, texture) 
    % find(texture == vi) to get all pos of vi minus one to get its neighbor
    % mod size(texture, 1) to get correct row_pos
    neighbors = texture( mod( find(texture == vi)-1,size(texture,1) )+1 ,:);  % get all row that represents the neighbor of vi
    vj = unique(neighbors); 
    vj( vj==vi ) = []; % get rid of duplicate vertices and vi itself
end