n = size(x, 1);
Bound = findBoundary(x, f); % in use for static anchors
L = sparse(n,n); % laplacian matrix

if uniform == true 
    % uniform weight
	for i=1:n
        vj = NeighborVex(i, f);
        di = size(vj, 1);
        % L = I - D^-1A
        L(i,i) = 1;
        L(i,vj) = -1/di;
	end
else
    % cot weight
    for i=1:n
        vj = NeighborVex(i, f);
        di = size(vj, 1); % same, finding neighbors

        if di < 3   % unable to calculate cot weight
            L(i,i) = 1;
            continue;
        end

        p = x(vj, :);    % p1...p_di
        p_i = x(i, :);
        angle = getAng([p(end, :); p(1:end - 1, :)], p_i, p);
        angle = 2*pi*angle / sum(angle); % unify
        angle = cumsum(angle);
        angle = angle - angle(1); % set angle(1) as basic line
        
        distance = vecnorm(p - p_i, 2, 2);
        p = [distance.*cos(angle), distance.*sin(angle)];
        p = p(:,1:2);       % in 2D
        p_i = [0, 0];
    
        L(i, vj) = cot(getAng(p_i, [p(end, :); p(1:end - 1, :)], p)) + cot(getAng(p_i, [p(2:end, :); p(1, :)], p));
        L(i,i) = - sum(L(i,:),2);
    end
end

%% utility
% find neighbors
% input vi = coordinate of vi
% input t = all coordinate of texture
% output vj = vi's neighbors
function vj = NeighborVex(vi, texture) 
    % find(texture == vi) to get all pos of vi minus one to get its neighbor
    % mod size(texture, 1) to get correct row_pos
    neighbors = texture( mod( find(texture == vi)-1,size(texture,1) )+1 ,:);  % get all row that represents the neighbor of vi
    vj = unique(neighbors); 
    vj( vj==vi ) = []; % get rid of duplicate vertices and vi itself
end
% output = ang value of p1-p-p2 
function ang = getAng(p1,p,p2)
    d1 = vecnorm(p1 - p, 2, 2);
    d2 = vecnorm(p2 - p, 2, 2);
    d3 = vecnorm(p1 - p2, 2, 2);
    ang = acos((d1.^2 + d2.^2 - d3.^2)./(2.*d1.*d2));
end
