function W = MVC(V, VC);

% number of model vertices
N = size(V,1);
% number of cage vertices
M = size(VC,1);
% initialse the weight matrix 
W = zeros(N,M);

for i = 1:N
    v = V(i,:);
    for j = 1:M
        if (j == 1)
            p1 = VC(end,:);
            p2 = VC(1,:);
            p3 = VC(2,:);
        elseif (j == M)
            p1 = VC(end-1,:);
            p2 = VC(end,:);
            p3 = VC(1,:);
        else 
            p1 = VC(j-1,:);
            p2 = VC(j,:);
            p3 = VC(j+1,:);
        end
        v_p2 = norm(v - p2);
        v_p1 = norm(v - p1);
        v_p3 = norm(v - p3);
        p1_p2 = norm(p1 - p2);
        p2_p3 = norm(p2 - p3);
        cos_alpha = (v_p1^2 + v_p2^2 - p1_p2^2) / (2 * v_p1 * v_p2);
        cos_beta = (v_p2^2 + v_p3^2 - p2_p3^2) / (2 * v_p2 * v_p3);
        alpha_2 = acos(cos_alpha) / 2;
        beta_2 = acos(cos_beta) / 2;
        weight = (tan(alpha_2) + tan(beta_2)) / v_p2;
        W(i,j) = weight;
    end
    % normalise
    W = W ./ sum(W,2);
end