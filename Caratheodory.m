function [S,w] = Caratheodory(P,u)
%% Explenation:
% Input: P- size(n,d), n points in R^d.
%        u- size(n,1), a weight function such that the sum, over all p in 
%        P, of u(p) is 1.
% Output: (S,w)- a Caratheodory set for (P,u), such that size(S)<=d+1, the
%        weighted mean is the same, and the sum of weights is 1.
%%
[n,d] = size(P);
if n <= d+1
    %fprintf('|P| is already small.\n')
    S = P;
    w = u;
else 
    A = zeros(d,n-1);
    for idx_p = 2:n
        A(:,idx_p-1) = P(idx_p,:)-P(1,:);
    end
    v_temp = null(A);                   
    v = v_temp(:,1);                       % size(n-1,1)
    if isempty(v)==1
        error('Error- the null space of A is empty.\n')
    end
    v1 = -sum(v);
    v_vec(1,1) = v1;
    v_vec(2:n,1) = v;
    idx_nnz1 = find(v_vec>0);
    [nz1,~] = size(idx_nnz1);
    alpha_vec = zeros(nz1,1);
    for i = 1:nz1
        idx = idx_nnz1(i,1);
        alpha_vec(i,1) = u(idx,1)/v_vec(idx,1);
    end
    alpha = min(alpha_vec);
    w = u - alpha*v_vec;
    idx_nnz2 = find(w>0);
    [nz2,~] = size(idx_nnz2);
    S = zeros(nz2,d);
    for i = 1:nz2
        idx = idx_nnz2(i,1);
        S(i,:) = P(idx,:);
    end
    w(w<=0) = [];
    if nz2 > d+1
        [S,w] = Caratheodory(S,w);   % recursive call.
    end
end

  