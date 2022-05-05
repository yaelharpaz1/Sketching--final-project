function [C,w] = FastCaratheodorySet(P,u,k)
%% Explenation:
% Input: P- size(n,d), n points in R^d.
%        u- size(n,1), a weight function such that the sum, over all p in
%        P, of u(p) is 1.
%        k- an integer, 1<k<n, this is the number of clusters.
% Output: (C,w)- a Caratheodory set for (P,u), such that |C|<=d+1, the
%        weighted mean is the same, and the sum of weights is 1.
%% Remove all points with zero weight:
idx_p = find(u==0);
if isempty(idx_p)~=1 
    P(idx_p,:) = [];
end
[n,d] = size(P);
%%
if n <= d+1
    fprintf('|P| is already small.\n')
    C = P;
    w = u;
else
    %% Step 1:compute a balanced partition of P into k disjoint subsets:
    if k > n
        k = d+2;
    end
    %fprintf(['Step 1:computing a balanced partition of P into '...
    %    'k=%d disjoint subsets\n'],k)
    [~,idx] = sort(rand(n,1));
    rem = mod(n,k);
    sz = (n-rem)/k;
    cl = zeros(n,1);
    for i = 1:k
        cl((1+(i-1)*sz):i*sz,1) = i;
    end
    for i = 1:rem
        cl(n-rem+i,1) = i;
    end
    cluster_vec = cl(idx,1);  
    %% Step 2: compute a sketch for each cluster:
    %fprintf('Step 2:computing a sketch for each cluster\n')
    mue = zeros(k,d);
    u_tag = zeros(k,1);
    for idx_k = 1:k
        idx_vec = find(cluster_vec==idx_k);
        [cluster_sz,~] = size(idx_vec);
        s1 = 0;
        s2 = zeros(1,d);
        for i = 1:cluster_sz
            idx = idx_vec(i,1);
            s1 = s1 + u(idx,1);             
            s2 = s2 + u(idx,1)*P(idx,:);   % size(1,d)
        end
        mue(idx_k,:) = (1/s1)*s2;         % The weighted mean of Pi
        u_tag(idx_k,1) = s1;              % The weight of the i'th cluster
    end
    %fprintf('Using function- Caratheodory\n')
    [mue_tild, w_tild] = Caratheodory(mue,u_tag);  
    %% Step 3:compute the coreset for the union of clusters correspond to 
    %  the selected sketches:
    %fprintf(['Step 3:computing the coreset for the union of clusters '...
    %'correspond to the selected sketches\n'])
    [sz_Cset,~] = size(mue_tild);
    C = zeros((sz+1)*sz_Cset,d);
    w = zeros((sz+1)*sz_Cset,1);
    for i = 1:sz_Cset
        [~,idx_Cset,~] = intersect(mue,mue_tild(i,:),'rows');
        idx_vec = find(cluster_vec==idx_Cset);
        [cluster_sz,~] = size(idx_vec);
        s1 = zeros(sz_Cset,1);
        for j = 1:cluster_sz
            idx = idx_vec(j,1);
            s1(i,1) = s1(i,1) + u(idx,:);
        end
        C((1+(i-1)*cluster_sz):i*cluster_sz,:) = P(idx_vec,:);
        % Assign weight for each point in C:
        w((1+(i-1)*cluster_sz):i*cluster_sz,1) = (w_tild(i,1)/s1(i,1))*u(idx_vec,1);
    end
    C( ~any(C,2),:) = [];
    w(w==0) = []; 
    [n_C,~] = size(C);
    %% Step 4: recursively compute a coreset C until a sufficiently small 
    %  coreset is obtaines:
    %fprintf(['Step 4:recursively computing a coreset C until a '...
    %'sufficiently small coreset is obtaines\n'])
    fprintf(['|C|=%d, applying recursive call in Fast-Caratheodory-Set '... 
    'function.\n'],n_C)
    [C,w] = FastCaratheodorySet(C,w,k);      % recursive call.
end
end

          
    