function [C,y] = LMS_Coreset(A,b,m,k)
%% Explenation:
% Input: A- size(n,d), n points in R^d.
%        b- size(n,1).
%        m-  an integer, of cross validation folds.
%        k- an integer, 1<k<n, for numerical accuracy/speed trade-off.
% Output: C- size(O(md^2)xd), whose rows are scaled rows from A.
%        y- size(d,1).
[n,d] = size(A);
A_tag = [A, b];
sz = fix(n/m);
a = (d+1)^2+1;
S = zeros(d+1,m*a);
for i = 1:m
    if i==m
        A_temp = A_tag((1+(i-1)*sz):n,:);
    else
        A_temp = A_tag((1+(i-1)*sz):i*sz,:);   %size of (n/m)x(d+1)
    end
    fprintf('Using function- CaratheodoryMatrix\n')
    S_temp = CaratheodoryMatrix(A_temp,k);     %size of (a)x(d+1)
    S(:,(1+(i-1)*a):i*a) = S_temp';
    fprintf('%d/m=%d is done.\n',i,m)
end 
S = S';                                       % size of (m*a)x(d+1)
C = S(:,1:d);
y = S(:,d+1);
end

    