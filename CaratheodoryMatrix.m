function S = CaratheodoryMatrix(A,k)
%% Explenation:
% Input: A- size(n,d), n points in R^d.
%        k- an integer, 1<k<n, for numerical accuracy/speed trade-off.
% Output: S- size(d^2+1,d), S rows are scaled rows from A, and A^TA=S^TS.
%% 
[n,d] = size(A);
P = zeros(n,d^2);
for i = 1:n
    temp = (A(i,:)).'*A(i,:);
    P(i,:) = (temp(:)).';
end
u = (1/n)*ones(n,1);
fprintf('Using function- FastCaratheodorySet\n')
[~,w] = FastCaratheodorySet(P,u,k);    %|c|=(d^2+1)xd, |w|=(d^2+1)x1.
S = zeros(d^2+1,d);
for i = 1:d^2+1
    S(i,:) = sqrt(n*w(i,1)).*A(i,:);
end
end
    