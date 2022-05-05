function test_1()
% for ridge regression solver
%% parameters:
d = [3, 5, 7];
m = 3;
n = [1, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6]*(5*10^5);
[~,sz] = size(n);
T = zeros(6,sz);
%% synthetic data:
for i = 1:3
    d_t = d(1,i);
    k = 2*(d_t+1)^2 +1;
    for idx = 1:sz
        A = randi([0 1000],n(1,idx),d_t);
        b = randi([0 1000],n(1,idx),1);
        
        t1 = tic;
        Lambda = logspace(-5,2,100);
        CVMdl = fitrlinear(A,b,'KFold',m,'Lambda',Lambda,...
        'Learner','leastsquares','Regularization','ridge'); 
        T(1+2*(i-1),idx) = toc(t1);
        
        t2 = tic;
        CVMdl = Ridgecv_Boost(A,b,m,k,Lambda);
        T(2*i,idx) = toc(t2);
    end
end    
figure
plot(n,T(1,:),'b-^',n,T(2,:),'b--',n,T(3,:),'g-^',n,T(4,:),'g--',...
    n,T(5,:),'m-^',n,T(6,:),'m--');
legend({'RidgeCV d=3','Ridge-Boost d=3','RidgeCV d=5','Ridge-Boost d=5',...
    'RidgeCV d=7','Ridge-Boost d=7'});
xlabel('Data size n')
ylabel('computation time (secounds)')
end
    