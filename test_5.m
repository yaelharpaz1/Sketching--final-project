function test_5()
%% parameters:
d = 7;
m = 3;
n = [1, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6]*(5*10^5);
[~,sz] = size(n);
T1 = zeros(2,sz);
T2 = zeros(2,sz);
T3 = zeros(2,sz);
k = 2*(d+1)^2 +1;
%% synthetic data:
l = 100;
Lambda = logspace(-5,2,l);
for idx = 1:sz
    
    A = randi([0 1000],n(1,idx),d);
    b = randi([0 1000],n(1,idx),1);
    
    t1 = tic;
    CVMdl = fitrlinear(A,b,'KFold',m,'Lambda',Lambda,...
    'Learner','leastsquares','Regularization','ridge'); 
    T1(1,idx) = toc(t1);
    t2 = tic;
    CVMdl = Ridgecv_Boost(A,b,m,k,Lambda);
    T1(2,idx) = toc(t2);
    
    t3 = tic;
    x = lasso(A,b,'Lambda',Lambda,'NumLambda',l,'cv',m);
    T2(1,idx) = toc(t3);
    t4 = tic;
    x = LassoCV_Boost(A,b,m,k,Lambda,l);
    T2(2,idx) = toc(t4);
    
    t5 = tic;
    x = lasso(A,b,'Lambda',Lambda,'NumLambda',l,'Alpha',0.5,'cv',m);
    T3(1,idx) = toc(t5);
    t6 = tic;
    x = ElasticCV_Boost(A,b,m,k,Lambda,l);
    T3(2,idx) = toc(t6);   
end
figure
plot(n,T1(1,:),'b-^',n,T1(2,:),'b--',n,T2(1,:),'g-^',n,T2(2,:),'g--',...
    n,T3(1,:),'m-^',n,T3(2,:),'m--');
legend({'RidgeCV','Ridge-Boost','Lasso-CV','Lasso-Boost',...
    'Elastic-CV','Elastic-Boost'})
xlabel('Data size n')
ylabel('computation time (seconds)')
end