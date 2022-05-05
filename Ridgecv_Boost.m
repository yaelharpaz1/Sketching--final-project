function CVMdl = Ridgecv_Boost(A,b,m,k,Lambda)

fprintf('Using function- LMS_Coreset\n')
[C,y] = LMS_Coreset(A,b,m,k);
CVMdl = fitrlinear(C,y,'KFold',m,'Lambda',Lambda,...
'Learner','leastsquares','Regularization','ridge'); 
end