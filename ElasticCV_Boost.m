function x = ElasticCV_Boost(A,b,m,k,Lambda,l)

fprintf('Using function- LMS_Coreset\n')
[C,y] = LMS_Coreset(A,b,m,k);
x = lasso(C,y,'Lambda',Lambda,'NumLambda',l,'Alpha',0.5,'cv',m);
end
