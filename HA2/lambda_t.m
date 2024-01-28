function lambda = lambda_t(l,b,T)

lambda = zeros(T,1);

b_expanded = [1,b,T];


for i = 1:length(b_expanded)-1
    lambda(b_expanded(i):b_expanded(i+1)) = l(i);
end