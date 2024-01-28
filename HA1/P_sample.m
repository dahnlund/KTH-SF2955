function newinds = P_sample(Z_inds)

newinds = zeros(size(Z_inds));

P = 1/20 * (ones(5,5) + 15*eye(5,5));


for i = 1:length(Z_inds)

    ind = randsample(5, 1, true, P(:, Z_inds(i)));
 
    newinds(i) = ind;


end

end