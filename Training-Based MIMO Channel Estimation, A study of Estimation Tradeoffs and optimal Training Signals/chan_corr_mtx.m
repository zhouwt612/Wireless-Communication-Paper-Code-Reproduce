function Rh = chan_corr_mtx(r,t,epsilon)
% Compute the channel correlation matrix
% [Rh]_n,m = r*epsilon^(abs(n-m))
% n and m: the indexes of the array sensors
% epsilon: the intensity of the correlated channels (0<=epsilon<1)

Rh = zeros(t,t);
for idx1 = 1:1:t
    for idx2 = 1:1:t
        Rh(idx1,idx2) = r*epsilon^(abs(idx1-idx2));
    end
end