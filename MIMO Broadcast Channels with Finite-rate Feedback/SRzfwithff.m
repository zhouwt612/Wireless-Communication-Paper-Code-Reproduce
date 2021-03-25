function Rp = SRzfwithff(P,B,M)

Rp = log2(1+P*2^(-B/(M-1)));
% Rp = log2(1+P*norm(h)^2*(1-2^(-B/(M-1))));
