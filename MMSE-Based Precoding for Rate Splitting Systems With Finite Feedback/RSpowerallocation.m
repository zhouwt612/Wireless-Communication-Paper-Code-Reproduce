function t = RSpowerallocation(M,K,B,pow)
% t: power allocation factor
% M: the number of transmit antenna
% K: the number of users
% B: the number of feedback bits
% pow: transmit power

delta = K/M;
rho = K/(pow*M^2);
g = 0.5 * (sqrt((1-delta)^2/rho^2 + 2*(1+delta)/rho + 1)+(1-delta)/rho-1);
epsilon = K/((1+g)^2);
B0 = log2(((2*(exp(1)-1)*(M-1))/(pow*M)-(M-1)*epsilon )^(-(M-1))-M/(2*(M-1)));
if B > B0
    t = 1;
else
    t = 1/(pow*M/2*(1/(M-1)*(2^B+M/(2*(M-1)))^(-1/(M-1))+epsilon) + 2 - exp(1) );
end