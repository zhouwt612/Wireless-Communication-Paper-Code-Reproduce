function [L,tau] = Ltaufinder(H,u)
tau = 2*(abs(3+1i*3)+1);
% tau = 0.5;
l_prime = 1000;
Ldata = zeros(1,l_prime);
for idx = 1:1:l_prime
    Hu = (u + tau*idx)'/(H*H')*(u + tau*idx);
    Ldata(1,idx) = Hu;
end
[~,L] = min(Ldata);