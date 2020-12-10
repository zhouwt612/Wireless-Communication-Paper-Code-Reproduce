function [L,tau] = LtaufinderM(H,u,alpha)
tau = 2*(abs(3+1i*3)+1);
% tau = 0.5;
l_prime = 1000;
Ldata = zeros(1,l_prime);
for idx = 1:1:l_prime
    Hu = norm(H'/(H*H'+alpha*eye(size(H*H')))*(u+tau*idx))^2;
    Ldata(1,idx) = Hu;
end
[~,L] = min(Ldata);