function precoder = RegZFprecoder(F,P)
% G: regularized zero-forcing precoder
% H: channel
% M: the number of transmit antennas
% P: the transmit power
[M,~,K] = size(F);
Fcom = zeros(M,K);
for idx = 1:1:K
    Fcom(:,idx) = F(:,1,idx);
end
H = Fcom';
P = H'/(H*H' + (M/P)*eye(K));
precoder = zeros(size(P));
for idx1 = 1:1:K
   precoder(:,idx1) = P(:,idx1)/norm(P(:,idx1));
end