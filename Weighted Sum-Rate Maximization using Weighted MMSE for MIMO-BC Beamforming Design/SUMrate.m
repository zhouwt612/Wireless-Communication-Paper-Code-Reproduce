function sumrate = SUMrate(H,B,Rn,u,Q,K)
% E: MMSE matrix
% H (QxP): Channel realization
% B (PxQ): Transmitter filter (Beamformer)
% Rn (QxQ): Noise covariance matrix
% u: parameter
% k: Number of users
R = zeros(1,K);
for i = 1:1:K
    Bk = B(:,((i-1)*Q+1):i*Q);
    Hk = H(((i-1)*Q+1:i*Q),:);
    Ek = eye(size(Hk*Hk')) + Bk'*Hk'/Rn(:,:,i)*Hk*Bk;
    Rk = log(det(Ek));
    R(i) = u*Rk;
end
sumrate = real(sum(R));