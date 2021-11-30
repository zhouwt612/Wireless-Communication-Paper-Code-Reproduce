function sumrate = SUMrateMISO(H,P,Rn,u,K)
R = zeros(1,K);
for idx = 1:1:K
    Pk = P(:,idx);
    Hk = H(idx,:);
    Ek = 1 + Pk'*Hk'/Rn(:,:,idx)*Hk*Pk;
    Rk = log2(det(Ek));
    R(idx) = u(idx)*Rk;
end
sumrate = real(sum(R));