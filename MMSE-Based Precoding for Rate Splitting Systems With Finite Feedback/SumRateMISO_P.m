function Rp = SumRateMISO_P(H,Pp,pow,t)
% Rp: sum rate of the private signal
% H: [h1,...,hK]^H channel realization
% Pp: the precoding matrix for private part
% pow: transmit power
% t: fraction of the transmitted power

[K,~] = size(H);
% Calculate the private sum rate
Rp = zeros(1,K);
for idx0 = 1:1:K
    Powk = abs(H(idx0,:)*Pp(:,idx0))^2;
    Powint = 0;
    for idx1 = 1:1:K
        if idx1 ~= idx0
            Powint = Powint + abs(H(idx0,:)*Pp(:,idx1))^2;
        end
    end
    gamma_k = (pow*t/K*Powk)/(1 + pow*t/K*Powint);
    Rp(idx0) = log2(1 + gamma_k);
end

