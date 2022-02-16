function R = SumRateMIMOforK2(H,P)
% H: channel realization (size: (M,N,K))
% P: precoder matrix (size: (M,N,K))
% R: sum rate of MIMO with K users
[M,N,K] = size(H);
sigpow = zeros(N);
interferencepow = zeros(N);
for idx1 = 1:1:K
    for idx2 = 1:1:K
        if idx1 == idx2
            sigpow = sigpow + H(:,:,idx1)'*P(:,:,idx1)*P(:,:,idx1)'*H(:,:,idx1);
        else
            interferencepow = interferencepow + H(:,:,idx1)'*P(:,:,idx2)*P(:,:,idx2)'*H(:,:,idx1);
        end
    end
end
sigpow = sigpow + eye(N);
interferencepow = interferencepow + eye(N);
% R = real(log2(det(sigpow/interferencepow)));
R = real(log2(det(sigpow)/det(interferencepow)));