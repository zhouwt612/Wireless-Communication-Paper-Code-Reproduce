function R = SumRateMIMOforK(H,P)
% H: channel realization (size: (M,N,K))
% P: precoder matrix (size: (M,N,K))
% R: sum rate of MIMO with K users

[M,N,K] = size(H);
interference = zeros(N,N,K);
for idx1 = 1:1:K
    for idx2 = 1:1:K
        if idx2 ~= idx1
            interference(:,:,idx1) = interference(:,:,idx1) + ...
                H(:,:,idx1)'*P(:,:,idx2)*P(:,:,idx2)'*H(:,:,idx1);
        else
            interference(:,:,idx1) = interference(:,:,idx1) + zeros(N);
        end
    end
end
interference = interference + eye(N);
T = zeros(1,K);
for idx3 = 1:1:K
    T(idx3) = log2(det(eye(N) + ...
        P(:,:,idx3)'*H(:,:,idx3)/interference(:,:,idx3)*H(:,:,idx3)'*P(:,:,idx3)));
end
R = sum(real(T));
