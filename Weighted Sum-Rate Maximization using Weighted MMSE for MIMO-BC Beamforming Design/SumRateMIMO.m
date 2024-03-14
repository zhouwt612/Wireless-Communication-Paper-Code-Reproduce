function R = SumRateMIMO(H,P)
[M,N,K] = size(H);
INTk = zeros(N,N,K);
for idx1 = 1:1:K
    for idx2 = 1:1:K
        if idx2 ~= idx1
            INTk(:,:,idx1) = INTk(:,:,idx1) + H(:,:,idx1)'*P(:,:,idx2)*P(:,:,idx2)'*H(:,:,idx1);
        end
    end
end
R = 0;
for idx3 = 1:1:K
    R = R + log2(det(eye(N) + P(:,:,idx3)'*H(:,:,idx3)*inv(INTk(:,:,idx3)+eye(N))*H(:,:,idx3)'*P(:,:,idx3)));
end
R = real(R);