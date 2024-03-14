function R = SumRateMIMOforK_RSMA(H,Pc,Pp)
% H: channel realization
% Pc: precoding matrix for common signal
% Pp: precoding matrix for private signal
[M,N,K] = size(H);
interference_c = zeros(N,N,K);
for idx1 = 1:1:K
    for idx2 = 1:1:K
        interference_c(:,:,idx1) = interference_c(:,:,idx1) + ...
            H(:,:,idx1)'*Pp(:,:,idx2)*Pp(:,:,idx2)'*H(:,:,idx1);
    end
end
interference_c = interference_c + eye(N);
Rc = zeros(1,K);
for idx3 = 1:1:K
    Rc(idx3) = log2(det(eye(N) +  Pc'*H(:,:,idx3)/interference_c(:,:,idx3)*H(:,:,idx3)'*Pc));
end
Rc = min(Rc);
interference_p = zeros(N,N,K);
for idx4 = 1:1:K
    for idx5 = 1:1:K
        if idx4 ~= idx5
            interference_p(:,:,idx4) = interference_p(:,:,idx4) + ...
                H(:,:,idx4)'*Pp(:,:,idx5)*Pp(:,:,idx5)'*H(:,:,idx4);
        end
    end
end
interference_p = interference_p + eye(N);
Rp = zeros(1,K);
for idx6 = 1:1:K
    Rp(idx6) = log2(det(eye(N) + Pp(:,:,idx6)'*H(:,:,idx6)/interference_p(:,:,idx6)*H(:,:,idx6)'*Pp(:,:,idx6)));
end
R = real(Rc + sum(Rp));