function [F,s_proj,gamma] = Algorithm_antenna_combining(H,C)
% H: channel realization (size: M,N,K)
% C: codebook (size: M, L=2^B, K)
% M: the number of transmit antennas
% N: the number of receive antennas
% K: the number of users
% L: the number of codewords
% F: quantized channel
% s_proj: the projection of quantization vector onto span H
% gamma: weighting vector

[M,N,K] = size(H);
[~,L,~] = size(C);
spanH = [];
for idx0 = 1:1:K
    spanH(:,:,idx0) = orth(H(:,:,idx0));
end

F = zeros(M,1,K);
for idx1 = 1:1:K
    sumk = zeros(1,L);
    for idx2 = 1:1:L
        for idx3 = 1:1:N
            sumk(idx2) = sumk(idx2) + abs(C(:,idx2,idx1)'*spanH(:,idx3,idx1));
        end
    end
    [~,Hindex] = max(sumk);
    F(:,1,idx1) = C(:,Hindex,idx1);
end

s_proj = zeros(M,1,K);
for idx4 = 1:1:K
     sk = spanH(:,:,idx4)*spanH(:,:,idx4)'*F(:,1,idx4);
     s_proj(:,1,idx4) = sk/norm(sk);
end

% s_proj = zeros(M,1,K);
% for idx4 = 1:1:K
%     sk = zeros(M,1);
%     for idx5 = 1:1:N
%         sk = sk + spanH(:,idx5,idx4)*(F(:,1,idx4)'*spanH(:,idx5,idx4));
%     end
%     s_proj(:,1,idx4) = sk/norm(sk);
% end

gamma = zeros(N,1,K);
for idx6 = 1:1:K
    gamma(:,1,idx6) = inv(H(:,:,idx6)'*H(:,:,idx6))*H(:,:,idx6)'*s_proj(:,1,idx6);
    gamma(:,1,idx6) = gamma(:,1,idx6)/norm(gamma(:,1,idx6));
end