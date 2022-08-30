function SR = SRofGPIP(H,F,w,sigma2,P)
% H: Channel (N x K)
% F: precoder (N x K)
[K,N] = size(H);
SumRate = zeros(1,K);
interference = zeros(1,K);
SINR = zeros(1,K);
for idx1 = 1:1:K
    Pk = abs(H(idx1,:)*F(:,idx1))^2;
%     Pk = F(:,idx1)'*H(:,idx1)*H(:,idx1)'*F(:,idx1);
    for idx2 = 1:1:K
        if idx2 ~= idx1
            interference(1,idx1) = interference(1,idx1) + abs(H(idx1,:)*F(:,idx2))^2;
%             interference(1,idx1) = interference(1,idx1) + H(:,idx1)'*F(:,idx2)*F(:,idx2)'*H(:,idx1);
        end
    end
    SINR(1,idx1) = Pk/(interference(idx1)+sigma2/P);
end

for idx3 = 1:1:K
    SumRate(idx3) = w(idx3)*log2(det(1+SINR(idx3)));
end
% SR = real(sum(SumRate));
SR = sum(SumRate);