function [Heff,F] = ASalgorithm(H,C)
% AS: antenna selection
% H: channel realization (size: M,N,K)
% C: codebook (size: M,2^B,K)
% M: the number of transmit antennas
% N: the number of receive antennas
% K: the number of UEs
% Heff: effective channel after the antenna selection process
% F: quantized channel (codeword)

[M,N,K] = size(H);
[~,L,~] = size(C);

F = zeros(M,1,K);
Heff = zeros(M,1,K);
for idx1 = 1:1:K
    hwmulti = zeros(N,L);
    for idx2 = 1:1:N
        for idx3 = 1:1:L
            hwmulti(idx2,idx3) = abs(H(:,idx2,idx1)'*C(:,idx3,idx1));
        end
    end
    val = zeros(1,L);
    index = zeros(1,L);
    for idx4 = 1:1:N
        [val(idx4),index(idx4)] = max(hwmulti(idx4,:));
    end
    [~,maxindex] = max(val);
    F(:,1,idx1) = C(:,index(maxindex),idx1);
    [heffindex,~] = find(hwmulti==max(max(hwmulti(:,:))));
    Heff(:,1,idx1) = H(:,heffindex,idx1);
end