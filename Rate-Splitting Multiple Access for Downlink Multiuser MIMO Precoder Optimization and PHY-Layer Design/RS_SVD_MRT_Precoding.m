function [Pc,Pp] = RS_SVD_MRT_Precoding(H,pow,sigmae)
% Pc, Pp: precoding matrix for the common and private data streams
% MMSE precoding based
% H: the quantized channel
% pow: transmit power
% Power allocation:
% Common: pow*(1-t)
% Private: pow*t/K
Pt_alpha = 1/sigmae;
t = min(Pt_alpha/pow,1);
[M,N,K] = size(H);
Hdedicated = [];
for idx2 = 1:1:K
    Hdedicated = [Hdedicated H(:,:,idx2)];
end
Pc_SVD = orth(Hdedicated);
Pc = Pc_SVD(:,1:N);
Pc = sqrt((pow*(1-t))/trace(Pc'*Pc))*Pc;
for idx3 = 1:1:K
    Pp(:,:,idx3) = sqrt(pow*t/K/trace(H(:,:,idx3)*H(:,:,idx3)'))*H(:,:,idx3);
end