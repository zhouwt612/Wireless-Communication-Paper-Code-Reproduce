function R = SumRateMIMOforK(H,P)
% Calculation for the sum rate
%This code refers to the following scientific article:
%
% Wentao Zhou, Di Zhang, Mérouane Debbah, and Inkyu Lee,
% "Robust Precoding Designs for Multiuser MIMO Systems with Limited Feedback,
%" IEEE Transactions on Wireless Communications, To appear.
% 
% This is version 1.0 (last edited: 2024-04-08)
% 
%License: This code is licensed under the GPLv2 license. If you in any way
%use this code for research that results in publications, please cite our
%original article listed above.
% Input:
% H: channel realization (size: (M,N,K))
% P: precoder matrix (size: (M,N,K))
% Output
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