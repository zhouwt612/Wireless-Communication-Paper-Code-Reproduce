function P = BD_MIMOforK(H,pow)
% Block diagonal precoding
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
% H: channel realizations
% pow: transmit power
% Output:
% P: precoding matrices

[M,N,K] = size(H);
P = zeros(M,N,K);
for idx1 = 1:1:K
    Hstack = [];
    for idx2 = 1:1:K
        if idx1 ~= idx2
            Hstack = [Hstack H(:,:,idx2)];
        end
    end
    precoder = null(Hstack');
%     [U S V] = svd(Hstack');
%     precoder = V(:,end-N+1:end);
    P(:,:,idx1) = sqrt(pow/K/trace(precoder*precoder'))*precoder;
end