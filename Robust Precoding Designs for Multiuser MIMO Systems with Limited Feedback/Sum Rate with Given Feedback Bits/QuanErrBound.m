function Dbar = QuanErrBound(M,N,B)
% Calculation for the quantization distortion
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
% M: the number of transmit antenna
% N: the number of receive antenna
% B: the number of feedback bits
% Output
% Dbar: quantization error bound
T = N*(M-N);
Cmn = 1;
for idx = 1:1:N
    Cmn = Cmn*(factorial(M-idx)/factorial(N-idx));
end
Cmn = 1/factorial(T)*Cmn;
Dbar = gamma(1/T)/T*(Cmn)^(-1/T)*2^(-B/T);