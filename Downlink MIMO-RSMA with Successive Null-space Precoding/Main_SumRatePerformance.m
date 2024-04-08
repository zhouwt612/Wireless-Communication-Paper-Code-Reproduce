clc
clear
close all

% M: # of tx antennas
% N: # of rx antennas
% K: # of users
K = 4;
N = 2;
M = N*K;

snr = 0:5:30;
sigma2 = 1;
Hnum = 1;

% Estimation error variance
sigmae = 0.05;
nu = 0.01;

% Allocate the storage of datas
SumRateSNS = zeros(1,length(snr));

for idx1 = 1:1:length(snr)
    disp(['SNR: ' num2str(snr(idx1))])
    pow = 10^(snr(idx1)/10);
    % Allocate the storage of datas
    RateSNS = 0;
    parfor idx2 = 1:1:Hnum
        [H,Hhat,E] = channel(M,N,K,sigmae);
        [H,Hhat] = Channel_Sort_Imp(H,Hhat,pow,sigma2);
        % SNS-Precoding
        [Pc_SNS,Pp_SNS] = SNS_Precoding_nu(Hhat,pow,sigma2,nu);
        RateSNS = RateSNS + SumRateMIMOforK_RSMA(H,Pc_SNS,Pp_SNS);
    end
    SumRateSNS(idx1) = RateSNS/Hnum;
end

SumRateSNS