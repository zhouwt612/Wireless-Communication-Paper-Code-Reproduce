clc
clear
close all

K = 4;
N = 2;
M = N*K;

snr = 0:5:35;
sigma2 = 1;
Hnum = 100;
% Number of samples
SampNum = 1000;

% Estimation error variance
sigmae = 0.05;
nu = 0.01;

% Allocate the storage of datas
SumRateSAA = zeros(1,length(snr));

for idx1 = 1:1:length(snr)
    disp(['SNR: ' num2str(snr(idx1))])
    pow = 10^(snr(idx1)/10);
    % Allocate the storage of datas
    RateSAA = 0;
    parfor idx2 = 1:1:Hnum
        [H,Hhat,E] = channel(M,N,K,sigmae);
        % SVD-MRT
        [Pc_MRT,Pp_MRT] = RS_SVD_MRT_Precoding(Hhat,pow,sigmae);
        % WMMSE-SAA
        warning off
        [Pc_SAA,Pp_SAA,~,~] = RWMMSE_RSMA_Precoding_SAA_nu(Hhat,Pc_MRT,Pp_MRT,pow,sigma2,sigmae,SampNum,nu);
        RateSAA = RateSAA + SumRateMIMOforK_RSMA(H,Pc_SAA,Pp_SAA);
    end
    SumRateSAA(idx1) = RateSAA/Hnum;
end

SumRateSAA
