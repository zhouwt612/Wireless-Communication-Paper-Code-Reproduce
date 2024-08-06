%This code can be used to generate the sum rate performance with given
%feedback bits.
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

clc
clear
% close all

K = 4;
N = 2;
M = K*N;

snr = -10:5:30;
sigma2 = 1;
B = [10];
hnum = 1000;

SumRateMRT = zeros(length(B),length(snr));
SumRateZF = zeros(length(B),length(snr));
SumRateMMSE = zeros(length(B),length(snr));
SumRateWMMSE = zeros(length(B),length(snr));
SumRateRMMSE = zeros(length(B),length(snr));
SumRateRWMMSE = zeros(length(B),length(snr));

for idx1 = 1:1:length(B)
    Dbar = QuanErrBound(M,N,B(idx1));
    gamma = Dbar/N;
    for idx2 = 1:1:length(snr)
        disp(['Feedback bits: ' num2str(B(idx1)) ' SNR: ' num2str(snr(idx2))])
        C = RVQ_MIMO_QRforK(M,N,B(idx1),K);
        pow = 10^(snr(idx2)/10);
        RateMRT = zeros(1,hnum);
        RateZF = zeros(1,hnum);
        RateMMSE = zeros(1,hnum);
        RateWMMSE = zeros(1,hnum);
        RateRMMSE = zeros(1,hnum);
        RateRWMMSE = zeros(1,hnum);
        parfor idx3 = 1:1:hnum
            H = channel(M,N,K);
            Htilde = Hbasis(H);
            F = quantizedchannel_MIMO(Htilde,C);
            % MRT
            PMRT = TransmitMF(F,pow);
            RateMRT(idx3) = SumRateMIMOforK(H,PMRT);
            % BD
            PZF = BD_MIMOforK(F,pow);
            RateZF(idx3) = SumRateMIMOforK(H,PZF);
            % MMSE and RMMSE
            PMMSE = MMSE_MIMOforK(F,pow);
            RateMMSE(idx3) = SumRateMIMOforK(H,PMMSE);
            PRMMSE = RMMSE_MIMOforK(F,pow,gamma);
            RateRMMSE(idx3) = SumRateMIMOforK(H,PRMMSE);
            % WMMSE and RWMMSE
            [PWMMSE,~,~] = WMMSE_Precoding(F,H,PMRT,ones(1,K),pow,sigma2,100);
            RateWMMSE(idx3) = SumRateMIMOforK(H,PWMMSE);
            [PRWMMSE,~,~] = RWMMSE_Precoding(F,H,PMRT,gamma,ones(1,K),pow,sigma2,100);
            RateRWMMSE(idx3) = SumRateMIMOforK(H,PRWMMSE);
        end
        SumRateMRT(idx1,idx2) = mean(RateMRT);
        SumRateZF(idx1,idx2) = mean(RateZF);
        SumRateMMSE(idx1,idx2) = mean(RateMMSE);
        SumRateRMMSE(idx1,idx2) = mean(RateRMMSE);
        SumRateWMMSE(idx1,idx2) = mean(RateWMMSE);
        SumRateRWMMSE(idx1,idx2) = mean(RateRWMMSE);
    end
end

SumRateMRT
SumRateZF
SumRateMMSE
SumRateWMMSE
SumRateRMMSE
SumRateRWMMSE


MyFigure_NewColor({SumRateMRT,SumRateZF,SumRateMMSE,SumRateWMMSE,SumRateRMMSE,SumRateRWMMSE},snr,1)