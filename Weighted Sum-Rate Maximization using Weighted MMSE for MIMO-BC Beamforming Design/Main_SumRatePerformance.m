clc
clear
close all

K = 4;
N = 2;
M = N*K;

snr = 0:5:30;
sigma2 = 1;
Hnum = 100;

itea_num = 100;

% Allocate the storage of datas
SumRateMRT = zeros(1,length(snr));
SumRateBD = zeros(1,length(snr));
SumRateMMSE = zeros(1,length(snr));
SumRateWMMSE = zeros(1,length(snr));

for idx1 = 1:1:length(snr)
    disp(['SNR: ' num2str(snr(idx1))])
    pow = 10^(snr(idx1)/10);
    % Allocate the storage of datas
    RateMRT = 0;
    RateBD = 0;
    RateMMSE = 0;
    RateWMMSE = 0;
    parfor idx2 = 1:1:Hnum
        H = channel(M,N,K);
        % MRT
        PMRT = TransmitMF(H,pow);
        RateMRT = RateMRT + SumRateMIMO(H,PMRT);
        % ZF
        PZF = BD_MIMOforK(H,pow);
        RateBD = RateBD + SumRateMIMO(H,PZF);
        % MMSE
        PMMSE = MMSE_MIMOforK(H,pow);
        RateMMSE = RateMMSE + SumRateMIMO(H,PMMSE);
        % RWMMSE
        [PWMMSE,~] = WMMSE_MIMOforK(H,PMRT,pow,itea_num);
        RateWMMSE = RateWMMSE + SumRateMIMO(H,PWMMSE);
    end
    SumRateMRT(idx1) = RateMRT/Hnum;
    SumRateBD(idx1) = RateBD/Hnum;
    SumRateMMSE(idx1) = RateMMSE/Hnum;
    SumRateWMMSE(idx1) = RateWMMSE/Hnum;
end

SumRateMRT
SumRateBD
SumRateMMSE
SumRateWMMSE

figure
hold on
plot(snr,SumRateMRT,'m-+','LineWidth',1)
plot(snr,SumRateBD,'k-*','LineWidth',1)
plot(snr,SumRateMMSE,'g-d','LineWidth',1)
plot(snr,SumRateWMMSE,'r-o','LineWidth',1)
xlabel('SNR (dB)')
ylabel('Sum Rate (bps/Hz)')
legend('MRT','BD','MMSE','WMMSE','Location','northwest')
set(gca,'FontSize',11);
set(gca,'LineWidth',1)
box on
grid on
set(gca,'GridLineStyle',':','GridColor','k','GridAlpha',0.5)