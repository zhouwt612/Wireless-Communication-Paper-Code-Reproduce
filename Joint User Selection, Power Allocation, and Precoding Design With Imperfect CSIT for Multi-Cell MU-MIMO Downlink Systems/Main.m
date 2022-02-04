clc
clear
close all

hnum = 10000;
K = 4;
N = 4;
w = ones(1,K);
Error_Pow = 0.1;
Phi = EstiErrCovMtx(N,K,Error_Pow*ones(K));
sigma2 = 1;
epsilon = 0.01;
snr = 0:5:30;

SR_GPIP = zeros(length(snr),hnum);
SR_GPIPnoCov = zeros(length(snr),hnum);
SR_RZF = zeros(length(snr),hnum);
SR_MRT = zeros(length(snr),hnum);
SR_RRZF = zeros(length(snr),hnum);
SR_WMMSE = zeros(length(snr),hnum);

for idx1 = 1:1:length(snr)
    disp(['SNR: ' num2str(snr(idx1))])
    P = 10^(snr(idx1)/10);
    for idx2 = 1:1:hnum
        H = sqrt(1/2)*(randn(K,N)+1i*randn(K,N));
        E = sqrt(Error_Pow/2)*(randn(K,N)+1i*randn(K,N));
        EH = H - E;
        
        MRTpre = MRTprecoder(EH,P);
        Rmrt = nointercov(H,MRTpre,N,1,K);
        SR_MRT(idx1,idx2) = SUMrateMISO(H,MRTpre,Rmrt,w,K);
        
        RZFpre = RegZFprecoder(EH,P);
        Rrzf = nointercov(H,RZFpre,N,1,K);
        SR_RZF(idx1,idx2) = SUMrateMISO(H,RZFpre,Rrzf,w,K);
        
        RRZFpre = RRZF(EH,Phi,sigma2,P);
        Rrrzf = nointercov(H,RRZFpre,N,1,K);
        SR_RRZF(idx1,idx2) = SUMrateMISO(H,RRZFpre,Rrrzf,w,K);
        
        finit = reshape(MRTpre,N*K,1);
        GPIPpre = GPIPalgorithm(EH,Phi,finit,w,sigma2,P,epsilon);
        SR_GPIP(idx1,idx2) = SRofGPIP(H,GPIPpre,w,sigma2,P);
        
        GPIPnocovpre = GPIPalgorithm(EH,zeros(N,N,K),finit,w,sigma2,P,epsilon);
        SR_GPIPnoCov(idx1,idx2) = SRofGPIP(H,GPIPnocovpre,w,sigma2,P);

        [Bwmmse,SR] = PREandSRwithImperfectCSI(H,EH,P,w,N,1,K,20);
        SR_WMMSE(idx1,idx2) = SR(end);
    end
end

SR_GPIPmean = mean(SR_GPIP,2);
SR_GPIPnoCovmean = mean(SR_GPIPnoCov,2);
SR_RZFmean = mean(SR_RZF,2);
SR_MRTmean = mean(SR_MRT,2);
SR_RRZFmean = mean(SR_RRZF,2);
SR_WMMSEmean = mean(SR_WMMSE,2);

figure(1)
plot(snr,SR_GPIPmean,'r-+','linewidth',1.5)
hold on
plot(snr,SR_GPIPnoCovmean,'y-<','linewidth',1.5)
plot(snr,SR_RZFmean,'g->','linewidth',1.5)
plot(snr,SR_MRTmean,'b-*','linewidth',1.5)
plot(snr,SR_WMMSEmean,'k-s','linewidth',1.5)
plot(snr,SR_RRZFmean,'m-d','linewidth',1.5)
xlabel('SNR (dB)')
ylabel('Sum Rate (bps/Hz)')
legend('GPIP','GPIP no Cov.','RZF','MRT','WMMSE','RRZF','Location','Northwest')
set(gca,'Fontsize',11)
grid on
set(gca,'ygrid','on','gridlinestyle',':','Gridalpha',0.5)