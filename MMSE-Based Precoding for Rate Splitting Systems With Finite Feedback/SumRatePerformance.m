clc
clear
close all

K = 2;
M = 4;
snr = 0:5:30;
FB = [4 10];
sig2 = 1;
u = ones(1,M);
hnum = 5000;
SumRateRS = zeros(length(FB),length(snr));
SumRate = zeros(length(FB),length(snr));
for idx1 = 1:1:length(FB)
    for idx2 = 1:1:length(snr)
        disp(['Feedback bits: ', num2str(FB(idx1)), ' SNR: ', num2str(snr(idx2))])
        pow = sig2 * 10^(snr(idx2)/10);
        t = RSpowerallocation(M,K,FB(idx1),pow);
        tau2 = imperfectCSIlevel(FB(idx1),M);
        sumrateRS = zeros(1,hnum);
        sumrate = zeros(1,hnum);
        for idx3 = 1:1:hnum
            H = sqrt(1/2)*(randn(K,M)+1i*randn(K,M));
            C = RVQ(M,FB(idx1));
            H_hat = quantizedchannel(H,C);
            % Rate splitting
            [Pc,Pp] = RCPrecoding(H_hat,pow,t,tau2);
            Rc = SumRateMISO_C(H,Pc,Pp,pow,t);
            Rp = SumRateMISO_P(H,Pp,pow,t);
            sumrateRS(idx3) = Rc + sum(Rp);
            % Linear precoding
            Pmmse = RobustMMSE(H_hat,pow,M,FB(idx1));
            Rn = nointercov(H,Pmmse,M,1,K);
            sumrate(idx3) = SUMrateMISO(H,Pmmse,Rn,ones(K),K);
        end
        SumRateRS(idx1,idx2) = mean(sumrateRS);
        SumRate(idx1,idx2) = mean(sumrate);
    end
end

SumRateRS
SumRate

figure
plot(snr,SumRateRS(1,:),'k--+','linewidth',2)
hold on
plot(snr,SumRateRS(2,:),'k-+','linewidth',2)
plot(snr,SumRate(1,:),'r--o','linewidth',2)
plot(snr,SumRate(2,:),'r-o','linewidth',2)
grid on
xlabel('SNR (dB)')
ylabel('Sum Rate (bps/Hz)')
legend('RS MMSE B = 5','RS MMSE B = 10','MMSE B = 5','MMSE B = 10','Location','Northwest')
set(gca,'FontSize',14);