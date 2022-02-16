clc
clear
close

M = 6;
N = 2;
K = 3;
snr = 0:5:30;
B = [5 10];
hnum = 1000;
SumRate = zeros(length(B),length(snr));
for idx1 = 1:1:length(B)
    for idx2 = 1:1:length(snr)
        disp(['Feedback bits: ' num2str(B(idx1)) ' SNR: ' num2str(snr(idx2))])
        C = RVQ_MIMOforK(M,N,B(idx1),K);
        pow = 10^(snr(idx2)/10);
        Rate = zeros(1,hnum);
        for idx3 = 1:1:hnum
            H = channel(M,N,K);
            Htilde = Hbasis(H);
            F = quantizedchannel_MIMO(Htilde,C);
            P = BD_MIMOforK(F,pow);
            Rate(idx3) = SumRateMIMOforK(H,P);
        end
        SumRate(idx1,idx2) = mean(Rate);
    end
end

SumRate

figure
plot(snr,SumRate(1,:),'r-+','linewidth',2)
hold on
plot(snr,SumRate(2,:),'k-*','linewidth',2)
xlabel('SNR (dB)')
ylabel('Sum Rate (bps/Hz)')
legend('B=5','B=10')
axis([0 30 0 10])