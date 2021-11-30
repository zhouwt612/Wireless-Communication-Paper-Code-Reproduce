clc
clear
% close
% A simulation of WMMSE procoder in MIMO broadcast channel with Limited-feedback
% Assuming that the number of the transmit antennas and the number of the
% UE with single antenna are the same.
% The power of noise is set to 1, and the precoder is normalized.

M = 4;
K = 4;
B = [5 10];
snr = 0:5:30;
Hnum = 5000;
sigma2 = 1;
sumratedata = zeros(length(B),length(snr));
for idx1 = 1:1:length(B)
    sumrate = zeros(1,length(snr));
    for idx2 = 1:1:length(snr)
        disp(['Feedback bits:' num2str(B(idx1)) ' SNR:' num2str(snr(idx2))])
        P = sigma2*10^(snr(idx2)/10);
        sumrate1 = zeros(1,Hnum);
        for idx3 = 1:1:Hnum
            C = RVQ(M,B(idx1));
            H = channel(K,M);
            F = quantizedchannel(H,C);
            Brmmse = RobustMMSE(F,P,M,B(idx1));
            sumrate1(idx3) = SRergodic(H,Brmmse);
%             Rn = nointercov(H,Brmmse,M,1,K);
%             sumrate1(idx3) = SUMrateMISO(H,Brmmse,Rn,ones(K),K);
        end
        sumrate(1,idx2) = K*mean(sumrate1);
    end
    sumratedata(idx1,:) = sumrate;
end

figure
plot(snr,sumratedata(1,:),'r-*','linewidth',2);
hold on
plot(snr,sumratedata(2,:),'g-+','linewidth',2)
grid on
xlabel('SNR(dB)');
ylabel('Throughput (bps/Hz)')
legend('B=5','B=10')
axis([0 30 0 16])
set(gca, 'YTick', 0:16/8:16);
set(gca, 'XTick', 0:30/5:30);
set(gca,'FontSize',14);