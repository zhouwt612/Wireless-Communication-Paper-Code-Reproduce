clc
clear
% close

M = 4;
K = 4;
B = [5 10];
snr = 0:5:30;
Hnum = 1000;
sigma2 = 1;
sumratedata = zeros(length(B),length(snr));
for idx1 = 1:1:length(B)
    for idx2 = 1:1:length(snr)
        disp(['Feedback bits:' num2str(B(idx1)) ' SNR:' num2str(snr(idx2))])
        P = sigma2*10^(snr(idx2)/10);
        sumrate1 = zeros(1,Hnum);
        for idx3 = 1:1:Hnum
            C = RVQcodebook(M,B(idx1));
            H = channel(K,M);
            F = quantizedchannel(H,C);
            precoder = RegZFprecoder(F,P);
            sumrate1(1,idx3) = SRviaSINR(H,precoder,P);
        end
        sumrate(1,idx2) = mean(sumrate1);
    end
    sumratedata(idx1,:) = sumrate;
end

figure
plot(snr,sumratedata(1,:),'r-*','linewidth',2);
hold on
plot(snr,sumratedata(2,:),'r-*','linewidth',2)
xlabel('SNR(dB)');
ylabel('Throughput (bps/Hz)')

