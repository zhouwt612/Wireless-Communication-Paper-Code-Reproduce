clc
clear
% close all

K = 4;
M = 4;
% B = [5,10,15,20];
B = [5 10];

snr = 0:5:30;
Hnum = 5000;

SRdata = zeros(length(B),length(snr));
for idx1 = 1:1:length(B)
    for idx2 = 1:1:length(snr)
        disp(['Feedback bits:' num2str(B(idx1)) ' SNR:' num2str(snr(idx2))])
        SRdata(idx1,idx2) = one_transmission(K,M,snr(idx2),B(idx1),Hnum);
    end
end

Ccsit = [2.20;3.63;5.16;6.81;8.46;10.11;11.78];

figure
plot(snr,Ccsit,'b-^','linewidth',2);
hold on
plot(snr,SRdata(1,:),'r-*','linewidth',2);
plot(snr,SRdata(2,:),'g-+','linewidth',2);
% plot(snr,SRdata(3,:),'b-<');
% plot(snr,SRdata(4,:),'k->');
xlabel('SNR(dB)');
ylabel('Throughput (bps/Hz)')
legend('CSIT','B=5','B=10')
axis([0 30 0 20])
