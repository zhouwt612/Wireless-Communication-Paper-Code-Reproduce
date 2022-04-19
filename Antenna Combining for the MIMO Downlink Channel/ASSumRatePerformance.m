% This program is for the reproduce of the paper:
% Antenna combining for the MIMO downlink siganl (Antenna selection)
% clc
clear
% close

M = 4;
N = 2;
K = 4;
snr = 0:5:30;
FB = [5 10];
sig2 = 1;
u = ones(1,M);
hnum = 5000;

ergodicC = zeros(length(FB),length(snr));
ergodicC2 = zeros(length(FB),length(snr));
for idx1 = 1:1:length(FB)
    for idx2 = 1:1:length(snr)
        disp(['Feedback bits: ', num2str(FB(idx1)), ' Signal-to-Noise Ratio: ', num2str(snr(idx2))])
        C = RVQforK(M,FB(idx1),K);
        pow = 10^(snr(idx2)/10);
        Capacity = zeros(1,hnum);
        Capacity2 = zeros(1,hnum);
        for idx3 = 1:1:hnum
            H = channel_forK(M,N,K);
            [Heff,F] = ASalgorithm(H,C);
            ZFpre = RegZFprecoder(F,pow);
            Capacity(idx3) = ACcapacity(Heff,ZFpre,pow);
        end
        ergodicC(idx1,idx2) = sum(Capacity)/hnum;
    end
end

ergodicC

figure
plot(snr,ergodicC(1,:),'r-+','linewidth',1.5)
hold on
plot(snr,ergodicC(2,:),'k-*','linewidth',1.5)
xlabel('SNR (dB)')
ylabel('Sum Rate (bps/Hz)')
axis([0 30 0 10])
legend('B = 5','B = 10')
set(gca,'FontSize',11);
grid on
set(gca,'GridLineStyle',':','GridColor','k','GridAlpha',1)