clc
clear
close all

M = 8;
N = 2;
K = 4;

snr = 20;
pow = 10^(snr/10);
Hnum = 1;
itea_num = 50;

Rate_Itea = zeros(1,itea_num+1);

for idx = 1:1:Hnum
    H = channel(M,N,K);
    P = TransmitMF(H,pow);
    [~,Ritea] = WMMSE_MIMOforK(H,P,pow,itea_num);
    Rate_Itea = Rate_Itea + Ritea;
end
Rate_Itea = Rate_Itea/Hnum

figure
hold on
plot(0:1:itea_num,Rate_Itea,'r-+','LineWidth',1)
xlabel('SNR (dB)')
ylabel('Sum Rate (bps/Hz)')
set(gca,'FontSize',11);
set(gca,'LineWidth',1)
box on
grid on
axis([0 50 0 50])
set(gca,'GridLineStyle',':','GridColor','k','GridAlpha',0.5)