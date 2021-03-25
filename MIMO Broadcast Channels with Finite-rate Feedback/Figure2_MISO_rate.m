clc
clear
close

snr = 0:5:30;
nopow = 1;
P_data = 10.^(snr./10).*nopow;

hnum = 10000;
M = 4;
K = 1;
B = M-1;

Ccsit_data = zeros(length(snr),hnum);
Cnocsit_data = zeros(length(snr),hnum);
Rfb_data = zeros(length(snr),hnum);

for idx1 = 1:1:length(snr)
    for idx2 = 1:1:hnum
        h = sqrt(1/2)*(randn(M,K)+1i*randn(M,K));
        Ccsit = log2(1+P_data(idx1)*norm(h)^2);
        Cnocsit = log2(1+P_data(idx1)/M*norm(h)^2);
        Rfb = log2(1+P_data(idx1)*norm(h)^2*(1-2^(-B/(M-1))));
        
        Ccsit_data(idx1,idx2) = Ccsit;
        Cnocsit_data(idx1,idx2) = Cnocsit;
        Rfb_data(idx1,idx2) = Rfb;
    end
end

Ccsit_mean = mean(Ccsit_data,2);
Cnocsit_mean = mean(Cnocsit_data,2);
Rfb_mean = mean(Rfb_data,2);

figure
plot(snr,Ccsit_mean,'r-*','linewidth',2)
hold on
plot(snr,Cnocsit_mean,'b-+','linewidth',2)
plot(snr,Rfb_mean,'k-<','linewidth',2)
xlabel('SNR (dB)');ylabel('Rate (bps/hz)')
legend('CSIT','No CSIT','Finite Rate Feedback')
grid on
