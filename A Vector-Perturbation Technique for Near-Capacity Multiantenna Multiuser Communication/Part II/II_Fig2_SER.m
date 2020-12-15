clc
clear
close

SNR = 5:5:35;
K = 10;
M = 10;
Modu = 4;
datanum = 100;
Lnum = 100;
tau = 2*(abs(3 + 1i*3)+2/2);

BER = zeros(1,length(SNR));
for idx = 1:1:length(SNR)
    snr = SNR(idx);
    disp(['SNR: ' num2str(snr)])
    BER(1,idx) = unreguperttrans(K,M,Modu,snr,datanum,Lnum,tau);
end
BER

figure
semilogy(SNR,BER,'r-<')
xlabel('rho (dB)');ylabel('Average Prob. (err)')