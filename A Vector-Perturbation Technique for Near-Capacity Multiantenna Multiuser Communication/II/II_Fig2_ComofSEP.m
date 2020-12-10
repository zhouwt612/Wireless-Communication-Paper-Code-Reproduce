%% Without perturbation
clc
clear
% close

snr = 5:5:35;
K = 10;
M = 16;
datanum = 10000;

CIberdata = zeros(1,length(snr));
RCIberdata = zeros(1,length(snr));
for ind = 1:1:length(snr)
    disp(['Signal to noise ratio: ' num2str(snr(ind)) ' dB']);
    CIber = invM(snr(ind),M,K,datanum);
    RCIber = reginvM(snr(ind),M,K,datanum);
    CIberdata(1,ind) = CIber;
    RCIberdata(1,ind) = RCIber;
end

figure
semilogy(snr,CIberdata,'r-x',snr,RCIberdata,'r-.d')
xlabel('rho (dB)');ylabel('Average Pro (err)');
legend('Channel inversion','Regularized Inversion');
axis([5 35 10^(-4) 10^0])

%% With perturbation
clc
clear
% close

snr = 5:5:35;
K = 10;
M = 4;
datanum = 10000;

CIpberdata = zeros(1,length(snr));
for ind = 1:1:length(snr)
    disp(['Signal to noise ratio: ' num2str(snr(ind)) ' dB']);
    alpha = K/snr(ind);
    ber = invMperM(snr(ind),M,K,datanum,alpha);
    CIpberdata(1,ind) = ber;    
end
CIpberdata
semilogy(snr,CIpberdata,'r->')
% 'Sphere Encoder'