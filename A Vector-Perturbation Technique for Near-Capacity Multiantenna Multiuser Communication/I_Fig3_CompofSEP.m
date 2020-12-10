clc
clear
% close

snr = -10:5:30;
datanum = 10000;
M = 4;
K = [4, 10];
invber = zeros(2,length(snr));
reginvber = zeros(2,length(snr));
for ind1 = 1:1:length(snr)
    disp(['Signal to noise ratio: ' num2str(snr(ind1)) ' dB'])
    invber4 = invtransmission(snr(ind1),K(1),datanum);
    reginvber4 = reginvtransmission(snr(ind1),K(1),datanum);
    invber10 = invtransmission(snr(ind1),K(2),datanum);
    reginvber10 = reginvtransmission(snr(ind1),K(2),datanum);
    invber(:,ind1) = [invber4; reginvber4];
    reginvber(:,ind1) = [invber10; reginvber10];
end

invber
reginvber

figure
semilogy(snr,invber(1,:),'r--x','linewidth',1.5)
hold on
semilogy(snr,invber(2,:),'r-d','linewidth',1.5)
semilogy(snr,reginvber(1,:),'b-+','linewidth',1.5)
semilogy(snr,reginvber(2,:),'b-*','linewidth',1.5)
xlabel('rho (dB)');ylabel('BER')
legend('4x4 Channal Inversion','4x4 Regularized Channal Inversion',...
    '10x10 Channal Inversion','10x10 Regularized Channal Inversion')
