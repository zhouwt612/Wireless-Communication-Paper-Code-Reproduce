SNR_Num=11;
SNR=zeros(1,SNR_Num);
for i2=1:SNR_Num,
    SNR(1,i2)=-15+5*(i2-1);
end

figure(1)
semilogy(SNR,BER(1,:),'g--');
axis([-15 35 10^(-5) 10^(0)])
xlabel('average SNR in dB');
ylabel('average BER');
hold on;
semilogy(SNR,BER(2,:),'b-*');
hold on;
semilogy(SNR,BER(3,:),'r-^');
hold on;
semilogy(SNR,BER(4,:),'k-o');

legend('CSIT','4 bits','8 bits','12 bits');
title('BER vs SNR');

figure(2)
plot(SNR,capacity_sum(1,:),'g--');
xlabel('average SNR in dB');
ylabel('sum capacity');
hold on;
plot(SNR,capacity_sum(2,:),'b-*');
hold on;
plot(SNR,capacity_sum(3,:),'r-^');
hold on;
plot(SNR,capacity_sum(4,:),'k-o');

legend('CSIT','4 bits','8 bits','12 bits');
title('sum capacity vs SNR');