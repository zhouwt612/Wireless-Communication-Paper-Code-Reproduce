% for i1=1:K_all_Num,
%     figure(i1*2-1)
%     semilogy(SNR,BER1(i1,:),'g-o');
%     axis([-15 35 10^(-6) 10^(0)])
%     xlabel('Average SNR (dB)');
%     ylabel('Average BER');
%     hold on;
%     semilogy(SNR,BER3(i1,:),'b-*');
%     hold on;
%     semilogy(SNR,BER2(i1,:),'r-^');
%     hold on;
%     semilogy(SNR,BER4(i1,:),'k-v');
%     hold on;
%     semilogy(SNR,BER5(i1,:),'c-x');
% %     hold on;
% %     semilogy(SNR,BER6(i1,:),'m-v');
%     legend('channel inversion','BD','MMSE','SLNR','SINR');
%     title('BER vs SNR');
% 
%     figure(i1*2)
%     plot(SNR,capacity_sum1(i1,:),'g-o');
%     grid on;
%     xlabel('Average SNR (dB)');
%     ylabel('Sum capacity (bps/Hz)');
%     plot(SNR,capacity_sum3(i1,:),'b-*');
%     hold on;
%     plot(SNR,capacity_sum2(i1,:),'r-^');
%     hold on;
%     plot(SNR,capacity_sum4(i1,:),'k-v');
%     hold on;
%     plot(SNR,capacity_sum5(i1,:),'c-x');
% %     hold on;
% %     plot(SNR,capacity_sum6(i1,:),'m-v');
%     legend('channel inversion','BD','MMSE','SLNR','SINR');
%     title('Sum capacity vs SNR');
%     
% end



% for i1=1:K_all_Num,
%     figure(i1*2-1)
%     semilogy(SNR,BER3(i1,:),'b-*');
%     axis([-15 35 10^(-6) 10^(0)])
%     xlabel('Average SNR (dB)');
%     ylabel('Average BER');
% %     hold on;
% %     hold on;
% %     semilogy(SNR,BER2(i1,:),'r-^');
%     hold on;
%     semilogy(SNR,BER4(i1,:),'k-v');
%     hold on;
%     semilogy(SNR,BER5(i1,:),'c-x');
% %     hold on;
% %     semilogy(SNR,BER6(i1,:),'m-v');
%     legend('ZF','SLNR','SINR');
%     title('BER vs SNR');
% 
%     figure(i1*2)
%     plot(SNR,capacity_sum3(i1,:),'b-*');
%     xlabel('Average SNR (dB)');
%     ylabel('Sum capacity (bps/Hz)');
% %     hold on;
% %     hold on;
% %     plot(SNR,capacity_sum2(i1,:),'r-^');
%     hold on;
%     plot(SNR,capacity_sum4(i1,:),'k-v');
%     hold on;
%     plot(SNR,capacity_sum5(i1,:),'c-x');
% %     hold on;
% %     plot(SNR,capacity_sum6(i1,:),'m-v');
%     legend('ZF','SLNR','SINR');
%     title('Sum capacity vs SNR');
%     
% end

% for i1=1:K_all_Num,
%     figure(i1*2-1)
%     semilogy(SNR,BER3(i1,:),'b-*');
%     axis([-15 35 10^(-6) 10^(0)])
%     xlabel('Average SNR (dB)');
%     ylabel('Average BER');
% %     hold on;
% %     hold on;
% %     semilogy(SNR,BER2(i1,:),'r-^');
%     hold on;
%     semilogy(SNR,BER4(i1,:),'k-v');
%     hold on;
%     semilogy(SNR,BER5(i1,:),'c-x');
% %     hold on;
% %     semilogy(SNR,BER6(i1,:),'m-v');
%     legend('SLNR','SLNR+PD','SINR');
%     title('BER vs SNR');
% 
%     figure(i1*2)
%     plot(SNR,capacity_sum3(i1,:),'b-*');
%     xlabel('Average SNR (dB)');
%     ylabel('Sum capacity (bps/Hz)');
% %     hold on;
% %     hold on;
% %     plot(SNR,capacity_sum2(i1,:),'r-^');
%     hold on;
%     plot(SNR,capacity_sum4(i1,:),'k-v');
%     hold on;
%     plot(SNR,capacity_sum5(i1,:),'c-x');
% %     hold on;
% %     plot(SNR,capacity_sum6(i1,:),'m-v');
%     legend('SLNR','SLNR+PD','SINR');
%     title('Sum capacity vs SNR');
%     
% end

% figure(i1*2-1)
% semilogy(SNR,BER1(i1,:),'b-o');
% grid on;
% axis([-15 35 10^(-6) 10^(0)])
% xlabel('Average SNR (dB)');
% ylabel('Average BER');
% hold on;
% semilogy(SNR,BER2(i1,:),'g-*');
% hold on;
% semilogy(SNR,BER3(i1,:),'k-^');
% hold on;
% semilogy(SNR,BER4(i1,:),'r-v');
% % hold on;
% % semilogy(SNR,BER5(i1,:),'c-x');
% %     hold on;
% %     semilogy(SNR,BER6(i1,:),'m-v');
% legend('SINR','SLNR','SLNR+CN','SLNR+SLNR');
% title('BER vs SNR');

figure(i1*2)
plot(SNR,capacity_sum1(i1,:),'b-o');
grid on
xlabel('Average SNR (dB)');
ylabel('Sum capacity (bps/Hz)');
hold on;
plot(SNR,capacity_sum2(i1,:),'g-*');
hold on;
plot(SNR,capacity_sum3(i1,:),'k-^');
hold on;
plot(SNR,capacity_sum4(i1,:),'r-v');
% hold on;
% plot(SNR,capacity_sum5(i1,:),'c-x');
%     hold on;
%     plot(SNR,capacity_sum6(i1,:),'m-v');
legend('SINR','SLNR','SLNR+CN','SLNR+SLNR');
title('Sum capacity vs SNR');