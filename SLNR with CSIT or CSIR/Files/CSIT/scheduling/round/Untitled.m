for i1=1:K_all_Num,
%     figure(i1*3-2)
%     semilogy(SNR,BER1(i1,:,1),'g--');
%     axis([-15 35 10^(-6) 10^(0)])
%     xlabel('Average SNR (dB)');
%     ylabel('Average BER');
%     hold on;
%     semilogy(SNR,BER1(i1,:,2),'k-*');
%     hold on;
%     semilogy(SNR,BER1(i1,:,3),'r-o');
%     hold on;
%     semilogy(SNR,BER1(i1,:,4),'b-+');
%     hold on;
%     semilogy(SNR,BER1(i1,:,5),'c-v');
% 
%     title('BER vs SNR');
%     legend('Q=0','Q=1','Q=2','Q=3','Q=4');
%     hold off;
% 
%     figure(i1*3-1)
%     plot(SNR,capacity_avg1(i1,:,1),'g--');
%     xlabel('Average SNR (dB)');
%     ylabel('Sum capacity');
%     hold on;
%     plot(SNR,capacity_avg1(i1,:,2),'k-*');
%     hold on;
%     plot(SNR,capacity_avg1(i1,:,3),'r-o');
%     hold on;
%     plot(SNR,capacity_avg1(i1,:,4),'b-+');
%     hold on;
%     plot(SNR,capacity_avg1(i1,:,5),'c-v');
% 
%     title('Sum capacity vs SNR');
%     legend('Q=0','Q=1','Q=2','Q=3','Q=4');
%     hold off;
    
    figure(i1*3)
    plot(SNR,eig_num_avg1(i1,:,1),'g--');
    axis([-15 35 0 260])
    xlabel('Average SNR (dB)');
    ylabel('Amount of calculating eigenvalue problem');
    hold on;
    plot(SNR,eig_num_avg1(i1,:,2),'k-*');
    hold on;
    plot(SNR,eig_num_avg1(i1,:,3),'r-o');
    hold on;
    plot(SNR,eig_num_avg1(i1,:,4),'b-+');
    hold on;
    plot(SNR,eig_num_avg1(i1,:,5),'c-v');

    title('Computational complexity');
    legend('Q=0','Q=1','Q=2','Q=3','Q=4');
    hold off;
    
end