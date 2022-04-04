for i1=1:K_all_Num,
    figure(i1*2-1)
    semilogy(SNR,BER1(i1,:),'g--');
    axis([-15 35 10^(-4) 10^(0)])
    xlabel('Average SNR (dB)');
    ylabel('Average BER');
    hold on;
    semilogy(SNR,BER4(i1,:),'b-*');
    hold on;
    semilogy(SNR,BER2(i1,:),'r-^');
    hold on;
    semilogy(SNR,BER3(i1,:),'k-x');
    legend('CSIT','EMMSE','DSLNR','ESLNR');
    title('BER vs SNR');

    figure(i1*2)
    plot(SNR,capacity_sum1(i1,:),'g--');
    xlabel('Average SNR (dB)');
    ylabel('Sum capacity (bps/Hz)');
    hold on;
    plot(SNR,capacity_sum4(i1,:),'b-*');
    hold on;
    plot(SNR,capacity_sum2(i1,:),'r-^');
    hold on;
    plot(SNR,capacity_sum3(i1,:),'k-x');
    legend('CSIT','EMMSE','DSLNR','ESLNR');
    title('Sum capacity vs SNR');
end

% for i1=1:K_all_Num,
%     figure(i1*2-1)
%     semilogy(SNR,BER1(i1,:),'g--');
%     axis([-15 35 10^(-4) 10^(0)])
%     xlabel('Average SNR (dB)');
%     ylabel('Average BER');
%     hold on;
%     semilogy(SNR,BER2(i1,:),'b-*');
%     hold on;
%     semilogy(SNR,BER3(i1,:),'r-^');
%     hold on;
%     semilogy(SNR,BER4(i1,:),'k-x');
%     legend('CSIT','B=4','B=8','B=12');
%     title('BER vs SNR');
% 
%     figure(i1*2)
%     plot(SNR,capacity_sum1(i1,:),'g--');
%     xlabel('Average SNR (dB)');
%     ylabel('Sum capacity (bps/Hz)');
%     hold on;
%     plot(SNR,capacity_sum2(i1,:),'b-*');
%     hold on;
%     plot(SNR,capacity_sum3(i1,:),'r-^');
%     hold on;
%     plot(SNR,capacity_sum4(i1,:),'k-x');
%     legend('CSIT','B=4','B=8','B=12');
%     title('Sum capacity vs SNR');
% end