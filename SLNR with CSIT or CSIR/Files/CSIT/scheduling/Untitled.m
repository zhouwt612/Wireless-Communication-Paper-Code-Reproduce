for i1=1:Kall_Num,
    figure(4*i1-3)
    semilogy(SNR,BER1(i1,:),'k--');
    axis([-15 35 10^(-6) 10^(0)])
    xlabel('Average SNR (dB)');
    ylabel('Average BER');
    hold on;
    semilogy(SNR,BER2(i1,:),'g-+');
%     hold on;
%     semilogy(SNR,BER3(i1,:),'r-o');
%     hold on;
%     semilogy(SNR,BER4(i1,:),'b-*');
%     hold on;
%     semilogy(SNR,BER5(i1,:),'c-v');
    hold on;
    semilogy(SNR,BER6(i1,:),'b-^');

    legend('Round Robin','Max SNR','MMSLNR');
    title('BER vs SNR');
    hold off;

    figure(4*i1-2)
    plot(SNR,capacity_sum1(i1,:),'k--');
    xlabel('Average SNR (dB)');
    ylabel('Sum capacity (bps/Hz)');
    hold on;
    plot(SNR,capacity_sum2(i1,:),'g-+');
%     hold on;
%     plot(SNR,capacity_sum3(i1,:),'r-o');
%     hold on;
%     plot(SNR,capacity_sum4(i1,:),'b-*');
%     hold on;
%     plot(SNR,capacity_sum5(i1,:),'c-v');
    hold on;
    plot(SNR,capacity_sum6(i1,:),'b-^');

    legend('Round Robin','Max SNR','MMSLNR');
    title('Sum capacity vs SNR');
    hold off;

    figure(4*i1-1)
    plot(SNR,round_num_sum4(i1,:)-1,'r-o');
    xlabel('Average SNR (dB)');
    ylabel('Amount of the replacement round');
    hold on;
    plot(SNR,round_num_sum6(i1,:)-1,'b-*');
    legend('Kall=20','Kall=50');
    title('Amount of the replacement round vs SNR');
    hold off;
    
    figure(4*i1)
    plot(SNR,cal_num_sum4(i1,:),'r-o');
    xlabel('average SNR in dB');
    ylabel('amount of the calculation of SLNR');
    hold on;
    plot(SNR,cal_num_sum6(i1,:),'b-*');
    legend('simplified MM-SLNR','MM-SLNR');
    title('computational complexity vs SNR');
    hold off;
    
end