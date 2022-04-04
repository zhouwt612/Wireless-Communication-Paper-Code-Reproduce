for i1=1:K_all_Num,
    figure(i1*2-1)
    semilogy(SNR,BER1(i1,:),'g--');
    axis([-15 35 10^(-3) 10^(0)])
    xlabel('Average SNR (dB)');
    ylabel('Average BER');
    hold on;
    semilogy(SNR,BER4(i1,:),'b-*');
    hold on;
    semilogy(SNR,BER3(i1,:),'r-^');
    hold on;
    semilogy(SNR,BER5(i1,:),'k-x');
    hold on;
    semilogy(SNR,BER6(i1,:),'c-v');
    legend('CSIT','EMMSE','DSLNR(Single)','ESLNR(Single)','ESLNR(MMSLNR)');
    title('BER vs SNR');

    figure(i1*2)
    plot(SNR,capacity_sum1(i1,:),'g--');
    xlabel('Average SNR (dB)');
    ylabel('Sum capacity (bps/Hz)');
    hold on;
    plot(SNR,capacity_sum4(i1,:),'b-*');
    hold on;
    plot(SNR,capacity_sum3(i1,:),'r-^');
    hold on;
    plot(SNR,capacity_sum5(i1,:),'k-x');
    hold on;
    plot(SNR,capacity_sum6(i1,:),'c-v');
    legend('CSIT','EMMSE','DSLNR(Single)','ESLNR(Single)','ESLNR(MMSLNR)');
    title('Sum capacity vs SNR');
end