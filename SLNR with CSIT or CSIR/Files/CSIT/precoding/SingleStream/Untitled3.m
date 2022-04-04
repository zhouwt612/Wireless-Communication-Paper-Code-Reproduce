for i1=1:K_all_Num,    
    figure(i1*2-1)
    semilogy(SNR,BER3(i1,:),'g-o');
    axis([-15 35 10^(-6) 10^(0)])
    xlabel('Average SNR (dB)');
    ylabel('Average BER');
    hold on;
    semilogy(SNR,BER4(i1,:),'b-*');
    hold on;
    semilogy(SNR,BER5(i1,:),'r-^');
%     hold on;
%     semilogy(SNR,BER4(i1,:),'k--');
%     hold on;
%     semilogy(SNR,BER5(i1,:),'c-x');
%     hold on;
%     semilogy(SNR,BER6(i1,:),'m-v');
    legend('channel norm','SINR','SLNR');
    title('BER vs SNR');

    figure(i1*2)
    plot(SNR,capacity_sum3(i1,:),'g-o');
    xlabel('Average SNR (dB)');
    ylabel('Sum capacity (bps/Hz)');
    hold on;
    plot(SNR,capacity_sum4(i1,:),'b-*');
    hold on;
    plot(SNR,capacity_sum5(i1,:),'r-^');
%     hold on;
%     plot(SNR,capacity_sum4(i1,:),'k--');
%     hold on;
%     plot(SNR,capacity_sum5(i1,:),'c-x');
%     hold on;
%     plot(SNR,capacity_sum6(i1,:),'m-v');
    legend('channel norm','SINR','SLNR');
    title('Sum capacity vs SNR');
end