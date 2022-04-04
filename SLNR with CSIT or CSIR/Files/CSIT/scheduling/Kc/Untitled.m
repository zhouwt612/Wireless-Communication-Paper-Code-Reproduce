for i1=1:K_all_Num,
    figure(4*i1-3)
    semilogy(SNR,BER1(i1,:),'g--');
    axis([-15 35 10^(-6) 10^(0)])
    xlabel('Average SNR (dB)');
    ylabel('Average BER');
    hold on;
    semilogy(SNR,BER2(i1,:),'k-*');
    hold on;
    semilogy(SNR,BER3(i1,:),'r-o');
    hold on;
    semilogy(SNR,BER4(i1,:),'b-+');
    legend('Kc=2K','Kc=3K','Kc=4K','Kc=5K');
    title('BER vs SNR');
    hold off;

    figure(4*i1-2)
    plot(SNR,capacity_avg1(i1,:),'g--');
    xlabel('Average SNR (dB)');
    ylabel('Sum capacity');
    hold on;
    plot(SNR,capacity_avg2(i1,:),'k-*');
    hold on;
    plot(SNR,capacity_avg3(i1,:),'r-o');
    hold on;
    plot(SNR,capacity_avg4(i1,:),'b-+');
    legend('Kc=2K','Kc=3K','Kc=4K','Kc=5K');
    title('Sum capacity vs SNR');
    hold off;

    figure(4*i1-1)
    plot(SNR,round_num_avg1(i1,:)-1,'g--');
    xlabel('Average SNR (dB)');
    ylabel('Amount of the replacement rounds');
    hold on;
    plot(SNR,round_num_avg2(i1,:)-1,'k-*');
    hold on;
    plot(SNR,round_num_avg3(i1,:)-1,'r-o');
    hold on;
    plot(SNR,round_num_avg4(i1,:)-1,'b-+');
    legend('Kc=2K','Kc=3K','Kc=4K','Kc=5K');
    title('Amount of the replacement rounds vs SNR');
    hold off;

    figure(4*i1)
    plot(SNR,eig_num_avg1(i1,:),'g--');
    xlabel('Average SNR (dB)');
    ylabel('Amount of the solution of the eigenvalue problem');
    hold on;
    plot(SNR,eig_num_avg2(i1,:),'k-*');
    hold on;
    plot(SNR,eig_num_avg3(i1,:),'r-o');
    hold on;
    plot(SNR,eig_num_avg4(i1,:),'b-+');
    legend('Kc=2K','Kc=3K','Kc=4K','Kc=5K');
    title('Computational complexity');
    hold off;
    
end