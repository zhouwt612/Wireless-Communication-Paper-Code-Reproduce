for i1=1:1,
    P=zeros(1,Kall);
    for i2=1:400000,
        P(1,scheduled_sum(i1,i2))=P(1,scheduled_sum(i1,i2))+1;
    end
end

stem([1:20],P/400000*100)
axis([0 20 0 8]);
xlabel('Average SNR (dB)');
ylabel('Percentage of time (%)');
title('Fairness');