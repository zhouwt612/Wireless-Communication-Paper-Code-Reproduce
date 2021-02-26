clc
clear
close all

snr = 1:1:19;
snr_linear = 10.^(snr/10);

t = 1:1:5;
r = 3;
epsilon = 0;
hnum = 5000;

SLSchanestimse = zeros(length(t),length(snr));
LSSLSchanestimse = zeros(length(t),length(snr));
LSSLSdata = zeros(length(t),length(snr),hnum);
for idx1 = 1:1:length(t)
    for idx2 = 1:1:length(snr)
        SLSchanestimse(idx1,idx2) = ...
            SLSmse(trace(chan_corr_mtx(r,t(idx1),epsilon)),t(idx1),r,1/snr_linear(idx2),1)/(t(idx1)*r);
        for idx3 = 1:1:hnum
            [~,Hls] = LStransmission(t(idx1),r,t(idx1),snr(idx2),1);
             LSSLSdata(idx1,idx2,idx3) = ...
                SLSmse(trace(Hls'*Hls),t(idx1),r,1/snr_linear(idx2),1)/(t(idx1)*r);
        end
    end
end

for idx4 = 1:1:length(t)
    for idx5 = 1:1:length(snr)
        LSSLSchanestimse(idx4,idx5) = mean(LSSLSdata(idx4,idx5,:));
    end
end



figure
semilogy(snr,SLSchanestimse(1,:),'-r','LineWidth',1)
hold on
semilogy(snr,SLSchanestimse(2,:),'-xg','LineWidth',1)
semilogy(snr,SLSchanestimse(3,:),'-vb','LineWidth',1)
semilogy(snr,SLSchanestimse(4,:),'-sk','LineWidth',1)
semilogy(snr,SLSchanestimse(5,:),'-cp','LineWidth',1)

semilogy(snr,LSSLSchanestimse(1,:),'--r','LineWidth',1)
semilogy(snr,LSSLSchanestimse(2,:),'--xg','LineWidth',1)
semilogy(snr,LSSLSchanestimse(3,:),'--vb','LineWidth',1)
semilogy(snr,LSSLSchanestimse(4,:),'--sk','LineWidth',1)
semilogy(snr,LSSLSchanestimse(5,:),'--cp','LineWidth',1)
xlabel('Pow/sigma^2');ylabel('Normalized MSE');
legend('t=1,r=3','t=2,r=3','t=3,r=3','t=4,r=3','t=5,r=3',...
    't=1,r=3','t=2,r=3','t=3,r=3','t=4,r=3','t=5,r=3')