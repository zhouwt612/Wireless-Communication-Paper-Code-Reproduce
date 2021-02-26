clc
clear
close all

r = 3;
t = [1 2 3 5];
N = t;
epsilon = 0.4;
pow = 1;

snr = 3:1:19;
snr_linear = 10.^(snr./10);

MMSEchanestimse = zeros(length(t),length(snr));
for idx1 = 1:1:length(t)
    for idx2 = 1:1:length(snr)
        P = orth_pilot(N(idx1),t(idx1),pow);
        Rh = chan_corr_mtx(r,t(idx1),epsilon);
        sigma2 = pow/snr_linear(idx2);
        [Jmmse,~] = MMSEmse(r,P,sigma2,Rh);
        MMSEchanestimse(idx1,idx2) = Jmmse/(t(idx1)*r);
    end
end

figure
semilogy(snr,MMSEchanestimse(1,:),'-.r','LineWidth',1)
hold on
semilogy(snr,MMSEchanestimse(2,:),'-xg','LineWidth',1)
semilogy(snr,MMSEchanestimse(3,:),'-vb','LineWidth',1)
semilogy(snr,MMSEchanestimse(4,:),'-sk','LineWidth',1)
xlabel('Pow/sigma^2');ylabel('Normalized MSE');
legend('t=1,r=3','t=2,r=3','t=3,r=3','t=5,r=3')