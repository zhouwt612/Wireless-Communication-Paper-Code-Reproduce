clc
clear
close all

% M: # of tx antennas
% N: # of rx antennas
% K: # of users
M = 8;
N = 2;
K = 4;

snr = 20;
pow = 10^(snr/10);
sigma2 = 1;
sigmae = 0;
Hnum = 1;
itea_num = 20;

R_total = zeros(1,itea_num);
for idx0 = 1:1:Hnum
    [H,Hhat,E] = channel(M,N,K,sigmae);
    [H,Hhat] = Channel_Sort_Imp(H,Hhat,pow,sigma2);
    [~,~,Pc_itea,Pp_itea] = SNS_Precoding(H,pow,sigma2,itea_num);
    for idx1 = 1:1:itea_num
        R_total(idx1) = R_total(idx1) + SumRateMIMOforK_RSMA(H,Pc_itea(:,:,idx1),Pp_itea(:,:,:,idx1));
    end
end
R_total = R_total/Hnum;
figure
plot(1:1:itea_num,R_total,'r-+','linewidth',1)
xlabel('Iteration')
ylabel('Ergodic Sum Rate (bps/Hz)')
set(gca,'FontSize',11);
set(gca,'LineWidth',1)
box on
grid on
axis([0 20 0 50])
set(gca,'GridLineStyle',':','GridColor','k','GridAlpha',0.5)