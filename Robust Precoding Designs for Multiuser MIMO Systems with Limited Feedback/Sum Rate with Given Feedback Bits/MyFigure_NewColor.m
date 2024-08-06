function MyFigure_NewColor(InputDataSet,snr,SetNum)

figure
hold on
plot(snr,InputDataSet{1}(SetNum,:),'Color',[0.68,0.04,0.68],'Marker','+','linewidth',1.5)
plot(snr,InputDataSet{2}(SetNum,:),'Color',[0.08,0.64,0.16],'Marker','^','linewidth',1.5)
plot(snr,InputDataSet{3}(SetNum,:),'Color',[0.10,0.07,0.75],'Marker','s','linewidth',1.5)
plot(snr,InputDataSet{4}(SetNum,:),'k->','linewidth',1.5)
plot(snr,InputDataSet{5}(SetNum,:),'Color',[0.10,0.07,0.75],'Marker','<','linewidth',1.5)
plot(snr,InputDataSet{6}(SetNum,:),'k-o','linewidth',1.5)
xlabel('SNR (dB)')
ylabel('Sum Rate (bps/Hz)')
legend('MRT','BD','MMSE','WMMSE','RMMSE','RWMMSE','Location','northwest')
set(gca,'FontSize',11);
set(gca,'LineWidth',1)
box on
grid on
set(gca,'GridLineStyle',':','GridColor','k','GridAlpha',0.5)
% title('M=8,N=2,K=4')