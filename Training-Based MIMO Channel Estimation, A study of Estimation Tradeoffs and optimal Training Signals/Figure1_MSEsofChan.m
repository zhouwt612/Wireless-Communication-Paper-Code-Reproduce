clc
clear
close all

snr = 1:1:19;
snr_linear = 10.^(snr/10);

t = 1:1:5;
r = 3;

chanestimse = zeros(length(t),length(snr));
for idx1 = 1:1:length(t)
    for idx2 = 1:1:length(snr)
        chanestimse(idx1,idx2) = LSmse(snr_linear(idx2),t(idx1),r)/(t(idx1)*r);
    end
end

figure
semilogy(snr,chanestimse(1,:),'-.r','LineWidth',1.5)
hold on
semilogy(snr,chanestimse(2,:),'-xg','LineWidth',1.5)
semilogy(snr,chanestimse(3,:),'-vb','LineWidth',1.5)
semilogy(snr,chanestimse(4,:),'-sk','LineWidth',1.5)
semilogy(snr,chanestimse(5,:),'-cp','LineWidth',1.5)
xlabel('pow/sigma^2');ylabel('Normalized MSE');
legend('t=1,r=3','t=2,r=3','t=3,r=3','t=4,r=3','t=5,r=3')
