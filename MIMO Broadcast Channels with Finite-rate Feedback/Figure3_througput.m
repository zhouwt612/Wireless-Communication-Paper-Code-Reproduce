% clc
% clear
% close

snr = 0:5:30;
nopow = 1;
P_data = 10.^(snr./10).*nopow;

for idx = 1:1:length(snr)
%     h = sqrt(1/2)*(randn(5)+1i*randn(5));
    Rp(idx) = SRzfwithff(P_data(idx),100,5);
end

Rp
