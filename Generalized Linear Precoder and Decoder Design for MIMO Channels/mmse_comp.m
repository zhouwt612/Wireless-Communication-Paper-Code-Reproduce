clc
clear
close

h_num = 5000;
snr = 0:2:20;
p0 = 1;
B = 2;

ue_data = randsrc(2,h_num,0:3);
trans_data = qammod(ue_data,4);

BER1 = zeros(1,length(snr));
BER2 = zeros(1,length(snr));
BER3 = zeros(1,length(snr));

for i = 1:1:length(snr)
    data_tilde1 = zeros(size(trans_data));
    data_tilde2 = zeros(size(trans_data));
    data_tilde3 = zeros(size(trans_data));
    snr_lin = 10^(snr(i)/10);
    for j = 1:1:h_num
        h1 = sqrt(1/2)*(randn(2,2)+1i*randn(2,2));
        h2 = sqrt(1/2)*(randn(2,4)+1i*randn(2,4));
        h3 = sqrt(1/2)*(randn(2,6)+1i*randn(2,6));
        
        [F1,G1] = uMMSE_design(h1,p0,B,1,snr_lin);
        [F2,G2] = uMMSE_design(h2,p0,B,1,snr_lin);
        [F3,G3] = uMMSE_design(h3,p0,B,1,snr_lin);
        
        afth_data1 = h1*F1*trans_data(:,j);
        afth_data2 = h2*F2*trans_data(:,j);
        afth_data3 = h3*F3*trans_data(:,j);
        
        noise = sqrt(1/snr_lin)*(randn(size(afth_data1))+1i*randn(size(afth_data1)));
        receive_data1 = afth_data1 + noise;
        receive_data2 = afth_data2 + noise;
        receive_data3 = afth_data3 + noise;
        
        s_tilde1 = G1*receive_data1;
        s_tilde2 = G2*receive_data2;
        s_tilde3 = G3*receive_data3;
        
        s_tilde_demod1 = qamdemod(s_tilde1,4);
        s_tilde_demod2 = qamdemod(s_tilde2,4);
        s_tilde_demod3 = qamdemod(s_tilde3,4);
        
        data_tilde1(:,j) = s_tilde_demod1;
        data_tilde2(:,j) = s_tilde_demod2;
        data_tilde3(:,j) = s_tilde_demod3;
        
    end
[num1,ber1] = biterr(data_tilde1(1,:),ue_data(1,:));
BER1(1,i) = ber1;
[num2,ber2] = biterr(data_tilde2(1,:),ue_data(1,:));
BER2(1,i) = ber2;
[num3,ber3] = biterr(data_tilde3(1,:),ue_data(1,:));
BER3(1,i) = ber3;
end
figure
semilogy(snr,BER1(1,:),'r-*')
hold on
semilogy(snr,BER2(1,:),'g-o')
semilogy(snr,BER3(1,:),'b-+')
legend('MT=2','MT=4','MT=6');
xlabel('SNR');ylabel('BER')
axis([0 20 10^(-4) 10^0])


