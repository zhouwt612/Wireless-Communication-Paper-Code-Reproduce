clc
clear
% close

h_num = 5000;
snr_dB = 0:2:20;
% UE1: video  UE2: audio
ue1_data = randsrc(1,h_num,0:3);
ue2_data = randsrc(1,h_num,0:3);
ue_data = [ue1_data;ue2_data];
D = diag([0.7597 0.2403]);
% D = diag([0.5 0.5]);
p0 = 1;

ue1_mod_data = pskmod(ue1_data,4,pi/4);
ue2_mod_data = pskmod(ue2_data,4,pi/4);
trans_data = [ue1_mod_data;ue2_mod_data];

BER = zeros(2,length(snr_dB));

for i = 1:1:length(snr_dB)
    data_tilde = zeros(size(trans_data));
    snr_lin = 10^(snr_dB(i)/10);
    for j = 1:1:h_num
        h = sqrt(1/2)*(randn(3)+1i*randn(3));
        Rnn = 1/snr_lin*eye(size(h));
        [v,lambda] = eig(h'/Rnn*h);
        [d,ind] = sort(diag(lambda),'descend');
        LAMBDA = lambda(ind,ind);
        V = v(:,ind);
        V3x2 = V(:,[1 2]);
        LAMBDA2x2 = LAMBDA([1 2],[1 2]);
        
        % Precoder
        gamma = p0/trace(D/LAMBDA2x2); % error
        phi_f = gamma^(1/2)*D^(1/2)/LAMBDA2x2^(1/2);
        F = V3x2*phi_f;
        
        % Decoder
        phi_g = gamma^(1/2)*D^(1/2)/LAMBDA2x2^(1/2)/(gamma*D+eye(size(D))); % error
        G = phi_g*V3x2'*h'/Rnn;
        
        % Through channel
        afth_data = h*F*trans_data(:,j);
        noise = sqrt(1/snr_lin)*(randn(size(afth_data))+1i*randn(size(afth_data)));
        receive_data = afth_data + noise;
        s_tilde = G*receive_data;
        s_tilde_demod = pskdemod(s_tilde,4,pi/4);
        data_tilde(:,j) = s_tilde_demod;
    end
[num1,ber1] = biterr(data_tilde(1,:),ue_data(1,:));
[num2,ber2] = biterr(data_tilde(2,:),ue_data(2,:));
BER(1,i) = ber1;
BER(2,i) = ber2;
end
figure
semilogy(snr_dB,BER(1,:),'r-*')
hold on
semilogy(snr_dB,BER(2,:),'g-o')
legend('UE1','UE2')
xlabel('SNR')
ylabel('BER')