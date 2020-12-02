clc
clear
close

h_num = 10000;
snr = 0:2:20;
p0 = 1;

ue_data = randsrc(2,h_num,0:3);
trans_data = qammod(ue_data,4);

BER = zeros(1,length(snr));

for i = 1:1:length(snr)
    data_tilde = zeros(size(trans_data));
    snr_lin = 10^(snr(i)/10);
    for j = 1:1:h_num
        h = sqrt(1/2)*(randn(2,2)+1i*randn(2,2));
        R = 1/snr_lin*eye(2);
        [v,lambda] = eig(h'/R*h);
        [d,ind] = sort(diag(lambda),'descend');
        LAMBDA = lambda(ind,ind);
        V = v(:,ind);
        V2x2 = V(:,[1 2]);
        LAMBDA2x2 = LAMBDA([1 2],[1 2]);
        % precoder
        w = eye(2);
        mu12 = trace(LAMBDA2x2^(-1/2)*w^(1/2))/(p0+trace(LAMBDA2x2^(-1)));
        phi_f = mu12^(-1)*LAMBDA2x2^(-1/2)*w^(1/2) - LAMBDA2x2^(-1);
                phi_f((phi_f<0)) = 0;
%         if phi_f(4) < 0
%             phi_f((phi_f<0)) = 0;
%             phi_f = phi_f/phi_f(1);
%         end

        phi_f = phi_f^(1/2);
        F = V2x2 * phi_f;
        
        
        %decoder
        phi_g = mu12*LAMBDA2x2^(-1/2)*w^(-1/2)-(mu12^2)*LAMBDA2x2^(-1)*w^(-1);
        phi_g((phi_g<0)) = 0;
        phi_g = phi_g^(1/2)*LAMBDA2x2^(-1/2);
        G = phi_g*V2x2'*h'/R;
        
        % through channel        
        afth_data = h*F*trans_data(:,j);
        noise = sqrt(1/snr_lin)*(randn(size(afth_data))+1i*randn(size(afth_data)));
        receive_data = afth_data + noise;
        s_tilde = G*receive_data;
        s_tilde_demod = qamdemod(s_tilde,4);
        data_tilde(:,j) = s_tilde_demod;
    end
[num,ber] = biterr(data_tilde(1,:),ue_data(1,:));
BER(1,i) = ber;
end

semilogy(snr,BER(1,:),'r-*')
legend('MT=2');
xlabel('SNR');ylabel('BER')
axis([0 20 10^(-4) 10^0])
        
