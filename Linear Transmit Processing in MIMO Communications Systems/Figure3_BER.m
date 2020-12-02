clc
clear
% close

sym_num = 100;
h_num = 10000;
Etr = 2;
% B: data streams  M: transmit antennas  N: receive antennas
B = 2;
M = 2;
N = 2;

% Data generation
sym_orig_data = randsrc(B,sym_num*h_num,0:3);
sym_mod = qammod(sym_orig_data,4);

SNRdB = -10:5:30;

Rs = eye(B);
P_orig = eye(B);
G_orig = eye(B);

BERmf = zeros(B,length(SNRdB));
BERzf = zeros(B,length(SNRdB));
BERwf = zeros(B,length(SNRdB));

for i = 1:1:length(SNRdB)
    disp(['SNR: ' num2str(SNRdB(i))])
    
    % Assign the storge
    sym_sto_rxmf = zeros(B,sym_num*h_num);
    sym_sto_rxzf = zeros(B,sym_num*h_num);
    sym_sto_rxwf = zeros(B,sym_num*h_num);
    sym_sto_txmf = zeros(B,sym_num*h_num);
    sym_sto_txzf = zeros(B,sym_num*h_num);
    sym_sto_txwf = zeros(B,sym_num*h_num);
    
    snr = 10^(SNRdB(i)/10);
    Rn = (1/snr)*eye(B);
    
    for j = 1:1:h_num
        H = sqrt(1/2)*(randn(N,M)+1i*randn(N,M));
        sym_trans = sym_mod(:,((j-1)*sym_num+1):(j*sym_num));
        
        % Receive Filters
        sym_HP = H*P_orig*sym_trans;
%         sym_noise = awgn(sym_HP,SNRdB(i));
        sym_noise = sym_HP + sqrt(1/(2*snr))*(randn(size(sym_HP))+1i*randn(size(sym_HP)));
        
        % MF ZF WF
        Gmf = RxMF(P_orig,H,Rs,Rn,1);
        Gzf = RxZF(P_orig,H,Rn);
        Gwf = RxWF(P_orig,H,Rs,Rn);
        
        % Receive signals
        sym_GHF_rxmf = Gmf*sym_noise;
        sym_GHF_rxzf = Gzf*sym_noise;
        sym_GHF_rxwf = Gwf*sym_noise;
        
        % Demodulation
        sym_demod_rxmf = qamdemod(sym_GHF_rxmf,4);
        sym_demod_rxzf = qamdemod(sym_GHF_rxzf,4);
        sym_demod_rxwf = qamdemod(sym_GHF_rxwf,4);
        
        % Save data
        sym_sto_rxmf(:,((j-1)*sym_num+1):(j*sym_num)) = sym_demod_rxmf;
        sym_sto_rxzf(:,((j-1)*sym_num+1):(j*sym_num)) = sym_demod_rxzf;
        sym_sto_rxwf(:,((j-1)*sym_num+1):(j*sym_num)) = sym_demod_rxwf;
        
        % Transmit Filters
        % MF ZF WF
        Pmf = TxMF(H,G_orig,Rs,Etr);
        [Pzf, beta_txzf] = TxZF(H,G_orig,Rs,Etr);
        Pwf = TxWF(H,G_orig,Rs,Rn,Etr);
        
        % Receive signals
        sym_receive_txmf = H*Pmf*sym_trans + sqrt(1/(2*snr))*(randn(size(sym_HP))+1i*randn(size(sym_HP)));
        sym_receive_txzf = H*Pzf*sym_trans + sqrt(1/(2*snr))*(randn(size(sym_HP))+1i*randn(size(sym_HP)));
        sym_receive_txwf = H*Pwf*sym_trans + sqrt(1/(2*snr))*(randn(size(sym_HP))+1i*randn(size(sym_HP)));
%         sym_receive_txmf = awgn(H*Pmf*sym_trans,SNRdB(i));
%         sym_receive_txzf = awgn(H*Pzf*sym_trans,SNRdB(i));
%         sym_receive_txwf = awgn(H*Pwf*sym_trans,SNRdB(i));
        
        sym_GHP_txmf = G_orig*sym_receive_txmf;
        sym_GHP_txzf = G_orig*sym_receive_txzf;
        sym_GHP_txwf = G_orig*sym_receive_txwf;
        
        % Demodulation
        sym_demod_txmf = qamdemod(sym_GHP_txmf,4);
        sym_demod_txzf = qamdemod(sym_GHP_txzf,4);
        sym_demod_txwf = qamdemod(sym_GHP_txwf,4);
        
        % Save data
        sym_sto_txmf(:,((j-1)*sym_num+1):(j*sym_num)) = sym_demod_txmf;
        sym_sto_txzf(:,((j-1)*sym_num+1):(j*sym_num)) = sym_demod_txzf;
        sym_sto_txwf(:,((j-1)*sym_num+1):(j*sym_num)) = sym_demod_txwf;
        
    end
    [~,berrxmf] = biterr(sym_sto_rxmf,sym_orig_data);
    [~,berrxzf] = biterr(sym_sto_rxzf,sym_orig_data);
    [~,berrxwf] = biterr(sym_sto_rxwf,sym_orig_data);
    BERmf(1,i) = berrxmf;
    BERzf(1,i) = berrxzf;
    BERwf(1,i) = berrxwf;
    
    [~,bertxmf] = biterr(sym_sto_txmf,sym_orig_data);
    [~,bertxzf] = biterr(sym_sto_txzf,sym_orig_data);
    [~,bertxwf] = biterr(sym_sto_txwf,sym_orig_data);
    BERmf(2,i) = bertxmf;
    BERzf(2,i) = bertxzf;
    BERwf(2,i) = bertxwf;
    
end
figure
semilogy(SNRdB,BERmf(1,:),'r-+');
hold on;
semilogy(SNRdB,BERzf(1,:),'g-o');
semilogy(SNRdB,BERwf(1,:),'b-*')
semilogy(SNRdB,BERmf(2,:),'r-^');
semilogy(SNRdB,BERzf(2,:),'g-v');
semilogy(SNRdB,BERwf(2,:),'b-d')
xlabel('Es/No in dB');ylabel('BER');
legend('RxMF','RxZF','RxWF','TxMF','TxZF','TxWF')