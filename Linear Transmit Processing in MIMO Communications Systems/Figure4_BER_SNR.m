clc
clear
close

sym_num = 100;
h_num = 10000;
Etr = 2;
% B: data streams  M: transmit antennas  N: receive antennas
B = 2;
M = 2;
N = 2;

% Data generation
sym_orig_data = randsrc(B,sym_num*h_num,0:3);
sym_mod = pskmod(sym_orig_data,4);

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
    sym_sto_txmf = zeros(B,sym_num*h_num);
    sym_sto_txzf = zeros(B,sym_num*h_num);
    sym_sto_txwf = zeros(B,sym_num*h_num);
    sym_sto_txcmmse1 = zeros(B,sym_num*h_num);
    sym_sto_txcmmse2 = zeros(B,sym_num*h_num);
    
    snr = 10^(SNRdB(i)/10);
    Rn = (1/snr)*eye(B);
    
    for j = 1:1:h_num
        H = sqrt(1/2)*(randn(N,M)+1i*randn(N,M));
        sym_trans = sym_mod(:,((j-1)*sym_num+1):(j*sym_num));
        
        % Transmit Filters
        % MF ZF WF
        Pmf = TxMF(H,G_orig,Rs,Etr);
        [Pzf, beta_txzf] = TxZF(H,G_orig,Rs,Etr);
        Pwf = TxWF(H,G_orig,Rs,Rn,Etr);
        
        % Receive signals
        sym_receive_txmf = H*Pmf*sym_trans + sqrt(1/snr)*(randn(2,sym_num)+1i*randn(2,sym_num));
        sym_receive_txzf = H*Pzf*sym_trans + sqrt(1/snr)*(randn(2,sym_num)+1i*randn(2,sym_num));
        sym_receive_txwf = H*Pwf*sym_trans + sqrt(1/snr)*(randn(2,sym_num)+1i*randn(2,sym_num));

        sym_GHP_txmf = G_orig*sym_receive_txmf;
        sym_GHP_txzf = G_orig*sym_receive_txzf;
        sym_GHP_txwf = G_orig*sym_receive_txwf;
        
        % Demodulation
        sym_demod_txmf = pskdemod(sym_GHP_txmf,4);
        sym_demod_txzf = pskdemod(sym_GHP_txzf,4);
        sym_demod_txwf = pskdemod(sym_GHP_txwf,4);
        
        % Save data
        sym_sto_txmf(:,((j-1)*sym_num+1):(j*sym_num)) = sym_demod_txmf;
        sym_sto_txzf(:,((j-1)*sym_num+1):(j*sym_num)) = sym_demod_txzf;
        sym_sto_txwf(:,((j-1)*sym_num+1):(j*sym_num)) = sym_demod_txwf;
        
        % CMMSE precoder
        Pcmmse9 = TxCMMSE(H,G_orig,Rs,2,1/9);
        sym_receive_txcmmse1 = H*Pcmmse9*sym_trans + sqrt(1/snr)*(randn(2,sym_num)+1i*randn(2,sym_num));
        sym_GHP_txcmmse1 = G_orig*sym_receive_txcmmse1;
        sym_demod_txcmmse1 = pskdemod(sym_GHP_txcmmse1,4);
        sym_sto_txcmmse1(:,((j-1)*sym_num+1):(j*sym_num)) = sym_demod_txcmmse1;
        
        Pcmmse205 = TxCMMSE(H,G_orig,Rs,2,1/20.5);
        sym_receive_txcmmse2 = H*Pcmmse205*sym_trans + sqrt(1/snr)*(randn(2,sym_num)+1i*randn(2,sym_num));
        sym_GHP_txcmmse2 = G_orig*sym_receive_txcmmse2;
        sym_demod_txcmmse2 = pskdemod(sym_GHP_txcmmse2,4);
        sym_sto_txcmmse2(:,((j-1)*sym_num+1):(j*sym_num)) = sym_demod_txcmmse2;
    end
    
    [~,bertxmf] = biterr(sym_sto_txmf,sym_orig_data);
    [~,bertxzf] = biterr(sym_sto_txzf,sym_orig_data);
    [~,bertxwf] = biterr(sym_sto_txwf,sym_orig_data);
    [~,bertxcmmse1] = biterr(sym_sto_txcmmse1,sym_orig_data);
    [~,bertxcmmse2] = biterr(sym_sto_txcmmse2,sym_orig_data);
    BERmf(2,i) = bertxmf;
    BERzf(2,i) = bertxzf;
    BERwf(2,i) = bertxwf;
    BERmf(1,i) = bertxcmmse1;
    BERzf(2,i) = bertxcmmse2;
end

semilogy(SNRdB,BERmf(2,:),'r-^');
hold on;
semilogy(SNRdB,BERzf(2,:),'g-v');
semilogy(SNRdB,BERwf(2,:),'b-d')
semilogy(SNRdB,BERmf(1,:),'r-+')
semilogy(SNRdB,BERzf(1,:),'g-<')
xlabel('Es/No in dB');ylabel('BER');
legend('TxMF','TxZF','TxWF','TxMMSE,Xi^{(-1)}=9dB','TxMMSE,Xi^{(-1)}=20.5dB')