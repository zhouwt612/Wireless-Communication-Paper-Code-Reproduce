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

SNRdB = -10:5:30;

Rs = eye(B);
P_orig = eye(B);
G_orig = eye(B);

MSERX = zeros(3,length(SNRdB));
MSETX = zeros(3,length(SNRdB));
for i = 1:1:length(SNRdB)
    disp(['SNR: ' num2str(SNRdB(i))])
    snr = 10^(SNRdB(i)/10);
    Rn = (1/snr)*eye(B);
    
    MSE_rxmf_data = zeros(1,h_num);
    MSE_rxzf_data = zeros(1,h_num);
    MSE_rxwf_data = zeros(1,h_num);
    MSE_txmf_data = zeros(1,h_num);
    MSE_txzf_data = zeros(1,h_num);
    MSE_txwf_data = zeros(1,h_num);
    for j = 1:1:h_num
        % Channel generation
        H = sqrt(1/2)*(randn(N,M)+1i*randn(N,M));
        
        % MSEs of the receive filters
        MSErxmf = REMFmse(H,P_orig,Rs,Rn);
        MSErxzf = REZFmse(H,P_orig,Rs,Rn);
        MSErxwf = REWFmse(H,P_orig,Rs,Rn);
        % MSEs of the transmit filters
        MSEtxmf = TXMFmse(H,G_orig,Rs,Rn,Etr);
        MSEtxzf = TXZFmse(H,G_orig,Rs,Rn,Etr);
        MSEtxwf = TXWFmse(H,G_orig,Rs,Rn,Etr);
        
        % Save data
        MSE_rxmf_data(1,j) = MSErxmf;
        MSE_rxzf_data(1,j) = MSErxzf;
        MSE_rxwf_data(1,j) = MSErxwf;
        MSE_txmf_data(1,j) = MSEtxmf;
        MSE_txzf_data(1,j) = MSEtxzf;
        MSE_txwf_data(1,j) = MSEtxwf;
    end
    % Save MSEs
    MSERX(1,i) = real(mean(MSE_rxmf_data));
    MSERX(2,i) = real(mean(MSE_rxzf_data));
    MSERX(3,i) = real(mean(MSE_rxwf_data));
    
    MSETX(1,i) = real(mean(MSE_txmf_data));
    MSERX(2,i) = real(mean(MSE_txzf_data));
    MSERX(3,i) = real(mean(MSE_txwf_data));
end

figure
semilogy(SNRdB,MSERX(1,:),'r-+');
hold on
semilogy(SNRdB,MSERX(2,:),'g-o');
semilogy(SNRdB,MSERX(3,:),'b-*');
semilogy(SNRdB,MSETX(1,:),'r-^');
semilogy(SNRdB,MSERX(2,:),'g-v');
semilogy(SNRdB,MSERX(3,:),'b-d');
xlabel('Es/No in dB');ylabel('BER');
legend('RxMF','RxZF','RxWF','TxMF','TxZF','TxWF')
