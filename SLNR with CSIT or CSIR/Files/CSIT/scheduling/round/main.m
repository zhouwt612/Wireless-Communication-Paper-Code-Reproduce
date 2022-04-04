% Simulation of the MMSLNR scheme with different Kall
% Author: Wang Xiaotian @ Sep-2009
% Revised @ Oct-2-2009 
% Email: wxtlovewlj520@126.com

% Reference:

% Parameters:
% N: amount of transmitting antennas
% K_all: amount of users
% K: amount of the scheduled users
% Mi: amount of receiving antennas of each user
% S: amount of data streams for one user
% Tc: the coherence time, amount of symbol periods per user 
% H: the channel matrix

% Archives depended on:
% None

clear;
% the parameters
N = 4;
K = 4;
Mi = 1;
S = 1;
Tc = 200;

K_all_Num = 1;
K_all=zeros(1,K_all_Num);
K_all=[20];

SNR_Num=11;
SNR=zeros(1,SNR_Num);
for i1=1:SNR_Num,
    SNR(1,i1)=-15+5*(i1-1);
end

BER1=zeros(K_all_Num,SNR_Num,K+1);

capacity_avg1=zeros(K_all_Num,SNR_Num,K+1);

LoopNum=50;
h = waitbar(0,'please wait...');
for i1=1:K_all_Num,
    for i2=1:SNR_Num,
        sigma2=1/(10.^(SNR(1,i2)/10));

        count1=zeros(1,K+1);
        capacity1=zeros(1,K+1);
        eig_num1=zeros(1,K+1);

        count_temp=0;
        capacity_temp=0;
        eig_num_temp=0;
        
        for i3=1:LoopNum,       
            s=zeros(K_all(1,i1)*S,2*Tc);
            x=zeros(K_all(1,i1)*S,Tc);
            [s,x] = source_QPSK(K_all(1,i1),S,Tc);

            H_real=randn(K_all(1,i1)*Mi,N);
            H_imag=randn(K_all(1,i1)*Mi,N);
            H=(1/sqrt(2))*complex(H_real,H_imag);

            w_real=sqrt(sigma2/2)*randn(K_all(1,i1)*Mi,Tc);
            w_imag=sqrt(sigma2/2)*randn(K_all(1,i1)*Mi,Tc);
            w=complex(w_real,w_imag);

            KK_old=zeros(1,K);
            u_old=zeros(N,K);
            SLNR_temp_old=zeros(1,K);
            
            over_old=0;
            for i4=1:(K+1),
                [KK1, u1, SLNR_temp, eig_num_temp, over_new] = MMSLNR_SLNR_SingleStream_OneRound(N,K_all(1,i1),K,Mi,H,sigma2,KK_old,u_old,SLNR_temp_old,over_old);
                eig_num1(1,i4)=eig_num1(1,i4)+eig_num_temp;
                
                [count_temp,capacity_temp] = coherence_receive_SingleStream(N,K_all(1,i1),K,Mi,KK1,Tc,s,x,H,sigma2,w,u1);
                count1(1,i4)=count1(1,i4)+count_temp;
                capacity1(1,i4)=capacity1(1,i4)+real(capacity_temp);
                
                KK_old=KK1;
                u_old=u1;
                SLNR_temp_old=SLNR_temp;
                over_old=over_new;
                
                waitbar(((i1-1)*SNR_Num*LoopNum*(K+1)+(i2-1)*LoopNum*(K+1)+i3*(K+1)+i4)/(K_all_Num*SNR_Num*LoopNum*(K+1)),h)
            end

        end
        
        for i5=1:(K+1),
            BER1(i1,i2,i5)=count1(1,i5)/(LoopNum*K*2*Tc);
            capacity_avg1(i1,i2,i5)=capacity1(1,i5)/(LoopNum);
            eig_num_avg1(i1,i2,i5)=sum(eig_num1(1,1:i5))/(LoopNum);
        end

    end
    
    figure(i1*3-2)
    semilogy(SNR,BER1(i1,:,1),'g--');
    axis([-15 35 10^(-6) 10^(0)])
    xlabel('Average SNR (dB)');
    ylabel('Average BER');
    hold on;
    semilogy(SNR,BER1(i1,:,2),'k-*');
    hold on;
    semilogy(SNR,BER1(i1,:,3),'r-o');
    hold on;
    semilogy(SNR,BER1(i1,:,4),'b-+');
    hold on;
    semilogy(SNR,BER1(i1,:,5),'c-v');
    title('BER vs SNR');
    legend('Q=1','Q=2','Q=3','Q=4','Q=5');
    hold off;

    figure(i1*3-1)
    plot(SNR,capacity_avg1(i1,:,1),'g--');
    xlabel('Average SNR (dB)');
    ylabel('Sum capacity (bps/Hz)');
    hold on;
    plot(SNR,capacity_avg1(i1,:,2),'k-*');
    hold on;
    plot(SNR,capacity_avg1(i1,:,3),'r-o');
    hold on;
    plot(SNR,capacity_avg1(i1,:,4),'b-+');
    hold on;
    plot(SNR,capacity_avg1(i1,:,5),'c-v');
    title('Sum capacity vs SNR');
    legend('Q=1','Q=2','Q=3','Q=4','Q=5');
    hold off;
    
    figure(i1*3)
    plot(SNR,eig_num_avg1(i1,:,1),'g--');
    xlabel('Average SNR (dB)');
    ylabel('');
    hold on;
    plot(SNR,eig_num_avg1(i1,:,2),'k-*');
    hold on;
    plot(SNR,eig_num_avg1(i1,:,3),'r-o');
    hold on;
    plot(SNR,eig_num_avg1(i1,:,4),'b-+');
    hold on;
    plot(SNR,eig_num_avg1(i1,:,5),'c-v');
    title('Computational complexity');
    legend('Q=1','Q=2','Q=3','Q=4','Q=5');
    hold off;
    
end