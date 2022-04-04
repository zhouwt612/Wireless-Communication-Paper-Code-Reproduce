% Simulation of precoding schemes for multi-stream per user with power distribution
% Author: Xiaotian Wang @ Nov-2009
% Revised: @ Jan-20-2010
% Email: wxtlovewlj520@126.com

% Reference:

% Parameters:
% N: amount of transmitting antennas
% K_all: amount of users waiting for being scheduled
% K: amount of the scheduled users
% Mi: amount of receiving antennas of each user
% S: amount of data streams for each user
% Tc: the coherence time, amount of symbol periods per user 
% H: the channel matrix

% Archives depended on:
% None

clear;
% the parameters
N = 9;
K = 4;
Mi = 2;
S = 2;
Tc = 200;

K_all_Num = 1;
K_all = zeros(1,K_all_Num);
K_all = [K];      % No scheduling

SNR_Num=11;
SNR=zeros(1,SNR_Num);
for i2=1:SNR_Num,
    SNR(1,i2)=-15+5*(i2-1);
end

BER3=zeros(K_all_Num,SNR_Num);
BER4=zeros(K_all_Num,SNR_Num);

capacity_sum3=zeros(K_all_Num,SNR_Num);
capacity_sum4=zeros(K_all_Num,SNR_Num);

h=waitbar(0,'please wait...');
LoopNum=100;
for i1=1:K_all_Num,
    for i2=1:SNR_Num,
        sigma2=1/(10.^(SNR(1,i2)/10));
        
        u3=zeros(N,K);
        u4=zeros(N,K);
        
        count3=0;
        count4=0;

        capacity3=0;
        capacity4=0;

        count_temp=0;
        capacity_temp=0;
        
        KK_NoScheduling=[1:K];
        
        candidate_num=K_all(1,i1);
        
        for i3=1:LoopNum,
            s=zeros(K_all(1,i1)*S,2*Tc);
            x=zeros(K_all(1,i1)*S,Tc);
            [s,x]=source_QPSK(K_all(1,i1),S,Tc);        % Attention: the modulation method influences not only the receiver design, but also the computation of BER

            H_real=randn(K_all(1,i1)*Mi,N);
            H_imag=randn(K_all(1,i1)*Mi,N);
            H=(1/sqrt(2))*complex(H_real,H_imag);

            w_real=sqrt(sigma2/2)*randn(K_all(1,i1)*Mi,Tc);
            w_imag=sqrt(sigma2/2)*randn(K_all(1,i1)*Mi,Tc);
            w=complex(w_real,w_imag);
            
            % BD
            [u3, r3, gama3] = BD_PD(N,Mi,S,K,H,sigma2);
            p3 = PowerDistribution_Reciprocal_MultiStream1(K,S,gama3);
            [count_temp,capacity_temp] = svd_receive(N,Mi,S,K,KK_NoScheduling,Tc,s,x,H,sigma2,w,u3,p3,r3);
            count3=count3+count_temp;
            capacity3=capacity3+real(capacity_temp);

            % SLNR orthogonal
            [u4, SLNR_val] = SLNR_orthogonal_PD(N,Mi,S,K,H,sigma2);
            p4 = PowerDistribution_Reciprocal_MultiStream1(K,S,SLNR_val);
            [count_temp,capacity_temp] = coherence_receive(N,Mi,S,K,KK_NoScheduling,Tc,s,x,H,sigma2,w,u4,p4);
            count4=count4+count_temp;
            capacity4=capacity4+real(capacity_temp);
            
            waitbar(((i1-1)*SNR_Num*LoopNum+(i2-1)*LoopNum+i3)/(K_all_Num*SNR_Num*LoopNum),h);
        end
        
        % If the modulation method change, the calculation of BER may change,too.
        BER3(i1,i2)=count3/(LoopNum*K*S*2*Tc);
        capacity_sum3(i1,i2)=capacity3/(LoopNum);
        
        BER4(i1,i2)=count4/(LoopNum*K*S*2*Tc);
        capacity_sum4(i1,i2)=capacity4/(LoopNum);

    end

    figure(i1*2-1)
    semilogy(SNR,BER3(i1,:),'b-*');
    axis([-15 35 10^(-6) 10^(0)])
    xlabel('Average SNR (dB)');
    ylabel('Average BER');
    hold on;
    semilogy(SNR,BER4(i1,:),'k--');
    legend('BD','SLNR');
    title('BER vs SNR');

    figure(i1*2)
    plot(SNR,capacity_sum3(i1,:),'b-*');
    xlabel('Average SNR (dB)');
    ylabel('Sum capacity (bps/Hz)');
    hold on;
    plot(SNR,capacity_sum4(i1,:),'k--');
    legend('BD','SLNR');
    title('Sum capacity vs SNR');
end