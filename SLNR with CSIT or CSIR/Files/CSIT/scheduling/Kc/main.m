% Simulation of the MMSLNR scheme with different Kc
% Author: Wang Xiaotian @ Sep-2009
% Revised @ Jan-20-2010 
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
K_all=[50];

SNR_Num=11;
SNR=zeros(1,SNR_Num);
for i1=1:SNR_Num,
    SNR(1,i1)=-15+5*(i1-1);
end

BER1=zeros(K_all_Num,SNR_Num);
BER2=zeros(K_all_Num,SNR_Num);
BER3=zeros(K_all_Num,SNR_Num);
BER4=zeros(K_all_Num,SNR_Num);

capacity_avg1=zeros(K_all_Num,SNR_Num);
capacity_avg2=zeros(K_all_Num,SNR_Num);
capacity_avg3=zeros(K_all_Num,SNR_Num);
capacity_avg4=zeros(K_all_Num,SNR_Num);

round_num_avg1=zeros(K_all_Num,SNR_Num);
round_num_avg2=zeros(K_all_Num,SNR_Num);
round_num_avg3=zeros(K_all_Num,SNR_Num);
round_num_avg4=zeros(K_all_Num,SNR_Num);

eig_num_avg1=zeros(K_all_Num,SNR_Num);
eig_num_avg2=zeros(K_all_Num,SNR_Num);
eig_num_avg3=zeros(K_all_Num,SNR_Num);
eig_num_avg4=zeros(K_all_Num,SNR_Num);

LoopNum=50;
h = waitbar(0,'please wait...');
for i1=1:K_all_Num,
    for i2=1:SNR_Num,
        sigma2=1/(10.^(SNR(1,i2)/10));
        
        round_num1=0;
        eig_num1=0;
        count1=0;
        capacity1=0;
        
        round_num2=0;
        eig_num2=0;
        count2=0;
        capacity2=0;
        
        round_num3=0;
        eig_num3=0;
        count3=0;
        capacity3=0;
        
        round_num4=0;
        eig_num4=0;
        count4=0;
        capacity4=0;
        
        round_num_temp=0;
        eig_num_temp=0;
        count_temp=0;
        capacity_temp=0;
        
        K_preselect=0;
        
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

            if(K_all(1,i1)>=2*K),
                K_preselect=2*K;
            else
                K_preselect=0;
            end
            if(K_preselect~=0)
                [round_num_temp, eig_num_temp, KK1, u1] = MMSLNR_preselect_SLNR_SingleStream(N,K_all(1,i1),K_preselect,K,Mi,H,sigma2);
                round_num1=round_num1+round_num_temp;
                eig_num1=eig_num1+eig_num_temp;

                [count_temp,capacity_temp] = coherence_receive_SingleStream(N,K_all(1,i1),K,Mi,KK1,Tc,s,x,H,sigma2,w,u1);
                count1=count1+count_temp;
                capacity1=capacity1+real(capacity_temp);
            end
            
            if(K_all(1,i1)>=3*K),
                K_preselect=3*K;
            else
                K_preselect=0;
            end
            if(K_preselect~=0)
                [round_num_temp, eig_num_temp, KK1, u1] = MMSLNR_preselect_SLNR_SingleStream(N,K_all(1,i1),K_preselect,K,Mi,H,sigma2);
                round_num2=round_num2+round_num_temp;
                eig_num2=eig_num2+eig_num_temp;

                [count_temp,capacity_temp] = coherence_receive_SingleStream(N,K_all(1,i1),K,Mi,KK1,Tc,s,x,H,sigma2,w,u1);
                count2=count2+count_temp;
                capacity2=capacity2+real(capacity_temp);
            end
            
            if(K_all(1,i1)>=4*K),
                K_preselect=4*K;
            else
                K_preselect=0;
            end
            if(K_preselect~=0)
                [round_num_temp, eig_num_temp, KK1, u1] = MMSLNR_preselect_SLNR_SingleStream(N,K_all(1,i1),K_preselect,K,Mi,H,sigma2);
                round_num3=round_num3+round_num_temp;
                eig_num3=eig_num3+eig_num_temp;

                [count_temp,capacity_temp] = coherence_receive_SingleStream(N,K_all(1,i1),K,Mi,KK1,Tc,s,x,H,sigma2,w,u1);
                count3=count3+count_temp;
                capacity3=capacity3+real(capacity_temp);
            end
            
            if(K_all(1,i1)>=5*K),
                K_preselect=5*K;
            else
                K_preselect=0;
            end
            if(K_preselect~=0)
                [round_num_temp, eig_num_temp, KK1, u1] = MMSLNR_preselect_SLNR_SingleStream(N,K_all(1,i1),K_preselect,K,Mi,H,sigma2);
                round_num4=round_num4+round_num_temp;
                eig_num4=eig_num4+eig_num_temp;

                [count_temp,capacity_temp] = coherence_receive_SingleStream(N,K_all(1,i1),K,Mi,KK1,Tc,s,x,H,sigma2,w,u1);
                count4=count4+count_temp;
                capacity4=capacity4+real(capacity_temp);
            end

            waitbar(((i1-1)*SNR_Num*LoopNum+(i2-1)*LoopNum+i3)/(K_all_Num*SNR_Num*LoopNum),h)
        end
        
        round_num_avg1(i1,i2)=round_num1/(LoopNum);
        eig_num_avg1(i1,i2)=eig_num1/(LoopNum);
        
        BER1(i1,i2)=count1/(LoopNum*K*2*Tc);
        capacity_avg1(i1,i2)=capacity1/(LoopNum);
        
        round_num_avg2(i1,i2)=round_num2/(LoopNum);
        eig_num_avg2(i1,i2)=eig_num2/(LoopNum);
        
        BER2(i1,i2)=count2/(LoopNum*K*2*Tc);
        capacity_avg2(i1,i2)=capacity2/(LoopNum);
        
        round_num_avg3(i1,i2)=round_num3/(LoopNum);
        eig_num_avg3(i1,i2)=eig_num3/(LoopNum);
        
        BER3(i1,i2)=count3/(LoopNum*K*2*Tc);
        capacity_avg3(i1,i2)=capacity3/(LoopNum);
        
        round_num_avg4(i1,i2)=round_num4/(LoopNum);
        eig_num_avg4(i1,i2)=eig_num4/(LoopNum);
        
        BER4(i1,i2)=count4/(LoopNum*K*2*Tc);
        capacity_avg4(i1,i2)=capacity4/(LoopNum);

    end
    
    figure(4*i1-3)
    semilogy(SNR,BER1(i1,:),'g--');
    axis([-15 35 10^(-6) 10^(0)])
    xlabel('Average SNR (dB)');
    ylabel('Average BER');
    hold on;
    semilogy(SNR,BER2(i1,:),'k-*');
    hold on;
    semilogy(SNR,BER3(i1,:),'r-o');
    hold on;
    semilogy(SNR,BER4(i1,:),'b-+');
    title('BER vs SNR');
    legend('Kc=2K','Kc=3K','Kc=4K','Kc=5K');
    hold off;

    figure(4*i1-2)
    plot(SNR,capacity_avg1(i1,:),'g--');
    xlabel('Average SNR (dB)');
    ylabel('Sum capacity (bps/Hz)');
    hold on;
    plot(SNR,capacity_avg2(i1,:),'k-*');
    hold on;
    plot(SNR,capacity_avg3(i1,:),'r-o');
    hold on;
    plot(SNR,capacity_avg4(i1,:),'b-+');
    title('Sum capacity vs SNR');
    legend('Kc=2K','Kc=3K','Kc=4K','Kc=5K');
    hold off;

    figure(4*i1-1)
    plot(SNR,round_num_avg1(i1,:),'g--');
    xlabel('Average SNR (dB)');
    ylabel('Iteration round');
    hold on;
    plot(SNR,round_num_avg2(i1,:),'k-*');
    hold on;
    plot(SNR,round_num_avg3(i1,:),'r-o');
    hold on;
    plot(SNR,round_num_avg4(i1,:),'b-+');
    title('Convergence');
    legend('Kc=2K','Kc=3K','Kc=4K','Kc=5K');
    hold off;

    figure(4*i1)
    plot(SNR,eig_num_avg1(i1,:),'g--');
    xlabel('Average SNR (dB)');
    ylabel('Amount of the solution of the eigenvalue problem');
    hold on;
    plot(SNR,eig_num_avg2(i1,:),'k-*');
    hold on;
    plot(SNR,eig_num_avg3(i1,:),'r-o');
    hold on;
    plot(SNR,eig_num_avg4(i1,:),'b-+');
    title('Computational complexity');
    legend('Kc=2K','Kc=3K','Kc=4K','Kc=5K');
    hold off;
    
end