% Simulation of Scheduling schemes for SLNR(SJNR) precoding
% Author: Wang Xiaotian @ Dec-2008
% Revised @ Jan-20-2010 
% Email: wxtlovewlj520@126.com

% Reference:

% Parameters:
% N: amount of transmitting antennas
% Kall: amount of users
% K: amount of the scheduled users
% Mi: amount of receiving antennas of each user
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

Kall_Num = 1;
Kall=zeros(1,Kall_Num);
Kall=[20];

SNR_Num=11;
SNR=zeros(1,SNR_Num);
for i1=1:SNR_Num,
    SNR(1,i1)=-15+5*(i1-1);
end

BER1=zeros(Kall_Num,SNR_Num);
BER2=zeros(Kall_Num,SNR_Num);
BER3=zeros(Kall_Num,SNR_Num);

capacity_sum1=zeros(Kall_Num,SNR_Num);
capacity_sum2=zeros(Kall_Num,SNR_Num);
capacity_sum3=zeros(Kall_Num,SNR_Num);

round_num_sum3=zeros(Kall_Num,SNR_Num);

cal_num_sum3=zeros(Kall_Num,SNR_Num);

LoopNum=100;

scheduled_sum=zeros(SNR_Num,LoopNum,Kall_Num);

h = waitbar(0,'please wait...');
for i1=1:Kall_Num,
    for i2=1:SNR_Num,
        sigma2=1/(10.^(SNR(1,i2)/10));

        count1=0;
        count2=0;
        count3=0;

        capacity1=0;
        capacity2=0;
        capacity3=0;
        
        round_num3=0;
        
        cal_num3=0;
        
        scheduled1_old=zeros(1,K);
        scheduled1=[1:K];
        scheduled2=zeros(1,K);
        scheduled3=zeros(1,K);

        for i3=1:LoopNum,
                s=zeros(Kall(1,i1)*S,2*Tc);
                x=zeros(Kall(1,i1)*S,Tc);
                [s,x]=source_QPSK(Kall(1,i1),S,Tc);      % attention: if the modulation method changes, the link between information bits and signals also changes, including something in function SLNR

                H_real=randn(Kall(1,i1)*Mi,N);
                H_imag=randn(Kall(1,i1)*Mi,N);
                H=(1/sqrt(2))*complex(H_real,H_imag);

                w_real=sqrt(sigma2/2)*randn(Kall(1,i1)*Mi,Tc);
                w_imag=sqrt(sigma2/2)*randn(Kall(1,i1)*Mi,Tc);
                w=complex(w_real,w_imag);
                
                % Round Robin
                [scheduled1] = RoundRobin(Kall(1,i1),K,scheduled1_old);
                scheduled1_old=scheduled1;

                H_now = zeros(K*Mi,N);
                for i4=1:K,
                    H_now((i4-1)*Mi+1:i4*Mi,:)=H((scheduled1(1,i4)-1)*Mi+1:scheduled1(1,i4)*Mi,:);
                end

                u1 = SLNR(N,Mi,S,K,H_now,sigma2);

                [count_temp,capacity_temp] = coherence_receive(N,Mi,S,K,scheduled1,Tc,s,x,H,sigma2,w,u1,ones(1,K));
                count1=count1+count_temp;
                capacity1=capacity1+real(capacity_temp);

                % maxH
                [scheduled2]=maxH(N,Mi,S,Kall(1,i1),K,H);
                H_now = zeros(K*Mi,N);
                for i4=1:K,
                    H_now((i4-1)*Mi+1:i4*Mi,:)=H((scheduled2(1,i4)-1)*Mi+1:scheduled2(1,i4)*Mi,:);
                end

                u2 = SLNR(N,Mi,S,K,H_now,sigma2);

                [count_temp,capacity_temp] = coherence_receive(N,Mi,S,K,scheduled2,Tc,s,x,H,sigma2,w,u2,ones(1,K));
                count2=count2+count_temp;
                capacity2=capacity2+real(capacity_temp);

                % MMSLNR
                candidate_num=Kall(1,i1);

                [scheduled4,u3,round_num_temp,cal_num_temp] = MMSLNR(N,Mi,S,Kall(1,i1),candidate_num,K,H,sigma2);

                cal_num3=cal_num3+cal_num_temp;
                round_num3=round_num3+round_num_temp;

                [count_temp,capacity_temp] = coherence_receive(N,Mi,S,K,scheduled4,Tc,s,x,H,sigma2,w,u3,ones(1,K));
                count3=count3+count_temp;
                capacity3=capacity3+real(capacity_temp);

            waitbar(((i1-1)*SNR_Num*LoopNum+(i2-1)*LoopNum+i3)/(Kall_Num*SNR_Num*LoopNum),h)
        end
        
        BER1(i1,i2)=count1/(LoopNum*K*2*Tc);
        capacity_sum1(i1,i2)=capacity1/(LoopNum);

        BER2(i1,i2)=count2/(LoopNum*K*2*Tc);
        capacity_sum2(i1,i2)=capacity2/(LoopNum);

        BER3(i1,i2)=count3/(LoopNum*K*2*Tc);
        capacity_sum3(i1,i2)=capacity3/(LoopNum);

        round_num_sum3(i1,i2)=round_num3/(LoopNum);
        cal_num_sum3(i1,i2)=cal_num3/(LoopNum);
        
    end

    figure(4*i1-3)
    semilogy(SNR,BER1(i1,:),'k--');
    axis([-15 25 10^(-6) 10^(0)])
    xlabel('average SNR in dB');
    ylabel('average BER');
    hold on;
    semilogy(SNR,BER2(i1,:),'g-+');
    hold on;
    semilogy(SNR,BER3(i1,:),'r-o');

    legend('RoundRobin','maxH','MMSLNR');
    title('BER vs SNR');
    hold off;

    figure(4*i1-2)
    plot(SNR,capacity_sum1(i1,:),'k--');
    xlabel('average SNR in dB');
    ylabel('sum capacity');
    hold on;
    plot(SNR,capacity_sum2(i1,:),'g-+');
    hold on;
    plot(SNR,capacity_sum3(i1,:),'r-o');

    legend('RoundRobin','maxH','MMSLNR');
    title('sum capacity vs SNR');
    hold off;

    figure(4*i1-1)
    plot(SNR,round_num_sum3(i1,:),'r-o');
    xlabel('average SNR in dB');
    ylabel('amount of the replacement round');
   
    legend('MMSLNR');
    title('amount of the replacement round vs SNR');
    hold off;

    figure(4*i1)
    plot(SNR,cal_num_sum3(i1,:),'r-o');
    xlabel('average SNR in dB');
    ylabel('amount of the calculation of SLNR');
    
    legend('MMSLNR');
    title('computational complexity vs SNR');
    hold off;

end