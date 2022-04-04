% Simulation of precoding schemes with RVQ for limited feedback systems
% Author: Xiaotian Wang @ May-2009
% Revised: @ Jan-20-2010
% Email: wxtlovewlj520@126.com

% Reference:

% Parameters:
% N: amount of transmitting antennas
% K_all: amount of users waiting for being scheduled
% K: amount of the scheduled users
% Mi: amount of receiving antennas of each user
% Tc: the coherence time, amount of symbol periods per user
% H: the channel matrix

% Archives depended on:
% None
clc
clear;
% the parameters
N = 8;
K = 4;
Mi = 2;
Tc = 1;

K_all_Num = 1;
K_all = zeros(1,K_all_Num);
K_all = [K];

SNR_Num=11;
SNR=zeros(1,SNR_Num);
for i2=1:SNR_Num
    SNR(1,i2)=-15+5*(i2-1);
end

BER1=zeros(K_all_Num,SNR_Num);

BER3=zeros(K_all_Num,SNR_Num);
BER4=zeros(K_all_Num,SNR_Num);
BER5=zeros(K_all_Num,SNR_Num);
BER6=zeros(K_all_Num,SNR_Num);

capacity_sum1=zeros(K_all_Num,SNR_Num);

capacity_sum3=zeros(K_all_Num,SNR_Num);
capacity_sum4=zeros(K_all_Num,SNR_Num);
capacity_sum5=zeros(K_all_Num,SNR_Num);
capacity_sum6=zeros(K_all_Num,SNR_Num);

number=zeros(K_all_Num,SNR_Num);        % observe the percent of the case that more than one receiving antenna is used

FB_bits = 5;
L = 2^FB_bits;;

h=waitbar(0,'please wait...');
LoopNum=1000;
for i1=1:K_all_Num,
    for i2=1:SNR_Num,
        sigma2=1/(10.^(SNR(1,i2)/10));

        u1=zeros(N,K);

        u3=zeros(N,K);
        u4=zeros(N,K);
        u5=zeros(N,K);
        u6=zeros(N,K);
        
        count1=0;

        count3=0;
        count4=0;
        count5=0;
        count6=0;

        capacity1=0;

        capacity3=0;
        capacity4=0;
        capacity5=0;
        capacity6=0;

        count_temp=0;
        capacity_temp=0;
        
%         KK_NoScheduling=[1:K];
%         KK_CSIT=zeros(1,K);
%         KK_CSIP1=zeros(1,K);     % the inaccurate CSI may result in different scheduling user's set
%         KK_CSIP2=zeros(1,K);
        
        KK_NoScheduling=[1:K];
        KK_CSIT=[1:K];
        KK_CSIP1=[1:K];
        KK_CSIP2=[1:K];
        
        candidate_num=K_all(1,i1);
        
        codebook=zeros(L,N,K_all(1,i1));
        for i4=1:K_all(1,i1),
            codebook(:,:,i4)=RVQ(N,1,L);
            % It is not advantageous for the SLNR scheme to use
            % more than one antenna to receive signals, so it is
            % reasonable to only consider the vector codebook.
        end
        
        for i3=1:LoopNum,
            s=zeros(K_all(1,i1),Tc);
            x=zeros(K_all(1,i1),Tc);
            [s,x]=source_BPSK(K_all(1,i1),1,Tc);        % Attention: the modulation method influences not only the receiver design, but also the computation of BER

            H_real=randn(K_all(1,i1)*Mi,N);
            H_imag=randn(K_all(1,i1)*Mi,N);
            H=(1/sqrt(2))*complex(H_real,H_imag);

            w_real=sqrt(sigma2/2)*randn(K_all(1,i1)*Mi,Tc);
            w_imag=sqrt(sigma2/2)*randn(K_all(1,i1)*Mi,Tc);
            w=complex(w_real,w_imag);

            FB_Info1=zeros(1,K_all(1,i1)*Mi);
            index1=zeros(1,K_all(1,i1)*Mi);
            % If more than one receiving antenna are used, all the
            % channel vectors will be quantized even though only one of
            % them may be used.

            % Attention: the quantized channel vectors are sorted
            % according to the quantization error from small to large.
            
            % quantize CDI
            for i4=1:K_all(1,i1),
                Hi=zeros(Mi,N);
                Hi=H((i4-1)*Mi+1:i4*Mi,:);
                for i5=1:Mi,
                    Hi_temp=zeros(1,N);
                    Hi_temp=Hi(i5,:);
                    Hi(i5,:)=Hi_temp/(sqrt(trace(Hi_temp'*Hi_temp)));
                end
                [FB_Info1(1,(i4-1)*Mi+1:i4*Mi), index1(1,(i4-1)*Mi+1:i4*Mi)]=feedback(N,Mi,L,Hi,codebook(:,:,i4));
            end

            H_LFB1=zeros(K_all(1,i1)*Mi,N);
            for i4=1:K_all(1,i1),
                for i5=1:Mi,
                    H_LFB1((i4-1)*Mi+i5,:)=codebook(FB_Info1(1,(i4-1)*Mi+i5),:,i4);
                end
            end
             
            % CSI-T
            u1 = SLNR_SingleAnt(N,Mi,K,KK_NoScheduling,H,sigma2);
            [count_temp,capacity_temp] = coherence_receive_SingleAnt(N,Mi,K_all(1,i1),K,KK_NoScheduling,Tc,s,x,H,sigma2,w,u1,ones(1,K*Mi));
            count1=count1+count_temp;
            capacity1=capacity1+real(capacity_temp);
          
            % DSLNR
            u3 = SLNR_SingleAnt(N,Mi,K,KK_NoScheduling,H_LFB1,sigma2);
            [count_temp,capacity_temp] = coherence_receive_SingleAnt(N,Mi,K_all(1,i1),K,KK_NoScheduling,Tc,s,x,H,sigma2,w,u3,index1);
            count3=count3+count_temp;
            capacity3=capacity3+real(capacity_temp);

            % RMMSE
            u4 = EMMSE(N,Mi,K,KK_NoScheduling,H_LFB1,sigma2,FB_bits);
            [count_temp,capacity_temp] = coherence_receive_SingleAnt(N,Mi,K_all(1,i1),K,KK_NoScheduling,Tc,s,x,H,sigma2,w,u4,index1);
            count4=count4+count_temp;
            capacity4=capacity4+real(capacity_temp);

            % ESLNR_SingleAnt
            u5 = SLNR_SingleAnt_modified(N,Mi,K,KK_NoScheduling,H_LFB1,sigma2,FB_bits);
            [count_temp,capacity_temp] = coherence_receive_SingleAnt(N,Mi,K_all(1,i1),K,KK_NoScheduling,Tc,s,x,H,sigma2,w,u5,index1);
            count5=count5+count_temp;
            capacity5=capacity5+real(capacity_temp);
            
            % ESLNR_MultiAnt
            [u6, ant_num] = SLNR_MultiAnt_modified(N,Mi,K,KK_NoScheduling,H_LFB1,sigma2,FB_bits);
            for i4=1:K,
                if ant_num(1,i4)~=1,
                    number(i1,i2)=number(i1,i2)+1;
                    break;
                end
            end
            
            [count_temp,capacity_temp] = coherence_receive_MultiAnt(N,Mi,K_all(1,i1),K,KK_NoScheduling,Tc,s,x,H,sigma2,w,u6,index1,ant_num);
            count6=count6+count_temp;
            capacity6=capacity6+real(capacity_temp);
            
            waitbar(((i1-1)*SNR_Num*LoopNum+(i2-1)*LoopNum+i3)/(K_all_Num*SNR_Num*LoopNum),h);
        end

        number(i1,i2)=number(i1,i2)/LoopNum;
        
        BER1(i1,i2)=count1/(LoopNum*K*Tc);
        capacity_sum1(i1,i2)=capacity1/(LoopNum);

        BER3(i1,i2)=count3/(LoopNum*K*Tc);
        capacity_sum3(i1,i2)=capacity3/(LoopNum);
        
        BER4(i1,i2)=count4/(LoopNum*K*Tc);
        capacity_sum4(i1,i2)=capacity4/(LoopNum);

        BER5(i1,i2)=count5/(LoopNum*K*Tc);
        capacity_sum5(i1,i2)=capacity5/(LoopNum);

        BER6(i1,i2)=count6/(LoopNum*K*Tc);
        capacity_sum6(i1,i2)=capacity6/(LoopNum);

    end

    figure(i1*2-1)
    semilogy(SNR,BER1(i1,:),'g--');
    axis([-15 35 10^(-5) 10^(0)])
    xlabel('Average SNR (dB)');
    ylabel('Average BER');
    hold on;
    semilogy(SNR,BER3(i1,:),'r-^');
    hold on;
    semilogy(SNR,BER4(i1,:),'k-o');
    hold on;
    semilogy(SNR,BER5(i1,:),'c-x');
    hold on;
    semilogy(SNR,BER6(i1,:),'b-v');
    legend('CSIT','DSLNR','EMMSE','ESLNR(Single)','ESLNR(MMSLNR)');
    title('BER vs SNR');

    figure(i1*2)
%     plot(SNR,capacity_sum1(i1,:),'g--');
    xlabel('Average SNR (dB)');
    ylabel('Sum capacity (bps/Hz)');
    hold on;
    plot(SNR,capacity_sum3(i1,:),'r-^','Linewidth',2);
    hold on;
    plot(SNR,capacity_sum4(i1,:),'k-o','Linewidth',2);
    hold on;
    plot(SNR,capacity_sum5(i1,:),'c-x','Linewidth',2);
    hold on;
    plot(SNR,capacity_sum6(i1,:),'b-v','Linewidth',2);
%     legend('CSIT','DSLNR','EMMSE','ESLNR(Single)','ESLNR(MMSLNR)');
    legend('DSLNR','EMMSE','ESLNR(Single)','ESLNR(MMSLNR)');
    title('Sum capacity vs SNR');
end