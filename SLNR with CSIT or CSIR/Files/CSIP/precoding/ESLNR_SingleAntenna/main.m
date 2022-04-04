% Simulation of precoding schemes by using RVQ for limited feedback systems
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

clear;
% the parameters
N = 4;
K = 4;
Mi = 1;
Tc = 200;

K_all_Num = 1;
K_all = zeros(1,K_all_Num);
K_all = K;

SNR_Num=11;
SNR=zeros(1,SNR_Num);
for i2=1:SNR_Num,
    SNR(1,i2)=-15+5*(i2-1);
end

BER1=zeros(K_all_Num,SNR_Num);
BER2=zeros(K_all_Num,SNR_Num);
BER3=zeros(K_all_Num,SNR_Num);
BER4=zeros(K_all_Num,SNR_Num);

capacity_sum1=zeros(K_all_Num,SNR_Num);
capacity_sum2=zeros(K_all_Num,SNR_Num);
capacity_sum3=zeros(K_all_Num,SNR_Num);
capacity_sum4=zeros(K_all_Num,SNR_Num);

FB_bits = 8;
L = 2^FB_bits;

h=waitbar(0,'please wait...');
LoopNum=50;
for i1=1:K_all_Num,
    for i2=1:SNR_Num,
        sigma2=1/(10.^(SNR(1,i2)/10));

        u1=zeros(N,K);
        u2=zeros(N,K);
        u3=zeros(N,K);
        u4=zeros(N,K);

        count1=0;
        count2=0;
        count3=0;
        count4=0;

        capacity1=0;
        capacity2=0;
        capacity3=0;
        capacity4=0;
        
        count_temp=0;
        capacity_temp=0;
        
        KK=zeros(1,K);
        KK=[1:K];

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

            FB_Info=zeros(1,K_all(1,i1)*Mi);
            index=zeros(1,K_all(1,i1)*Mi);
            % If more than one receiving antenna are used, all the
            % channel vectors will be quantized even though only one of
            % them may be used.

            % Attention: the quantized channel vectors are sorted
            % according to the quantization error from small to large.

            for i4=1:K_all(1,i1),
                [FB_Info(1,(i4-1)*Mi+1:i4*Mi), index(1,(i4-1)*Mi+1:i4*Mi)]=feedback(N,Mi,L,H((i4-1)*Mi+1:i4*Mi,:),codebook(:,:,i4));
            end

            H_LFB=zeros(K_all(1,i1)*Mi,N);
            for i4=1:K_all(1,i1),
                for i5=1:Mi,
                    H_LFB((i4-1)*Mi+i5,:)=codebook(FB_Info(1,(i4-1)*Mi+i5),:,i4);
                end
            end

            u1 = SLNR_SingleAnt(N,Mi,K,KK,H,sigma2);
            [count_temp,capacity_temp] = coherence_receive_SingleAnt(N,Mi,K_all(1,i1),K,KK,Tc,s,x,H,sigma2,w,u1);
            count1=count1+count_temp;
            capacity1=capacity1+real(capacity_temp);

            u2 = SLNR_SingleAnt(N,Mi,K,KK,H_LFB,sigma2);
            [count_temp,capacity_temp] = coherence_receive_SingleAnt(N,Mi,K_all(1,i1),K,KK,Tc,s,x,H,sigma2,w,u2);
            count2=count2+count_temp;
            capacity2=capacity2+real(capacity_temp);

            u3 = SLNR_SingleAnt_modified(N,Mi,K,KK,H_LFB,sigma2,FB_bits);
            [count_temp,capacity_temp] = coherence_receive_SingleAnt(N,Mi,K_all(1,i1),K,KK,Tc,s,x,H,sigma2,w,u3);
            count3=count3+count_temp;
            capacity3=capacity3+real(capacity_temp);

            u4 = EMMSE(N,Mi,K,KK,H_LFB,sigma2,FB_bits);
            [count_temp,capacity_temp] = coherence_receive_SingleAnt(N,Mi,K_all(1,i1),K,KK,Tc,s,x,H,sigma2,w,u4);
            count4=count4+count_temp;
            capacity4=capacity4+real(capacity_temp);

            waitbar(((i1-1)*SNR_Num*LoopNum+(i2-1)*LoopNum+i3)/(K_all_Num*SNR_Num*LoopNum),h);
        end

        BER1(i1,i2)=count1/(LoopNum*K*Tc);
        capacity_sum1(i1,i2)=capacity1/(LoopNum);

        BER2(i1,i2)=count2/(LoopNum*K*Tc);
        capacity_sum2(i1,i2)=capacity2/(LoopNum);

        BER3(i1,i2)=count3/(LoopNum*K*Tc);
        capacity_sum3(i1,i2)=capacity3/(LoopNum);

        BER4(i1,i2)=count4/(LoopNum*K*Tc);
        capacity_sum4(i1,i2)=capacity4/(LoopNum);
    end
    
    figure(i1*2-1)
    semilogy(SNR,BER1(i1,:),'g--');
    axis([-15 35 10^(-5) 10^(0)])
    xlabel('average SNR in dB');
    ylabel('average BER');
    hold on;
    semilogy(SNR,BER2(i1,:),'b-*');
    hold on;
    semilogy(SNR,BER3(i1,:),'r-^');
    hold on;
    semilogy(SNR,BER4(i1,:),'k-x');
    legend('CSIT','DSLNR','ESLNR','EMMSE');
    title('BER vs SNR');

    figure(i1*2)
    plot(SNR,capacity_sum1(i1,:),'g--');
    xlabel('average SNR in dB');
    ylabel('sum capacity');
    hold on;
    plot(SNR,capacity_sum2(i1,:),'b-*');
    hold on;
    plot(SNR,capacity_sum3(i1,:),'r-^');
    hold on;
    plot(SNR,capacity_sum4(i1,:),'k-x');
    legend('CSIT','DSLNR','ESLNR','EMMSE');
    title('sum capacity vs SNR');
end
