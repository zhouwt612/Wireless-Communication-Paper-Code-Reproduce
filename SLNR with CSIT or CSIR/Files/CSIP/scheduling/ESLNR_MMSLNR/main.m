% Simulation of scheduling schemes with RVQ and ESLNR precoding for limited feedback systems
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
S = 1;
Tc = 200;

K_all_Num = 1;
K_all = zeros(1,K_all_Num);
K_all = [50];

SNR_Num=11;
SNR=zeros(1,SNR_Num);
for i2=1:SNR_Num,
    SNR(1,i2)=-15+5*(i2-1);
end

BER1=zeros(K_all_Num,SNR_Num);
BER2=zeros(K_all_Num,SNR_Num);

BER5=zeros(K_all_Num,SNR_Num);
BER6=zeros(K_all_Num,SNR_Num);

capacity_sum1=zeros(K_all_Num,SNR_Num);
capacity_sum2=zeros(K_all_Num,SNR_Num);

capacity_sum5=zeros(K_all_Num,SNR_Num);
capacity_sum6=zeros(K_all_Num,SNR_Num);

round_num_sum5=zeros(K_all_Num,SNR_Num);
round_num_sum6=zeros(K_all_Num,SNR_Num);

cal_num_sum5=zeros(K_all_Num,SNR_Num);
cal_num_sum6=zeros(K_all_Num,SNR_Num);

FB_bits = 8;
L = 2^FB_bits;;

h=waitbar(0,'please wait...');
LoopNum=50;
for i1=1:K_all_Num,
    for i2=1:SNR_Num,
        sigma2=1/(10.^(SNR(1,i2)/10));

        u1=zeros(N,K);
        u2=zeros(N,K);
        
        u5=zeros(N,K);
        u6=zeros(N,K);
        
        count11=0;
        count22=0;

        count55=0;
        count66=0;

        capacity11=0;
        capacity22=0;

        capacity55=0;
        capacity66=0;
        
        round_num5=0;
        round_num6=0;
        
        cal_num5=0;
        cal_num6=0;
        
%         KK_NoScheduling=[1:K];
%         KK_CSIT=zeros(1,K);
%         KK_CSIP1=zeros(1,K);     % the inaccurate CSI may result in different scheduling user's set
%         KK_CSIP2=zeros(1,K);
        
        KK_NoScheduling=[1:K];
        KK_CSIT=[1:K];
        KK_CSIP1=[1:K];
        KK_CSIP2=[1:K];
        KK1=[1:K];
        KK1_old=[1:K];
        KK2=[1:K];
        
        candidate_num=K_all(1,i1);
        
        codebook=zeros(L,N,K_all(1,i1));
        for i4=1:K_all(1,i1),
            codebook(:,:,i4)=RVQ(N,1,L);
            % It is not advantageous for the SLNR scheme to use
            % more than one antenna to receive signals, so it is
            % reasonable to only consider the vector codebook.
        end
        
        for i3=1:LoopNum,
            s=zeros(K_all(1,i1)*S,2*Tc);
            x=zeros(K_all(1,i1)*S,Tc);
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
            
            [KK1] = RRS(K_all(1,i1),K,KK1_old);    % inaccurate CSI is applied to the scheduling scheme
            u1 = SLNR_modified(N,Mi,K,KK1,H_LFB,sigma2,FB_bits);
            [count_temp,capacity_temp] = coherence_receive(N,Mi,K_all(1,i1),K,KK1,Tc,s,x,H,sigma2,w,u1);
            count11=count11+count_temp;
            capacity11=capacity11+real(capacity_temp);
            
            [KK2] = maxH(N,Mi,S,K_all(1,i1),K,H);    % inaccurate CSI is applied to the scheduling scheme
            u2 = SLNR_modified(N,Mi,K,KK2,H_LFB,sigma2,FB_bits);
            [count_temp,capacity_temp] = coherence_receive(N,Mi,K_all(1,i1),K,KK2,Tc,s,x,H,sigma2,w,u2);
            count22=count22+count_temp;
            capacity22=capacity22+real(capacity_temp);

            [KK_CSIP1,round_num,cal_num] = MMSLNR_SingleAnt(N,K_all(1,i1),candidate_num,K,Mi,H_LFB,sigma2);    % inaccurate CSI is applied to the scheduling scheme
            round_num5=round_num5+round_num;
            cal_num5=cal_num5+cal_num;
            u5 = SLNR(N,Mi,K,KK_CSIP1,H_LFB,sigma2);
            [count_temp,capacity_temp] = coherence_receive(N,Mi,K_all(1,i1),K,KK_CSIP1,Tc,s,x,H,sigma2,w,u5);
            count55=count55+count_temp;
            capacity55=capacity55+real(capacity_temp);
            
            [KK_CSIP2,round_num,cal_num] = MMSLNR_SingleAnt_modified(N,K_all(1,i1),candidate_num,K,Mi,H_LFB,sigma2,FB_bits);
            round_num6=round_num6+round_num;
            cal_num6=cal_num6+cal_num;
            u6 = SLNR_modified(N,Mi,K,KK_CSIP2,H_LFB,sigma2,FB_bits);
            [count_temp,capacity_temp] = coherence_receive(N,Mi,K_all(1,i1),K,KK_CSIP2,Tc,s,x,H,sigma2,w,u6);
            count66=count66+count_temp;
            capacity66=capacity66+real(capacity_temp);
 
            waitbar(((i1-1)*SNR_Num*LoopNum+(i2-1)*LoopNum+i3)/(K_all_Num*SNR_Num*LoopNum),h);
        end

        BER1(i1,i2)=count11/(LoopNum*K*Tc);
        capacity_sum1(i1,i2)=capacity11/(LoopNum);

        BER2(i1,i2)=count22/(LoopNum*K*Tc);
        capacity_sum2(i1,i2)=capacity22/(LoopNum);

        BER5(i1,i2)=count55/(LoopNum*K*Tc);
        capacity_sum5(i1,i2)=capacity55/(LoopNum);
        
        round_num_sum5(i1,i2)=round_num5/(LoopNum);
        cal_num_sum5(i1,i2)=cal_num5/(LoopNum);

        BER6(i1,i2)=count66/(LoopNum*K*Tc);
        capacity_sum6(i1,i2)=capacity66/(LoopNum);
        
        round_num_sum6(i1,i2)=round_num6/(LoopNum);
        cal_num_sum6(i1,i2)=cal_num6/(LoopNum);

    end

    figure(i1*2-1)
    semilogy(SNR,BER1(i1,:),'g--');
    axis([-15 35 10^(-4) 10^(0)])
    xlabel('Average SNR (dB)');
    ylabel('Average BER');
    hold on;
    semilogy(SNR,BER2(i1,:),'b-*');

    hold on;
    semilogy(SNR,BER5(i1,:),'c-x');
    hold on;
    semilogy(SNR,BER6(i1,:),'m-v');
    legend('Round Robin','Max SNR','DMMSLNR','EMMSLNR');
    title('BER vs SNR');

    figure(i1*2)
    plot(SNR,capacity_sum1(i1,:),'g--');
    xlabel('Average SNR (dB)');
    ylabel('Sum capacity (bps/Hz)');
    hold on;
    plot(SNR,capacity_sum2(i1,:),'b-*');

    hold on;
    plot(SNR,capacity_sum5(i1,:),'c-x');
    hold on;
    plot(SNR,capacity_sum6(i1,:),'m-v');
    legend('Round Robin','Max SNR','DMMSLNR','EMMSLNR');
    title('Sum capacity vs SNR');
end