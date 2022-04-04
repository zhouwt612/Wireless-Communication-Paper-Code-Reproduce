% Simulation of limited feedback by using modified quantization version
% Author: Xiaotian Wang @ May-2009
% Revised: @ May-19-2009
% Email: wxtlovewlj520@126.com

% Reference:

% Parameters:
% N: amount of transmitting antennas
% K_all: amount of users waiting for being scheduled
% K: amount of the scheduled users
% Mi: amount of receiving antennas of each user
% Tc: amount of bits transmitted in a coherence time
% H: the channel matrix

% Archives depended on:
% None

clear;
% the parameters
N = 4;
K = 4;
Mi = 2;
Tc = 200;

K_all_Num = 1;
K_all = zeros(1,K_all_Num);
K_all = [1*N];

SNR_Num=11;
SNR=zeros(1,SNR_Num);
for i2=1:SNR_Num,
    SNR(1,i2)=-15+5*(i2-1);
end

BER1=zeros(K_all_Num,SNR_Num);
BER2=zeros(K_all_Num,SNR_Num);
BER3=zeros(K_all_Num,SNR_Num);
BER4=zeros(K_all_Num,SNR_Num);
BER5=zeros(K_all_Num,SNR_Num);
BER6=zeros(K_all_Num,SNR_Num);

capacity_sum1=zeros(K_all_Num,SNR_Num);
capacity_sum2=zeros(K_all_Num,SNR_Num);
capacity_sum3=zeros(K_all_Num,SNR_Num);
capacity_sum4=zeros(K_all_Num,SNR_Num);
capacity_sum5=zeros(K_all_Num,SNR_Num);
capacity_sum6=zeros(K_all_Num,SNR_Num);

number1=zeros(K_all_Num,SNR_Num);        % observe the percent of the case that more than one receiving antenna is used
number2=zeros(K_all_Num,SNR_Num);

FB_bits = 6;
L = 2^FB_bits;;

h=waitbar(0,'please wait...');
LoopNum=50000;
for i1=1:K_all_Num,
    codebook=zeros(L,N,K_all(1,i1));
    for i4=1:K_all(1,i1),
        codebook(:,:,i4)=RVQ(N,1,L);
        % It is not advantageous for the SLNR scheme to use
        % more than one antenna to receive signals, so it is
        % reasonable to only consider the vector codebook.
    end
    
    for i2=1:SNR_Num,
        sigma2=1/(10.^(SNR(1,i2)/10));

        u1=zeros(N,K);
        u2=zeros(N,K);
        u3=zeros(N,K);
        u4=zeros(N,K);
        u5=zeros(N,K);
        u6=zeros(N,K);
        
        count1=0;
        count2=0;
        count3=0;
        count4=0;
        count5=0;
        count6=0;

        capacity1=0;
        capacity2=0;
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

            % These two parameters are for the quantization of CSI
%             FB_Info2=zeros(1,K_all(1,i1)*Mi);
%             index2=zeros(1,K_all(1,i1)*Mi);
            
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
            
%             % quantize CSI
%             for i4=1:K_all(1,i1),
%                 Hi=zeros(Mi,N);
%                 Hi=H((i4-1)*Mi+1:i4*Mi,:);
%                 [FB_Info2(1,(i4-1)*Mi+1:i4*Mi), index2(1,(i4-1)*Mi+1:i4*Mi)]=feedback(N,Mi,L,Hi,codebook(:,:,i4));
%             end
% 
%             H_LFB2=zeros(K_all(1,i1)*Mi,N);
%             for i4=1:K_all(1,i1),
%                 for i5=1:Mi,
%                     H_LFB2((i4-1)*Mi+i5,:)=codebook(FB_Info2(1,(i4-1)*Mi+i5),:,i4);
%                 end
%             end

%             % CSI-T
%             u1 = SLNR_SingleAnt(N,Mi,K,KK_NoScheduling,H,sigma2);
%             [count_temp,capacity_temp] = coherence_receive_SingleAnt(N,Mi,K_all(1,i1),K,KK_NoScheduling,Tc,s,x,H,sigma2,w,u1,ones(1,K*Mi));
%             count1=count1+count_temp;
%             capacity1=capacity1+real(capacity_temp);

            % directly quantize CSI, then MMSE
            u2 = CMMSE(N,Mi,K,KK_NoScheduling,H_LFB1,sigma2);
            [count_temp,capacity_temp] = coherence_receive_SingleAnt(N,Mi,K_all(1,i1),K,KK_NoScheduling,Tc,s,x,H,sigma2,w,u2,index1);
            count2=count2+count_temp;
            capacity2=capacity2+real(capacity_temp);

            % RMMSE
            u3 = RMMSE(N,Mi,K,KK_NoScheduling,H_LFB1,sigma2,FB_bits);
            [count_temp,capacity_temp] = coherence_receive_SingleAnt(N,Mi,K_all(1,i1),K,KK_NoScheduling,Tc,s,x,H,sigma2,w,u3,index1);
            count3=count3+count_temp;
            capacity3=capacity3+real(capacity_temp);
            
            % CSLNR
            [u4, ant_num] = SLNR_MultiAnt(N,Mi,K,KK_NoScheduling,H_LFB1,sigma2);
            for i4=1:K,
                if ant_num(1,i4)~=1,
                    number1(i1,i2)=number1(i1,i2)+1;
                    break;
                end
            end
            
            [count_temp,capacity_temp] = coherence_receive_MultiAnt(N,Mi,K_all(1,i1),K,KK_NoScheduling,Tc,s,x,H,sigma2,w,u4,index1,ant_num);
            count4=count4+count_temp;
            capacity4=capacity4+real(capacity_temp);

%             % ESLNR_SingleAnt
%             u5 = SLNR_SingleAnt_modified(N,Mi,K,KK_NoScheduling,H_LFB1,sigma2,FB_bits);
%             [count_temp,capacity_temp] = coherence_receive_SingleAnt(N,Mi,K_all(1,i1),K,KK_NoScheduling,Tc,s,x,H,sigma2,w,u5,index1);
%             count5=count5+count_temp;
%             capacity5=capacity5+real(capacity_temp);
            
            % ESLNR_MultiAnt
            [u5, ant_num] = SLNR_MultiAnt_modified(N,Mi,K,KK_NoScheduling,H_LFB1,sigma2,FB_bits);
            for i4=1:K,
                if ant_num(1,i4)~=1,
                    number2(i1,i2)=number2(i1,i2)+1;
                    break;
                end
            end
            
            [count_temp,capacity_temp] = coherence_receive_MultiAnt(N,Mi,K_all(1,i1),K,KK_NoScheduling,Tc,s,x,H,sigma2,w,u5,index1,ant_num);
            count5=count5+count_temp;
            capacity5=capacity5+real(capacity_temp);
            
            waitbar(((i1-1)*SNR_Num*LoopNum+(i2-1)*LoopNum+i3)/(K_all_Num*SNR_Num*LoopNum),h);
        end

        number1(i1,i2)=number1(i1,i2)/LoopNum;
        number2(i1,i2)=number2(i1,i2)/LoopNum;
        
        BER1(i1,i2)=count1/(LoopNum*K*Tc);
        capacity_sum1(i1,i2)=capacity1/(LoopNum);

        BER2(i1,i2)=count2/(LoopNum*K*Tc);
        capacity_sum2(i1,i2)=capacity2/(LoopNum);

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
    semilogy(SNR,BER1(i1,:),'k--');
    axis([-15 35 10^(-5) 10^(0)])
    xlabel('average SNR in dB');
    ylabel('average BER');
    hold on;
    semilogy(SNR,BER2(i1,:),'g-*');
    hold on;
    semilogy(SNR,BER3(i1,:),'b-^');
    hold on;
    semilogy(SNR,BER4(i1,:),'r-o');
    hold on;
    semilogy(SNR,BER5(i1,:),'m-x');
%     hold on;
%     semilogy(SNR,BER6(i1,:),'m-v');
    legend('CSIT','CMMSE','RMMSE','CSLNR','ESLNR');
    title('BER vs SNR');

    figure(i1*2)
    plot(SNR,capacity_sum1(i1,:),'k--');
    axis([-15 35 0 15])
    xlabel('average SNR in dB');
    ylabel('sum capacity');
    hold on;
    plot(SNR,capacity_sum2(i1,:),'g-*');
    hold on;
    plot(SNR,capacity_sum3(i1,:),'b-^');
    hold on;
    plot(SNR,capacity_sum4(i1,:),'r-o');
    hold on;
    plot(SNR,capacity_sum5(i1,:),'m-x');
%     hold on;
%     plot(SNR,capacity_sum6(i1,:),'m-v');
    legend('CSIT','CMMSE','RMMSE','CSLNR','ESLNR');
    title('sum capacity vs SNR');
end