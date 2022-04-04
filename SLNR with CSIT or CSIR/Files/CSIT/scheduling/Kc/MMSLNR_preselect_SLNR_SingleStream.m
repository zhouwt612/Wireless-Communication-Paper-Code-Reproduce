function [round_num, eig_num, KK, u] = MMSLNR_preselect_SLNR_SingleStream(N,K_all,K_preselect,K,Mi,H,sigma2)

% Scheduling
KK=zeros(1,K);
u=zeros(N,K);       % Single data stream
round_num=0;
eig_num=0;

channel_gain=zeros(1,K_all);
HHH=zeros(Mi,N);
for i1=1:K_all,
    HHH=H((i1-1)*Mi+1:i1*Mi,:);
    channel_gain(1,i1)=trace(HHH'*HHH);
end

[ch_gain,index]=another_sort(channel_gain);

KK_temp=zeros(1,K);

% the first round
KK_temp=index(1,1:K);
H_now=zeros(K*Mi,N);
for i1=1:K,
    H_now((i1-1)*Mi+1:i1*Mi,:)=H((KK_temp(1,i1)-1)*Mi+1:KK_temp(1,i1)*Mi,:);
end

SLNR_temp=zeros(1,K);
for i1=1:K,
    Hi=zeros(Mi,N);
    Hk=zeros(Mi,N);

    Hi=H_now((i1-1)*Mi+1:i1*Mi,:);
    C=Hi'*Hi;

    D=Mi*sigma2*eye(N);
    for i2=1:K,
        if i2~=i1,
            Hk=H_now((i2-1)*Mi+1:i2*Mi,:);
            D=D+Hk'*Hk;
        end
    end

    [eig_vec, eig_val]=eig(C,D);
    temp1=zeros(1,N);
    for i2=1:N,
        temp1(1,i2)=eig_val(i2,i2);
    end
    [max_val,max_pos]=max(temp1);
    v=eig_vec(:,max_pos);
    u(:,i1)=v/(sqrt((v')*v));
    
    SLNR_temp(1,i1)=real(max_val);
end

round_num=round_num+1;
eig_num=eig_num+K;

% the following rounds
over=1;
temp1=zeros(1,K);
temp2=zeros(1,K);
u_temp=zeros(N,K);
while over~=0,
    [min_val, min_pos]=min(SLNR_temp);
    changed=0;
    for i1=1:K_preselect,
        scheduled=0;
        % check whether this user has been scheduled
        for i2=1:K,
            if index(1,i1)==KK_temp(1,i2);
                scheduled=1;
                break;
            end
        end

        if scheduled==0,
            temp1=KK_temp;      % KK_temp is used to store the scheduling result of last round
            temp1(1,min_pos)=index(1,i1);       % try to replace，注意index中的元素才代表真正的用户序号

            for i2=1:K,
                H_now((i2-1)*Mi+1:i2*Mi,:)=H((temp1(1,i2)-1)*Mi+1:temp1(1,i2)*Mi,:);
            end

            Hi=zeros(Mi,N);
            Hk=zeros(Mi,N);
            for i2=1:K,
                Hi=H_now((i2-1)*Mi+1:i2*Mi,:);
                C=Hi'*Hi;

                D=Mi*sigma2*eye(N);
                for i3=1:K,
                    if i3~=i2,
                        Hk=H_now((i3-1)*Mi+1:i3*Mi,:);
                        D=D+Hk'*Hk;
                    end
                end

                [eig_vec, eig_val]=eig(C,D);
                temp3=zeros(1,N);
                for i3=1:N,
                    temp3(1,i3)=eig_val(i3,i3);
                end
                [max_val, max_pos]=max(temp3);
                v=eig_vec(:,max_pos);
                u_temp(:,i2)=v/(sqrt((v')*v));

                temp2(1,i2)=real(max_val);
            end
            
            eig_num=eig_num+K;

            if min(temp2)>min_val,
                KK_temp=temp1;
                SLNR_temp=temp2;
                u=u_temp;
                changed=1;
                break;
            end
        end
    end
    
    round_num=round_num+1;
    
    if changed==0,
        over=0;
    end
end
KK=KK_temp;