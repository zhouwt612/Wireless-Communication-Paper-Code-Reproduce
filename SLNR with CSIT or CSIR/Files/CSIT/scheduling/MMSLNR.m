function [scheduled,u,round_num,cal_num] = MMSLNR(N,Mi,S,Kall,candidate_num,K,H,sigma2)

scheduled=zeros(1,K);
u=zeros(N,K);
round_num=0;
cal_num=0;

channel_gain=zeros(1,Kall);
HHH=zeros(Mi,N);
for i1=1:Kall,
    HHH=H((i1-1)*Mi+1:i1*Mi,:);
    channel_gain(1,i1)=trace(HHH'*HHH);
end

[ch_gain,index]=another_sort(channel_gain);

scheduled_temp=zeros(1,K);

% the first round
scheduled_temp=index(1,1:K);

H_now=zeros(K*Mi,N);
for i1=1:K,
    H_now((i1-1)*Mi+1:i1*Mi,:)=H((scheduled_temp(1,i1)-1)*Mi+1:scheduled_temp(1,i1)*Mi,:);
end

SLNR_temp=zeros(1,K);
u_temp=zeros(N,K);

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
    
    [eig_vec,eig_val]=eig(C,D);
    for i3=1:N,
        temp(1,i3)=eig_val(i3,i3);
    end
    [max_val,max_pos]=max(temp);
    SLNR_temp(1,i1)=real(max_val);
    v=eig_vec(:,max_pos);
    u_temp(:,i1)=v/(sqrt((v')*v));

end

u=u_temp;
cal_num=cal_num+K;
round_num=round_num+1;

% the following rounds
over_all=1;
temp1=zeros(1,K);
temp2=zeros(1,K);
while over_all~=0,
    [min_val, min_pos]=min(SLNR_temp);
    changed=0;
    for i1=1:candidate_num,
        over=0;
        % check whether this user has been scheduled
        for i2=1:K,
            if index(1,i1)==scheduled_temp(1,i2),
                over=1;
                break;
            end
        end

        if over==0,
            temp1=scheduled_temp;      % scheduled_temp is used to store the scheduling result of last round
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
                
                [eig_vec,eig_val]=eig(C,D);
                for i3=1:N,
                    temp(1,i3)=eig_val(i3,i3);
                end
                [max_val,max_pos]=max(temp);
                temp2(1,i2)=real(max_val);
                v=eig_vec(:,max_pos);
                u_temp(:,i2)=v/(sqrt((v')*v));

            end
            cal_num=cal_num+K;

            if min(temp2)>min_val,
                scheduled_temp=temp1;
                SLNR_temp=temp2;
                u=u_temp;
                changed=1;
                round_num=round_num+1;
                break;
            end
        end
    end
    
    if changed==0,
        over_all=0;
        round_num=round_num+1;
    end
end

scheduled=scheduled_temp;