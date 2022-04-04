function [KK_new,round_num,cal_num] = MMSLNR_SingleAnt_modified(MT,K_all,candidate_num,K,R,H,sigma2,FB_bits)

KK_new=zeros(1,K);
round_num=0;
cal_num=0;

err_UpperBound=2^(-FB_bits/(MT-1));

channel_gain=zeros(1,K_all);
HHH=zeros(1,MT);
for i1=1:K_all,
    HHH=H(i1,:);
    channel_gain(1,i1)=trace(HHH'*HHH);
end

[ch_gain,index]=another_sort(channel_gain);

KK_temp=zeros(1,K);

% the first round
KK_temp=index(1,1:K);
H_now=zeros(K,MT);
for i1=1:K,
    H_now(i1,:)=H((KK_temp(1,i1)-1)*R+1,:);
end

SLNR_temp=zeros(1,K);
for i1=1:K,
    Hi=zeros(1,MT);
    Hk=zeros(1,MT);

    Hi=H_now(i1,:);
    C=(MT-MT*err_UpperBound)*Hi'*Hi+err_UpperBound*MT*eye(MT);

    D=sigma2*eye(MT);
    for i2=1:K,
        if i2~=i1,
            Hk=H_now(i2,:);
            D=D+(MT-MT*err_UpperBound)*Hk'*Hk+err_UpperBound*MT*eye(MT);
        end
    end

    SLNR_temp(1,i1)=real(max(eig(inv(D)*C)));
end
cal_num=cal_num+K;

% the following rounds
over=1;
temp1=zeros(1,K);
temp2=zeros(1,K);
while over~=0,
    [min_val, min_pos]=min(SLNR_temp);
    changed=0;
    for i1=1:candidate_num,
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
                H_now(i2,:)=H((temp1(1,i2)-1)*R+1,:);
            end

            Hi=zeros(1,MT);
            Hk=zeros(1,MT);
            for i2=1:K,
                Hi=H_now(i2,:);
                C=(MT-MT*err_UpperBound)*Hi'*Hi+err_UpperBound*MT*eye(MT);

                D=sigma2*eye(MT);
                for i3=1:K,
                    if i3~=i2,
                        Hk=H_now(i3,:);
                        D=D+(MT-MT*err_UpperBound)*Hk'*Hk+err_UpperBound*MT*eye(MT);
                    end
                end

                temp2(1,i2)=real(max(eig(inv(D)*C)));
            end
            cal_num=cal_num+K;

            if min(temp2)>min_val,
                KK_temp=temp1;
                SLNR_temp=temp2;
                changed=1;
                round_num=round_num+1;
                break;
            end
        end
    end
    
    if changed==0,
        over=0;
        round_num=round_num+1;
    end
end
KK_new=KK_temp;
    
    % 暂不考虑有用户退出的情况
%     otherwise
%         for i1=1:K,
%             KK_temp=zeros(1,K);
%             KK_temp=KK_old;
% 
%             initial=0;
%             if KK_temp(1,i1)==0,
%                 need=need-1;
%                 changed=0;
% 
%                 sub_max=0;
%                 sub_index=0;
%                 for i2=1:K_all,
%                     % check whether this user has been scheduled
%                     scheduled=0;
%                     for i3=1:K,
%                         if KK_temp(1,i3)==index(1,i2),
%                             scheduled=1;
%                         end
%                     end
% 
%                     if scheduled==0,
%                         KK_temp(1,i1)=index(1,i2);
% 
%                         % get the total number of the scheduled user
%                         amount=0;
%                         for i3=1:K,
%                             if KK_temp(1,i3)~=0,
%                                 amount=amount+1;
%                             end
%                         end
% 
%                         temp5=zeros(1,amount);
%                         H_temp=zeros(amount*R,MT);
%                         for i3=1:K,
%                             temp6=1;
%                             if KK_temp(1,i3)~=0,
%                                 H_temp((temp6-1)*R+1:temp6*R,:)=H((KK_temp(1,i3)-1)*R+1:KK_temp(1,i3)*R,:);
%                                 temp6=temp6+1;
%                             end
%                         end
% 
%                         Hi=zeros(R,MT);
%                         Hk=zeros(R,MT);
%                         for i3=1:amount,
%                             Hi=H_temp((i3-1)*R+1:i3*R,:);
%                             C=Hi'*Hi;
% 
%                             D=R*sigma2*eye(MT);
%                             for i4=1:K,
%                                 if i4~=i3,
%                                     Hk=HHH((i3-1)*R+1:i3*R,:);
%                                     D=D+Hk'*Hk;
%                                 end
%                             end
% 
%                             temp5(1,i3)=max(eig(inv(D)*C));
%                         end
% 
%                         [min_val, min_pos]=min(SLNR_old);
%                         if min(temp5)>min_val,
%                             changed=1;
%                             if need==0,
%                                 SLNR_new=temp5;
%                             end
%                             break;
%                         else
%                             if initial==0,
%                                 sub_max=min(temp5);
%                                 sub_index=index(1,i2);
%                                 initial=1;
%                                 if need==0,
%                                     SLNR_new=temp5;
%                                 end
%                             else
%                                 if min(temp5)>sub_max,
%                                     sub_max=min(temp5);
%                                     sub_index=index(1,i2);
%                                     if need==0,
%                                         SLNR_new=temp5;
%                                     end
%                                 end
%                             end
%                             KK_temp(1,i1)=0;
%                         end
%                     end
%                 end
%                 if changed==0,
%                     KK_temp(1,i1)=sub_index;
%                 end
%             end
%         end
%         KK_new=KK_temp;