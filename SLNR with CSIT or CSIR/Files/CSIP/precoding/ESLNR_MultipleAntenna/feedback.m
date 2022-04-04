function [FB_Info, index] = feedback(MT,R,L,H,codebook)

FB_Info=zeros(1,R);
index=zeros(1,R);

error=zeros(1,R);
temp1=zeros(1,R);

for i1=1:R,
    HHH=zeros(1,MT);
    HHH=H(i1,:);

    temp2=zeros(1,L);
    for i2=1:L,
        code=zeros(1,MT);
        code=codebook(i2,:);

        temp2(1,i2)=sqrt(abs(trace(HHH'*code)));
    end
    
    [max_val,max_pos]=max(temp2);

    error(1,i1)=max_val;
    temp1(1,i1)=max_pos;
end

value=zeros(1,R);
[value, index]=another_sort(error);

for i1=1:R,
    FB_Info(1,i1)=temp1(1,index(1,i1));
end