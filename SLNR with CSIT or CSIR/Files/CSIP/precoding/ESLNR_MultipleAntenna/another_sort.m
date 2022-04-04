function [value,index] = another_sort(X)

% For a row vector X, this function will sort the elements of X in descending order and show the index of the element in the old order.
% If there are some elements with the same value, they will be ordered from the left to the right 

[row_num, col_num]=size(X);
Y=zeros(1,col_num);
index=Y;
value=Y;

for i1=1:col_num,
    Y(1,i1)=i1;
end

for i1=1:col_num,
    [max_val, max_pos]=max(X);
    value(1,i1)=max_val;
    index(1,i1)=Y(1,max_pos);
    X(:,max_pos)=[];
    Y(:,max_pos)=[];
end