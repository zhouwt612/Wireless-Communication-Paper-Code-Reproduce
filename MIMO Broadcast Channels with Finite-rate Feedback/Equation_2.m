clc
clear
close

x = -10:0.1:10;
y1 = sin(x);
y2 = cos(x);
y3 = sin(x).^2;
y4 = cos(x).^2;

figure
plot(x,y1,'r-+',x,y2,'g-x',x,y3,'b-*',x,y4,'k-d')
legend('Sin','Cos','Sin^2','Cos^2')