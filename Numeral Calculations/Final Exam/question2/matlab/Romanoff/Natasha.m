%Powered By Walter
%期末试题第二题
%Forever black widow, Natasha Romanoff.
close all
clear
clc
%%
h=0.02;
f=@(x) exp(2.*x.^2);
F=@(x) 4+16.*x.^2;
Fx=@(x) 18.*exp(6.*x-12)/(1+exp(3.*x-6))^2;
Gx=@(x) -9.*exp(3.*x-6)/(1+exp(3.*x-6))^2;
x=0:0.02:5;
%%
y0=f(x);
y1=zeros(1,251);
y1(1)=f(0.0);
y1(2)=f(0.02);
%%
for i=2:250
    y1(i+1)=((2+5.0*h^2/6.0.*F(x(i))).*y1(i)-(1-h^2/12.0*F(x(i-1))).*y1(i-1))/(1-h^2/12.0*F(x(i+1)));
end
%%
figure(1)
plot(x,y0)
hold on
plot(x,y1)
%%
y=zeros(1,251);
y(1)=0.99752738;
y(2)=0.99737488;
%%
for i=2:250
    y(i+1)=((2.0+5.0*h^2/6.0.*Fx(x(i))).*y(i)-(1.0-h^2/12.0*Fx(x(i-1))).*y(i-1)+h^2/12.0*(Gx(x(i-1))+10.0.*Gx(x(i))+Gx(x(i+1))))/(1-h^2/12.0*Fx(x(i+1)));
end
%%
figure(2)
plot(x,y)