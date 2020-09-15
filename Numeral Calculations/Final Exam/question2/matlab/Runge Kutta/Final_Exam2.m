%Powered By Walter
%期末试题第二题
%include: "Final_Exam2_Fun1.m","Final_Exam2_Fun2.m"
close all
clear
clc
%%
%测试部分
f = @(x) exp(-x).*(4 - 2.*x);
[x0,y0] = ode45(@Final_Exam2_Fun1,[0.0 5.0],[4.0 -2.0]);
y1 = f(x0);
figure(1)
plot(x0,y0(:,1))
hold on
plot(x0,y1)
%%
%计算部分
[x,y] = ode45(@Final_Exam2_Fun2,[0.0 5.0],[0.99752738 -0.0073995279]);
figure(2)
plot(x,y(:,1))