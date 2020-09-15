%Powered By Walter
%期末试题第二题的画图部分
%infile: "FinalExam2-1.txt","FinalExam2-2.txt"
close all
clear
clc
%%
%读取数据
A = importdata('FinalExam2-1.txt');
B = importdata('FinalExam2-2.txt');
x = A(1,:);
y = A(2,:);
z = A(3,:);
x1 = B(1,:);
y1 = B(2,:);
%%
%画图
figure(1);
plot(x,y);
hold on;
plot(x,z);
figure(2);
plot(x1,y1);