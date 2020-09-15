%Powered By Walter
%期末试题第一题
%infile: "Final_Exam1.txt"
close all;
clear
clc
%%
%测试部分
f = @(x) x.^3 + x.^2 + x + 1;
F = @(x) (x.^4)./4.0 + (x.^3)./3.0 + (x.^2)./2.0 + x;
B = zeros(2,19);
B(1,:) = 1.2:0.1:3.0;
B(2,:) = f(B(1,:));
sum1 = 0;
for i = 1:3:16
    x0 = B(1, i : i + 3);
    y0 = B(2, i : i + 3);
    x0i = x0(1) : 0.00005 : x0(4);
    y0i = interp1(x0, y0, x0i, 'spline');
    sum1 = sum1 + trapz(x0i, y0i);
end
sum2 = F(3.0)-F(1.2);
%%
%计算部分
A = importdata('Final_Exam1.txt');
sum3 = 0;
for i = 1:3:16
    x = A(1, i : i + 3);
    y = A(2, i : i + 3);
    xi = x(1) : 0.00005 : x(4);
    yi = interp1(x, y, xi, 'spline');
    sum3 = sum3 + trapz(xi, yi);
end
sum1
sum2
dif = abs(sum1-sum2)
sum3