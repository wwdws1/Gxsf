%Powered By Walter
%期末试题第二题的题目函数
function dydx = Final_Exam2_Fun2(x, y)
dydx = zeros(2,1);
dydx(1) = y(2);
dydx(2) = ((18.0.*exp(6.0.*x-12.0))/((1.0+exp(3.0.*x-6.0))^2)).*y(1)-(9.0.*exp(3.0.*x-6.0))/((1.0+exp(3.0.*x-6.0))^2);
end