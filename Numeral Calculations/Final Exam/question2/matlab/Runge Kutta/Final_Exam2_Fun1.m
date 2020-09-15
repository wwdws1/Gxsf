%Powered By Walter
%期末试题第二题的测试用函数
function dydx = Final_Exam2_Fun1(x,y)
dydx = zeros(2,1);
dydx(1) = y(2);
dydx(2) = -y(1)-2.*y(2);
end