%Powered By Walter
%��ĩ����ڶ���Ĳ����ú���
function dydx = Final_Exam2_Fun1(x,y)
dydx = zeros(2,1);
dydx(1) = y(2);
dydx(2) = -y(1)-2.*y(2);
end