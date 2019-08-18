function y=testFun(x,index)
%x代表参数，index代表测试的函数的选择
%该测试函数为通用测试函数，可以移植
%目录
%  函数名            位置                   最优值
%1.Sphere             0                       0
%2.Camel             多个      
%3.Rosenbrock
switch index
    case 1 %Sphere函数
        y=sum(x.^2);
    case 2 %Camel函数,Dim只能取2
        if length(x)>2
            error('x的维度超出了2');
        end
        xx=x(1);yy=x(2);y=(4-2.1*xx^2+xx^4/3)*xx^2+xx*yy+(-4+4*yy^2)*yy^2;
    case 3 %Rosenbrock函数
        y=0;
        for i=2:length(x)
        	y=y+100*(x(i)-x(i-1)^2)^2+(x(i-1)-1)^2;
        end
    otherwise
        disp('no such function, please choose another');
end