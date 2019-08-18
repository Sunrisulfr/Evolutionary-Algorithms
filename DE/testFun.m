function y=testFun(x,index)
%x���������index������Եĺ�����ѡ��
%�ò��Ժ���Ϊͨ�ò��Ժ�����������ֲ
%Ŀ¼
%  ������            λ��                   ����ֵ
%1.Sphere             0                       0
%2.Camel             ���      
%3.Rosenbrock
switch index
    case 1 %Sphere����
        y=sum(x.^2);
    case 2 %Camel����,Dimֻ��ȡ2
        if length(x)>2
            error('x��ά�ȳ�����2');
        end
        xx=x(1);yy=x(2);y=(4-2.1*xx^2+xx^4/3)*xx^2+xx*yy+(-4+4*yy^2)*yy^2;
    case 3 %Rosenbrock����
        y=0;
        for i=2:length(x)
        	y=y+100*(x(i)-x(i-1)^2)^2+(x(i-1)-1)^2;
        end
    otherwise
        disp('no such function, please choose another');
end