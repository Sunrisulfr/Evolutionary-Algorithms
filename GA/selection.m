function [newpop]=selection(pop,fitvalue) 
totalfit=sum(fitvalue);                   %����Ӧֵ֮��
fitvalue=fitvalue/totalfit;                %�������屻ѡ��ĸ���
fitvalue=cumsum(fitvalue);            %�� fitvalue=[1 2 3 4]���� cumsum(fitvalue)=[1 3 6 10]��������ΪʲôҪ�ۼ� 
[px,py]=size(pop);                       %20*10
ms=sort(rand(px,1));                   %��С��������
fitin=1;
newin=1;
while newin<=px                          %ѡ��20���¸��壬���ظ��������������ܵķ�����̫һ��
        if(ms(newin))<fitvalue(fitin)
                newpop(newin,:)=pop(fitin,:);
                newin=newin+1;
        else
                fitin=fitin+1;
        end
end
