%�������¡�Differential Evolution Algorithm With Strategy Adaptation for Global Numerical Optimization�����㷨��ALGORITHMIC DESCRIPTION OF DE
%@written by Zhan Qian,2015-5-24
%���Ժ�����ֵ�ú���testFun(x,FunIndex)
%���������ú���mutation(X,bestX,F,mutationStrategy)
%���������ú���crossover(X,V,CR,crossStrategy)
%mutation
%mutationStrategy=1��DE/rand/1,
%mutationStrategy=2��DE/best/1,
%mutationStrategy=3��DE/rand-to-best/1,
%mutationStrategy=4��DE/best/2,
%mutationStrategy=5��DE/rand/2.
%crossover
%crossStrategy=1:binomial crossover
%crossStrategy=2:Exponential crossover
clear
maxIteration=1000;%����������
Generation=0;%�������������ߵ�ǰ��������
Xmax=30;%�����Ͻ磬���Ը�����Ҫ��Ϊ������ʽ
Xmin=-30;%�����½�
Dim=2;%����ά��
NP=200;%population size,��Ⱥ��ģ
F=0.5;%scaling factor ��������
CR=0.3;%crossover rate �������
FunIndex=1;%���Է�����������ֵͬ��Ӧ��ͬ�Ĳ��Ժ���
mutationStrategy=1;%�������
crossStrategy=2;%�������
%%
%step1 ��ʼ��
%X represent population
%Generation=0;
X=(Xmax-Xmin)*rand(NP,Dim)+Xmin;%X�д������i���д������i��ά��j
 
%%
%step2 mutation,crossover,selection
while Generation<maxIteration
%��bestX
    for i=1:NP
        fitnessX(i)=testFun(X(i,:),FunIndex);%fitnessX��ʾX����Ӧֵ
    end
    [fitnessbestX,indexbestX]=min(fitnessX);
    bestX=X(indexbestX,:);%bestX��ʾ����ֵ��Ӧ��λ��
%%
%step2.1 mutation
%mutationStrategy=1��DE/rand/1,
%mutationStrategy=2��DE/best/1,
%mutationStrategy=3��DE/rand-to-best/1,
%mutationStrategy=4��DE/best/2,
%mutationStrategy=5��DE/rand/2,
%����Ϊÿһ������Xi,G ����һ����������Vi,G�� G�����������
    V=mutation(X,bestX,F,mutationStrategy);
 %%   
%step2.2 crossover
%crossStrategy=1:binomial crossover
%crossStrategy=2:Exponential crossover
%����Ϊÿһ������Xi,G ����һ����������Ui,G�� G�����������
    U=crossover(X,V,CR,crossStrategy);
%%    
%step2.3 selection
    for i=1:NP
        fitnessU(i)=testFun(U(i,:),FunIndex);
        if fitnessU(i)<=fitnessX(i)
            X(i,:)=U(i,:);
            fitnessX(i)=fitnessU(i);
            if fitnessU(i)<fitnessbestX
                bestX=U(i,:);
                fitnessbestX=fitnessU(i);
            end
        end
    end
%%
    Generation=Generation+1;
    bestfitnessG(Generation)=fitnessbestX;
end
 
%%
%��ͼ
plot(bestfitnessG);
optValue=num2str(fitnessbestX);
Location=num2str(bestX);
disp(strcat('the optimal value','=',optValue));
disp(strcat('the best location','=',Location));
