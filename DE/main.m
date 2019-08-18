%根据文章《Differential Evolution Algorithm With Strategy Adaptation for Global Numerical Optimization》的算法：ALGORITHMIC DESCRIPTION OF DE
%@written by Zhan Qian,2015-5-24
%测试函数求值用函数testFun(x,FunIndex)
%变异向量用函数mutation(X,bestX,F,mutationStrategy)
%交叉向量用函数crossover(X,V,CR,crossStrategy)
%mutation
%mutationStrategy=1：DE/rand/1,
%mutationStrategy=2：DE/best/1,
%mutationStrategy=3：DE/rand-to-best/1,
%mutationStrategy=4：DE/best/2,
%mutationStrategy=5：DE/rand/2.
%crossover
%crossStrategy=1:binomial crossover
%crossStrategy=2:Exponential crossover
clear
maxIteration=1000;%最大迭代次数
Generation=0;%进化代数，或者当前迭代代数
Xmax=30;%搜索上界，可以根据需要改为向量形式
Xmin=-30;%搜索下界
Dim=2;%个体维数
NP=200;%population size,种群规模
F=0.5;%scaling factor 缩放因子
CR=0.3;%crossover rate 交叉概率
FunIndex=1;%测试方程索引，不同值对应不同的测试函数
mutationStrategy=1;%变异策略
crossStrategy=2;%交叉策略
%%
%step1 初始化
%X represent population
%Generation=0;
X=(Xmax-Xmin)*rand(NP,Dim)+Xmin;%X行代表个体i，列代表个体i的维度j
 
%%
%step2 mutation,crossover,selection
while Generation<maxIteration
%求bestX
    for i=1:NP
        fitnessX(i)=testFun(X(i,:),FunIndex);%fitnessX表示X的适应值
    end
    [fitnessbestX,indexbestX]=min(fitnessX);
    bestX=X(indexbestX,:);%bestX表示最优值对应的位置
%%
%step2.1 mutation
%mutationStrategy=1：DE/rand/1,
%mutationStrategy=2：DE/best/1,
%mutationStrategy=3：DE/rand-to-best/1,
%mutationStrategy=4：DE/best/2,
%mutationStrategy=5：DE/rand/2,
%产生为每一个个体Xi,G 产生一个变异向量Vi,G。 G代表进化代数
    V=mutation(X,bestX,F,mutationStrategy);
 %%   
%step2.2 crossover
%crossStrategy=1:binomial crossover
%crossStrategy=2:Exponential crossover
%产生为每一个个体Xi,G 产生一个交叉向量Ui,G。 G代表进化代数
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
%画图
plot(bestfitnessG);
optValue=num2str(fitnessbestX);
Location=num2str(bestX);
disp(strcat('the optimal value','=',optValue));
disp(strcat('the best location','=',Location));
