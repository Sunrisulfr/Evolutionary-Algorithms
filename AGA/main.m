% ���Ժ���ͼ��
% ���Ժ���ͼ��
% �Ľ�������Ӧ�Ŵ��㷨��
% �ο����ף�[7] M. Srinivas and L. M. Patnaik, "Adaptive probabilities of crossover and mutation in genetic algorithms," 
%              in IEEE Transactions on Systems, Man, and Cybernetics, vol. 24, no. 4, pp. 656-667, April 1994.
%              doi: 10.1109/21.286385
clc;
clear all;
mode = 'Schaffer';
% mode = 'self_define';
if strcmp(mode, 'Schaffer')
    figure(1)
    x = -4:0.1:4;
    y = -4:0.1:4;
    [X,Y] = meshgrid(x,y);
    % Z = 3*cos(X.*Y)+X+Y.^2;
    Z = 0.5-((sin(sqrt(X.^2+Y.^2)).^2)-0.5)./(1+0.001.*(X.^2+Y.^2)).^2;
    surf(X,Y,Z);
    title('Schaffer Function');
    xlabel('X-��');
    ylabel('Y-��');
    zlabel('Z-��');
    
    figure(2);
    contour(X, Y, Z, 8);
    title('Schaffer�����ȸ���');
    xlabel('X-��');
    ylabel('Y-��');
end
 
if strcmp(mode, 'self_define')
    figure(1);
    x = -4:0.1:4;
    y = -4:0.1:4;
    [X,Y] = meshgrid(x,y);
    % Z = 100.*(Y-X.^2).^2+(1-X).^2;
    Z = (cos(X.^2+Y.^2)-0.1)./(1+0.3*(X.^2+Y.^2).^2)+3;
    surf(X,Y,Z);
    %title('Rosen Brock valley Function');
    title('Self define Function');
    xlabel('X-��');
    ylabel('Y-��');
    zlabel('Z-��');
end
 
clc;
clearvars -except mode;
 
r = 0.2;
b = 3;
NP=100;
% Pc=0.65;   % ��Pc,Pm�����Ľ�Ϊ����Ӧ����
% Pm=0.20;
G=100;    %        �ǵø�
D=2;    % ��������
 
k1 = 1;
k3 = 1;
 
k2 = 0.5;
k4 = 0.5;
 
X_min=-4;
X_max=4;
Y_min=-4;
Y_max=4;
% optimization_trace = [];    % ��ά����, �У��У�Ҷ
for count_1=1:NP  % ������ʼ��
    temp1 = X_min+rand()*(X_max-X_min);
    temp2 = Y_min+rand()*(Y_max-Y_min);
    x(count_1,:) = [temp1,temp2];
end
 
save_pic_cnt = 1;
A = figure(3);
 
for gen=1:G
    pause(0.2);
    if rem(gen, 2)==1
        scatter(x(:,1), x(:, 2));
        axis([-4, 4, -4, 4]);
        title(['��', num2str(gen), '�ε���']);
        xlabel('����X');
        ylabel('����Y');
        base_path = 'C:\Users\18811\Desktop\graph\';
        cnt = num2str(save_pic_cnt);
        tail_path = '.jpg';
        frame = getframe(A);
        im=frame2im(frame);
        path_img = [base_path, cnt, tail_path];
        % imwrite(im, path_img);
        % save_x(:, :, save_pic_cnt) = x;
        save_pic_cnt = save_pic_cnt + 1;
    end
 
    % scatter(0, 0, 'o', 'r');
    for count_2=1:NP
        fitness(count_2)=func(x(count_2,:), mode);
    end
    fitness_ = fitness;
    %[fitness_min,index0] = min(fitness);
    %fitness_max = max(fitness);
    [fitness_max,index0] = max(fitness);
    fitness_average = sum(fitness)/(length(fitness));  % ��Ⱥ��ƽ��ֵ
    collect_fit_average(gen) = fitness_average;   % ������Ӧ�ȵ�ƽ��ֵ
    collect_fitmax_subtract_fit_average(gen) = fitness_max - fitness_average;  % ����f_max-f_average ;
    fitness_min = min(fitness);
    best_indiv = x(index0,:);  % ���ŵĸ���
    % optimization_trace(gen,: , global_count) = best_indiv;
    % best_solution(gen) = fitness_min;
    best_solution(gen) = fitness_max;
    % �����һ������Ӧ��ֵ
    fitness = (fitness - fitness_min)/(fitness_max - fitness_min);
    fitness_sum = sum(fitness);
    fitness = fitness./fitness_sum;
    fitness = cumsum(fitness);  
    
    % ���̶�ѡ��
    newi = 1;
    while newi<=NP
        random_num = rand();   % ���������
        if random_num<fitness(1)
            clone_x(newi, :) = x(1, :);
            newi = newi+1;
        else
            for ct=1:NP-1
                if random_num>fitness(ct) && random_num<fitness(ct+1)
                    clone_x(newi,:) = x(ct,:);
                    newi = newi+1;
                    break;
                end
            end
        end
    end
    % disp(clone_x - x);
    % ���н��棬�������
    % count=0;
    for count=1:2:NP
        % ����Ӧ����Pc.
        % ѡ����������ĸ���Ľϴ����Ӧ��ֵ
        if fitness_(count)>=fitness_(count+1)
            fitness_selected = fitness_(count);
        else
            fitness_selected = fitness_(count+1);
        end
        % ����Pc
        if fitness_selected >= fitness_average
            Pc = k1*(fitness_max-fitness_selected)/(fitness_max-fitness_average);
        else
            Pc = k3;
        end
        collect_Pc(gen, count) = Pc;   % ����Pc��ֵ
        temp_cross = rand();
        if temp_cross < Pc
            % ��������   ע�����ֽ�������Ч������
            temp_alpha = 0.6;  
            cross_x(count,:) = temp_alpha*clone_x(count,:)+(1-temp_alpha)*clone_x(count+1,:);
            cross_x(count+1,:) = temp_alpha*clone_x(count+1,:)+(1-temp_alpha)*clone_x(count,:);
            % �Ľ��Ľ�������  �ο����ף���С��. ʵ���������Ŵ��㷨�ĸĽ�����Ӧ��[D].�����ѧ,2012.   ע�������ֽ�������ʵ�ʵ�Ч��������
            % temp_gama = rand();
            % temp_alpha = 0.98;
            % cross_x(count,:) = temp_alpha*clone_x(count,:)+(1-temp_alpha)*clone_x(count+1,:)+temp_gama*(clone_x(count,:)-clone_x(count+1,:));
            % cross_x(count+1,:) = temp_alpha*clone_x(count+1,:)+(1-temp_alpha)*clone_x(count,:)+temp_gama*(clone_x(count,:)-clone_x(count+1,:));
        else
            cross_x(count,:)=clone_x(count,:);
            cross_x(count+1,:)=clone_x(count+1,:);
        end
        % �߽��������
        if cross_x(count,1)>X_max || cross_x(count,1)<X_min || cross_x(count,2)>Y_max || cross_x(count,2)<Y_min
                temp1 = X_min+rand()*(X_max-X_min);
                temp2 = Y_min+rand()*(Y_max-Y_min);
                cross_x(count,:) = [temp1,temp2];
        end
    end
    cross_x = cross_x(1:NP,:);
    %   cross_xΪ��ɽ���ĸ��壻
    % �������
    for count=1:1:NP
        % ����Pm
        if fitness_(count)>=fitness_average
            Pm = k2*(fitness_max-fitness_(count))/(fitness_max-fitness_average);
        else
            Pm = k4;
        end
        collect_Pm(gen,count) = Pm;     %  ����Pm��ֵ
        temp_mutation=rand();
        if temp_mutation<Pm
            %mutation_x(count,:) = (1+0.01).*cross_x(count,:);       %���ֱ�������Ч��������
            % ��������   �ο����ף���С��. ʵ���������Ŵ��㷨�ĸĽ�����Ӧ��[D].�����ѧ,2012       
            mutation_pos = randi(D);
            if mutation_pos==1
                low = X_min;
                high = X_max;
            else
                low = Y_min;
                high = Y_max;
            end
            s_t(gen) = 1-r^((1-gen/G)^b);
            new_low = cross_x(count, mutation_pos)-s_t(gen)*(cross_x(count, mutation_pos)-low);
            new_high = cross_x(count, mutation_pos)+s_t(gen)*(high-cross_x(count, mutation_pos));
            mutation_x(count, :) = cross_x(count, :);
            mutation_x(count, mutation_pos) = new_low+rand()*(new_high-new_low);
            if mutation_x(count,1)>X_max || mutation_x(count,1)<X_min || mutation_x(count,2)>Y_max || mutation_x(count,2)<Y_min
                temp1 = X_min+rand()*(X_max-X_min);
                temp2 = Y_min+rand()*(Y_max-Y_min);
                mutation_x(count,:) = [temp1,temp2];
            end
        else
            mutation_x(count,:) = cross_x(count,:);
        end        
    end
    %�߽���������
    x=mutation_x(1:NP, :);
    x(1,:)= best_indiv;
end
%% ��ͼ
figure(4)
plot(best_solution);
%hold on;
xlabel('��������');
ylabel('��Ӧ��ֵ');
title('��Ӧ�Ƚ�������');
 
figure(5);
plot(collect_fitmax_subtract_fit_average);
title('f_{max}-f_{average}����');
xlabel('��������');
ylabel('f_{max}-f_{average}');
 
% function f=func(buf)
%     f=0.5-((sin(sqrt(buf(1).^2+buf(2).^2)).^2)-0.5)./(1+0.001.*(buf(1).^2+buf(2).^2)).^2;
% end
 
function f=func(buf, md)
    if strcmp(md, 'Schaffer')
        f=0.5-((sin(sqrt(buf(1).^2+buf(2).^2)).^2)-0.5)./(1+0.001.*(buf(1).^2+buf(2).^2)).^2;
    end
    
    if strcmp(md,'self_define')
        % f = 100*(buf(2)-buf(1).^2).^2+(1-buf(1)).^2;
        f = (cos(buf(1).^2+buf(2).^2)-0.1)./(1+0.3*(buf(1).^2+buf(2).^2).^2)+3;
    end
end