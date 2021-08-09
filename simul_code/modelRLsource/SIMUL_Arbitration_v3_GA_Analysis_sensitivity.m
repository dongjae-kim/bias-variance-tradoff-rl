clear all
close all

% option
effect_range=[30 80]; % trial range to see effect
effect_pt2=2*effect_range; %actual column

% file read option
R_cond=[1:1:3];
R_str{1,1}='minimum Risk';
R_str{1,2}='medium Risk';
R_str{1,3}='maximum Risk';
D_cond=[1:1:3];
D_str{1,1}='Devalued';
D_str{1,2}='Rewarded';
D_str{1,3}='Degraded';

%% read all files
folder_work='F:\0-Program\MATLABwork\work\ModelRL\result_simul\';

invFano_mat_fwd=cell(length(R_cond),length(D_cond));
invFano_mat_sarsa=cell(length(R_cond),length(D_cond));
for i=1:1:length(R_cond)
    for j=1:1:length(D_cond)
        file_name=sprintf('R%dD%d.mat',R_cond(i),D_cond(j));
        eval(['load ' folder_work file_name]);
        ind_array=trial_array;
        invFano_mat_fwd{i,j}=X_in_fwd;
        invFano_mat_sarsa{i,j}=X_in_sarsa;
    end
end


%% display (all together)
% inv Fano
str_fig=sprintf('Average inv Fano (Sbjs0 use only) - ALL');
figure('Name',str_fig);
for i=1:1:length(R_cond)
    for j=1:1:length(D_cond)
        
        subplot(3,3,(3*(i-1)+j));
        var_invFano_fwd=var(invFano_mat_fwd{R_cond(i),D_cond(j)});  var_invFano_sarsa=var(invFano_mat_sarsa{R_cond(i),D_cond(j)});
        mean_invFano_fwd=mean(invFano_mat_fwd{R_cond(i),D_cond(j)});  mean_invFano_sarsa=mean(invFano_mat_sarsa{R_cond(i),D_cond(j)});
        shade_mat_invFano_fwd=reshape(var_invFano_fwd',[size(var_invFano_fwd',1) 1 size(var_invFano_fwd',2)]);
        shade_mat_invFano_sarsa=reshape(var_invFano_sarsa',[size(var_invFano_sarsa',1) 1 size(var_invFano_sarsa',2)]);
        boundedline(trial_array,mean_invFano_fwd,shade_mat_invFano_fwd,'r',trial_array,mean_invFano_sarsa,shade_mat_invFano_sarsa,'b');
        axis([min(trial_array) max(trial_array) 0 1]);
        
%         str_t=sprintf('Risk%d-Dev%d',R_cond(i),D_cond(j));
        str_t=[R_str{1,i} '-' D_str{1,j}];
        title(str_t);
    end
end



%% box plot and t-test2
Risk_condition_A_vs_B=[1 3];
Dev_condition_A_vs_B=[3 1];

% display of inv Fano
for i=1:1:length(Risk_condition_A_vs_B)
    
    str_fig=sprintf('Average inv Fano (Sbjs0 use only) - R%dD%d',Risk_condition_A_vs_B(i),Dev_condition_A_vs_B(i));
    figure('Name',str_fig);
    var_invFano_fwd=var(invFano_mat_fwd{Risk_condition_A_vs_B(i),Dev_condition_A_vs_B(i)});  var_invFano_sarsa=var(invFano_mat_sarsa{Risk_condition_A_vs_B(i),Dev_condition_A_vs_B(i)});
    mean_invFano_fwd=mean(invFano_mat_fwd{Risk_condition_A_vs_B(i),Dev_condition_A_vs_B(i)});  mean_invFano_sarsa=mean(invFano_mat_sarsa{Risk_condition_A_vs_B(i),Dev_condition_A_vs_B(i)});
    shade_mat_invFano_fwd=reshape(var_invFano_fwd',[size(var_invFano_fwd',1) 1 size(var_invFano_fwd',2)]);
    shade_mat_invFano_sarsa=reshape(var_invFano_sarsa',[size(var_invFano_sarsa',1) 1 size(var_invFano_sarsa',2)]);
    boundedline(trial_array,mean_invFano_fwd,shade_mat_invFano_fwd,'r',trial_array,mean_invFano_sarsa,shade_mat_invFano_sarsa,'b');
    axis([min(trial_array) max(trial_array) 0 1]);
    
end


% display of boxplot and t-test2
if((max(length(Risk_condition_A_vs_B),length(Dev_condition_A_vs_B))==2)&&(min(length(Risk_condition_A_vs_B),length(Dev_condition_A_vs_B))==2))
    figure('Name',file_name);
    
    % condition A
    i=1; j=1;
    show_mat_fwd=invFano_mat_fwd{Risk_condition_A_vs_B(i),Dev_condition_A_vs_B(j)}(:,effect_pt2(1):effect_pt2(2));
    show_mat_fwd0_A=reshape(show_mat_fwd,[size(show_mat_fwd,1)*size(show_mat_fwd,2) 1]);
    show_mat_sarsa=invFano_mat_sarsa{Risk_condition_A_vs_B(i),Dev_condition_A_vs_B(j)}(:,effect_pt2(1):effect_pt2(2));
    show_mat_sarsa0_A=reshape(show_mat_sarsa,[size(show_mat_sarsa,1)*size(show_mat_sarsa,2) 1]);
    
    % condition B
    i=2; j=2;    
    show_mat_fwd=invFano_mat_fwd{Risk_condition_A_vs_B(i),Dev_condition_A_vs_B(j)}(:,effect_pt2(1):effect_pt2(2));
    show_mat_fwd0_B=reshape(show_mat_fwd,[size(show_mat_fwd,1)*size(show_mat_fwd,2) 1]);
    show_mat_sarsa=invFano_mat_sarsa{Risk_condition_A_vs_B(i),Dev_condition_A_vs_B(j)}(:,effect_pt2(1):effect_pt2(2));
    show_mat_sarsa0_B=reshape(show_mat_sarsa,[size(show_mat_sarsa,1)*size(show_mat_sarsa,2) 1]);
    
    % fwd under condition A & B
    subplot(1,2,1);
    boxplot([show_mat_fwd0_A show_mat_fwd0_B])
    [h,p] = ttest2(show_mat_fwd0_A,show_mat_fwd0_B);
    str=sprintf('fwd (t-test p=%e)',p);
    title(str)
    
    % sarsa under condition A & B
    subplot(1,2,2);
    boxplot([show_mat_sarsa0_A show_mat_sarsa0_B])
    [h,p] = ttest2(show_mat_sarsa0_A,show_mat_sarsa0_B);
    str=sprintf('sarsa (t-test p=%e)',p);
    title(str)
end

disp('done');