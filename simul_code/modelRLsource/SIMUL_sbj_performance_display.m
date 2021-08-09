clear all
close all

LIST_SBJ={'david', 'DeDe', 'rosemary', 'Boyu', 'melissa', 'Rehevolew', 'joel', 'clarke', 'angela', 'william', 'josephine'};
pop_all=[1:1:size(LIST_SBJ,2)];
pop_non_caltech=[1 2 4:8 11]; pop_caltech=[3 9 10];
pop_male=[1 3 7 8 10]; pop_female=[2 4 5 6 9 11];
pop_asian=[2 3 6 9 10 11];  pop_caucasian=[1 4 5 7 8];
pop_worst=[3 5];    pop_exclude_worst=[1 2 4 6:11];
num_sbj_total=size(LIST_SBJ,2);

%% option
ind_sbj_included=pop_exclude_worst; 


%% process files
num_sbj_included=length(ind_sbj_included);
for k=1:1:num_sbj_included
    LIST_SBJ_included{1,k}=LIST_SBJ{1,ind_sbj_included(k)};
end

mat_in_all=cell(num_sbj_included,4); % amount record for each subject, for each condition
mat_suc_all=cell(num_sbj_included,4); % binary hit record for each subject, for each condition
total_reward=zeros(num_sbj_included,1);
for i=1:1:num_sbj_included
    
    % file load
    file_name=[LIST_SBJ{1,ind_sbj_included(i)} '_fmri_info.mat'];
    file_name_full=[pwd '\result_save\' file_name];
    load(file_name_full);
    
    mat_in=zeros(5,4);
    for session_ind=1:1:5
        
        tot_block=size(HIST_block_condition{1,session_ind},2);
        
        for block_ind=1:1:tot_block % reward, hit rate
            block_condition=HIST_block_condition{1,session_ind}(2,block_ind); % G:1,G':2,H:3,H':4
            mat_in(session_ind,block_condition)=mat_in(session_ind,block_condition)+HIST_behavior_info{1,session_ind}(block_ind,16);
            mat_in_all{i,block_condition}=[mat_in_all{i,block_condition} HIST_behavior_info{1,session_ind}(block_ind,16)];
            
            if(block_ind>=3) % G condition
                mat_suc_all{i,block_condition}=[mat_suc_all{i,block_condition} (HIST_behavior_info{1,session_ind}(block_ind,16)>0)];
            else
                mat_suc_all{i,block_condition}=[mat_suc_all{i,block_condition} (HIST_behavior_info{1,session_ind}(block_ind,16)>0)];
            end
        end
        
    end
    
    % sanity check
    total_earn=0;
    for j=1:1:5
        total_earn=total_earn+sum(HIST_behavior_info{1,j}(:,16));
    end
    total_reward(i,1)=total_earn;
    if(sum(sum(mat_in))~=total_earn)
        disp('### WARNING: number of success mismatch ###');
    end
    
end

%% box plot (individual)
XX1_all=cell(1,4); 
XX2_all=cell(1,4);
new_order_condition_ind=[1 2 4 3];
new_order_condition_name={'G_l','G_h','H_l','H_h'};
[score_tot sub_ind2]=sort(total_reward,'descend');
for k=1:1:num_sbj_included % subject
    % indentify rank in terms of total reward
    ranking=find(sub_ind2==k);    
    figure('Name',sprintf('Subject ''%s'' (rank#%d, total rwd=%d)',LIST_SBJ_included{1,k},ranking,total_reward(k)));
    XX=[];  ORG=[];
    XX1=[]; STD1=[];
    XX2=[]; STD2=[];
    for cc=1:1:4 % condition
        %         XX=[XX mat_in_all{k,cc}];
        %         ORG=[ORG cc*ones(1,size(mat_in_all{k,cc},2))];
        XX1=[XX1 mean(mat_in_all{k,cc})];
        STD1=[STD1 std(mat_in_all{k,cc})];
        XX2=[XX2 mean(mat_suc_all{k,cc})];
        STD2=[STD2 std(mat_suc_all{k,cc})];
        
        % collect every single entity
%         XX1_all{1,cc}=[XX1_all{1,cc} mat_in_all{k,cc}];
%         XX2_all{1,cc}=[XX2_all{1,cc} mat_suc_all{k,cc}];
        
        % collect mean for each subject
        XX1_all{1,cc}=[XX1_all{1,cc} mean(mat_in_all{k,cc})];
        XX2_all{1,cc}=[XX2_all{1,cc} mean(mat_suc_all{k,cc})];
    end    
    %     boxplot(XX',ORG');    
    subplot(1,2,1)
    errorbar(XX1(new_order_condition_ind),STD1(new_order_condition_ind),'or')
    set(gca,'XtickL',new_order_condition_name);
    title('reward - G_l|G_h|H_l|H_h');
    subplot(1,2,2)
    errorbar(XX2(new_order_condition_ind),STD2(new_order_condition_ind),'or')
    set(gca,'XtickL',new_order_condition_name);
    title('hit rate - G_l|G_h|H_l|H_h')
end

%% blox plot (total)
figure('Name','all subjects');
XX1_all0=[];    XX2_all0=[];
STD1_all0=[];    STD2_all0=[];
for cc=1:1:4
    XX1_all0=[XX1_all0 mean(XX1_all{1,cc})];
    XX2_all0=[XX2_all0 mean(XX2_all{1,cc})];
    STD1_all0=[STD1_all0 std(XX1_all{1,cc})];
    STD2_all0=[STD2_all0 std(XX2_all{1,cc})];
end

% (1) error bar plot
subplot(1,2,1)
errorbar(XX1_all0(new_order_condition_ind),STD1_all0(new_order_condition_ind),'or')
set(gca,'XtickL',new_order_condition_name);
title('reward - G_l|G_h|H_l|H_h');
subplot(1,2,2)
errorbar(XX2_all0(new_order_condition_ind),STD2_all0(new_order_condition_ind),'or')
set(gca,'XtickL',new_order_condition_name);
title('hit rate - G_l|G_h|H_l|H_h')

% (2) bar chart with error bar
colormap(summer)

h1=subplot(1,2,1);
y_plot1=[XX1_all0(new_order_condition_ind(1:2)); XX1_all0(new_order_condition_ind(3:4))]; % (2 groups of 2 parameters) 
errY1 = zeros(2,2,2);
lower_bar_bound_size1=[XX1_all0(new_order_condition_ind(1:2)); XX1_all0(new_order_condition_ind(3:4))]-0;
upper_bar_bound_size1=40-[XX1_all0(new_order_condition_ind(1:2)); XX1_all0(new_order_condition_ind(3:4))];
errY1(:,:,1)=min(lower_bar_bound_size1,[STD1_all0(new_order_condition_ind(1:2)); STD1_all0(new_order_condition_ind(3:4))]); % lower error
errY1(:,:,2)=min(upper_bar_bound_size1,[STD1_all0(new_order_condition_ind(1:2)); STD1_all0(new_order_condition_ind(3:4))]); % upper error
barwitherr(errY1, y_plot1);    % Plot with errorbars
set(gca,'XTickLabel',{'Goal-directed','Habitual'})
legend('low uncertainty','high uncertainty')
ylabel('reward value')
current_axis=axis;
axis([current_axis(1:2) 0 41]);

h2=subplot(1,2,2);
y_plot2=[XX2_all0(new_order_condition_ind(1:2)); XX2_all0(new_order_condition_ind(3:4))]; % (2 groups of 2 parameters) 
errY2 = zeros(2,2,2);
lower_bar_bound_size2=[XX1_all0(new_order_condition_ind(1:2)); XX1_all0(new_order_condition_ind(3:4))]-0;
upper_bar_bound_size2=40-[XX1_all0(new_order_condition_ind(1:2)); XX1_all0(new_order_condition_ind(3:4))];
errY2(:,:,1)=min(lower_bar_bound_size2,[STD2_all0(new_order_condition_ind(1:2)); STD2_all0(new_order_condition_ind(3:4))]); % lower error
errY2(:,:,2)=min(upper_bar_bound_size2,[STD2_all0(new_order_condition_ind(1:2)); STD2_all0(new_order_condition_ind(3:4))]); % upper error
barwitherr(errY2, y_plot2);    % Plot with errorbars
set(gca,'XTickLabel',{'Goal-directed','Habitual'})
legend('low uncertainty','high uncertainty')
ylabel('hit rate')
current_axis=axis;
axis([current_axis(1:2) 0 1.05]);

%   y = randn(3,4);         % random y values (3 groups of 4 parameters) 
%   errY = 0.1.*y;          % 10% error
%   barwitherr(errY, y);    % Plot with errorbars
%
%   set(gca,'XTickLabel',{'Group A','Group B','Group C'})
%   legend('Parameter 1','Parameter 2','Parameter 3','Parameter 4')
%   ylabel('Y Value')

disp('done')
