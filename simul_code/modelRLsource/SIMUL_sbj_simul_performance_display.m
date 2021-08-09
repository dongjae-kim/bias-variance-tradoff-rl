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


%% file load
file_name_full=[pwd '\result_save\swsw_total_amount_earned.mat'];
load(file_name_full)
tot_amount_main=swsw_amount_main;

% (0) arbitration: fwd+sarsa
% from save swsw_score_NegLogLik_without_Bwdupdate swsw_score swsw_score_each_mode
file_name_full=[pwd '\result_save\swsw_score_NegLogLik_without_Bwdupdate.mat'];
load(file_name_full)
neg_log_lik_0=swsw_score;
avg_neg_log_lik_each_mode_0=swsw_score_each_mode;

% (1) arbitration: fwdbwd+sarsa
% from save swsw_score_NegLogLik_with_Bwdupdate swsw_score swsw_score_each_mode
file_name_full=[pwd '\result_save\swsw_score_NegLogLik_with_Bwdupdate.mat'];
load(file_name_full)
neg_log_lik_1=swsw_score;
avg_neg_log_lik_each_mode_1=swsw_score_each_mode;

% (2) sarsa
% from save swsw_score_NegLogLik_sarsa_only swsw_score swsw_score_each_mode
file_name_full=[pwd '\result_save\swsw_score_NegLogLik_sarsa_only.mat'];
load(file_name_full)
neg_log_lik_2=swsw_score;
avg_neg_log_lik_each_mode_2=swsw_score_each_mode;

% (3) fwd
% from save swsw_score_NegLogLik_fwd_only swsw_score swsw_score_each_mode
file_name_full=[pwd '\result_save\swsw_score_NegLogLik_fwd_only.mat'];
load(file_name_full)
neg_log_lik_3=swsw_score;
avg_neg_log_lik_each_mode_3=swsw_score_each_mode;

% (4) fwdbwd
% from save swsw_score_NegLogLik_fwdbwd_only swsw_score swsw_score_each_mode
file_name_full=[pwd '\result_save\swsw_score_NegLogLik_fwdbwd_only.mat'];
load(file_name_full)
neg_log_lik_4=swsw_score;
avg_neg_log_lik_each_mode_4=swsw_score_each_mode;


%
XX1_all0=[];    STD1_all0=[];
for jjj=ind_sbj_included
    XX1_all0=[XX1_all0; [mean(neg_log_lik_3(:,jjj)) mean(neg_log_lik_2(:,jjj)) mean(neg_log_lik_1(:,jjj))]];
    STD1_all0=[STD1_all0; [std(neg_log_lik_3(:,jjj)) std(neg_log_lik_2(:,jjj)) std(neg_log_lik_1(:,jjj))]];
end


%% total amount plot (total)
figure('Name','total point earned');
% (1) error bar plot
errorbar(mean(tot_amount_main),std(tot_amount_main),'or')
xlabel('subject #');
ylabel('points')


% (2) bar chart with error bar
figure('Name','Negative log likelihood');
colormap(summer)
y_plot1=XX1_all0; % (9 groups of 3 parameters) 
errY1 = zeros(length(ind_sbj_included),3,2);
errY1(:,:,1)=STD1_all0; % lower error
errY1(:,:,2)=STD1_all0; % upper error
barwitherr(errY1, y_plot1);    % Plot with errorbars
set(gca,'XTickLabel',{'1','2','4','6','7','8','9','10','11'})
legend('fwd', 'sarsa', 'arbitrator(fwd-bwd,sarsa)');
xlabel('subject #');
ylabel('NegLogLik')
current_axis=axis;
axis([current_axis(1:2) min(min(XX1_all0))*0.9 max(max(XX1_all0))*1.1]);


%% Avg neg log-likelihood for each mode
Avg_Neg_LogLik_each_all0=[];    Avg_Neg_LogLik_each_all1=[];    Avg_Neg_LogLik_each_all2=[];    Avg_Neg_LogLik_each_all3=[];    Avg_Neg_LogLik_each_all4=[];
Std_Neg_LogLik_each_all0=[];    Std_Neg_LogLik_each_all1=[];    Std_Neg_LogLik_each_all2=[];    Std_Neg_LogLik_each_all3=[];    Std_Neg_LogLik_each_all4=[];
reject_null_hypo=zeros(length(ind_sbj_included),4); % column: each mode (G,G',H,H')
reject_null_hypo_pval=zeros(length(ind_sbj_included),4); % column: each mode (G,G',H,H')
ii=0;
for jjj=ind_sbj_included
    ii=ii+1;
    Avg_Neg_LogLik_each_all0=[Avg_Neg_LogLik_each_all0; mean(avg_neg_log_lik_each_mode_0{1,jjj})];
    Avg_Neg_LogLik_each_all1=[Avg_Neg_LogLik_each_all1; mean(avg_neg_log_lik_each_mode_1{1,jjj})];
    Avg_Neg_LogLik_each_all2=[Avg_Neg_LogLik_each_all2; mean(avg_neg_log_lik_each_mode_2{1,jjj})];
    Avg_Neg_LogLik_each_all3=[Avg_Neg_LogLik_each_all3; mean(avg_neg_log_lik_each_mode_3{1,jjj})];
    Avg_Neg_LogLik_each_all4=[Avg_Neg_LogLik_each_all4; mean(avg_neg_log_lik_each_mode_4{1,jjj})];
    
    Std_Neg_LogLik_each_all0=[Std_Neg_LogLik_each_all0; std(avg_neg_log_lik_each_mode_0{1,jjj})];
    Std_Neg_LogLik_each_all1=[Std_Neg_LogLik_each_all1; std(avg_neg_log_lik_each_mode_1{1,jjj})];
    Std_Neg_LogLik_each_all2=[Std_Neg_LogLik_each_all2; std(avg_neg_log_lik_each_mode_2{1,jjj})];
    Std_Neg_LogLik_each_all3=[Std_Neg_LogLik_each_all3; std(avg_neg_log_lik_each_mode_3{1,jjj})];
    Std_Neg_LogLik_each_all4=[Std_Neg_LogLik_each_all4; std(avg_neg_log_lik_each_mode_4{1,jjj})];
    
    [reject_null_hypo(ii,:),reject_null_hypo_pval(ii,:)] = ttest2(avg_neg_log_lik_each_mode_4{1,jjj},avg_neg_log_lik_each_mode_2{1,jjj},0.1,[],'unequal');
end

% block condition of each trial
figure_name={'G_{low} block', 'G_{high} block', 'H_{high} block', 'H_{low} block'};
for ind_mode=1:1:4 % 1:G(low T uncertainty), 2:G'(high T uncertainty), 3:H(high T uncertainty), 4:H'(low T uncertainty)
    figure('Name',['Average negative log-likelihood of ' figure_name{1,ind_mode}]);
    colormap(summer)
    y_plot1=[Avg_Neg_LogLik_each_all4(:,ind_mode) Avg_Neg_LogLik_each_all2(:,ind_mode)]; % (9 groups of 2 models (fwdbwd,sarsa))
    errY1 = zeros(length(ind_sbj_included),2,2);
    errY1(:,:,1)=[Std_Neg_LogLik_each_all4(:,ind_mode) Std_Neg_LogLik_each_all2(:,ind_mode)]; % lower error
    errY1(:,:,2)=errY1(:,:,1); % upper error
    barwitherr(errY1, y_plot1);    % Plot with errorbars
    set(gca,'XTickLabel',{'1','2','4','6','7','8','9','10','11'})
    legend('fwd-bwd', 'sarsa');
    xlabel('subject #');
    ylabel('Average NegLogLik')
    current_axis=axis;
    axis([current_axis(1:2) min(min(y_plot1))*0.9 max(max(y_plot1))*1.1]);
end

% avgNegLogLik of models in goal-directed mode across all subjects
figure('Name','Contribution of fwdbwd and sarsa to goal-directed blocks (multiplicative inverse of average negative log-likelihood)');
colormap(summer)
y_plot1=[mean(Avg_Neg_LogLik_each_all4(:,1:2))' mean(Avg_Neg_LogLik_each_all2(:,1:2))']; % (2 groups of 2 models (fwdbwd,sarsa))
errY1 = zeros(2,2,2);
errY1(:,:,1)=[mean(Std_Neg_LogLik_each_all4(:,1:2))' mean(Std_Neg_LogLik_each_all2(:,1:2))']; % lower error
errY1(:,:,2)=errY1(:,:,1); % upper error
barwitherr(errY1, y_plot1);    % Plot with errorbars
set(gca,'XTickLabel',{'low uncertainty','high uncertainty'})
legend('fwd-bwd', 'sarsa');
xlabel('state-transition');
ylabel('multiplicative inverse of average NegLogLik')
current_axis=axis;
axis([current_axis(1:2) min(min(y_plot1))*0.9 max(max(y_plot1))*1.1]);

% avgNegLogLik of models in habitual mode across all subjects
figure('Name','Contribution of fwdbwd and sarsa to habitual blocks (inverse Average negative log-likelihood)');
colormap(summer)
y_plot1=[mean(Avg_Neg_LogLik_each_all4(:,[4 3]))' mean(Avg_Neg_LogLik_each_all2(:,[4 3]))']; % (2 groups of 2 models (fwdbwd,sarsa))
errY1 = zeros(2,2,2);
errY1(:,:,1)=[mean(Std_Neg_LogLik_each_all4(:,[4 3]))' mean(Std_Neg_LogLik_each_all2(:,[4 3]))']; % lower error
errY1(:,:,2)=errY1(:,:,1); % upper error
barwitherr(errY1, y_plot1);    % Plot with errorbars
set(gca,'XTickLabel',{'low uncertainty','high uncertainty'})
legend('fwd-bwd', 'sarsa');
xlabel('state-transition');
ylabel('multiplicative inverse average NegLogLik')
current_axis=axis;
axis([current_axis(1:2) min(min(y_plot1))*0.9 max(max(y_plot1))*1.1]);

disp('done')
