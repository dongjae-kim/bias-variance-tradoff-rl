clear all
close all

RT1_column=7;
RT2_column=8;
session_half=80;

load('F:\0-Program\MATLABwork\work\ModelRL\behavior_data_glascher\data.mat')

total_epoch=size(data,2);
data0=cell(1,total_epoch-2);
ind_included=[1:3,5:12,14:total_epoch];
for j=1:1:length(ind_included)
    data0{1,j}=data{1,ind_included(j)};
end

total_epoch=size(data0,2);
total_trial=size(data0{1,1},1);

RT1=zeros(total_epoch,total_trial);
RT2=zeros(total_epoch,total_trial);

for i=1:1:total_epoch
    RT1(i,:)=data0{1,i}(:,RT1_column)';
    RT2(i,:)=data0{1,i}(:,RT2_column)';
end


RT1_r_1st=reshape(RT1(:,1:session_half),[total_epoch*session_half,1]);
RT1_r_2nd=reshape(RT1(:,session_half+1:end),[total_epoch*(total_trial-session_half),1]);
RT2_r_1st=reshape(RT2(:,1:session_half),[total_epoch*session_half,1]);
RT2_r_2nd=reshape(RT2(:,session_half+1:end),[total_epoch*(total_trial-session_half),1]);

%% SPE correlation analysis
[myMap N_state N_action N_transition]=Model_Map_Init('glascher2010');
SPE_1=zeros(total_epoch,total_trial-session_half);
SPE_2=zeros(total_epoch,total_trial-session_half);
for i=1:1:total_epoch
    for j=1:1:(total_trial-session_half)
    end
end


%%
figure('Name','all trials')
subplot(2,1,1)
boxplot(RT1)
subplot(2,1,2)
boxplot(RT2)

%%
figure('Name','session 1st vs 2nd / state 1 vs 2')


subplot(2,2,1)
boxplot([RT1_r_1st RT1_r_2nd])
[h,p]=ttest(RT1_r_1st,RT1_r_2nd);
str=sprintf('State1: session1 vs session2 (t-test p=%e)',p);
title(str)

subplot(2,2,2)
boxplot([RT2_r_1st RT2_r_2nd])
[h,p]=ttest(RT2_r_1st,RT2_r_2nd);
str=sprintf('State2: session1 vs session2 (t-test p=%e)',p);
title(str)

subplot(2,2,3)
boxplot([RT1_r_1st RT2_r_1st])
[h,p]=ttest(RT1_r_1st,RT2_r_1st);
str=sprintf('Session1: state1 vs state2 (t-test p=%e)',p);
title(str)

subplot(2,2,4)
boxplot([RT1_r_2nd RT2_r_2nd])
[h,p]=ttest(RT1_r_2nd,RT2_r_2nd);
str=sprintf('Session2: state1 vs state2 (t-test p=%e)',p);
title(str)


%%
figure('Name','RT in state1 vs RT in state 2')
RT1_r=reshape(RT1,[total_epoch*total_trial,1]);
RT2_r=reshape(RT2,[total_epoch*total_trial,1]);
boxplot([RT1_r RT2_r])
[h,p]=ttest(RT1_r,RT2_r);
str=sprintf('State1 vs State 2 (t-test p=%e)',p);
title(str)

%%
figure('Name','mean RT trace: state1 vs state2')
ind=[1:1:total_trial-session_half];
plot(ind,mean(RT1(:,session_half+1:end)),'b',ind,mean(RT2(:,session_half+1:end)),'r')
legend('distal','proximal')
grid on