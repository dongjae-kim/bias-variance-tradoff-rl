clear all
close all

state1_column=2;
state2_column=3;
action1_column=4;
action2_column=5;
RT1_column=7;
RT2_column=8;
RWD_column=14;
session_half=80;

load('F:\0-Program\MATLABwork\work\ModelRL\behavior_data_glascher\data.mat')

total_epoch=size(data,2);
ind_included=[1:3,5:12,14:total_epoch];
data0=cell(1,length(ind_included));
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


if(0)
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
end




%% SPE correlation analysis
[map N_state N_action N_transition]=Model_Map_Init('glascher2010');

SPE1=zeros(total_epoch,total_trial-session_half);
SPE2=zeros(total_epoch,total_trial-session_half);
for i=1:1:total_epoch
    for j=1:1:(total_trial-session_half)
        
        action_s1=data0{1,i}(session_half+j,action1_column);
        prob_mat=map.action(1,action_s1).prob(1,:); %(ex) [... 0 0.3 0 0.7]
        [tmp state_expected]=max(prob_mat);
        state_actual=data0{1,i}(session_half+j,state1_column);
        if(state_expected~=state_actual)
            SPE1(i,j)=1;
        end
        
        action_s2=data0{1,i}(session_half+j,action2_column);
        prob_mat=map.action(1,action_s2).prob(state_actual,:); %(ex) [... 0 0.3 0 0.7]
        [tmp state_expected]=max(prob_mat);
        state_actual=data0{1,i}(session_half+j,state2_column);
        if(state_expected~=state_actual)
            SPE2(i,j)=1;
        end
    end
end

RT1_normalized=RT1(:,session_half+1:end);
% RT1_normalized=(RT1_normalized-min(min(RT1_normalized)))/(max(max(RT1_normalized))-min(max(RT1_normalized)));
RT2_normalized=RT2(:,session_half+1:end);
% RT2_normalized=(RT2_normalized-min(min(RT2_normalized)))/(max(max(RT2_normalized))-min(max(RT2_normalized)));

%% RT analysis
% corr_RT11=RT1_normalized*SPE1';
% corr_RT12=RT1_normalized*SPE2';
% [corr_RT11 corr_RT11_pval]=corrcoef(RT1_normalized(1,:)',SPE1(1,:)');
div=8;
p_history=[];
range_index=[];
rr_history=[];      pp_history=[];
for kk=1:1:div
    if(kk==1)
        figure('Name','RT(distal)-RT(proximal) -1');
    end
    if(kk==div/2+1)
        figure('Name','RT(distal)-RT(proximal) -2');
    end
    if(kk<=div/2)
        subplot(div/2,1,kk);
    end
    if(kk>div/2)
        subplot(div/2,1,kk-div/2);
    end
    range=[(1+(kk-1)*total_trial/(2*div)):1:kk*total_trial/(2*div)];
    range_index=[range_index session_half+kk*total_trial/(2*div)];
    rtt1=reshape(RT1_normalized(:,range),[1 size(RT1_normalized(:,range),1)*size(RT1_normalized(:,range),2)]);
    rtt2=reshape(RT2_normalized(:,range),[1 size(RT2_normalized(:,range),1)*size(RT2_normalized(:,range),2)]);
    hold on
    scatter(rtt1,rtt2,'*');
    p = polyfit(rtt1,rtt2,1);    
    p_history=[p_history; p];
    line([0 1000],[p(1)*0+p(2) p(1)*1000+p(2)],[1 1],'LineStyle','-','Color','r')
    
    [RR PP] = corrcoef(rtt1,rtt2);
    rr_history=[rr_history RR(1,2)];    pp_history=[pp_history PP(1,2)];
    
    str=sprintf('[%d-%d trial] RT2/RT1 slope=%f, Corr.Coeff=%f(p-val=%f)',session_half+min(range),session_half+max(range),p(1),RR(1,2),PP(1,2));
    title(str)
    hold off
        
end
figure('Name','RT2-RT1 trace');
subplot(2,1,1)
plot(range_index,p_history(:,1),'-o')
title('1st order A (RT2=A*RT1+B)');
subplot(2,1,2)
plot(range_index,p_history(:,2),'-o')
title('0th order B (RT2=A*RT1+B)');

figure('Name','RT2-RT1 correlation');
errorbar(range_index,rr_history,pp_history,'r')
title('correlation coefficient (with p-val)');



%% RT2|SPE1

SPE_condition=1;
RT_condition=2;
epoch_index=[1:1:total_epoch];

if(SPE_condition==1)
    SPE_in=SPE1;
else    
    SPE_in=SPE2;
end
if(RT_condition==1)
    RT_in=RT1_normalized;
else    
    RT_in=RT2_normalized;
end
 
div=8;
mean_diff=zeros(div,1);
mean1=zeros(div,1); mean2=zeros(div,1);
std1=zeros(div,1); std2=zeros(div,1);
range_ind=[];
RT2_SPE1_history=cell(1,div);
RT2_SPE0_history=cell(1,div);
for kk=1:1:div    
    if(kk==1)
        str=sprintf('SPE%d-RT%d (i)',SPE_condition,RT_condition);
        figure('Name',str);
    end
    if(kk==div/2+1)
        str=sprintf('SPE%d-RT%d (ii)',SPE_condition,RT_condition);
        figure('Name',str);
    end
    if(kk<=div/2)
        subplot(div/2,1,kk);
    end
    if(kk>div/2)
        subplot(div/2,1,kk-div/2);
    end
    range_ind=[range_ind kk*total_trial/(2*div)+session_half];
    range=[(1+(kk-1)*total_trial/(2*div)):1:kk*total_trial/(2*div)];
    [row_spe0 col_spe0]=find(SPE_in(epoch_index,range)==0);
    rtt0=[];
    for m=1:1:length(row_spe0)
        rtt0=[rtt0 RT_in(row_spe0(m),col_spe0(m))];
    end
    %     rtt0=reshape(rtt0,[1 size(rtt0,1)*size(rtt0,2)]);
    [row_spe1 col_spe1]=find(SPE_in(epoch_index,range)==1);
    rtt2=[];
    for m=1:1:length(row_spe1)
        rtt2=[rtt2 RT_in(row_spe1(m),col_spe1(m))];
    end
    %     rtt2=reshape(rtt2,[1 size(rtt2,1)*size(rtt2,2)]);
    if(RT_condition==2)
        RT2_SPE0_history{1,kk}=rtt0;
        RT2_SPE1_history{1,kk}=rtt2;        
    end
        
    origin=[repmat('RT|SPE=0',[length(rtt0) 1]); repmat('RT|SPE=1',[length(rtt2) 1])];
    boxplot([rtt0'; rtt2'],origin);
    
    [h p]=ttest2(rtt0,rtt2); % h=1: means are different.    
    mean1(kk)=mean(rtt0);       mean2(kk)=mean(rtt2);    
    std1(kk)=std(rtt0);       std2(kk)=std(rtt2);    
    mean_diff(kk)=mean2(kk)-mean1(kk);
    str=sprintf('[%d-%d trial] m_R-m_L=%f, t-test (h,p)=(%d,%e)',session_half+min(range),session_half+max(range),mean_diff(kk),h,p);
    title(str)    
end

str=sprintf('RT%d|SPE%d -Summary',RT_condition,SPE_condition);
figure('Name',str);
hold on
errorbar(range_ind,mean1,std1/2,'r')
errorbar(range_ind,mean2,std2,'b')
title('RT|SPE=0 (red), RT|SPE=1 (blue)')
hold off



%% RT|RPE (reaction time 1/2 after RPE=1/0
RT_condition=2;
epoch_index=[1:1:total_epoch];

RT1_rpe0=zeros(total_epoch,1);  RT2_rpe0=zeros(total_epoch,1);
RT1_rpe1=zeros(total_epoch,1);  RT2_rpe1=zeros(total_epoch,1);
div=8;
mean_diff=zeros(div,1);
mean1=zeros(div,1); mean2=zeros(div,1);
std1=zeros(div,1); std2=zeros(div,1);
range_ind=[];
RT2_RPE1_history=cell(1,div);
RT2_RPE0_history=cell(1,div);

for kk=1:1:div    
    if(kk==1)
        str=sprintf('RPE-RT%d (i)',RT_condition);
        figure('Name',str);        
    end
    if(kk==div/2+1)
        str=sprintf('RPE-RT%d (ii)',RT_condition);
        figure('Name',str);
    end
    if(kk<=div/2)
        subplot(div/2,1,kk);
    end
    if(kk>div/2)
        subplot(div/2,1,kk-div/2);
    end
    
    range_ind=[range_ind kk*total_trial/(2*div)+session_half];
    range=[(1+(kk-1)*total_trial/(2*div)):1:kk*total_trial/(2*div)];
    
    rpe0_index=cell(1,length(epoch_index));     rpe1_index=cell(1,length(epoch_index));
    RT_in_rpe0=[];      RT_in_rpe1=[];
    for n=1:1:length(epoch_index) % for each individual
        
        rpe0_index=find(data0{1,epoch_index(n)}(session_half+range,RWD_column)~=0);
        rpe1_index=find(data0{1,epoch_index(n)}(session_half+range,RWD_column)==0);
        
        next_trial0=session_half+rpe0_index+1;
        next_trial_ind0=find(next_trial0<=total_trial);
        next_trial0=next_trial0(next_trial_ind0);
        next_trial1=session_half+rpe1_index+1;
        next_trial_ind1=find(next_trial1<=total_trial);
        next_trial1=next_trial1(next_trial_ind1);
          
        % collection of all RT in this range of all individuals
        if(RT_condition==1)
            RT_in_rpe0=[RT_in_rpe0 RT1(epoch_index(n),next_trial0)];
            RT_in_rpe1=[RT_in_rpe1 RT1(epoch_index(n),next_trial1)];
        else
            RT_in_rpe0=[RT_in_rpe0 RT2(epoch_index(n),next_trial0)];
            RT_in_rpe1=[RT_in_rpe1 RT2(epoch_index(n),next_trial1)];
        end
        
    end
    if(RT_condition==2)
        RT2_RPE0_history{1,kk}=RT_in_rpe0;
        RT2_RPE1_history{1,kk}=RT_in_rpe1;
    end
    str=sprintf('- ## Range[%d-%d]: %d(RPE=0), %d(RPE=1)',session_half+min(range),session_half+max(range),length(rpe0_index),length(rpe1_index));
    disp(str)
    
    str0=sprintf('RT%d|RPE=0',RT_condition);    str1=sprintf('RT%d|RPE=1',RT_condition);
    origin=[repmat(str0,[length(RT_in_rpe0) 1]); repmat(str1,[length(RT_in_rpe1) 1])];
    boxplot([RT_in_rpe0'; RT_in_rpe1'],origin);
    
    [h p]=ttest2(RT_in_rpe0, RT_in_rpe1); % h=1: means are different.
    mean1(kk)=mean(RT_in_rpe0);       mean2(kk)=mean(RT_in_rpe1);
    std1(kk)=std(RT_in_rpe0);       std2(kk)=std(RT_in_rpe1);
    mean_diff(kk)=mean2(kk)-mean1(kk);
    str=sprintf('[%d-%d trial] m_R-m_L=%f, t-test (h,p)=(%d,%e)',session_half+min(range),session_half+max(range),mean_diff(kk),h,p);
    title(str)
end
[mean(RT_in_rpe1) mean(RT_in_rpe0)]

str=sprintf('RT%d|RPE -Summary',RT_condition);
figure('Name',str);
hold on
errorbar(range_ind,mean1,std1,'r')
errorbar(range_ind,mean2,std2,'b')
str=sprintf('RT|RPE=0 (red), RT|RPE=1 (blue)',RT_condition);
title(str)
hold off

%% RT2|SPE=1 vs RT2|RPE=1
% for this, index size and RT_condition should be the same.
figure('Name','RT2|SPE=1 vs RT2|RPE=1');
for j=1:1:length(range_ind)
    subplot(div,1,j);
    str0=sprintf('RT%d|SPE=1',RT_condition);    str1=sprintf('RT%d|RPE=1',RT_condition);
    origin=[repmat(str0,[length(RT2_SPE1_history{1,j}) 1]); repmat(str1,[length(RT2_RPE1_history{1,j}) 1])];
    boxplot([RT2_SPE1_history{1,j}'; RT2_RPE1_history{1,j}'],origin);
    [h p]=ttest2(RT2_SPE1_history{1,j}', RT2_RPE1_history{1,j}'); % h=1: means are different.
    mean_diff=mean(RT2_RPE1_history{1,j})-mean(RT2_SPE1_history{1,j});
    str=sprintf('[%d-%d trial] m_R-m_L=%f, t-test (h,p)=(%d,%e)',range_index(j)-total_trial/(2*div)+1,range_index(j),mean_diff,h,p);
    title(str)
%     [RR2 PP2]=corrcoef(RT2_SPE1_history,RT2_RPE1_history);
end
figure('Name','RT2|SPE=0 vs RT2|RPE=0');
for j=1:1:length(range_ind)
    subplot(div,1,j);
    str0=sprintf('RT%d|SPE=0',RT_condition);    str1=sprintf('RT%d|RPE=0',RT_condition);
    origin=[repmat(str0,[length(RT2_SPE0_history{1,j}) 1]); repmat(str1,[length(RT2_RPE0_history{1,j}) 1])];
    boxplot([RT2_SPE0_history{1,j}'; RT2_RPE0_history{1,j}'],origin);
    [h p]=ttest2(RT2_SPE0_history{1,j}', RT2_RPE0_history{1,j}'); % h=1: means are different.
    mean_diff=mean(RT2_RPE0_history{1,j})-mean(RT2_SPE0_history{1,j});
    str=sprintf('[%d-%d trial] m_R-m_L=%f, t-test (h,p)=(%d,%e)',range_index(j)-total_trial/(2*div)+1,range_index(j),mean_diff,h,p);
    title(str)
%     [RR2 PP2]=corrcoef(RT2_SPE0_history,RT2_RPE0_history);
end
figure('Name','RT2|SPE=1 vs RT2|RPE=0');
for j=1:1:length(range_ind)
    subplot(div,1,j);
    str0=sprintf('RT%d|SPE=1',RT_condition);    str1=sprintf('RT%d|RPE=0',RT_condition);
    origin=[repmat(str0,[length(RT2_SPE1_history{1,j}) 1]); repmat(str1,[length(RT2_RPE0_history{1,j}) 1])];
    boxplot([RT2_SPE1_history{1,j}'; RT2_RPE0_history{1,j}'],origin);
    [h p]=ttest2(RT2_SPE1_history{1,j}', RT2_RPE0_history{1,j}'); % h=1: means are different.
    mean_diff=mean(RT2_RPE0_history{1,j})-mean(RT2_SPE1_history{1,j});
    str=sprintf('[%d-%d trial] m_R-m_L=%f, t-test (h,p)=(%d,%e)',range_index(j)-total_trial/(2*div)+1,range_index(j),mean_diff,h,p);
    title(str)
%     [RR2 PP2]=corrcoef(RT2_SPE1_history,RT2_RPE0_history);
end
figure('Name','RT2|SPE=0 vs RT2|RPE=1');
for j=1:1:length(range_ind)
    subplot(div,1,j);
    str0=sprintf('RT%d|SPE=0',RT_condition);    str1=sprintf('RT%d|RPE=1',RT_condition);
    origin=[repmat(str0,[length(RT2_SPE0_history{1,j}) 1]); repmat(str1,[length(RT2_RPE1_history{1,j}) 1])];
    boxplot([RT2_SPE0_history{1,j}'; RT2_RPE1_history{1,j}'],origin);
    [h p]=ttest2(RT2_SPE0_history{1,j}', RT2_RPE1_history{1,j}'); % h=1: means are different.
    mean_diff=mean(RT2_RPE1_history{1,j})-mean(RT2_SPE0_history{1,j});
    str=sprintf('[%d-%d trial] m_R-m_L=%f, t-test (h,p)=(%d,%e)',range_index(j)-total_trial/(2*div)+1,range_index(j),mean_diff,h,p);
    title(str)
%     [RR2 PP2]=corrcoef(RT2_SPE0_history,RT2_RPE1_history);
end

%% MISC

% corr_RT22=RT2_normalized*SPE2';
% corr_RT21=RT2_normalized*SPE1';

RT1_diff=RT1_normalized(:,2:end)-RT1_normalized(:,1:end-1);
RT1_diff_normalized=(RT1_diff-min(min(RT1_diff)))/(max(max(RT1_diff)-min(min(RT1_diff))));
corr_RT1_diff=RT1_diff_normalized*SPE1(:,2:end)';

RT2_diff=RT2_normalized(:,2:end)-RT2_normalized(:,1:end-1);
RT2_diff_normalized=(RT2_diff-min(min(RT2_diff)))/(max(max(RT2_diff)-min(min(RT2_diff))));
corr_RT2_diff=RT2_diff_normalized*SPE2(:,2:end)';

% figure('Name','SPE - RT');
% ind=[1:1:(total_trial-session_half)];
% plot(ind,SPE1(1,1:end),'k',ind,RT1_normalized,'b');


disp('-done.');