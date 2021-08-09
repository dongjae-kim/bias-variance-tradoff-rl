IS_PARAMETER_PLUG_IN=1;



if(IS_PARAMETER_PLUG_IN==1)
    [aaa pop_id]=max(fitness_val);
    close all
    % pop_id=1; the best one at the previous iteration
else
    clear all
    close all
end
   





%% create MAP
% [myMap N_state N_action N_transition]=Model_Map_Init('daw2005');
[myMap N_state N_action N_transition]=Model_Map_Init('glascher2010');

%% create my arbitrator
myArbitrator=Bayesian_Arbitration_Init(N_state,N_action,N_transition);

%% create my RL
myState=Model_RL_Init(N_state,N_action,N_transition);



%% Load behaviour data
state_column=[1 2 3];   action_column=[4 5];
load('F:\0-Program\MATLABwork\work\ModelRL\behavior_data_glascher\data.mat')
% and attach it to the map
ind_included=[1:3,5:12,14:size(data,2)];
data0=cell(1,length(ind_included));
for j=1:1:length(ind_included)
    data0{1,j}=data{1,ind_included(j)};
end
myMap.data=data0;
Sum_NegLogLik=0.0;




%% model parameter - for the functional mode of RL
% SARSA model
param_sarsa.gamma=1.0; % fixed - not actual parameter
param_sarsa.alpha=0.15; % learning rate (0.1~0.2)
param_sarsa.tau=0.5; % decision focus parameter
% FWD model (for the latent stage)
param_fwd0.alpha=0.1; % learning rate
param_fwd0.tau=param_sarsa.tau; % decision focus parameter
% FWD model
param_fwd.alpha=0.15; % learning rate
param_fwd.tau=param_sarsa.tau; % decision focus parameter





%% [SW] Parameter insertion/change
Sum_NegLogLik=0.0;
if(IS_PARAMETER_PLUG_IN==1)
    if(iscell(parameter_converted))
        parameter_converted=cell2mat(parameter_converted);
    end
    % arbitrator
    myArbitrator.PE_tolerance=parameter_converted(pop_id,1); % defines threshold for zero PE
    myArbitrator.Time_Step=parameter_converted(pop_id,2); % the smaller, the slower
    myArbitrator.A_12=parameter_converted(pop_id,3);  myArbitrator.B_12=parameter_converted(pop_id,4);
    myArbitrator.A_21=parameter_converted(pop_id,5);  myArbitrator.B_21=parameter_converted(pop_id,6);
    % SARSA
    param_sarsa.alpha=parameter_converted(pop_id,7); % learning rate (0.1~0.2)
    % FWD
    param_fwd.alpha=parameter_converted(pop_id,8); % learning rate
    % Softmax parameter for all models
    myArbitrator.tau_softmax=parameter_converted(pop_id,9); % use the same value as sarsa/fwd
    param_sarsa.tau=parameter_converted(pop_id,9);
    param_fwd.tau=parameter_converted(pop_id,9);
    param_fwd0.tau=parameter_converted(pop_id,9);
    % arbitrator start
    myArbitrator.ind_active_model=parameter_converted(pop_id,10);    
end

%% no-Devaluation/Devaluation test
pt_devaluation=[1000]; % if devaluation point > num_max_trial, then no devaluation
pt_ref=pt_devaluation(1);

total_simul=1; 
num_max_trial0=size(data0{1,1},1)/2; % for latent learning
num_max_trial=size(data0{1,1},1)/2;

response_rate_distal_sarsa=zeros(total_simul,2); response_rate_proximal_sarsa=zeros(total_simul,2);
response_rate_distal_fwd=zeros(total_simul,2); response_rate_proximal_fwd=zeros(total_simul,2);
HIST_TRANSITION_RATE_SARSA=zeros(total_simul,num_max_trial);
HIST_PROB_CHOICE_SARSA=zeros(total_simul,num_max_trial);
HIST_PE_MEAN_SARSA=zeros(myArbitrator.K,num_max_trial);
HIST_PE_VAR_SARSA=zeros(myArbitrator.K,num_max_trial);
HIST_invFano_SARSA=zeros(myArbitrator.K,num_max_trial);
HIST_TRANSITION_RATE_FWD=zeros(total_simul,num_max_trial);
HIST_PROB_CHOICE_FWD=zeros(total_simul,num_max_trial);
HIST_PE_MEAN_FWD=zeros(myArbitrator.K,num_max_trial);
HIST_PE_VAR_FWD=zeros(myArbitrator.K,num_max_trial);
HIST_invFano_FWD=zeros(myArbitrator.K,num_max_trial);
HIST_MODEL_CHOICE=zeros(total_simul,num_max_trial);
ZZZ=[];
ZZZ1=[];
ZZZ2=[];
YYY1=[];
YYY2=[];
for ind_dev=[1:1:length(pt_devaluation)]
    
    %% Multiple simulations    
    win_fwd=0; win_sarsa=0;        
    for kk=1:1:total_simul
        
        
%         if(0)
%             %% SARSA model
%             state_sarsa=myState;
%             map=myMap;  map0=myMap; % no reward
%             map0.reward=zeros(N_state,1);
%             map0.epoch=kk;   map.epoch=kk;
%             myArbitrator_sarsa=myArbitrator;
%             
%             i=0;  cond=1;
%             time_start=tic;
%             HIST_RWD_SARSA=zeros(1,num_max_trial);
%             HIST_Q_SARSA=zeros(N_state,2,num_max_trial);
%             while ((i<num_max_trial)&&(cond))
%                 i=i+1;
%                 map.trial=i;
%                 %     disp(sprintf('- [%d/%d] session...',i,num_max_trial));
%                 % initializing the state
%                 [state_sarsa map]=StateClear(state_sarsa,map);
%                 while (~state_sarsa.JobComplete)
%                     % devaluation
%                     if(i>=pt_devaluation(ind_dev))
%                         map.reward(5)=0;
%                     end
%                     % current action selection : (s,a)
%                     state_sarsa=Model_RL(state_sarsa, param_sarsa, map, 'decision');
%                     % state update (get reward and next state) : (r,s')
%                     [state_sarsa map]=StateSpace(state_sarsa,map);  % state index ++
%                     % next action selection : (s',a') - if s' is terminal, then no decision
%                     state_sarsa=Model_RL(state_sarsa, param_sarsa, map, 'decision_hypo');
%                     % model upate
%                     state_sarsa=Model_RL(state_sarsa, param_sarsa, map, 'sarsa_update');
%                     % ARBITRATOR: transition rate of the sarsa only
%                     [myArbitrator_sarsa state_sarsa state_sarsa]=Bayesian_Arbitration_v2(myArbitrator_sarsa, state_sarsa, state_sarsa, map, 2);
%                 end
%                 HIST_TRANSITION_RATE_SARSA(kk,i)=myArbitrator_sarsa.transition_rate21;                
%                 HIST_PE_MEAN_SARSA(:,i)=HIST_PE_MEAN_SARSA(:,i)+myArbitrator_sarsa.m2_mean/total_simul;
%                 HIST_PE_VAR_SARSA(:,i)=HIST_PE_VAR_SARSA(:,i)+myArbitrator_sarsa.m2_var/total_simul;
%                 HIST_invFano_SARSA(:,i)=HIST_invFano_SARSA(:,i)+myArbitrator_sarsa.m2_inv_Fano/total_simul;
%                 HIST_RWD_SARSA(i)=state_sarsa.reward_history(state_sarsa.index);
%                 HIST_Q_SARSA(:,:,i)=state_sarsa.Q(1:N_state,1:N_action);
%                 % response rate
%                 if(i>=pt_ref)
%                     if(state_sarsa.action_history(1)==1)
%                         response_rate_distal_sarsa(kk,ind_dev)=response_rate_distal_sarsa(kk,ind_dev)+1;
%                     end
%                     if(state_sarsa.action_history(2)==2)
%                         response_rate_proximal_sarsa(kk,ind_dev)=response_rate_proximal_sarsa(kk,ind_dev)+1;
%                     end
%                 end
%                 %             %  step check
%                 %             reply = input('- continue? [y/n]: ', 's');
%                 %             if(reply=='n')
%                 %                 cond=0;
%                 %             end
%             end
%             processing_time=toc(time_start);
%         end
        
        
        %% Simulation
        % (1) phase 1 - random action, no reward   
        state_fwd=myState;
        state_sarsa=myState;               
        map=myMap;  map0=myMap; % no reward
        map0.reward=zeros(N_state,1);       
        map0.epoch=kk;   map.epoch=kk;        
        
        i=0;  cond=1;        
        % time_start=tic;        
        HIST_T_FWD_1=zeros(N_state,N_action,N_state,num_max_trial0);        
        if(1) % latent learning
            while ((i<num_max_trial0)&&(cond))
                i=i+1;
                map0.trial=i;
                %             disp(sprintf('- [%d/%d] session...',i,num_max_trial0));
                % initializing the state
                [state_fwd map0]=StateClear(state_fwd,map0);
                while (~state_fwd.JobComplete)
                    % current action selection : (s,a)
%                     state_fwd=Model_RL(state_fwd, param_fwd0, map0, 'decision_behavior_data');        
                    state_fwd=Model_RL(state_fwd, param_fwd0, map0, 'decision_random');
                    % state update (get reward and next state) : (r,s')
                    [state_fwd map0]=StateSpace(state_fwd,map0);  % state index ++                    
                    % model upate
                    state_fwd=Model_RL(state_fwd, param_fwd0, map0, 'fwd_update');                    
                end
                HIST_T_FWD_1(:,:,:,i)=state_fwd.T;
                %  step check
                %     reply = input('- continue? [y/n]: ', 's');
                %     if(reply=='n')
                %         cond=0;
                %     end
            end            
        end
        
        % (2) phase 2 - intact action, rewards given                
%         state_fwd=myState;
        
        i=0;  cond=1;
        time_start=tic;
        myArbitrator_top=myArbitrator;
        
        HIST_RWD_SARSA=zeros(1,num_max_trial);
        HIST_Q_SARSA=zeros(N_state,2,num_max_trial);                
        HIST_RWD_FWD=zeros(1,num_max_trial);
        HIST_Q_FWD=zeros(N_state,N_action,num_max_trial);
        HIST_T_FWD_2=zeros(N_state,N_action,N_state,num_max_trial);      
        HIST_MODEL_CHOICE2=[];
        ZZZ3=[];
        
        while ((i<num_max_trial)&&(cond))
            i=i+1;
            map.trial=num_max_trial0+i;            
            %     disp(sprintf('- [%d/%d] session...',i,num_max_trial));
            % initializing the state            
            [state_fwd map]=StateClear(state_fwd,map);
            [state_sarsa map]=StateClear(state_sarsa,map);
            
            HIST_MODEL_CHOICE(kk,i)=myArbitrator_top.ind_active_model;
            while (((myArbitrator_top.ind_active_model==1)&&(~map.JobComplete))||((myArbitrator_top.ind_active_model==2)&&(~map.JobComplete)))
                            
                HIST_MODEL_CHOICE2=[HIST_MODEL_CHOICE2 myArbitrator_top.ind_active_model];
                % devaluation
                if(i>=pt_devaluation(ind_dev))
                    map.reward(:)=0;
                end
                
                %% Compute negative log-likelihood
                for ll=1:1:length(ind_included) % compute map.trial's log-lik for all subjects
                    state_data=data0{1,ll}(map.trial,state_column(map.index));
                    action_data=data0{1,ll}(map.trial,action_column(map.index));
                    var_exp=exp(myArbitrator_top.tau_softmax*myArbitrator_top.Q(state_data,:)); % (N_actionx1)                    
                    Sum_NegLogLik=Sum_NegLogLik-log(var_exp(action_data)/sum(var_exp))/total_simul;
                end
                %%
                
                % index synchronization
                state_fwd.index=map.index;                state_sarsa.index=map.index;
                                
                if(myArbitrator_top.ind_active_model==1) % fwd                    
                    state_fwd=Model_RL(state_fwd, param_fwd, map, 'decision_arbitrator', myArbitrator_top);
                    % state update (get reward and next state) : (r,s')
                    [state_fwd map]=StateSpace(state_fwd,map);  % map&state index ++
                    % model upate
                    state_fwd=Model_RL(state_fwd, param_fwd, map, 'fwd_update');
                    % state synchronization
                    state_sarsa.state_history(state_fwd.index)=state_fwd.state_history(state_fwd.index);
                end
                
                if(myArbitrator_top.ind_active_model==2) % sarsa
                    % current action selection : (s,a) - using arbitrator's Q-value
                    state_sarsa=Model_RL(state_sarsa, param_sarsa, map, 'decision_arbitrator', myArbitrator_top);
                    % state update (get reward and next state) : (r,s')
                    [state_sarsa map]=StateSpace(state_sarsa,map);  % map&state index ++
                    % next action selection : (s',a') - if s' is terminal, then no decision
                    state_sarsa=Model_RL(state_sarsa, param_sarsa, map, 'decision_hypo');
                    % model upate
                    state_sarsa=Model_RL(state_sarsa, param_sarsa, map, 'sarsa_update');            
                    % state synchronization
                    state_fwd.state_history(state_sarsa.index)=state_sarsa.state_history(state_sarsa.index);
                end
                
                % ARBITRATOR: transition rate of the sarsa only
                [myArbitrator_top, state_fwd, state_sarsa]=Bayesian_Arbitration_v2(myArbitrator_top, state_fwd, state_sarsa, map);
                               
                                
                ZZZ=[ZZZ [map.index; state_fwd.index; state_sarsa.index; myArbitrator_top.index; myArbitrator_top.temp; myArbitrator_top.ind_active_model]];                
                    ZZZ1=[ZZZ1 [i; myArbitrator_top.m1_wgt'*myArbitrator_top.m1_inv_Fano; sum(myArbitrator_top.m1_inv_Fano)]];
                    ZZZ2=[ZZZ2 [i; myArbitrator_top.m2_wgt'*myArbitrator_top.m2_inv_Fano; sum(myArbitrator_top.m2_inv_Fano)]];
                    YYY1=[YYY1 [i; myArbitrator_top.transition_rate12]];
                    YYY2=[YYY2 [i; myArbitrator_top.transition_rate21]];
                
               
            end
               
            ZZZ3=[ZZZ3 state_sarsa.RPE_history];
            HIST_TRANSITION_RATE_SARSA(kk,i)=myArbitrator_top.transition_rate21;
            HIST_PROB_CHOICE_SARSA(kk,i)=myArbitrator_top.m2_prob;
            HIST_PE_MEAN_SARSA(:,i)=HIST_PE_MEAN_SARSA(:,i)+myArbitrator_top.m2_mean/total_simul;
            HIST_PE_VAR_SARSA(:,i)=HIST_PE_VAR_SARSA(:,i)+myArbitrator_top.m2_var/total_simul;
            HIST_invFano_SARSA(:,i)=HIST_invFano_SARSA(:,i)+myArbitrator_top.m2_inv_Fano/total_simul;
            HIST_RWD_SARSA(i)=state_sarsa.reward_history(state_sarsa.index);
            HIST_Q_SARSA(:,:,i)=state_sarsa.Q(1:N_state,1:N_action);
            
            HIST_TRANSITION_RATE_FWD(kk,i)=myArbitrator_top.transition_rate12;
            HIST_PROB_CHOICE_FWD(kk,i)=myArbitrator_top.m1_prob;
            HIST_PE_MEAN_FWD(:,i)=HIST_PE_MEAN_FWD(:,i)+myArbitrator_top.m1_mean/total_simul;
            HIST_PE_VAR_FWD(:,i)=HIST_PE_VAR_FWD(:,i)+myArbitrator_top.m1_var/total_simul;
            HIST_invFano_FWD(:,i)=HIST_invFano_FWD(:,i)+myArbitrator_top.m1_inv_Fano/total_simul;
            HIST_RWD_FWD(i)=state_fwd.reward_history(state_fwd.index);
            HIST_Q_FWD(:,:,i)=state_fwd.Q;
            HIST_T_FWD_2(:,:,:,i)=state_fwd.T;
            % response rate
            if(i>=pt_ref)
                if(state_fwd.action_history(1)==1)
                    response_rate_distal_fwd(kk,ind_dev)=response_rate_distal_fwd(kk,ind_dev)+1;
                end
                if(state_fwd.action_history(2)==2)
                    response_rate_proximal_fwd(kk,ind_dev)=response_rate_proximal_fwd(kk,ind_dev)+1;
                end
            end
            %  step check
            %     reply = input('- continue? [y/n]: ', 's');
            %     if(reply=='n')
            %         cond=0;
            %     end
        end
        processing_time=toc(time_start);
        
        % determine who learns fast (sum of rwd in early stage)
        pt1=1;
        pt2=round(num_max_trial/3);
        if(sum(HIST_RWD_FWD(pt1:pt2))>sum(HIST_RWD_SARSA(pt1:pt2)))
            win_fwd=win_fwd+1;
        else
            win_sarsa=win_sarsa+1;
        end
        disp(sprintf('- condition[%d],simulation[%d/%d] - (fwd_win|sarsa_win)=(%d|%d)',ind_dev,kk,total_simul,win_fwd,win_sarsa))
        
    end
    
    %% figure
    
    str=sprintf('[Condition-%d] Q(S1,L(line)/R(dotted)) - SARSA(blue) vs. FWD(red)',ind_dev);
    figure('Name',str)
    ind=[1:1:length(squeeze(HIST_Q_SARSA(1,1,:)))];
    plot(ind,squeeze(HIST_Q_SARSA(1,1,:)),'b-',ind,squeeze(HIST_Q_SARSA(1,2,:)),'b-.',ind,squeeze(HIST_Q_FWD(1,1,:)),'r-',ind,squeeze(HIST_Q_FWD(1,2,:)),'r-.')
    
    str=sprintf('[Condition-%d] Q(S4,L(line)/R(dotted)) - SARSA(blue) vs. FWD(red)',ind_dev);
    figure('Name',str)
    ind=[1:1:length(squeeze(HIST_Q_SARSA(4,2,:)))];
    plot(ind,squeeze(HIST_Q_SARSA(4,1,:)),'b-',ind,squeeze(HIST_Q_SARSA(4,2,:)),'b-.',ind,squeeze(HIST_Q_FWD(4,1,:)),'r-',ind,squeeze(HIST_Q_FWD(4,2,:)),'r-.')
    
    str=sprintf('[Condition-%d] T(S1,R,S4) - FWD',ind_dev);
    figure('Name',str)
    plot([squeeze(HIST_T_FWD_1(1,2,4,:)); squeeze(HIST_T_FWD_2(1,2,4,:))])

end
disp('-done');


%% figure - alpha/beta
if(1)
    
    figure('Name','transition rate and winner')
    ind_line=[1:1:num_max_trial];
    
    subplot(3,1,1)
    plot(ind_line,mean(HIST_TRANSITION_RATE_SARSA,1),'b-',ind_line,mean(HIST_TRANSITION_RATE_FWD,1),'r-')
    title('transition rate');
    legend('SARSA->FWD (sarsa disbelief)','FWD->SARSA (fwd disbelief)')
    
    subplot(3,1,2)
    winner_fwd=mean(HIST_TRANSITION_RATE_SARSA,1)-mean(HIST_TRANSITION_RATE_FWD,1);
    winner_prob=mean(HIST_TRANSITION_RATE_SARSA-HIST_TRANSITION_RATE_FWD,1);
    plot(ind_line,winner_prob);
    title('transition rate difference (FWD>SARSA)');
    
    subplot(3,1,3)
    plot(ind_line,mean(HIST_PROB_CHOICE_FWD,1),'b-',ind_line,mean(HIST_PROB_CHOICE_SARSA,1),'r-')
    title('prob of choosing fwd');
    legend('p(FWD choice)','p(SARSA choice)')
    
    figure('Name','means (line) with variances (shade)')
    str_title{1}='negative PE';     str_title{2}='zero PE';     str_title{3}='positive PE';
    for j=1:1:myArbitrator.K
        subplot(myArbitrator.K,1,j);
        shade_mat_sarsa=reshape(HIST_PE_VAR_SARSA(j,:)',[size(HIST_PE_VAR_SARSA(j,:)',1) 1 size(HIST_PE_VAR_SARSA(j,:)',2)]);
        shade_mat_fwd=reshape(HIST_PE_VAR_FWD(j,:)',[size(HIST_PE_VAR_FWD(j,:)',1) 1 size(HIST_PE_VAR_FWD(j,:)',2)]);
        boundedline(ind_line,HIST_PE_MEAN_SARSA(j,:),shade_mat_sarsa,'alpha','b',ind_line,HIST_PE_MEAN_FWD(j,:),shade_mat_fwd,'alpha','r');
        legend('E_{sarsa}(\theta)','E_{fwd}(\theta)');
        axis([min(ind_line) max(ind_line) 0 1])
        title(str_title{j});
    end
    
    figure('Name','inverse Fano')
    str_title{1}='negative PE';     str_title{2}='zero PE';     str_title{3}='positive PE';
    for j=1:1:myArbitrator.K
        subplot(myArbitrator.K,1,j);
        plot(ind_line,HIST_invFano_SARSA(j,:),'b',ind_line,HIST_invFano_FWD(j,:),'r');
        legend('invFano(SARSA)','invFano(FWD)');
        title(str_title{j});
    end
    
end

%% figure

% figure('Name','SARSA model reward')
% rwd_acc=[];
% for i=100:10:length(HIST_RWD_SARSA)
%     rwd_acc=[rwd_acc sum(HIST_RWD_SARSA(i-99:i))];
% end
% plot(rwd_acc);
% 
% figure('Name','FWD model reward')
% rwd_acc2=[];
% for i=100:10:length(HIST_RWD_SARSA)
%     rwd_acc2=[rwd_acc2 sum(HIST_RWD_FWD(i-99:i))];
% end
% plot(rwd_acc2);

if(0)
    
    %% t-test
    
    figure('Name','average actions - non-devaluation vs. devaluation (overtraining)')
    
    subplot(2,2,1)
    boxplot(response_rate_distal_sarsa)
    [h,p] = ttest(response_rate_distal_sarsa(:,2),mean(response_rate_distal_sarsa(:,1)));
    str=sprintf('sarsa - distal (t-test p=%e)',p);
    title(str)
    
    subplot(2,2,2)
    boxplot(response_rate_proximal_sarsa)
    [h,p] = ttest(response_rate_proximal_sarsa(:,2),mean(response_rate_proximal_sarsa(:,1)));
    str=sprintf('sarsa - proximal (t-test p=%e)',p);
    title(str)
    
    subplot(2,2,3)
    boxplot(response_rate_distal_fwd)
    [h,p] = ttest(response_rate_distal_fwd(:,2),mean(response_rate_distal_fwd(:,1)));
    str=sprintf('fwd - distal (t-test p=%e)',p);
    title(str)
    
    subplot(2,2,4)
    boxplot(response_rate_proximal_fwd)
    [h,p] = ttest(response_rate_proximal_fwd(:,2),mean(response_rate_proximal_fwd(:,1)));
    str=sprintf('fwd - proximal (t-test p=%e)',p);
    title(str)
    
    %% one way anova table display (similar to t-test result in terms of p-val)
    p=anova1(response_rate_distal_sarsa);
    p=anova1(response_rate_proximal_sarsa);
    p=anova1(response_rate_distal_fwd);
    p=anova1(response_rate_proximal_fwd);
    
end
