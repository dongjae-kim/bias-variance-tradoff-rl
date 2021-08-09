clear all
close all

%% Option
N_state=5;
N_action=2; % left/right
N_transition=3; % depth:  inital, choice1, choice2

%% my state
myState.N_state=N_state;
myState.N_action=N_action;
myState.N_transition=N_transition;
myState.index=1; % current index
myState.state_history=zeros(N_transition,1);    myState.state_history(1,1)=1;
myState.action_history=zeros(N_transition,1);
myState.action_prob=zeros(N_transition,1);
myState.reward_history=zeros(N_transition,1);
myState.RPE_history=zeros(N_transition,1);
myState.SPE_history=zeros(N_transition,1);
myState.Q=zeros(N_state,N_action);
% myState.IsTerminal=zeros(N_state,1);    myState.IsTerminal(6:21)=1;
myState.JobComplete=0; % if the current state is final, then 1.
% for SARSA model
myState.SARSA=zeros(1,6); % (s,a,r,s',a',complete_index)
% for FWD model
myState.T=zeros(N_state,2,N_state); % (s,a,s')
for j=1:1:N_state % initialized to uniform distributions
    for k=1:1:2
        myState.T(j,k,:)=1/N_state;
    end
end


%% model parameter
% SARSA model
param_sarsa.gamma=1.0; % fixed - not actual parameter
param_sarsa.alpha=0.1; % learning rate
param_sarsa.tau=1.0; % decision focus parameter
% FWD model (for the latent stage)
param_fwd0.alpha=0.1; % learning rate
param_fwd0.tau=param_sarsa.tau; % decision focus parameter
% FWD model
param_fwd.alpha=0.1; % learning rate
param_fwd.tau=param_sarsa.tau; % decision focus parameter


%% MAP: transition prob. for the left(1,1)/right(1,2) action
prob_seed=1.0;
prob_seed_mat=[prob_seed 1-prob_seed];
prob_seed_mat0=[1-prob_seed prob_seed];
myMap.action(1,1).prob=zeros(N_state,N_state);
myMap.action(1,2).prob=zeros(N_state,N_state);
% 1. transition for the left action
myMap.action(1,1).prob(1,2:3)=prob_seed_mat;
myMap.action(1,1).prob(2,4:5)=prob_seed_mat;
% 2. transition for the right action
myMap.action(1,2).prob(1,2:3)=prob_seed_mat0;
myMap.action(1,2).prob(2,4:5)=prob_seed_mat0;
% 3. get a connection matrix
myMap.action(1,1).connection=double(myMap.action(1,1).prob&ones(N_state,N_state));
myMap.action(1,2).connection=double(myMap.action(1,2).prob&ones(N_state,N_state));
% 4. set terminal states
myMap.IsTerminal=zeros(N_state,1);
myMap.IsTerminal(3:5)=[1 1 1];
% 5. reward matrix
myMap.reward=zeros(N_state,1);
myMap.reward(3:5)=[0 0 1]; % [0 10 10 0 10 0 25 0 0 25 0 10 10 0 25 0];




%% no-Devaluation/Devaluation test
pt_devaluation=[1000 120]; % [no-devaluation devaluation]
pt_ref=pt_devaluation(2);

total_simul=250;
num_max_trial0=150; % for latent learning
num_max_trial=150;

response_rate_distal_sarsa=zeros(total_simul,2); response_rate_proximal_sarsa=zeros(total_simul,2);
response_rate_distal_fwd=zeros(total_simul,2); response_rate_proximal_fwd=zeros(total_simul,2);
for ind_dev=[1:1:length(pt_devaluation)]
    
    %% Multiple simulations    
    win_fwd=0; win_sarsa=0;        
    for kk=1:1:total_simul
               
        
        %% SARSA model        
        state_sarsa=myState;        
        map=myMap;  map0=myMap; % no reward
        map0.reward=zeros(N_state,1);
        
        i=0;  cond=1;
        time_start=tic;
        HIST_RWD_SARSA=zeros(1,num_max_trial);
        HIST_Q_SARSA=zeros(5,2,num_max_trial);        
        while ((i<num_max_trial)&&(cond))
            i=i+1;
            %     disp(sprintf('- [%d/%d] session...',i,num_max_trial));
            % initializing the state
            state_sarsa=StateClear(state_sarsa);
            while (~state_sarsa.JobComplete)
                % devaluation
                if(i>=pt_devaluation(ind_dev))
                    map.reward(5)=0;
                end
                % current action selection : (s,a)
                state_sarsa=Model_RL(state_sarsa, param_sarsa, map, 'decision');
                % state update (get reward and next state) : (r,s')
                state_sarsa=StateSpace(state_sarsa,map);  % state index ++
                % next action selection : (s',a') - if s' is terminal, then no decision
                state_sarsa=Model_RL(state_sarsa, param_sarsa, map, 'decision_hypo');
                % model upate
                state_sarsa=Model_RL(state_sarsa, param_sarsa, map, 'sarsa_update');
            end
            HIST_RWD_SARSA(i)=state_sarsa.reward_history(state_sarsa.index);
            HIST_Q_SARSA(:,:,i)=state_sarsa.Q(1:N_state,1:N_action);
            % response rate
            if(i>=pt_ref)
                if(state_sarsa.action_history(1)==1)
                    response_rate_distal_sarsa(kk,ind_dev)=response_rate_distal_sarsa(kk,ind_dev)+1;
                end
                if(state_sarsa.action_history(2)==2)
                    response_rate_proximal_sarsa(kk,ind_dev)=response_rate_proximal_sarsa(kk,ind_dev)+1;
                end
            end
            %             %  step check
            %             reply = input('- continue? [y/n]: ', 's');
            %             if(reply=='n')
            %                 cond=0;
            %             end
        end
        processing_time=toc(time_start);
        
        %% FWD model
        % (1) phase 1 - random action, no reward   
        state_fwd=myState;
        map=myMap;  map0=myMap; % no reward
        map0.reward=zeros(N_state,1);
        i=0;  cond=1;        
        % time_start=tic;
        HIST_T_FWD_1=zeros(N_state,N_action,N_state,num_max_trial0);
        if(0)
            while ((i<num_max_trial0)&&(cond))
                i=i+1;
                %             disp(sprintf('- [%d/%d] session...',i,num_max_trial0));
                % initializing the state
                state_fwd=StateClear(state_fwd);
                while (~state_fwd.JobComplete)
                    % current action selection : (s,a)
                    state_fwd=Model_RL(state_fwd, param_fwd0, map0, 'decision_random');
                    % state update (get reward and next state) : (r,s')
                    state_fwd=StateSpace(state_fwd,map0);  % state index ++
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
            % processing_time=toc(time_start)
            %
            % % sanity check
            % mat_check=[];
            % for j=1:1:N_state
            %     for k=1:1:N_state
            %         for a=1:1:2
            %             if(state_fwd.T(j,a,k)>0.1)
            %                 tmp=[j a k state_fwd.T(j,a,k)];
            %                 mat_check=[mat_check; tmp];
            %             end
            %         end
            %     end
            % end
        end
        
        % (2) phase 2 - intact action, rewards given                
%         state_fwd=myState;
        
        i=0;  cond=1;
        time_start=tic;
        HIST_RWD_FWD=zeros(1,num_max_trial);
        HIST_Q_FWD=zeros(N_state,N_action,num_max_trial);
        HIST_T_FWD_2=zeros(N_state,N_action,N_state,num_max_trial);
        while ((i<num_max_trial)&&(cond))
            i=i+1;
            %     disp(sprintf('- [%d/%d] session...',i,num_max_trial));
            % initializing the state
            state_fwd=StateClear(state_fwd);
            while (~state_fwd.JobComplete)                
                % devaluation
                if(i>=pt_devaluation(ind_dev))
                    map.reward(5)=0;
                end
                % current action selection : (s,a)
                state_fwd=Model_RL(state_fwd, param_fwd, map, 'decision');
                % state update (get reward and next state) : (r,s')
                state_fwd=StateSpace(state_fwd,map);  % state index ++
                % model upate
                state_fwd=Model_RL(state_fwd, param_fwd, map, 'fwd_update');
            end
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
    
    str=sprintf('[Condition-%d] Q(S2,L(line)/R(dotted)) - SARSA(blue) vs. FWD(red)',ind_dev);
    figure('Name',str)
    ind=[1:1:length(squeeze(HIST_Q_SARSA(2,2,:)))];
    plot(ind,squeeze(HIST_Q_SARSA(2,1,:)),'b-',ind,squeeze(HIST_Q_SARSA(2,2,:)),'b-.',ind,squeeze(HIST_Q_FWD(2,1,:)),'r-',ind,squeeze(HIST_Q_FWD(2,2,:)),'r-.')
    
    str=sprintf('[Condition-%d] T(S2,R,S5) - FWD',ind_dev);
    figure('Name',str)
    plot([squeeze(HIST_T_FWD_1(2,2,5,:)); squeeze(HIST_T_FWD_2(2,2,5,:))])

end
disp('-done');

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
