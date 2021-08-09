clear all
close all

%% Option
N_state=(1+4+4*4);
N_action=2; % left/right
N_transition=3; % inital, choice1, choice2

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
myState.IsTerminal=zeros(N_state,1);    myState.IsTerminal(6:21)=1;
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
param_sarsa.alpha=0.2; % learning rate
param_sarsa.tau=0.1; % decision focus parameter
% FWD model
param_fwd.alpha=0.2; % learning rate
param_fwd.tau=param_sarsa.tau; % decision focus parameter


%% MAP: transition prob. for the left(1,1)/right(1,2) action
prob_seed=0.7;
prob_seed_mat=[prob_seed 1-prob_seed];
myMap.action(1,1).prob=zeros(N_state,N_state);
myMap.action(1,2).prob=zeros(N_state,N_state);
% 1. transition for the left action
myMap.action(1,1).prob(1,2:3)=prob_seed_mat;
myMap.action(1,1).prob(2,6:7)=prob_seed_mat;
myMap.action(1,1).prob(3,10:11)=prob_seed_mat;
myMap.action(1,1).prob(4,14:15)=prob_seed_mat;
myMap.action(1,1).prob(5,18:19)=prob_seed_mat;
% 2. transition for the right action
myMap.action(1,2).prob(1,4:5)=prob_seed_mat;
myMap.action(1,2).prob(2,8:9)=prob_seed_mat;
myMap.action(1,2).prob(3,12:13)=prob_seed_mat;
myMap.action(1,2).prob(4,16:17)=prob_seed_mat;
myMap.action(1,2).prob(5,20:21)=prob_seed_mat;
% 3. get a connection matrix
myMap.action(1,1).connection=double(myMap.action(1,1).prob&ones(N_state,N_state));
myMap.action(1,2).connection=double(myMap.action(1,2).prob&ones(N_state,N_state));

%% MAP: reward matrix
myMap.reward=zeros(N_state,1);
myMap.reward(6:21)=[10 0 0 10 0 10 0 25 25 0 10 0 0 10 0 25]; % [0 10 10 0 10 0 25 0 0 25 0 10 10 0 25 0];

%% sanity check - confirmed
% num_test_trial=50000;
% HIST_state=zeros(N_transition,num_test_trial);
% HIST_action=zeros(N_transition,num_test_trial);
% HIST_reward=zeros(N_transition,num_test_trial);
% tic;
% for j=1:1:num_test_trial    
%     action_seq=double(rand(1,2)>=0.5)+1;
%     for i=1:1:length(action_seq)
%         % make a choice
%         myState.action_history(i)=action_seq(i);
%         % state update
%         myState=StateSpace(myState,myMap);
%     end
%     % save the history
%     HIST_state(:,j)=myState.state_history;
%     HIST_action(:,j)=myState.action_history;
%     HIST_reward(:,j)=myState.reward_history;
%     % initializing the state
%     myState=StateClear(myState);
% end
% computing_time=toc
% % sanity check
%  Z=sort(HIST_state(3,:));
%  Z_prob=zeros(1,N_state);
%  for k=1:1:N_state
%      [tmp dd]=find(HIST_state(3,:)==k);
%      Z_prob(k)=sum(length(dd))/num_test_trial;
%  end
 
%% simulation
state_sarsa=myState;
state_fwd=myState;
map=myMap;  map0=myMap; % no reward
map0.reward=zeros(N_state,1);

%% SARSA model
num_max_trial=80;
i=0;  cond=1;
time_start=tic;
HIST_RWD_SARSA=zeros(1,num_max_trial);
HIST_Q_SARSA=zeros(5,2,num_max_trial);
while ((i<num_max_trial)&&(cond))
    i=i+1;
    disp(sprintf('- [%d/%d] session...',i,num_max_trial));
    % initializing the state
    state_sarsa=StateClear(state_sarsa);    
    for j=1:1:(N_transition-1)                
        % current action selection : (s,a)
        state_sarsa=Model_RL(state_sarsa, param_sarsa, map, 'decision');
        % state update (get reward and next state) : (r,s')
        state_sarsa=StateSpace(state_sarsa,map);  % state index ++     
        % next action selection : (s',a') - if s' is terminal, then no decision
        state_sarsa=Model_RL(state_sarsa, param_sarsa, map, 'decision_hypo');
        % model upate
        state_sarsa=Model_RL(state_sarsa, param_sarsa, map, 'sarsa_update');        
    end
    HIST_RWD_SARSA(i)=state_sarsa.reward_history(N_transition);
    HIST_Q_SARSA(:,:,i)=state_sarsa.Q(1:5,1:2);
%  step check
%     reply = input('- continue? [y/n]: ', 's');    
%     if(reply=='n')
%         cond=0;
%     end    
end
processing_time=toc(time_start)

%% FWD model
% (1) phase 1 - random action, no reward
num_max_trial=80;
i=0;  cond=1;
time_start=tic;
HIST_T_FWD_1=zeros(5,2,21,num_max_trial);
while ((i<num_max_trial)&&(cond))
    i=i+1;
    disp(sprintf('- [%d/%d] session...',i,num_max_trial));
    % initializing the state    
    state_fwd=StateClear(state_fwd);
    for j=1:1:(N_transition-1)                        
        % current action selection : (s,a)
        state_fwd=Model_RL(state_fwd, param_fwd, map0, 'decision_random');
        % state update (get reward and next state) : (r,s')
        state_fwd=StateSpace(state_fwd,map0);  % state index ++     
        % model upate
        state_fwd=Model_RL(state_fwd, param_fwd, map0, 'fwd_update');
    end        
    HIST_T_FWD_1(:,:,:,i)=state_fwd.T(1:5,1:2,1:21);
%  step check
%     reply = input('- continue? [y/n]: ', 's');    
%     if(reply=='n')
%         cond=0;
%     end    
end
processing_time=toc(time_start)

% sanity check
mat_check=[];
for j=1:1:N_state
    for k=1:1:N_state
        for a=1:1:2
            if(state_fwd.T(j,a,k)>0.1)
                tmp=[j a k state_fwd.T(j,a,k)];
                mat_check=[mat_check; tmp];
            end
        end
    end
end

% (2) phase 2 - intact action, rewards given
num_max_trial=80;
i=0;  cond=1;
time_start=tic;
HIST_RWD_FWD=zeros(1,num_max_trial);
HIST_Q_FWD=zeros(5,2,num_max_trial);
HIST_T_FWD_2=zeros(5,2,21,num_max_trial);
while ((i<num_max_trial)&&(cond))
    i=i+1;
    disp(sprintf('- [%d/%d] session...',i,num_max_trial));
    % initializing the state    
    state_fwd=StateClear(state_fwd);
    for j=1:1:(N_transition-1)                        
        % current action selection : (s,a)
        state_fwd=Model_RL(state_fwd, param_fwd, map, 'decision');
        % state update (get reward and next state) : (r,s')
        state_fwd=StateSpace(state_fwd,map);  % state index ++     
        % model upate
        state_fwd=Model_RL(state_fwd, param_fwd, map, 'fwd_update');
    end
    HIST_RWD_FWD(i)=state_fwd.reward_history(N_transition);
    HIST_Q_FWD(:,:,i)=state_fwd.Q(1:5,1:2);
    HIST_T_FWD_2(:,:,:,i)=state_fwd.T(1:5,1:2,1:21);
%  step check
%     reply = input('- continue? [y/n]: ', 's');    
%     if(reply=='n')
%         cond=0;
%     end    
end
processing_time=toc(time_start)

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

figure('Name','Q(S1,R) - SARSA vs. FWD')
ind=[1:1:length(squeeze(HIST_Q_SARSA(1,2,:)))];
plot(ind,squeeze(HIST_Q_SARSA(1,2,:)),'g',ind,squeeze(HIST_Q_FWD(1,2,:)),'r')

figure('Name','Q(S4,L) - SARSA vs. FWD')
ind=[1:1:length(squeeze(HIST_Q_SARSA(4,2,:)))];
plot(ind,squeeze(HIST_Q_SARSA(4,1,:)),'g',ind,squeeze(HIST_Q_FWD(4,1,:)),'r')

figure('Name','T(S1,R,S4) - FWD')
plot([squeeze(HIST_T_FWD_1(1,2,4,:)); squeeze(HIST_T_FWD_2(1,2,4,:))])
