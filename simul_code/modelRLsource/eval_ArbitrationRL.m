function [data_out] = eval_ArbitrationRL(param_in, data_in, mode)
% Load the parameter files from the optimization folder, and then evaluate the model.



%% create MAP
map_opt.transition_prob_seed=[0.5 0.5];
map_opt.reward_seed=[40 20 10 0];
[myMap N_state N_action N_transition]=Model_Map_Init2('sangwan2012b',map_opt);


%% create my arbitrator
myArbitrator=Bayesian_Arbitration_Init(N_state,N_action,N_transition);

%% create my RL
myState=Model_RL_Init(N_state,N_action,N_transition);



%% model parameter - for the functional mode of RL
% SARSA model
param_sarsa.gamma=1.0; % fixed - not actual parameter
%     param_sarsa.alpha=0.15; % learning rate (0.1~0.2)
%     param_sarsa.tau=0.5; % decision focus parameter
% FWD model (for the latent stage)
% param_fwd0.alpha=0.1; % learning rate
%     param_fwd0.tau=param_sarsa.tau; % decision focus parameter
% FWD model
%     param_fwd.alpha=0.15; % learning rate
%     param_fwd.tau=param_sarsa.tau; % decision focus parameter





%% SIMULATION


%% parameter plug-in
% param_in(1): myArbitrator.PE_tolerance_m1 (m1's threshold for zero PE)
% param_in(2): myArbitrator.PE_tolerance_m2 (m2's threshold for zero PE)
% param_in(3): myArbitrator.A_12
% param_in(4): myArbitrator.B_12
% param_in(5): myArbitrator.A_21
% param_in(6): myArbitrator.B_21
% param_in(7): myArbitrator.tau_softmax/param_sarsa.tau/param_fwd.tau : better to fix at 0,2. This should be determined in a way that maintains softmax values in a reasonable scale. Otherwise, this will drive the fitness value!
% param_in(8): % param_sarsa.alpha/param_fwd.alpha 0.01~0.2 to ensure a good "state_fwd.T" in phase 1
param_fixed(1)=1; % 1: fwd-start, 2:sarsa-start
param_fixed(2)=1; % myArbitrator.p
param_fixed(3)=1e-1; % myArbitrator.Time_Step : time constant (1e0 (fast) ~ 1e-2 (slow))

pop_id=1; % fixed
Sum_NegLogLik=0.0;
% arbitrator
myArbitrator.PE_tolerance_m1=param_in(pop_id,1); % defines threshold for zero PE
myArbitrator.PE_tolerance_m2=param_in(pop_id,2); % defines threshold for zero PE
myArbitrator.Time_Step=param_fixed(3); % the smaller, the slower
myArbitrator.A_12=param_in(pop_id,3);  myArbitrator.B_12=param_in(pop_id,4);
myArbitrator.A_21=param_in(pop_id,5);  myArbitrator.B_21=param_in(pop_id,6);
% SARSA
param_sarsa.alpha=param_in(pop_id,8)*2; % learning rate (0.1~0.2)
% FWD
param_fwd.alpha=param_in(pop_id,8);
% pop_id=8 takes no effect.
% Softmax parameter for all models
myArbitrator.tau_softmax=param_in(pop_id,7); % use the same value as sarsa/fwd
param_sarsa.tau=param_in(pop_id,7);
param_fwd.tau=param_in(pop_id,7);
% arbitrator start
myArbitrator.ind_active_model=param_fixed(1);
if(myArbitrator.ind_active_model==1)
    myArbitrator.m1_prob_prev=0.7; % do not use 0.5 which causes oscillation.
    myArbitrator.m2_prob_prev=1-myArbitrator.m1_prob_prev;
else
    myArbitrator.m1_prob_prev=0.3;
    myArbitrator.m2_prob_prev=1-myArbitrator.m1_prob_prev;
end
% non-linear weight : p
myArbitrator.p=param_fixed(2);


%% Devaluation/coningency degradation/risk schedule
% (1) devaluation schedule
pt_devaluation=[1000];%[41:1:80];%[1:1:40]; % second phase 1~80. if devaluation point > num_max_trial, then no devaluation
pt_ref=pt_devaluation(1);
% (2) degradation schedule
pt_degradation=[1000];
% (3) transition prob schedule
pt_prob_chg_max_risk=[1000];%[1:1:30 61:1:90 121:1:150];%[11:1:40]; % second phase 1~80.
pt_prob_chg_min_risk=[1000];%[31:1:60 91:1:120 151:1:180];%[41:1:80]; %[41:1:80]; % second phase 1~80.
pt_prob_chg_reversed=[1000]; % second phase 1~80.
% (4) reward prob schedule
pt_prob_rwd_chg=[1000];
% (5) goal-directed mode (keep changing rwd values)
pt_devaluation_goal=[1000];%[31:1:60 91:1:120 151:1:180];
rwd_aversion_factor=0; % 1: no change, 0: zero reward (ok), -1:opposite (best)

% scenario 3
ss0=[1:1:30]; ss1=[31:1:60]; ss2=[61:1:90]; ss3=[91:1:120];   ss4=[121:1:150];   ss5=[151:1:180];
pt_prob_chg_max_risk=[ss0 ss2 ss4];
pt_prob_chg_min_risk=[ss1 ss3 ss5];
pt_devaluation_goal=[ss1 ss3 ss5];



for ind_dev=[1:1:length(pt_devaluation)]
    
    %% Multiple simulations
    win_fwd=0; win_sarsa=0;
    
    tot_num_sbj=size(data_in,2);
    for ll=1:1:tot_num_sbj % each subject
        
        
        n_tmp=0;
        for kk=1:1:mode.total_simul
            
            Sum_NegLogLik_each_simul=0;
            Sum_NegLogLik_each_simul_eachMode=zeros(1,4);
            Num_occur_eachMode=zeros(1,4); % # of each mode (G,G',H,H')
            OBS=[];
            
            
            
            %% Simulation
            % (1) phase 1 - random action, no reward
            state_fwd=myState;  state_fwd.name='fwd';
            state_sarsa=myState;    state_sarsa.name='sarsa';
            
            map=myMap;  map0=myMap; map0_s=myMap;
            
            
            %% (1) phase 1 - 'pre' training
            
            num_max_trial0=size(data_in{1,ll}.HIST_behavior_info_pre{1,1},1);
            map0.epoch=kk;                map0_s.epoch=kk;
            map0.data=data_in{1,ll}.HIST_behavior_info_pre{1,1};  map0_s.data=data_in{1,ll}.HIST_behavior_info_pre{1,1};
            
            opt_state_space.use_data=mode.experience_sbj_events(1); % use subject data for state-transition            
            if(1)
                % fwd learning
                i=0;  cond=1;
                while ((i<num_max_trial0)&&(cond))
                    i=i+1;
                    map0.trial=i;
                    %             disp(sprintf('- [%d/%d] session...',i,num_max_trial0));
                    % initializing the state
                    [state_fwd map0]=StateClear(state_fwd,map0);
                    while (~state_fwd.JobComplete)
                        % decision
                        if(mode.experience_sbj_events(1)==1)
                            state_fwd=Model_RL(state_fwd, param_fwd, map0, 'decision_behavior_data_save');
                        else
                            state_fwd=Model_RL(state_fwd, param_fwd, map0, 'decision_random');
                        end
                        % state transition
                        [state_fwd map0]=StateSpace_v1(state_fwd,map0,opt_state_space);  % map&state index ++
                        % 1. fwd model update
                        state_fwd=Model_RL(state_fwd, param_fwd, map0, 'fwd_update');
                    end
                end
                % sarsa learning
                i=0;  cond=1;
%                 if(mode.experience_sbj_events(1)==0)                    num_max_trial0=num_max_trial0/2;                end
                while ((i<num_max_trial0)&&(cond))
                    i=i+1;
                    map0_s.trial=i;
                    %             disp(sprintf('- [%d/%d] session...',i,num_max_trial0));
                    % initializing the state
                    [state_sarsa map0_s]=StateClear(state_sarsa,map0_s);
                    while (~state_sarsa.JobComplete)
                        % 0. current action selection : (s,a) - using arbitrator's Q-value
                        if(mode.experience_sbj_events(1)==1)
                            state_sarsa=Model_RL(state_sarsa, param_sarsa, map0_s, 'decision_behavior_data_save');
                        else
                            state_sarsa=Model_RL(state_sarsa, param_sarsa, map0_s, 'decision_random');
                        end
                        % 1. sarsa state update (get reward and next state) : (r,s')
                        [state_sarsa map0_s]=StateSpace_v1(state_sarsa,map0_s,opt_state_space); % map&state index ++
                        % 1. sarsa next action selection : (s',a') - if s' is terminal, then no decision
                        state_sarsa=Model_RL(state_sarsa, param_sarsa, map0_s, 'decision_hypo');
                        % 1. sarsa model upate
                        state_sarsa=Model_RL(state_sarsa, param_sarsa, map0_s, 'sarsa_update');
                    end
                end
                % Q-value synchronization
%                 state_sarsa.Q=state_fwd.Q;
                
            end
            
            %                 disp('- pretraining completed.')
            
            %% (2) phase 2 - intact action, rewards given
            
%             param_sarsa.alpha=param_sarsa.alpha/0.5;
%             param_fwd.alpha=param_fwd.alpha/0.5;
            
            num_max_session=size(data_in{1,ll}.HIST_behavior_info,2);
            
            cond=1;
            myArbitrator_top=myArbitrator;
            mode_data_prev=6;
            block_cond_prev=-1; on_h=-1; on_g=-1;
            
            for ind_sess=1:1:num_max_session
                
                % enter each session data into map
                i=0;
                num_max_trial=size(data_in{1,ll}.HIST_behavior_info{1,ind_sess},1);
                map.epoch=kk;
                map.data=data_in{1,ll}.HIST_behavior_info{1,ind_sess};
                
                %                     disp(sprintf('- Sbj[%d], Simul[%d/%d], Session [%d/%d]...',###,kk,mode.total_simul,ind_sess,num_max_session));
                
                while ((i<num_max_trial)&&(cond))
                    i=i+1;
                    map.trial=i;
                    %                         disp(sprintf('- Simul[%d/%d], Session [%d/%d], Trial [%d/%d]...',kk,mode.total_simul,ind_sess,num_max_session,i,num_max_trial));
                    
                    % initializing the state
                    [state_fwd map]=StateClear(state_fwd,map);
                    [state_sarsa map]=StateClear(state_sarsa,map);
                    
                    % read a mode from data
                    mode_data=map.data(map.trial,18);
                    
                    % mode change whenever there is a visual cue change
                    case_number=0;
                    myArbitrator_top.backward_flag=0;
                    if(mode_data_prev~=mode_data) % if there is any visual cue change, then goal mode
                        myArbitrator_top.backward_flag=1; % after backward update, should set it to 0.
                        myArbitrator_top.ind_active_model=1; % switching the mode
                        myArbitrator_top.m1_prob_prev=0.9; % changing the choice prob accordingly
                        myArbitrator_top.m1_prob=myArbitrator_top.m1_prob_prev;
                        if((mode_data_prev~=-1)&&(mode_data==-1)) % goal mode -> habitual mode
                            myArbitrator_top.ind_active_model=2; % switching the mode
                            myArbitrator_top.m1_prob_prev=0.1; % changing the choice prob accordingly
                            myArbitrator_top.m1_prob=myArbitrator_top.m1_prob_prev;
                        end
                    end
                    if(mode.USE_FWDSARSA_ONLY==1) % forward only
                        myArbitrator_top.ind_active_model=1; % switching the mode                        
                        myArbitrator_top.m1_prob_prev=0.99; % changing the choice prob accordingly
                        myArbitrator_top.m1_prob=myArbitrator_top.m1_prob_prev;
                        myArbitrator_top.Time_Step=1e-2; % extremely slow, so do not switch to the other learner
                    end
                    if(mode.USE_FWDSARSA_ONLY==2) % sarsa only
                        myArbitrator_top.ind_active_model=2; % switching the mode
                        myArbitrator_top.m1_prob_prev=0.01; % changing the choice prob accordingly
                        myArbitrator_top.m1_prob=myArbitrator_top.m1_prob_prev;
                        myArbitrator_top.Time_Step=1e-2; % extremely slow, so do not switch to the other learner
                    end
                    
                    if(mode.out==0) % debug
                        OBS=[OBS [myArbitrator_top.ind_active_model; mode_data; myArbitrator_top.m1_prob; state_fwd.Q(1,1); state_fwd.Q(1,2); state_sarsa.Q(1,1); state_sarsa.Q(1,2)]];
                    end
                    %
                    
                    opt_state_space.use_data=mode.experience_sbj_events(2); % use subject data for state-transition
                    
                    while (((myArbitrator_top.ind_active_model==1)&&(~map.JobComplete))||((myArbitrator_top.ind_active_model==2)&&(~map.JobComplete)))
                        
                        
                        % index synchronization
                        state_fwd.index=map.index;                state_sarsa.index=map.index;
                        
                        
                        %% fwd mode: backward update of fwd model
                        if(myArbitrator_top.backward_flag==1) % if the agent detects context change
                            % (1) revaluation
                            if(mode_data==-1) % reevaluation for habitual mode
                                mode_data_mat=[6 7 8 9];
                            else % reevaluation for goal mode
                                mode_data_mat=mode_data;
                            end
                            map.reward=zeros(N_state,1);    map.reward(mode_data_mat)=map.reward_save(mode_data_mat);
                            % (2) backward update of the fwd model
                            if(mode.USE_BWDupdate_of_FWDmodel==1)
                                state_fwd=Model_RL(state_fwd, param_fwd, map, 'bwd_update');
                            end
                        else
                            % preserve reward values
                            map.reward=map.reward_save;
                        end
                        
                        %% Compute negative log-likelihood : evaluate the arbitrator softmax using sbj's state,action
                        state_data=map.data(map.trial,3+map.index);
                        action_data=map.data(map.trial,6+map.index);
                        % compute real Q-value by merging two Q(fwd,sarsa)
                        myArbitrator_top.Q=...
                            ((myArbitrator_top.m1_prob*state_fwd.Q).^myArbitrator_top.p+...
                            (myArbitrator_top.m2_prob*state_sarsa.Q).^myArbitrator_top.p).^(1/myArbitrator_top.p);
                        var_exp=exp(myArbitrator_top.tau_softmax*myArbitrator_top.Q(state_data,:)); % (N_actionx1)
                        eval_num=log(var_exp(action_data)/sum(var_exp));
                        Sum_NegLogLik=Sum_NegLogLik-eval_num/(mode.total_simul*tot_num_sbj);
                        Sum_NegLogLik_each_simul=Sum_NegLogLik_each_simul-eval_num/(1*1);
                        % block condition of each trial: 1:G(low T uncertainty), 2:G'(high T uncertainty), 3:H(high T uncertainty), 4:H'(low T uncertainty)
                        block_cond=data_in{1,ll}.HIST_block_condition{1,ind_sess}(2,map.trial);
                        Sum_NegLogLik_each_simul_eachMode(block_cond)=Sum_NegLogLik_each_simul_eachMode(block_cond)-eval_num/(mode.total_simul*tot_num_sbj);
                        Num_occur_eachMode(block_cond)=Num_occur_eachMode(block_cond)+1;
                        
                        
                        
                        
                        %% main computation
                        QQ_prev=state_fwd.Q;
                        if(myArbitrator_top.ind_active_model==1) % fwd
                            % 0. current action selection : (s,a) - using arbitrator's Q-value
                            if(mode.experience_sbj_events(2)==1)
                                state_fwd=Model_RL(state_fwd, param_fwd, map, 'decision_behavior_data_save');
                            else
                                state_fwd=Model_RL(state_fwd, param_fwd, map, 'decision_arbitrator', myArbitrator_top, state_sarsa);
                            end
                            
                            % 1. fwd state update (get reward and next state) : (r,s')
                            [state_fwd map]=StateSpace_v1(state_fwd,map,opt_state_space); % map&state index ++
                            % 1. fwd model update
                            state_fwd=Model_RL(state_fwd, param_fwd, map, 'fwd_update');
                            % 2. state synchronization
                            state_sarsa.state_history(state_fwd.index)=state_fwd.state_history(state_fwd.index);
                            state_sarsa.SARSA=state_fwd.SARSA;
                            % 3. sarsa next action selection : (s',a') - if s' is terminal, then no decision
                            state_sarsa=Model_RL(state_sarsa, param_sarsa, map, 'decision_hypo');
                            % 3. sarsa model upate
                            state_sarsa=Model_RL(state_sarsa, param_sarsa, map, 'sarsa_update');
                            % history synchronization
                            myArbitrator_top.state_history=state_fwd.state_history;
                            myArbitrator_top.reward_history=state_fwd.reward_history;
                        end
                        QQ_current=state_fwd.Q;
                        
                        if(myArbitrator_top.ind_active_model==2) % sarsa
                            % 0. current action selection : (s,a) - using arbitrator's Q-value                            
                            if(mode.experience_sbj_events(2)==1)
                                state_sarsa=Model_RL(state_sarsa, param_sarsa, map, 'decision_behavior_data_save');
                            else
                                state_sarsa=Model_RL(state_sarsa, param_sarsa, map, 'decision_arbitrator', myArbitrator_top, state_fwd);
                            end
                            % 1. sarsa state update (get reward and next state) : (r,s')
                            [state_sarsa map]=StateSpace_v1(state_sarsa,map,opt_state_space); % map&state index ++
                            % 1. sarsa next action selection : (s',a') - if s' is terminal, then no decision
                            state_sarsa=Model_RL(state_sarsa, param_sarsa, map, 'decision_hypo');
                            % 1. sarsa model upate
                            state_sarsa=Model_RL(state_sarsa, param_sarsa, map, 'sarsa_update');
                            % 2. state synchronization
                            state_fwd.state_history(state_sarsa.index)=state_sarsa.state_history(state_sarsa.index);
                            state_fwd.SARSA=state_sarsa.SARSA;
                            % 3. fwd model update
                            state_fwd=Model_RL(state_fwd, param_fwd, map, 'fwd_update');
                            % history synchronization
                            myArbitrator_top.state_history=state_sarsa.state_history;
                            myArbitrator_top.reward_history=state_sarsa.reward_history;
                        end
                        
                        % ARBITRATOR: transition rate of the sarsa only
                        [myArbitrator_top, state_fwd, state_sarsa]=Bayesian_Arbitration_v3(myArbitrator_top, state_fwd, state_sarsa, map);
                        
                        
                        % DISPLAY Q-VAL change according to goal change
                        if(mode.DEBUG_Q_VALUE_CHG==1)
                            disp(sprintf('goal: %d->%d, Q_debug=[state|| fwdQ(L)|fwdQ(R) -> fwdQ(L)|fwdQ(R) | arbQ(L)|arbQ(R)]',mode_data_prev, mode_data))
                            %                                 Q_debug=[[1:1:9]' QQ_prev QQ_current myArbitrator_top.Q(:,:)]
                            %                                 myArbitrator_top.m1_prob
                            %                                 state_fwd.SARSA
                            [state_data, action_data, var_exp, (-1)*eval_num]
                            input('pause...')
                        end
                        
                    end
                    
                    % save to previous mode data
                    mode_data_prev=mode_data;
                    
                end
                
                
                
            end
            
            if(mode.simul_process_display==1)
            % display
            %             msg=[sprintf('- processing subject(%d/%d)',ll,tot_num_sbj) repmat('.',1,kk)];
            msg=[sprintf('       - processing subject(%d/%d),simulation(%d/%d),NegLogLik=%04.2f',ll,tot_num_sbj,kk,mode.total_simul,Sum_NegLogLik_each_simul)];
            fprintf(repmat('\b',1,n_tmp));            fprintf(msg);            n_tmp=numel(msg);
            end
        end
        if(mode.simul_process_display==1)        fprintf('\n');     end
    end
    
    
end


%% fitness values (to be maximized)
fitness_val(pop_id)=1e6*1.0/Sum_NegLogLik;
%     fitness_val(pop_id)=(-1.0)*Sum_NegLogLik;
Sum_NegLogLik_val(pop_id)=Sum_NegLogLik;
%     disp(sprintf('    = fitness value: %04.4f (NegLogLik: %f)',max(fitness_val),min(Sum_NegLogLik)));


%% function returns
if(mode.out==0)
    data_out=OBS;
end
if(mode.out==1)
    data_out=Sum_NegLogLik_val(pop_id,1);
end


% function ends
end








