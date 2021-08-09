function [fitness_val]=test_fit_v5(pop,opt)




% opt.bits : coding bits for each param
% opt.minmax : each column=[min; max] for each param
% opt.generation_id : ?-th generation

%% [SW] gene --> parameter conversion
epsilon=1e-300; % never go below 1e-300
[pop_size total_gene_length]=size(pop);
parameter_converted=zeros(pop_size,length(opt.bits));
for pop_id=1:1:pop_size
    for gene_id=1:1:length(opt.bits)
        % breaking into genes
        binary_vector=pop(pop_id,(1+sum(opt.bits(1:gene_id-1))):sum(opt.bits(1:gene_id)));
        % binary 2 decimal
        binary_vector2=[];
        for k=1:1:opt.bits(gene_id)
            binary_vector2=[binary_vector2 sprintf('%d',binary_vector(k))]; % make a binary vector by concatenation
        end
        parameter_decimal2=bin2dec(binary_vector2)/(2^opt.bits(gene_id)-1); % [0,1]
        if(opt.IsCodingLog(gene_id)==0) % uniform coding
            parameter_converted(pop_id,gene_id)=parameter_decimal2*(opt.minmax(2,gene_id)-opt.minmax(1,gene_id))+opt.minmax(1,gene_id); %[min,max]
        else % log coding
            parameter_decimal3=parameter_decimal2*(1.0-epsilon)+epsilon; % [epsilon, 1]
            parameter_decimal4=(-1.0)*log(parameter_decimal3); % [0, -log(epsilon)] - log scaling
            parameter_decimal4=parameter_decimal4/(-1.0*log(epsilon)); % [0,1]
            parameter_converted(pop_id,gene_id)=parameter_decimal4*(opt.minmax(2,gene_id)-opt.minmax(1,gene_id))+opt.minmax(1,gene_id); %[min,max]
        end
    end
end
disp(sprintf('- [%s][%d-th generation] gene -> parameter conversion done.',opt.EXP_NAME,opt.generation_id));

fitness_val=zeros(1,pop_size);
Sum_NegLogLik_val=zeros(pop_size,1);








%% create MAP
map_opt.transition_prob_seed=[0.5 0.5];
map_opt.reward_seed=[40 20 10 0];
[myMap N_state N_action N_transition]=Model_Map_Init2('sangwan2012',map_opt);


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


% %% [SW] Load behaviour data (pseudo data)
% state_column=[1 2 3];   action_column=[4 5];
% % eval(['load ',[opt.path_ext opt.str_path_div 'behavior_data_glascher' opt.str_path_div 'data.mat']]);
% eval(['load ',[opt.path_ext opt.str_path_div 'behavior_data_pseudo' opt.str_path_div 'data_save.mat']]);
% % and attach it to the map
% ind_included=opt.ind_included;
% if(opt.ind_included(1)==0) % for glascher dataset only
%     ind_included=[1:1:20];
% end
% data0=cell(1,length(ind_included));
% for j=1:1:length(ind_included)
%     data0{1,j}=data_behavior{1,ind_included(j)}.data;
% end
% myMap.data=data0;



%note: Time_Step: 1e-1 ~ 1e-2
parameter_converted=[0.7,1e-2,2.17142857142857,2.77777777777778,1.88888888888889,27.3015873015873,0.12,1,0.100000000000000,2,10];
parameter_converted(7)=0.05; % 0.01~0.2 to ensure a good "state_fwd.T" in phase 1
parameter_converted(9)=0.5;  % better to fix at 0,5. This should be determined in a way that maintains softmax values in a reasonable scale. Otherwise, this will drive the fitness value!





%% SIMULATION
for pop_id=1:1:pop_size
    
    
    disp(sprintf('  :: processing %d-th generation, %d-th population...',opt.generation_id,pop_id));
    
    
    %% [SW] Parameter insertion/change
    Sum_NegLogLik=0.0;
    % arbitrator
    myArbitrator.PE_tolerance_m1=parameter_converted(pop_id,1); % defines threshold for zero PE
    myArbitrator.PE_tolerance_m2=parameter_converted(pop_id,11); % defines threshold for zero PE
    myArbitrator.Time_Step=parameter_converted(pop_id,2); % the smaller, the slower
    myArbitrator.A_12=parameter_converted(pop_id,3);  myArbitrator.B_12=parameter_converted(pop_id,4);
    myArbitrator.A_21=parameter_converted(pop_id,5);  myArbitrator.B_21=parameter_converted(pop_id,6);
    myArbitrator.p=parameter_converted(pop_id,8);
    % SARSA
    param_sarsa.alpha=parameter_converted(pop_id,7); % learning rate (0.1~0.2)
    % FWD
    param_fwd.alpha=parameter_converted(pop_id,7);
    % pop_id=8 takes no effect.
    % Softmax parameter for all models
    myArbitrator.tau_softmax=parameter_converted(pop_id,9); % use the same value as sarsa/fwd
    param_sarsa.tau=parameter_converted(pop_id,9);
    param_fwd.tau=parameter_converted(pop_id,9);
    % arbitrator start
    myArbitrator.ind_active_model=parameter_converted(pop_id,10);
    if(myArbitrator.ind_active_model==1)
        myArbitrator.m1_prob_prev=0.7; % do not use 0.5 which causes oscillation.
        myArbitrator.m2_prob_prev=1-myArbitrator.m1_prob_prev;
    else
        myArbitrator.m1_prob_prev=0.3;
        myArbitrator.m2_prob_prev=1-myArbitrator.m1_prob_prev;
    end
    
    
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
        
        for ll=1:1:length(ind_included) % each subject
            
            for kk=1:1:total_simul
                
                Sum_NegLogLik_each_simul=0;
%                 OBS=[];

                
                
                %% Simulation
                % (1) phase 1 - random action, no reward
                state_fwd=myState;
                state_sarsa=myState;
                
                map=myMap;  map0=myMap; % no reward                
                                
                
                %% (1) phase 1 - 'pre' training
                
                i=0;  cond=1;                                
                num_max_trial0=size(HIST_behavior_info_pre{1,1},1);
                map0.epoch=kk;                
                map0.data=HIST_behavior_info_pre{1,1};
                
                opt_state_space.use_data=0; % do not use subject data
                if(1) 
                    while ((i<num_max_trial0)&&(cond))
                        i=i+1;
                        map0.trial=i;
                        %             disp(sprintf('- [%d/%d] session...',i,num_max_trial0));
                        % initializing the state
                        [state_fwd map0]=StateClear(state_fwd,map0);
                        while (~state_fwd.JobComplete)                            
                            % decision
%                             state_fwd=Model_RL(state_fwd, param_fwd, map0, 'decision_behavior_data_save');
                            state_fwd=Model_RL(state_fwd, param_fwd, map0, 'decision_random');
                            % state transition
                            [state_fwd map0]=StateSpace_v1(state_fwd,map0,opt_state_space);  % map&state index ++                            
                            % 1. fwd model update
                            state_fwd=Model_RL(state_fwd, param_fwd, map0, 'fwd_update');                           
                        end                        
                    end
                    % Q-value synchronization
                    state_sarsa.Q=state_fwd.Q;
                end
         
                disp('- pretraining completed.')
                
                %% (2) phase 2 - intact action, rewards given
                %         state_fwd=myState;
                
                num_max_session=size(HIST_behavior_info,2);
                
                cond=1;
                myArbitrator_top=myArbitrator;                
                mode_data_prev=6;
                
                for ind_sess=1:1:num_max_session
                    
                    % enter each session data into map
                    i=0;
                    num_max_trial=size(HIST_behavior_info{1,ind_sess},1);
                    map.epoch=kk;
                    map.data=HIST_behavior_info{1,ind_sess};
                    
                    disp(sprintf('- Simul[%d/%d], Session [%d/%d]...',kk,total_simul,ind_sess,num_max_session));
                                        
                    while ((i<num_max_trial)&&(cond))
                        i=i+1;
                        map.trial=i;
%                         disp(sprintf('- Simul[%d/%d], Session [%d/%d], Trial [%d/%d]...',kk,total_simul,ind_sess,num_max_session,i,num_max_trial));
                        
                        % initializing the state
                        [state_fwd map]=StateClear(state_fwd,map);
                        [state_sarsa map]=StateClear(state_sarsa,map);
                        
                        % read a mode from data
                        mode_data=map.data(map.trial,18);
                        
                        % mode change whenever there is a visual cue change
                        case_number=0;
                        if(mode_data_prev~=mode_data) % if there is any visual cue change, then goal mode
                            myArbitrator_top.backward_flag=1; % after backward update, should set it to 0.
                            myArbitrator_top.ind_active_model=1; % switching the mode
                            myArbitrator_top.m1_prob_prev=0.9; % changing the choice prob accordingly
                            if((mode_data_prev~=-1)&&(mode_data==-1)) % goal mode -> habitual mode
                                myArbitrator_top.ind_active_model=2; % switching the mode
                                myArbitrator_top.m1_prob_prev=0.1; % changing the choice prob accordingly
                            end
                        end
                        if(USE_SARSA_ONLY==1)
                            myArbitrator_top.ind_active_model=2; % switching the mode
                            myArbitrator_top.m1_prob_prev=0.01; % changing the choice prob accordingly
                        end
                        
%                         OBS=[OBS [myArbitrator_top.ind_active_model; mode_data; myArbitrator_top.m1_prob]];
                        %
                        
                        opt_state_space.use_data=0; % do not use subject data for state-transition
                                                                        
                        while (((myArbitrator_top.ind_active_model==1)&&(~map.JobComplete))||((myArbitrator_top.ind_active_model==2)&&(~map.JobComplete)))
                            
                            
                            % Compute negative log-likelihood                                                        
                            state_data=map.data(map.trial,3+map.index);
                            action_data=map.data(map.trial,6+map.index);                                                        
                            var_exp=exp(myArbitrator_top.tau_softmax*myArbitrator_top.Q(state_data,:)); % (N_actionx1)
                            eval_num=log(var_exp(action_data)/sum(var_exp));
                            Sum_NegLogLik=Sum_NegLogLik-eval_num/total_simul;                            
                            Sum_NegLogLik_each_simul=Sum_NegLogLik_each_simul-eval_num;                            
                            %
%                             ZZZ=[ZZZ; [state_data, action_data, var_exp, (-1)*eval_num]];
                                                        
                            % index synchronization
                            state_fwd.index=map.index;                state_sarsa.index=map.index;

                            QQ_prev=state_fwd.Q;
                            if(myArbitrator_top.ind_active_model==1) % fwd
                                % backward update
                                if(myArbitrator_top.backward_flag==1) % if the agent detects context change
                                    % (1) revaluation
                                    map.reward=zeros(N_state,1); 
                                    map.reward(mode_data)=map.reward_save(mode_data);
                                    % (2) backward update                                    
                                    if(USE_BWDupdate_of_FWDmodel==1)
                                        state_fwd=Model_RL(state_fwd, param_fwd, map, 'bwd_update');
                                    end
                                else
                                    % preserve reward values
                                    map.reward=map.reward_save;
                                end
                                % 0. current action selection : (s,a) - using arbitrator's Q-value
                                state_fwd=Model_RL(state_fwd, param_fwd, map, 'decision_arbitrator', myArbitrator_top);
%                                 state_fwd=Model_RL(state_fwd, param_fwd, map, 'decision_behavior_data_save');
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
                                state_sarsa=Model_RL(state_sarsa, param_sarsa, map, 'decision_arbitrator', myArbitrator_top);
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
                            if(DEBUG_Q_VALUE_CHG==1)
                                disp(sprintf('goal: %d->%d, Q_debug=[state|| fwdQ(L)|fwdQ(R) -> fwdQ(L)|fwdQ(R) | arbQ(L)|arbQ(R)]',mode_data_prev, mode_data))
%                                 Q_debug=[[1:1:9]' QQ_prev QQ_current myArbitrator_top.Q(:,:)]
%                                 myArbitrator_top.m1_prob
%                                 state_fwd.SARSA
%                                 [state_data, action_data, var_exp, (-1)*eval_num]
                                input('pause...')                                
                            end
                            
                        end
                        
                        % save to previous mode data
                        mode_data_prev=mode_data;
                        
                    end
                    
                    
                    
                end
                
               
            end
        end
        
        
    end
    
    
    % fitness values (to be maximized)
    fitness_val(pop_id)=1e6*1.0/Sum_NegLogLik;
    %     fitness_val(pop_id)=(-1.0)*Sum_NegLogLik;
    Sum_NegLogLik_val(pop_id)=Sum_NegLogLik;
    disp(sprintf('    = fitness value: %04.4f (NegLogLik: %f)',max(fitness_val),min(Sum_NegLogLik)));
    
end

%% file save
file_name=sprintf('generation%03d_bestfit(%04.4f).mat',opt.generation_id,max(fitness_val));
eval(['save ',[opt.path_ext opt.str_path_div opt.TARGET_FOLDER_NAME opt.str_path_div opt.EXP_NAME opt.str_path_div file_name ' fitness_val Sum_NegLogLik_val parameter_converted opt' ]]);


end
