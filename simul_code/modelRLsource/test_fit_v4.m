function [fitness_val]=test_fit_v4(pop,opt)




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
% [myMap N_state N_action N_transition]=Model_Map_Init('daw2005');
% [myMap N_state N_action N_transition]=Model_Map_Init('glascher2010');
[myMap N_state N_action N_transition]=Model_Map_Init('sangwan2012');

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
param_fwd0.alpha=0.1; % learning rate
%     param_fwd0.tau=param_sarsa.tau; % decision focus parameter
% FWD model
%     param_fwd.alpha=0.15; % learning rate
%     param_fwd.tau=param_sarsa.tau; % decision focus parameter


%% [SW] Load behaviour data (pseudo data)
state_column=[1 2 3];   action_column=[4 5];
% eval(['load ',[opt.path_ext opt.str_path_div 'behavior_data_glascher' opt.str_path_div 'data.mat']]);
eval(['load ',[opt.path_ext opt.str_path_div 'behavior_data_pseudo' opt.str_path_div 'data_save.mat']]);
% and attach it to the map
ind_included=opt.ind_included;
if(opt.ind_included(1)==0) % for glascher dataset only
    ind_included=[1:1:20];
end
data0=cell(1,length(ind_included));
for j=1:1:length(ind_included)
    data0{1,j}=data_behavior{1,ind_included(j)}.data;
end
myMap.data=data0;



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
    param_fwd.alpha=parameter_converted(pop_id,7); % learning rate set to be the same as sarsa.
    % pop_id=8 takes no effect.
    % Softmax parameter for all models
    myArbitrator.tau_softmax=parameter_converted(pop_id,9); % use the same value as sarsa/fwd
    param_sarsa.tau=parameter_converted(pop_id,9);
    param_fwd.tau=parameter_converted(pop_id,9);
    param_fwd0.tau=parameter_converted(pop_id,9);
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
    num_max_trial=180;
    % (1) devaluation schedule
    pt_devaluation=[1000];%[41:1:80];%[1:1:40]; % second phase 1~80. if devaluation point > num_max_trial, then no devaluation
    pt_ref=pt_devaluation(1);
    % (2) degradation schedule
    pt_degradation=[1000];
    % (3) transition prob schedule
    pt_prob_chg_max_risk=[1:1:30 61:1:90 121:1:150];%[11:1:40]; % second phase 1~80.
    pt_prob_chg_min_risk=[31:1:60 91:1:120 151:1:180];%[41:1:80]; %[41:1:80]; % second phase 1~80.
    pt_prob_chg_reversed=[1000]; % second phase 1~80.
    % (4) reward prob schedule
    pt_prob_rwd_chg=[1000];
    % (5) goal-directed mode (keep changing rwd values)
    pt_devaluation_goal=[31:1:60 91:1:120 151:1:180];
    rwd_aversion_factor=0; % 1: no change, 0: zero reward (ok), -1:opposite (best)
    
    % scenario 3
    ss0=[1:1:30]; ss1=[31:1:60]; ss2=[61:1:90]; ss3=[91:1:120];   ss4=[121:1:150];   ss5=[151:1:180];
    pt_prob_chg_max_risk=[ss0 ss2 ss4];
    pt_prob_chg_min_risk=[ss1 ss3 ss5];
    pt_devaluation_goal=[ss1 ss3 ss5];
    
    total_simul=1;
    num_max_trial0=size(data0{1,1},1)/2; % for latent learning
%     num_max_trial=size(data0{1,1},1)/2;
    num_max_trial=180;
    
    %     response_rate_distal_sarsa=zeros(total_simul,2); response_rate_proximal_sarsa=zeros(total_simul,2);
    %     response_rate_distal_fwd=zeros(total_simul,2); response_rate_proximal_fwd=zeros(total_simul,2);
    % HIST_TRANSITION_RATE_SARSA=zeros(total_simul,num_max_trial);
    % HIST_PROB_CHOICE_SARSA=zeros(total_simul,num_max_trial);
    % HIST_PE_MEAN_SARSA=zeros(myArbitrator.K,num_max_trial);
    % HIST_PE_VAR_SARSA=zeros(myArbitrator.K,num_max_trial);
    % HIST_invFano_SARSA=zeros(myArbitrator.K,num_max_trial);
    % HIST_TRANSITION_RATE_FWD=zeros(total_simul,num_max_trial);
    % HIST_PROB_CHOICE_FWD=zeros(total_simul,num_max_trial);
    % HIST_PE_MEAN_FWD=zeros(myArbitrator.K,num_max_trial);
    % HIST_PE_VAR_FWD=zeros(myArbitrator.K,num_max_trial);
    % HIST_invFano_FWD=zeros(myArbitrator.K,num_max_trial);
    % HIST_MODEL_CHOICE=zeros(total_simul,num_max_trial);
    % ZZZ=[];
    % ZZZ1=[];
    % ZZZ2=[];
    % YYY1=[];
    % YYY2=[];
    for ind_dev=[1:1:length(pt_devaluation)]
        
        %% Multiple simulations
        win_fwd=0; win_sarsa=0;
        
        for ll=1:1:length(ind_included)
            for kk=1:1:total_simul
                
                
                
                
                %% Simulation
                % (1) phase 1 - random action, no reward
                state_fwd=myState;
                state_sarsa=myState;
                
                map=myMap;  map0=myMap; % no reward
                map0.reward=zeros(N_state,1);
                map0.epoch=kk;   map.epoch=kk;
                
                i=0;  cond=1;
                % time_start=tic;
                %         HIST_T_FWD_1=zeros(N_state,N_action,N_state,num_max_trial0);
                if(0) % latent learning
                    while ((i<num_max_trial0)&&(cond))
                        i=i+1;
                        map0.trial=i;
                        %             disp(sprintf('- [%d/%d] session...',i,num_max_trial0));
                        % initializing the state
                        [state_fwd map0]=StateClear(state_fwd,map0);
                        while (~state_fwd.JobComplete)
                            % current action selection : (s,a)
                            %                         state_fwd=Model_RL(state_fwd, param_fwd0, map0, 'decision_behavior_data');
                            state_fwd=Model_RL(state_fwd, param_fwd0, map0, 'decision_random');
                            % state update (get reward and next state) : (r,s')
                            [state_fwd map0]=StateSpace(state_fwd,map0);  % state index ++
                            % model upate
                            state_fwd=Model_RL(state_fwd, param_fwd0, map0, 'fwd_update');
                        end
                        %                 HIST_T_FWD_1(:,:,:,i)=state_fwd.T;
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
                %             time_start=tic;
                myArbitrator_top=myArbitrator;
                
                %         HIST_RWD_SARSA=zeros(1,num_max_trial);
                %         HIST_Q_SARSA=zeros(N_state,2,num_max_trial);
                %         HIST_RWD_FWD=zeros(1,num_max_trial);
                %         HIST_Q_FWD=zeros(N_state,N_action,num_max_trial);
                %         HIST_T_FWD_2=zeros(N_state,N_action,N_state,num_max_trial);
                %         HIST_MODEL_CHOICE2=[];
                %         ZZZ3=[];
                
                while ((i<num_max_trial)&&(cond))
                    i=i+1;
                    map.trial=num_max_trial0+i;
                    %     disp(sprintf('- [%d/%d] session...',i,num_max_trial));
                    % initializing the state
                    [state_fwd map]=StateClear(state_fwd,map);
                    [state_sarsa map]=StateClear(state_sarsa,map);
                    
                    %             HIST_MODEL_CHOICE(kk,i)=myArbitrator_top.ind_active_model;
                    while (((myArbitrator_top.ind_active_model==1)&&(~map.JobComplete))||((myArbitrator_top.ind_active_model==2)&&(~map.JobComplete)))
                        
                        
                       
                        
                        % map save
                        if(map.index==1) % distal
                            map=data_behavior{1,ll}.map_1{1,i};
                        end
                        if(map.index==2) % proximal
                            map=data_behavior{1,ll}.map_2{1,i};
                        end
                        
                        
                        %% Compute negative log-likelihood
%                         for ll=1:1:length(ind_included) % compute map.trial's log-lik for all subjects
                            state_data=data0{1,ll}(map.trial,state_column(map.index));
                            action_data=data0{1,ll}(map.trial,action_column(map.index));
                            var_exp=exp(myArbitrator_top.tau_softmax*myArbitrator_top.Q(state_data,:)); % (N_actionx1)
                            Sum_NegLogLik=Sum_NegLogLik-log(var_exp(action_data)/sum(var_exp))/total_simul;
%                         end
                        %%
                        
                        % index synchronization
                        state_fwd.index=map.index;                state_sarsa.index=map.index;
                        
                        if(myArbitrator_top.ind_active_model==1) % fwd                            
                            % 0. current action selection : (s,a) - using arbitrator's Q-value
                            state_fwd=Model_RL(state_fwd, param_fwd, map, 'decision_arbitrator', myArbitrator_top);
                            % 1. fwd state update (get reward and next state) : (r,s')
                            [state_fwd map]=StateSpace(state_fwd,map);  % map&state index ++                            
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
                        
                        if(myArbitrator_top.ind_active_model==2) % sarsa
                            % 0. current action selection : (s,a) - using arbitrator's Q-value
                            state_sarsa=Model_RL(state_sarsa, param_sarsa, map, 'decision_arbitrator', myArbitrator_top);
                            % 1. sarsa state update (get reward and next state) : (r,s')
                            [state_sarsa map]=StateSpace(state_sarsa,map);  % map&state index ++
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
                        
                        
                        %                 ZZZ=[ZZZ [map.index; state_fwd.index; state_sarsa.index; myArbitrator_top.index; myArbitrator_top.temp]];
                        %                     ZZZ1=[ZZZ1 [i; myArbitrator_top.m1_wgt'*myArbitrator_top.m1_inv_Fano; sum(myArbitrator_top.m1_inv_Fano)]];
                        %                     ZZZ2=[ZZZ2 [i; myArbitrator_top.m2_wgt'*myArbitrator_top.m2_inv_Fano; sum(myArbitrator_top.m2_inv_Fano)]];
                        %                     YYY1=[YYY1 [i; myArbitrator_top.transition_rate12]];
                        %                     YYY2=[YYY2 [i; myArbitrator_top.transition_rate21]];
                        
                        
                    end
                    
                    %             ZZZ3=[ZZZ3 state_sarsa.RPE_history];
                    %             HIST_TRANSITION_RATE_SARSA(kk,i)=myArbitrator_top.transition_rate21;
                    %             HIST_PROB_CHOICE_SARSA(kk,i)=myArbitrator_top.m2_prob;
                    %             HIST_PE_MEAN_SARSA(:,i)=HIST_PE_MEAN_SARSA(:,i)+myArbitrator_top.m2_mean/total_simul;
                    %             HIST_PE_VAR_SARSA(:,i)=HIST_PE_VAR_SARSA(:,i)+myArbitrator_top.m2_var/total_simul;
                    %             HIST_invFano_SARSA(:,i)=HIST_invFano_SARSA(:,i)+myArbitrator_top.m2_inv_Fano/total_simul;
                    %             HIST_RWD_SARSA(i)=state_sarsa.reward_history(state_sarsa.index);
                    %             HIST_Q_SARSA(:,:,i)=state_sarsa.Q(1:N_state,1:N_action);
                    %
                    %             HIST_TRANSITION_RATE_FWD(kk,i)=myArbitrator_top.transition_rate12;
                    %             HIST_PROB_CHOICE_FWD(kk,i)=myArbitrator_top.m1_prob;
                    %             HIST_PE_MEAN_FWD(:,i)=HIST_PE_MEAN_FWD(:,i)+myArbitrator_top.m1_mean/total_simul;
                    %             HIST_PE_VAR_FWD(:,i)=HIST_PE_VAR_FWD(:,i)+myArbitrator_top.m1_var/total_simul;
                    %             HIST_invFano_FWD(:,i)=HIST_invFano_FWD(:,i)+myArbitrator_top.m1_inv_Fano/total_simul;
                    %             HIST_RWD_FWD(i)=state_fwd.reward_history(state_fwd.index);
                    %             HIST_Q_FWD(:,:,i)=state_fwd.Q;
                    %             HIST_T_FWD_2(:,:,:,i)=state_fwd.T;
                    % response rate
                    %             if(i>=pt_ref)
                    %                 if(state_fwd.action_history(1)==1)
                    %                     response_rate_distal_fwd(kk,ind_dev)=response_rate_distal_fwd(kk,ind_dev)+1;
                    %                 end
                    %                 if(state_fwd.action_history(2)==2)
                    %                     response_rate_proximal_fwd(kk,ind_dev)=response_rate_proximal_fwd(kk,ind_dev)+1;
                    %                 end
                    %             end
                    %  step check
                    %     reply = input('- continue? [y/n]: ', 's');
                    %     if(reply=='n')
                    %         cond=0;
                    %     end
                end
                %             processing_time=toc(time_start);
                
                %         % determine who learns fast (sum of rwd in early stage)
                %         pt1=1;
                %         pt2=round(num_max_trial/3);
                %         if(sum(HIST_RWD_FWD(pt1:pt2))>sum(HIST_RWD_SARSA(pt1:pt2)))
                %             win_fwd=win_fwd+1;
                %         else
                %             win_sarsa=win_sarsa+1;
                %         end
                %         disp(sprintf('- condition[%d],simulation[%d/%d] - (fwd_win|sarsa_win)=(%d|%d)',ind_dev,kk,total_simul,win_fwd,win_sarsa))
                
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
