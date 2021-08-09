


clear all
close all

%% subject parameter
EXP_NAME0='ANY_start_18_Sbj';
SUB_ARRAY=[0];
NUM_FILE_TO_PROCESS_FOR_EACH_SBJ=10; %  # of models to simulate for each subject.
model_id_to_process=2;
num_repeat_simulation=1; % # of repetition for each model.

% model start assumption
Force_FWD_start=0;
use_FWD_start_model=1; % 0: read all, 1: read fwd-start model only, 2: read sarsa-start model only


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


% scenario 1
% pt_prob_chg_max_risk=[1000];
% pt_prob_chg_min_risk=[91:1:180];
% pt_devaluation_goal=[91:1:180];

% scenario 2
% pt_prob_chg_max_risk=[91:1:180];
% pt_prob_chg_min_risk=[1000];
% pt_devaluation_goal=[1000];

% scenario 3
ss0=[1:1:30]; ss1=[31:1:60]; ss2=[61:1:90]; ss3=[91:1:120];   ss4=[121:1:150];   ss5=[151:1:180];     
pt_prob_chg_max_risk=[ss0 ss2 ss4];
pt_prob_chg_min_risk=[ss1 ss3 ss5];
pt_devaluation_goal=[ss1 ss3 ss5];

% scenario 4
% ss0=[1:1:30]; ss1=[31:1:60]; ss2=[61:1:90]; ss3=[91:1:120];   ss4=[121:1:150];   ss5=[151:1:180];     
% pt_prob_chg_max_risk=[ss1 ss3 ss5];
% pt_prob_chg_min_risk=[ss0 ss2 ss4];
% pt_devaluation_goal=[ss0 ss2 ss4];


%% variables
X_all=[]; % behavior profile vectors
Y_all=[]; % [parameter vector, subject index]
X_in_fwd=[];  X_in_sarsa=[];
X_in_fwd_mean=[];  X_in_sarsa_mean=[];
X_in_fwd_var=[];  X_in_sarsa_var=[];
HIST_best_SumNegLogLik=[];
% collection_behavior_profile=[];

%% data variable to which all the bevahior are to be saved
% variables:: data_behavior, data_model
data_behavior=cell(1,NUM_FILE_TO_PROCESS_FOR_EACH_SBJ);
data_model.id=model_id_to_process;
data_model.pt_prob_chg_max_risk=pt_prob_chg_max_risk;
data_model.pt_prob_chg_min_risk=pt_prob_chg_min_risk;
data_model.pt_devaluation_goal=pt_devaluation_goal;



sss_ind=0;
for subj_ind=SUB_ARRAY
    sss_ind=sss_ind+1;

    for repeat_ind=1:1:num_repeat_simulation
        disp(sprintf('[simulation repeatition %d/%d] start .........................................................',repeat_ind,num_repeat_simulation));
        
        
    %% Parameters
    % 0. experiment
    Opt.TARGET_FOLDER_NAME='result_simul';
    Opt.EXP_NAME=[EXP_NAME0 sprintf('%02d',subj_ind)];
    
    % 1. file
    Opt.SIMUL_PC=1; % 1:local   2:gpu
    if(Opt.SIMUL_PC==1)    Opt.str_path_div='\';   else     Opt.str_path_div='/'; end
    
    % get file names
%     Opt.path_ext=['F:' Opt.str_path_div '0-Program' Opt.str_path_div 'MATLABwork' Opt.str_path_div 'work' Opt.str_path_div 'ModelRL'];
    Opt.path_ext=pwd;
    base_folder=[Opt.path_ext Opt.str_path_div Opt.TARGET_FOLDER_NAME Opt.str_path_div Opt.EXP_NAME];
    list_file = what(base_folder);
    num_file=size(list_file.mat,1);
    
    disp(['#### PROCESSING ' Opt.EXP_NAME '##########']);
    
    %% collect all
    ALL.parameter_converted=cell(1,num_file);
    ALL.fitness_val=cell(1,num_file);
    ALL.Sum_NegLogLik_val=cell(1,num_file);
    all_fitness_mat=[];
    all_sumnegloglik_mat=[];
    for mmm=1:1:num_file
        
        % population file read
        file_name=list_file.mat{mmm,1};
        %     disp(sprintf('$$$$$ Reading {%s}...(%d/%d) $$$$$$$',file_name,mmm,num_file));
        eval(['load ',[Opt.path_ext Opt.str_path_div Opt.TARGET_FOLDER_NAME Opt.str_path_div Opt.EXP_NAME Opt.str_path_div file_name ]]);
        
        %% read all/fwd/sarsa start model
        if(use_FWD_start_model~=0)
            [ind_conform tmp]=find(parameter_converted(:,10)==use_FWD_start_model);
        else
            ind_conform=[1:1:size(parameter_converted,1)]; % select all
        end
        parameter_converted=parameter_converted(ind_conform,:);
        fitness_val=fitness_val(ind_conform);
        Sum_NegLogLik_val=Sum_NegLogLik_val(ind_conform);
        
        %% parameter read
        ALL.parameter_converted{1,mmm}=parameter_converted;
        ALL.fitness_val{1,mmm}=fitness_val;
        ALL.Sum_NegLogLik_val{1,mmm}=Sum_NegLogLik_val;
        
        all_fitness_mat=[all_fitness_mat [mmm*ones(1,length(fitness_val)); [1:1:length(fitness_val)]; fitness_val]];
        all_sumnegloglik_mat=[all_sumnegloglik_mat [mmm*ones(1,length(fitness_val)); [1:1:length(fitness_val)]; Sum_NegLogLik_val']];
        
        % if(IS_PARAMETER_PLUG_IN==1)
        %     [aaa pop_id]=max(fitness_val);
        %     close all
        % %      pop_id=1; %the best one at the previous iteration
        % else
        %     clear all
        %     close all
        %     IS_PARAMETER_PLUG_IN=0;
        % end
        
    end
    
    %% select the top 'NUM_FILE_TO_PROCESS_FOR_EACH_SBJ' sets :    
    all_fitness_mat_sorted=sortrows(all_fitness_mat',3); % ascending sort by fitness value
    all_sumnegloglik_mat_sorted=sortrows(all_sumnegloglik_mat',3); % ascending sort by fitness value
    HIST_best_SumNegLogLik=[HIST_best_SumNegLogLik; [sss_ind min(all_sumnegloglik_mat_sorted(:,3))]];
    
    
    %% Load behaviour data
    state_column=[1 2 3];   action_column=[4 5];
%     load('F:\0-Program\MATLABwork\work\ModelRL\behavior_data_glascher\data.mat')
    load([pwd '\behavior_data_glascher\data.mat'])
    % and attach it to the map
    ind_included=[1];%[1:3,5:12,14:size(data,2)];
    data0=cell(1,length(ind_included));
    for j=1:1:length(ind_included)
        data0{1,j}=data{1,ind_included(j)};
    end
    
    
    X_in=[];    X_in0=[];   
    Y_bestSbj=[];
    Y_param=[];
    
    %% Simulating each case
    % for nnn=1:1:20
%     for nnn=1:1:NUM_FILE_TO_PROCESS_FOR_EACH_SBJ
    for nnn0=1:1:NUM_FILE_TO_PROCESS_FOR_EACH_SBJ
        
        
        
        nnn=model_id_to_process;
        
        % data save
        data_behavior{1,nnn0}.map_1=cell(1,num_max_trial);
        data_behavior{1,nnn0}.map_2=cell(1,num_max_trial);
        
        
        index_nnn=size(all_fitness_mat_sorted,1)-nnn+1;
        index_gen=all_fitness_mat_sorted(index_nnn,1);
        index_pop=all_fitness_mat_sorted(index_nnn,2);
        index_fitness=all_fitness_mat_sorted(index_nnn,3);
        
        if(index_fitness==ALL.fitness_val{1,index_gen}(1,index_pop))
            pop_id=1;
            parameter_converted=ALL.parameter_converted{1,index_gen}(index_pop,:);
        else
            error('-somethings wrong in file reading. fitness in the list should match.')
            break;
        end
        
        disp(sprintf('$$$$$ Processing (%d/%d) $$$$$$$',nnn0,NUM_FILE_TO_PROCESS_FOR_EACH_SBJ));
        
        
        %% create MAP
        % [myMap N_state N_action N_transition]=Model_Map_Init('daw2005');
%         [myMap N_state N_action N_transition]=Model_Map_Init('glascher2010');
        [myMap N_state N_action N_transition]=Model_Map_Init('sangwan2012');
        
        %% create my arbitrator
        myArbitrator=Bayesian_Arbitration_Init(N_state,N_action,N_transition);
        
        %% create my RL
        myState=Model_RL_Init(N_state,N_action,N_transition);
        
        
        
        
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
        Sum_NegLogLik_for_eachSbj=zeros(1,length(ind_included));
        
        if(iscell(parameter_converted))
            parameter_converted=cell2mat(parameter_converted);
        end
        % 1.parameters
        % [PE_tolerance_m1, Time_Step,...
        % A_12, B_12, A_21, B_21,...
        % param_sarsa.alpha, p, .tau(softmax),...
        % myArbitrator.ind_active_model, PE_tolerance_m2]
        
        % [NOTE] the codingbit for 'param_fwd.alpha' is disabled. Instead,
        % 'param_fwd.alpha'='param_sarsa.alpha'
        
        
        data_model.parameter_converted=parameter_converted;
        
        
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
        param_fwd.alpha=parameter_converted(pop_id,7); % learning rate
        % Softmax parameter for all models
        myArbitrator.tau_softmax=parameter_converted(pop_id,9); % use the same value as sarsa/fwd
        param_sarsa.tau=parameter_converted(pop_id,9);
        param_fwd.tau=parameter_converted(pop_id,9);
        param_fwd0.tau=parameter_converted(pop_id,9);
        if(Force_FWD_start==0)
            % arbitrator start
            myArbitrator.ind_active_model=parameter_converted(pop_id,10);
            if(myArbitrator.ind_active_model==1)
                myArbitrator.m1_prob_prev=0.7; % do not use 0.5 which causes oscillation.
                myArbitrator.m2_prob_prev=1-myArbitrator.m1_prob_prev;
            else
                myArbitrator.m1_prob_prev=0.3;
                myArbitrator.m2_prob_prev=1-myArbitrator.m1_prob_prev;
            end
        end
        
        
        %% simulation
        
        total_simul=1;
        num_max_trial0=size(data0{1,1},1)/2; % for latent learning
%         num_max_trial=size(data0{1,1},1)/2;
        
        
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
        
        for ind_dev=1:1:1%[1:1:length(pt_devaluation)]
            
            %% Multiple simulations
            win_fwd=0; win_sarsa=0;
            for kk=1:1:total_simul
                
                ZZZ=[];
                ZZZ_R=[];
                ZZZ1=[];
                ZZZ2=[];
                YYY1=[];
                YYY2=[];
                VVV=[];
                
                
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
                if(0) % latent learning
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
                    map.trial=i;
                    %     disp(sprintf('- [%d/%d] session...',i,num_max_trial));
                    % initializing the state
                    [state_fwd map]=StateClear(state_fwd,map);
                    [state_sarsa map]=StateClear(state_sarsa,map);
                    
                    HIST_MODEL_CHOICE(kk,i)=myArbitrator_top.ind_active_model;
                    while (((myArbitrator_top.ind_active_model==1)&&(~map.JobComplete))||((myArbitrator_top.ind_active_model==2)&&(~map.JobComplete)))
                        
                        HIST_MODEL_CHOICE2=[HIST_MODEL_CHOICE2 myArbitrator_top.ind_active_model];
                        
                        Is_map_changed=0;
                        Is_rwd_prob_changed=0;
                        Is_rwd_changed=0;
                        % devaluation
                        if(length(find(i==pt_devaluation))==1)
                            map.reward=myMap.reward;
                            map.reward(6:8)=0;%[0 10 25];
                            Is_map_changed=1;                            
                        end
                        % contingency degradation
                        if(length(find(i==pt_degradation))==1)
                            map.reward=myMap.reward;
                            [row_nzp col_nzp]=find(map.IsTerminal==1);
                            max_rwd=max(map.reward);
                            map.reward(row_nzp)=max_rwd;
                            Is_map_changed=1;
                        end
                        
                        %change of state transition prob
                        if(length(find(i==pt_prob_chg_max_risk))==1)
                            map.action=myMap.action;
                            % maximum risk
                            for j_act=1:1:size(map.action,2)
                                [row_nzp col_nzp]=find(map.action(1,j_act).prob~=0);
                                for k_state=1:1:length(row_nzp)
%                                     if(row_nzp(k_state)==1) % distal only
                                        map.action(1,j_act).prob(row_nzp(k_state),col_nzp(k_state))=0.5; % all prob to .5
%                                     end
                                end
                            end
                            Is_map_changed=1;
                        end
                        
                        if(length(find(i==pt_prob_chg_min_risk))==1)
                            map.action=myMap.action;                                                        
                            for j_act=1:1:size(map.action,2)
                                [row_nzp col_nzp]=find((map.action(1,j_act).prob>0.69)&(map.action(1,j_act).prob<0.71));
                                [row_nzp2 col_nzp2]=find((map.action(1,j_act).prob>0.29)&(map.action(1,j_act).prob<0.31));
                                for k_state=1:1:length(row_nzp)
                                    map.action(1,j_act).prob(row_nzp(k_state),col_nzp(k_state))=0.9;
                                end
                                for k_state=1:1:length(row_nzp2)
                                    map.action(1,j_act).prob(row_nzp2(k_state),col_nzp2(k_state))=0.1;
                                end
                            end
                            Is_map_changed=1;
                        end
                        
                        if(length(find(i==pt_prob_chg_reversed))==1)
                            map.action=myMap.action;
                            % prob reversed
                            for j_act=1:1:size(map.action,2)
                                [row_nzp col_nzp]=find((map.action(1,j_act).prob>0.69)&(map.action(1,j_act).prob<0.71));
                                [row_nzp2 col_nzp2]=find((map.action(1,j_act).prob>0.29)&(map.action(1,j_act).prob<0.31));
                                for k_state=1:1:length(row_nzp)
                                    map.action(1,j_act).prob(row_nzp(k_state),col_nzp(k_state))=0.3;
                                end
                                for k_state=1:1:length(row_nzp2)
                                    map.action(1,j_act).prob(row_nzp2(k_state),col_nzp2(k_state))=0.7;
                                end
                            end
                            Is_map_changed=1;
                        end
                        
                        % change of reward delivery prob
                        if(length(find(i==pt_prob_rwd_chg))==1)
                            map.reward_prob=myMap.reward_prob;
                            [row_nzrp col_nzrp]=find(map.reward_prob>0.05);
                            for k_state=1:1:length(row_nzrp)
                                map.reward_prob(row_nzrp(k_state))=0.5;
                            end
                            Is_rwd_prob_changed=1;
                        end
                        
                        % keep changing reward delivery prob (goal-directed mode)
                        if(length(find(i==pt_devaluation_goal))==1)
                            map.reward=myMap.reward;
                            if(map.index==1) % distal
                                goal_cadidate=[6 8];
                                goal_ind=goal_cadidate(ceil(2*rand(1)));   goal_val=map.reward(goal_ind);
                                map.reward(6:9)=rwd_aversion_factor*map.reward(6:9);      map.reward(goal_ind)=goal_val;
                                Is_rwd_changed=1;
                            end
                            if(map.index==2) % proximal
                                if(myArbitrator_top.ind_active_model==1) % fwd
                                    current_state_2=state_fwd.state_history(2);
                                else
                                    current_state_2=state_sarsa.state_history(2);
                                end
                                if((current_state_2==2)||(current_state_2==3))
                                    goal_cadidate=[7 8];
                                end
                                if((current_state_2==4)||(current_state_2==5))
                                    goal_cadidate=[6 7];
                                end
                                goal_ind=goal_cadidate(ceil(2*rand(1)));   goal_val=map.reward(goal_ind);
                                map.reward(6:9)=rwd_aversion_factor*map.reward(6:9);      map.reward(goal_ind)=goal_val;
                                Is_rwd_changed=1;
                            end
                            % fwd model update
                            state_fwd=Model_RL(state_fwd, param_fwd, map, 'bwd_update');
                        end
                        
                        % none of the above case, then keep the original plan
                        if(Is_map_changed==0)
                            map.action=myMap.action;                            
                        end
                        if(Is_rwd_prob_changed==0)
                            map.reward_prob=myMap.reward_prob;
                        end
                        if(Is_rwd_changed==0)
                            map.reward=myMap.reward;
                        end
                        
                        % map save
                        if(map.index==1) % distal
                            data_behavior{1,nnn0}.map_1{1,i}=map;
                        end
                        if(map.index==2) % proximal
                            data_behavior{1,nnn0}.map_2{1,i}=map;
                        end
                        
                        
                        %% Compute negative log-likelihood
%                         for ll=1:1:length(ind_included) % compute map.trial's log-lik for ll-th subjects
%                             state_data=data0{1,ll}(map.trial,state_column(map.index));
%                             action_data=data0{1,ll}(map.trial,action_column(map.index));
%                             var_exp=exp(myArbitrator_top.tau_softmax*myArbitrator_top.Q(state_data,:)); % (N_actionx1)
%                             Sum_NegLogLik=Sum_NegLogLik-log(var_exp(action_data)/sum(var_exp))/total_simul;
%                             Sum_NegLogLik_for_eachSbj(1,ll)=Sum_NegLogLik_for_eachSbj(1,ll)-log(var_exp(action_data)/sum(var_exp))/total_simul;
%                         end
                        %%
                        
                        % index synchronization
                        state_fwd.index=map.index;                state_sarsa.index=map.index;
                        
                        ZZZ=[ZZZ [map.index; state_fwd.index; state_sarsa.index; [0;0;0]; myArbitrator_top.ind_active_model; myArbitrator_top.m1_prob_prev]];
                        ZZZ_R=[ZZZ_R [map.index; map.reward(6:9)]];
                        
                
                        
                        
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
                
                                                
                        %                     ZZZ=[ZZZ [map.index; state_fwd.index; state_sarsa.index; myArbitrator_top.index; myArbitrator_top.temp; myArbitrator_top.ind_active_model; myArbitrator_top.m1_prob]];
                        ZZZ1=[ZZZ1 [i; myArbitrator_top.m1_wgt'*myArbitrator_top.m1_inv_Fano; sum(myArbitrator_top.m1_inv_Fano); myArbitrator_top.m1_mean(2); myArbitrator_top.m1_var(2)]];
                        ZZZ2=[ZZZ2 [i; myArbitrator_top.m2_wgt'*myArbitrator_top.m2_inv_Fano; sum(myArbitrator_top.m2_inv_Fano); myArbitrator_top.m2_mean(2); myArbitrator_top.m2_var(2)]];
                        YYY1=[YYY1 [i; myArbitrator_top.transition_rate12]];
                        YYY2=[YYY2 [i; myArbitrator_top.transition_rate21]];                        
                        VVV=[VVV myArbitrator_top.action_history];                        
                        
                        
                        % data save
                        if(map.index==3)
                            if(i>1)
                                data_behavior{1,nnn0}.data(i,:)=[myArbitrator_top.state_history' myArbitrator_top.action_history(1:2)' 2 zeros(1,7) myArbitrator_top.reward_history(3) sum(data_behavior{1,nnn0}.data(1:i-1,14))+myArbitrator_top.reward_history(3)];
                            else % if i=1
                                data_behavior{1,nnn0}.data(i,:)=[myArbitrator_top.state_history' myArbitrator_top.action_history(1:2)' 2 zeros(1,7) myArbitrator_top.reward_history(3) myArbitrator_top.reward_history(3)];
                            end
                        end
                        
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
                %             disp(sprintf('- condition[%d],simulation[%d/%d] - (fwd_win|sarsa_win)=(%d|%d)',ind_dev,kk,total_simul,win_fwd,win_sarsa))
                
            end
            
            
            %         %% figure
            %
            %         str=sprintf('[Condition-%d] Q(S1,L(line)/R(dotted)) - SARSA(blue) vs. FWD(red)',ind_dev);
            %         figure('Name',str)
            %         ind=[1:1:length(squeeze(HIST_Q_SARSA(1,1,:)))];
            %         plot(ind,squeeze(HIST_Q_SARSA(1,1,:)),'b-',ind,squeeze(HIST_Q_SARSA(1,2,:)),'b-.',ind,squeeze(HIST_Q_FWD(1,1,:)),'r-',ind,squeeze(HIST_Q_FWD(1,2,:)),'r-.')
            %
            %         str=sprintf('[Condition-%d] Q(S4,L(line)/R(dotted)) - SARSA(blue) vs. FWD(red)',ind_dev);
            %         figure('Name',str)
            %         ind=[1:1:length(squeeze(HIST_Q_SARSA(4,2,:)))];
            %         plot(ind,squeeze(HIST_Q_SARSA(4,1,:)),'b-',ind,squeeze(HIST_Q_SARSA(4,2,:)),'b-.',ind,squeeze(HIST_Q_FWD(4,1,:)),'r-',ind,squeeze(HIST_Q_FWD(4,2,:)),'r-.')
            %
            %         str=sprintf('[Condition-%d] T(S1,R,S4) - FWD',ind_dev);
            %         figure('Name',str)
            %         plot([squeeze(HIST_T_FWD_1(1,2,4,:)); squeeze(HIST_T_FWD_2(1,2,4,:))])
            
        end
        %     disp('-done');
        
        %% save to matrix
        X_in_fwd=[X_in_fwd; ZZZ1(2,:)./ZZZ1(3,:)];     X_in_sarsa=[X_in_sarsa; ZZZ2(2,:)./ZZZ2(3,:)];
        X_in_fwd_mean=[X_in_fwd_mean; ZZZ1(4,:)];   X_in_sarsa_mean=[X_in_sarsa_mean; ZZZ2(4,:)];
        X_in_fwd_var=[X_in_fwd_var; ZZZ1(5,:)];   X_in_sarsa_var=[X_in_sarsa_var; ZZZ2(5,:)];
        X_in0=[X_in0; ZZZ(7,:)]; % ind_active_model
        X_in=[X_in; ZZZ(8,:)]; % p(m1)
        [tmp ind_best]=max(Sum_NegLogLik_for_eachSbj);
        Y_bestSbj=[Y_bestSbj ind_included(ind_best)];
        Y_param=[Y_param parameter_converted'];
        
    end
    
    X_all=[X_all; X_in];
    Y_all=[Y_all [Y_param; subj_ind*ones(1,NUM_FILE_TO_PROCESS_FOR_EACH_SBJ)]];
    
    %% save
    eval(['save ',[Opt.path_ext Opt.str_path_div Opt.TARGET_FOLDER_NAME Opt.str_path_div Opt.EXP_NAME ' X_in Y_bestSbj Y_param' ]]);

    
end


% collection_behavior_profile=[collection_behavior_profile; X_all];


end % end of repeat_ind

%% figure - alpha/beta
if(0)
    
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
        boundedline(ind_line,HIST_PE_MEAN_SARSA(j,:),shade_mat_sarsa,'alpha','b',ind_line,HIST_PE_MEAN_FWD(j,:),shade_mat_fwd,'alpha','r'); % with 'alpha' option you cannot make a proper screen capture.
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



%% data save
save data_save data_behavior data_model


%% Display - behaviour profile
% 1.P(fwd)
figure('Name','P of choosing fwd')
[X,Y] = meshgrid([1:1:(length(SUB_ARRAY)*NUM_FILE_TO_PROCESS_FOR_EACH_SBJ*num_repeat_simulation+1)], [1:1:size(X_all,2)]);
half_range=max(max(max(X_all))-0.5,0.5-min(min(X_all))); % maximum dist from 0.5
STATE_mat=(ZZZ(1,:)-1)*half_range*2+0.5-half_range;
surf(X,Y,[STATE_mat' X_all'])
axis([1 size(X,2) 1 size(X,1) 0.5-half_range 0.5+half_range])
title('P(fwd choice)');
view(90,90) % view point
colormap jet
caxis([0.5-half_range 0.5+half_range])
colorbar



%% Run "mani"

% [Note] use 'X_in' for input matrix: p(m1)

% [Note] 1. Y_all(j,:)
% [PE_tolerance, Time_Step,...
% A_12, B_12, A_21, B_21,...
% param_sarsa.alpha, param_fwd.alpha, .tau(softmax),...
% myArbitrator.ind_active_model subj_ind]
% [Note] 2. Y_bestSbj: subject index that the model fits best.

mani

%% [IMPORTANT] and run this in workspace!!!!
% SIMUL_Arbitration_v3_GA_Analysis_level2(maniX,maniY,Y_all)


disp('### after mani, run     SIMUL_Arbitration_v3_GA_Analysis_level2(maniX,maniY,Y_all,X_in_fwd,X_in_sarsa,X_in_fwd_mean,X_in_fwd_var,X_in_sarsa_mean,X_in_sarsa_var)')