
clear all
close all
warning('off')

path0='F:\0-Program\MATLABwork\work\ModelRL';
seed_path_result=[path0 '\result_save\'];
save_path_result=[path0 '\result_simul\'];
save_for_SPM=['C:\DATA_fmri\uncertainty_arbitration\regressors_contrasts\'];



% 1. Behavioral data
% LIST_SBJ={'david', 'DeDe', 'rosemary', 'Boyu', 'melissa', 'Rehevolew', 'joel', 'clarke', 'angela', 'william', 'josephine'}; % (good in pre which is mostly habitual - rosemary, melissa)
% mode.map_type=?;

% 2. behavioral + fmri data (for map config, see SIMUL_arbitraion_fmri2.m)
% [note] 'Oliver' uses an old map. the rest of them use a new map.
LIST_SBJ={'Oliver', 'Hao', 'Breanna', 'Derek', 'Timothy', 'Teagan', 'Jeffrey', 'Seung', 'Carole', 'Tony', 'Surendra', 'Lark'};
mode.map_type=1;

% regressor list
LIST_REGRESSOR={'SPE', 'RPE', 'uncertaintyM1', 'uncertaintyM2', 'meanM1', 'meanM2', 'invFanoM1', 'invFanoM2', 'weigtM1', 'weigtM2'};
row_mat=[7 7 8 8 8 8 8 8 7 7]; % from which row in the SBJ{}.regressor matrix the signal needs to be extracted. e.g., uncertainty of 0 prediction error


%% subject parameter
list_sbj_included=[2];%[2:1:12];

%% model options
option_optimizing_model=3; % 0: optimizing the model for each sbj, 1: for all sbj, 2: do not optimize; load saved model
mode.USE_FWDSARSA_ONLY=0; % 0: arbitration, 1: use fwd only, 2: use sarsa only
mode.USE_BWDupdate_of_FWDmodel=1; % 1: use the backward update for goal-directed model (fwd model), 0: do not use
mode.DEBUG_Q_VALUE_CHG=0; % Debug option 1: show Q-value before/after whenever there is a goal change.
mode.path_ext=path0;
mode.total_simul=1; % # of total simulation repetition per subject
mode.simul_process_display=0; % 1: display model's process, 0: no diplay
mode.experience_sbj_events=[0 1]; % [pre main] - 1: experience exactly the same events(decision,state) as subjects. 0: model's own experience
mode.out=1; % 1: normal evaluation mode, 99: regressor added to the SBJ, 0: debug mode
mode.max_iter=200; % maximum iteration for optimization

%% Regressor options
% 'SPE', 'RPE', 'uncertaintyM1', 'uncertaintyM2', 'meanM1', 'meanM2', 'invFanoM1', 'invFanoM2', 'weigtM1', 'weigtM2'
% should add the regressors in the order of importance
param_regressor_type_cue={'SPE', 'RPE'};


%% initialization
use_model_regressor_cue=0;  ind_regressor=[];
for ii=1:1:size(param_regressor_type_cue,2)
    if(strcmp(param_regressor_type_cue{1,ii},'SPE')==1)    use_model_regressor_cue=1;  ind_regressor=[ind_regressor 1];    end
    if(strcmp(param_regressor_type_cue{1,ii},'RPE')==1)    use_model_regressor_cue=1;  ind_regressor=[ind_regressor 2];    end    
    if(strcmp(param_regressor_type_cue{1,ii},'uncertaintyM1')==1)    use_model_regressor_cue=1;  ind_regressor=[ind_regressor 3];    end
    if(strcmp(param_regressor_type_cue{1,ii},'uncertaintyM2')==1)    use_model_regressor_cue=1;  ind_regressor=[ind_regressor 4];    end
    if(strcmp(param_regressor_type_cue{1,ii},'meanM1')==1)    use_model_regressor_cue=1;  ind_regressor=[ind_regressor 5];    end
    if(strcmp(param_regressor_type_cue{1,ii},'meanM2')==1)    use_model_regressor_cue=1;  ind_regressor=[ind_regressor 6];    end
    if(strcmp(param_regressor_type_cue{1,ii},'invFanoM1')==1)    use_model_regressor_cue=1;  ind_regressor=[ind_regressor 7];    end
    if(strcmp(param_regressor_type_cue{1,ii},'invFanoM2')==1)    use_model_regressor_cue=1;  ind_regressor=[ind_regressor 8];    end
    if(strcmp(param_regressor_type_cue{1,ii},'weigtM1')==1)    use_model_regressor_cue=1;  ind_regressor=[ind_regressor 9];    end
    if(strcmp(param_regressor_type_cue{1,ii},'weigtM2')==1)    use_model_regressor_cue=1;  ind_regressor=[ind_regressor 10];    end
end


%% subject data loading 

% which subject to be included
% ### READ ONLY ONE SBJ BECAUSE EACH MODEL WILL LEARN EACH SBJ BEHAVIOR.
ind_sbj_included=list_sbj_included;      SUB_ARRAY=list_sbj_included;
num_sbj_included=length(ind_sbj_included);
ind_included=ind_sbj_included;

for k=1:1:num_sbj_included
    LIST_SBJ_included{1,k}=LIST_SBJ{1,ind_sbj_included(k)};
end
for i=1:1:num_sbj_included %=1. process only 1 subject
    
    SBJ{1,i}.name=LIST_SBJ{1,ind_sbj_included(i)};
    
    % 'pre' file load: HIST_block_condition{1,session_ind}, HIST_behavior_info{1,session_ind}
    file_name=[LIST_SBJ{1,ind_sbj_included(i)} '_pre_info.mat'];
    file_name_full=[mode.path_ext '\result_save\' file_name];
    load(file_name_full);
    SBJ{1,i}.HIST_block_condition_pre=HIST_block_condition;
    SBJ{1,i}.HIST_behavior_info_pre=HIST_behavior_info;
    
    % 'fmri' file load: HIST_block_condition{1,session_ind}, HIST_behavior_info{1,session_ind}
    file_name=[LIST_SBJ{1,ind_sbj_included(i)} '_fmri_info.mat'];
    file_name_full=[mode.path_ext '\result_save\' file_name];
    load(file_name_full);
    SBJ{1,i}.HIST_behavior_info=HIST_behavior_info;
    SBJ{1,i}.HIST_behavior_info_Tag=HIST_behavior_info_Tag;
    SBJ{1,i}.HIST_event_info=HIST_event_info;
    SBJ{1,i}.HIST_event_info_Tag=HIST_event_info_Tag;
    SBJ{1,i}.HIST_block_condition=HIST_block_condition;
    SBJ{1,i}.HIST_block_condition_Tag=HIST_block_condition_Tag;
    num_tot_session=size(SBJ{1,i}.HIST_behavior_info,2);
    
    % [fixing part!!! - for Oliver]
    if(strcmp(SBJ{1,i}.name,'Oliver'))
        for mm=1:1:size(SBJ{1,i}.HIST_event_info,2) % each session
            mat_fixing=SBJ{1,i}.HIST_event_info{1,mm};
            index_delete=zeros(1,size(mat_fixing,2));
            [r_fix, c_fix]= find(mat_fixing(7,:)==9);
            for nn=1:1:length(c_fix)
                % check the previous event
                if(mat_fixing(7, c_fix(nn)-1)~=0.5)
                    index_delete(c_fix(nn))=1;
                end
            end
            [tmp c_keep]=find(index_delete==0);
            mat_fixed=mat_fixing(:,c_keep);
            SBJ{1,i}.HIST_event_info{1,mm}=mat_fixed;
        end
    end
     
    
    % [NOTE] now we have 4 variables: mode.HIST_block_condition_pre, mode.HIST_block_condition, mode.HIST_behavior_info_pre, mode.HIST_behavior_info
    % to read a block condition, use "block_condition=mode.HIST_block_condition{1,session_ind}(2,block_ind); % G:1,G':2,H:3,H':4"
    
    %     swsw_amount_pre = [swsw_amount_pre mode.HIST_behavior_info_pre{1,1}(end,17)];
    tot_amount_earned_main_each_sbj =[];
    for jk=1:1:size(SBJ{1,i}.HIST_behavior_info,2)    tot_amount_earned_main_each_sbj = [tot_amount_earned_main_each_sbj; SBJ{1,i}.HIST_behavior_info{1,jk}(end,17)]; end
    %     swsw_amount_main=[swsw_amount_main tot_amount_earned_main_each_sbj];
end







%% model optimization

param_init=[0.5, 10.0, 4, 2.67, 11.7, 10.0, 0.15, 0.1];
param_BoundL=[0.01, 3, 0.5*param_init(3:1:6), 0.01, 0.01]; %[0.01, 3, 1.0, 2.0, 0.5, 5, 0.01, 0.01];
param_BoundU=[1, 20, 2*param_init(3:1:6), 0.5, 0.2]; %[1, 20, 6, 10, 4, 30, 0.5, 0.2];

% ## (way1-each) optimizing for *each* subject and plug the result into each SBJ structure
if(option_optimizing_model==0)    
    for ind_sbj=1:1:size(SBJ,2)
        clear SBJ_test;
        SBJ_test{1,1}=SBJ{1,ind_sbj};
        disp('############################################')
        disp(['#### optimizing RL-arbitrator for ' sprintf('SBJ#%02d...',ind_sbj)]);
        disp('############################################')
        % [1] model optimization
        mode.out=1;
        myFunc_bu = @(x) eval_ArbitrationRL2(x, SBJ_test, mode); % define a new anonymous function
        [model_BayesArb.param, model_BayesArb.val]=fminsearchbnd(myFunc_bu, param_init, param_BoundL, param_BoundU, optimset('Display','iter','MaxIter',mode.max_iter));   % X0,LB,UB
        % [2-1] add regressor vector to SBJ
        mode.out=99;
        SBJ_test=eval_ArbitrationRL2(model_BayesArb.param,SBJ_test,mode);
        % [3] Save
        model_BayesArb.mode=mode;
        SBJ_test{1,1}.model_BayesArb=model_BayesArb;
        SBJ{1,ind_sbj}=SBJ_test{1,1};        
        save_file_name=['SBJ_structure.mat'];
        eval(['save ' save_path_result save_file_name ' SBJ'])
    end
end
if(option_optimizing_model==1)
    % ## (way2-batch) optimizing for *all* subjects and plug the result into each SBJ structure
    % [1] model optimization
    mode.out=1;
    myFunc_bu = @(x) eval_ArbitrationRL2(x, SBJ, mode); % define a new anonymous function
    [model_BayesArb.param, model_BayesArb.val]=fminsearchbnd(myFunc_bu, param_init, param_BoundL, param_BoundU, optimset('Display','iter'));   % X0,LB,UB
    % [2-1] add regressor vector to SBJ
    mode.out=99;
    SBJ=eval_ArbitrationRL2(model_BayesArb.param,SBJ,mode);
    % save
%     model_BayesArb.mode=mode;
%     SBJ_test{1,1}.model_BayesArb=model_BayesArb;
end
if(option_optimizing_model==2)
    load_file_name=['SBJ_structure.mat'];
    eval(['load ' save_path_result load_file_name])
end



%test mode
if(option_optimizing_model==3)
    model_BayesArb.param=param_init;
    mode.out=99;
    SBJ=eval_ArbitrationRL2(model_BayesArb.param,SBJ,mode);
end



%% Create regressors
% state. 0.5: fixation mark on, 1: S1, 2: S2, 3: S3, 4: S4, 5: S5,
% 6(+/-)0.1: O1(with win/lost msg), 7(+/-)0.1: O2(with win/lost msg), 8(+/-)0.1: O3(with win/lost msg), 9: O4,
% 10:A1, 11:A2, 20: a short blank page display, -99:fail to choose in time limit, (-) when display off

for jj2=1:1:size(SBJ,2)        % each subject
    for kk2=1:1:size(SBJ{1,jj2}.HIST_behavior_info,2)  % each main session
        
        %% regressor generation for each main session and save it to a single file that is compatible with SPM
        ind_reg=0;  % corresponds to the size of the structure
        ind_reg_abs=0; % actual number of regressors (including parametric)
        
        mat_work=SBJ{1,jj2}.HIST_event_info{1,kk2};
        
        durations={};
        onsets={};
        names={};
        pmod=struct('name',{},'param',{},'poly',{});
        
        
        %% Regressor cue presentation
        
        ind_reg=ind_reg+1;  use_model_regressor_cue=1;
        
        % (1) durations, name, onset        
        [tmp col_on]=find((mat_work(7,:)==1)|(mat_work(7,:)==2)|(mat_work(7,:)==3)|(mat_work(7,:)==4)|(mat_work(7,:)==5)...
            |(mat_work(7,:)==5.9)|(mat_work(7,:)==6.1)|(mat_work(7,:)==6.9)|(mat_work(7,:)==7.1)|(mat_work(7,:)==7.9)|(mat_work(7,:)==8.1));
        
        RT_mat=[];  onset_mat=[];   param_mat=[];
        prev_trial=0;  show_n_th_times_t=0;        
        param_mat=zeros(length(ind_regressor),length(col_on));
        for ll2=1:1:length(col_on)
            
            if(ll2<length(col_on)) % usual case
                pt_on=mat_work(4,col_on(ll2));
                col_off=col_on(ll2)-1+find(mat_work(7,[col_on(ll2):1:(col_on(ll2)+2)])==0.5); % find the next fixation mark presentation
                pt_off=mat_work(4,col_off);
                RT=pt_off-pt_on;
            else % last event in the session is the outcome presentation
                RT=2.0;
            end
            RT_mat=[RT_mat RT];
            onset_t=mat_work(4,col_on(ll2));
            onset_mat=[onset_mat onset_t];
            
            if(use_model_regressor_cue==1)
                
                for nn=1:1:length(ind_regressor)
                    mysession=kk2;
                    myblock=mat_work(1,col_on(ll2)); % block in session
                    mytrial=mat_work(2,col_on(ll2)); % trial in block
                    mytrial_s=mat_work(3,col_on(ll2))-1; % trial_s in trial (arbitration index: 1 at the second stage, 2 at the third stage)
                    
                    
                    if(mytrial_s==0) % the first decision stage
                        param_mat(nn,ll2)=0.0; % [NOTE] NO signal at the first stage
                    else % the second and the third decision stage
                        mat_work_reg=SBJ{1,jj2}.regressor{1,ind_regressor(nn)}.value(1:4,:);
                        identity_tmp=sum(abs(mat_work_reg-repmat([mysession myblock mytrial mytrial_s]',1,size(mat_work_reg,2))));
                        col_event=find(identity_tmp==0);
                        param_mat(nn,ll2)=SBJ{1,jj2}.regressor{1,ind_regressor(nn)}.value(row_mat(ind_regressor(nn)),col_event);
                    end
                    
                end
            end
            
        end        
        
        durations{1,ind_reg}=RT_mat;
        onsets{1,ind_reg}=onset_mat;
        names{1,ind_reg}=['Cue'];
        ind_reg_abs=ind_reg_abs+1;      list_name_for_contrast{1,ind_reg_abs}=names{1,ind_reg}; % add to the global list of regresssors
        
        % (2) pmod: how many times each cue presented
        if(use_model_regressor_cue==1)
            for nn=1:1:length(ind_regressor)
                pmod(1,ind_reg).name{1,nn}=[SBJ{1,jj2}.regressor{1,ind_regressor(nn)}.name];
                pmod(1,ind_reg).poly{1,nn}=1;
                pmod(1,ind_reg).param{1,nn}=param_mat(nn,:);
                ind_reg_abs=ind_reg_abs+1;      list_name_for_contrast{1,ind_reg_abs}=pmod(1,ind_reg).name{1,nn};
            end
        end
        
        
        
        %% Saving regressor file
        tot_num_myregressor=length(list_name_for_contrast);
        save_file_name=['Regressor--' SBJ{1,jj2}.name '_sess' sprintf('%02d.mat',kk2)];
        eval(['save ' save_path_result save_file_name ' durations names onsets pmod'])
        eval(['save ' save_for_SPM save_file_name ' durations names onsets pmod'])
        
        
    end
end






%% Saving Contrast file
% [index of my regressors for contrast vector] : total main regressor=6, total regressors=7
clear contrast_spm

total_number_regressor=tot_num_myregressor+6; % # + 6 movements (+ 1 constant)
ind_contrast_vec=0;


% individual : (ex) [0 1 0 0 0 0 0 0]
for ii=1:1:tot_num_myregressor 
    ind_contrast_vec=ind_contrast_vec+1;
    contrast=zeros(1,total_number_regressor);
    contrast(1,ii)=1;
    contrast_spm{1,ind_contrast_vec}.name=list_name_for_contrast{1,ii};
    contrast_spm{1,ind_contrast_vec}.vec=contrast;    
end

% difference : (ex) [0 0 0 0 1 -1 0 0]
for ii=[2:1:3]
    for jj=[(ii+1):1:3]
        % A-B
        ind_contrast_vec=ind_contrast_vec+1;
        contrast=zeros(1,total_number_regressor);
        contrast(1,ii)=1;    contrast(1,jj)=-1;
        contrast_spm{1,ind_contrast_vec}.name=[list_name_for_contrast{1,ii} '>' list_name_for_contrast{1,jj}];
        contrast_spm{1,ind_contrast_vec}.vec=contrast;
        % B-A
        ind_contrast_vec=ind_contrast_vec+1;
        contrast=zeros(1,total_number_regressor);
        contrast(1,jj)=1;    contrast(1,ii)=-1;
        contrast_spm{1,ind_contrast_vec}.name=[list_name_for_contrast{1,ii} '<' list_name_for_contrast{1,jj}];
        contrast_spm{1,ind_contrast_vec}.vec=contrast;
    end
end


eval(['save ' save_path_result 'contrast_spm.mat' ' contrast_spm'])
eval(['save ' save_for_SPM 'contrast_spm.mat' ' contrast_spm'])










%% measure the degree of habit in habitual conditions
% block condition - 1: G(with low uncertainty), 2: G''(with high uncertainty), 3:H(with high uncertainty), 4:H''(with low uncertainty)';
% ### [note]: using "state_action_vec_ref" might be stupid idea. using the
% actual action taken by the model would make more sense!!!
if(1)
    
    state_action_vec_ref=[1 2; 2 1; 3 2; 4 2; 5 1]; % col1:state, col2:corresponding action
    for i=1:1:num_sbj_included
        num_tot_sess=size(SBJ{1,i}.HIST_behavior_info,2);
        sum_mat_percentage=zeros(4,2);
        for i_sess=1:1:num_tot_sess
            condi_to_check=[1 2 3 4]; %for all conditions
            mat_percentage=[];
            for kk=1:1:length(condi_to_check)
                % (1) haibual condition
                row_condi=find(SBJ{1,i}.HIST_behavior_info{1,i_sess}(:,3)==condi_to_check(kk));
                num_total_trial=size(SBJ{1,i}.HIST_behavior_info{1,i_sess},1);
                val_score=0;
                for j=1:1:length(row_condi)
                    state_vec=SBJ{1,i}.HIST_behavior_info{1,i_sess}(row_condi(j),[4:5]);
                    state_action_vec0=state_action_vec_ref(state_vec,2)'; % strong habitual action
                    state_action_vec1=SBJ{1,i}.HIST_behavior_info{1,i_sess}(row_condi(j),[7:8]); % subject's action
                    eval_vec=abs(state_action_vec0-state_action_vec1); % all zero = same actions
                    val_score=val_score+length(find(eval_vec==0));
                end
                mat_percentage=[mat_percentage; [condi_to_check(kk) 100*val_score/(2*length(row_condi))]];
            end
            SBJ{1,i}.HIST_block_condition_habit_score{1,i_sess}=mat_percentage;
            SBJ{1,i}.HIST_block_condition_habit_score_Tag{1,i_sess}='col1: block condition, col2: percentage of habitual action';
            sum_mat_percentage=sum_mat_percentage+mat_percentage;
        end
        SBJ{1,i}.HIST_block_condition_habit_score_mean=sum_mat_percentage/num_tot_sess;
        
    end
    
end






disp('- all done.')




















