%%
%% run in MATLAB NOTE R2011b 32bit !!!!
%%
function func_SIMUL_arbitration_regressor_generator_v6_BMS_kdj_2(ori)

warning('off')

% For saving figures, 'savefigs'

path0='\\143.248.30.94\bmlsamba\swlee\repro_fmri\fmri_arbitration\modelRLsource';
seed_path_result=[path0 '\result_save\'];
save_path_result=[path0 '\result_simul\'];
save_for_SPM=['\\143.248.30.94\bmlsamba\swlee\repro_fmri\fmri_arbitration\regressors_contrasts\'];
save_path_neuroecon='\\143.248.30.94\bmlsamba\swlee\repro_fmri\fmri_arbitration\regressors_contrasts\';

% 1. Behavioral data
% LIST_SBJ={'david', 'DeDe', 'rosemary', 'Boyu', 'melissa', 'Rehevolew', 'joel', 'clarke', 'angela', 'william', 'josephine'}; % (good in pre which is mostly habitual - rosemary, melissa)
% mode.map_type=?;

% 2. behavioral + fmri data (for map config, see SIMUL_arbitraion_fmri2.m)
% [note] 'Oliver' uses an old map. the rest of them use a new map.
LIST_SBJ={'Oliver', 'Hao', 'Breanna', 'Derek', 'Timothy', 'Teagan', 'Jeffrey', 'Seung', 'Carole', 'Tony', 'Surendra', 'Lark',...
    'Joaquin', 'DavidB', 'Christopher', 'Gjergji', 'Charles', 'Erin', 'Connor', 'Domenick', 'Thao', 'Arin', 'Pauline', 'Tho'};
LIST_sbj_map_type=[1*ones(1,12) 2*ones(1,50)]; %1:'sangwan2012b', 2:'sangwan2012c'

% regressor list
% [CAUTION] DO NOT change the order!!!
% [NOTE] if "TYPE_REGRESSOR" changed, change "param_regressor_type_cue_abs_pos_in_design_mat" accordingly!!!
LIST_REGRESSOR={'SPE', 'RPE', 'uncertaintyM1', 'uncertaintyM2', 'meanM1', 'meanM2', 'invFanoM1', 'invFanoM2', 'weigtM1', 'weigtM2', 'Qfwd', 'Qsarsa', 'Qarb', 'dQbwdEnergy', 'dQbwdMean',...
    'duncertaintyM1', 'dinvFanoM1','TR_alpha','TR_beta','dinvFano12','ABSdinvFano12','MAXinvFano12','CONFLICTinvFano12','PMB',...
    'invFanoM1_meancorrected','invFanoM2_meancorrected'};
LIST_REGRESSOR={'SPE', 'RPE', 'uncertaintyM1', 'uncertaintyM2', 'meanM1', 'meanM2', 'invFanoM1', 'invFanoM2', 'weigtM1', 'weigtM2', 'Qfwd', 'Qsarsa', 'Qarb', 'dQbwdEnergy', 'dQbwdMean', 'duncertaintyM1', 'dinvFanoM1', 'TR_alpha', 'TR_beta', 'dinvFano12', 'ABSdinvFano12','MAXinvFano12', 'CONFLICTinvFano12',  'PMB', 'invFanoM1_meancorrected', 'invFanoM2_meancorrected', 'Prob_actionR', 'PMF','SPE_BL', 'RPE_BL', 'SPE_z', 'RPE_z'};

TYPE_REGRESSOR=[1 1, 1.5 1.5, 1.5 1.5, 1.5 1.5, 2 2, 1.5 1.5 1.5 1.5 1.5 1.5 1.5 1.5 1.5 1.5 1.5 1.5 1.5 1.5 1.5 1.5]; % 1: parametric modulation (0-duration), 1.5:parmetric modulation (non-zero duration), 1.7:parametric modulation (with decision onset)  2: extra continuous parametric modulation (TR-fixed) - this will be used by "dummy" regressor.
TYPE_REGRESSOR=[1 1, 1.5 1.5, 1.5 1.5, 1.5 1.5, 2 2, 1.5 1.5 1.5 1.5 1.5 1.5 1.5, 3 3 3 3 1.5 3 3 3 3 3 3 1.5 1.5 1.5 1.5]; % 1: parametric modulation (0-duration), 1.5:parmetric modulation (non-zero duration), 1.7:parametric modulation (with decision onset)  2: extra continuous parametric modulation (TR-fixed) - this will be used by "dummy" regressor. 3 : IDK

row_mat=[7 7 7 7 7 7 7 7 7 7 7 7 7 7 7 7 7 7 7 7 7 7 7 7 7 7 7 7 7 7 7 7 7 7 7 7]; % from which row in the SBJ{}.regressor matrix the signal needs to be extracted. e.g., uncertainty of 0 prediction error




%% OPTION - subject
% [note] DO NOT USE sbj#[20] - he pressed wrong buttons in session1,2, so need to shrink all the SBJ matrix size by deleting the session#1,2
list_sbj_included=[2:1:19 21:1:24]; % 15,22: badly fitted sbjs - [2:1:14 16:1:19 21 23:1:24]; 

%% OPTION - model optimization
option_optimizing_model=2; % 0: optimizing the model for each sbj, 1: for all sbj, 2: do not optimize; load saved model
% BEST full Bayes: name_paramfile_to_test='SBJ_structure_each_exp_BESTcollection_FullBayes.mat', BEST simpleBayes: 'SBJ_structure_each_exp_BESTcollectionJan29.mat'
% name_paramfile_to_test='SBJ_structure_each_exp_BESTcollectionJan29.mat'; % valid only if option_optimizing_model=2.
% name_paramfile_to_test='SBJ_structure_each_exp_BESTcollection_Daw.mat'; % Daw's model

if ori == 2
name_paramfile_to_test='SBJ_structure_new2.mat'; 
else
name_paramfile_to_test='SBJ_structure_ori2.mat'; 
end

update_SBJ_structure=0; % 0: no update/just read and use, 1: update the changes to the saved SBJ file
mode.opt_ArbModel=1; % 0: full arbitrator, 1: invF-based, 2: mean-based, 3: uncertainty-based arbitrator
mode.USE_FWDSARSA_ONLY=0; % 0: arbitration, 1: use fwd only, 2: use sarsa only
mode.USE_BWDupdate_of_FWDmodel=1; % 1: use the backward update for goal-directed model (fwd model), 0: do not use
mode.DEBUG_Q_VALUE_CHG=0; % Debug option 1: show Q-value before/after whenever there is a goal change.
mode.path_ext=path0;
mode.total_simul=1; % # of total simulation repetition per subject
mode.simul_process_display=0; % 1: display model's process, 0: no diplay
mode.experience_sbj_events=[1 1]; % [pre main]  +1: experience exactly the same events(decision,state) as subjects. 0: model's own experience -1: use saved setting
mode.max_iter=100; % maximum iteration for optimization
% mode.out=1; % 1: normal evaluation mode, 99: regressor added to the SBJ, 0: debug mode

%% OPTION - Regressor arrangement
% original : {'SPE', 'RPE', 'uncertaintyM1', 'invFanoM1', 'Qsarsa','Qfwd','Qarb','dQbwdEnergy', 'weigtM1'}; reg_type_go_first=[1 1.5 2];
% Q-test : {'uncertaintyM1', 'invFanoM1', 'Qsarsa','Qfwd','Qarb','dQbwdEnergy','SPE', 'RPE', 'weigtM1'};   reg_type_go_first=[1.5 1 2];
% should add the regressors in the order of importance
% param_regressor_type_cue={'SPE', 'RPE', 'invFanoM2','invFanoM1','Qsarsa','Qfwd','Qarb','dQbwdEnergy', 'weigtM1'};

% Daw's version test
% param_regressor_type_cue={'RPE', 'SPE', 'MAXinvFano12', 'Qfwd','Qsarsa','Qarb', 'dQbwdEnergy', 'weigtM1'};

% full test
% param_regressor_type_cue={'SPE', 'RPE', 'uncertaintyM1', 'MAXinvFano12', 'Qarb', 'dQbwdEnergy', 'weigtM1'};
param_regressor_type_cue={'SPE', 'RPE', 'uncertaintyM1','MAXinvFano12','Qarb','Qsarsa','Qfwd'}; % 7 regressor
param_regressor_type_cue={'SPE', 'RPE', 'MAXinvFano12','Qarb','Qfwd','Qsarsa', 'weigtM1'}; % value
param_regressor_type_cue={'SPE', 'RPE', 'invFanoM2','invFanoM1','Qsarsa','Qfwd','Qarb','dQbwdEnergy', 'weigtM1'};
param_regressor_type_cue={'SPE', 'RPE','uncertaintyM1','MAXinvFano12','Qarb'}; % 5 regressor
param_regressor_type_cue={'SPE', 'RPE','MAXinvFano12'}; % 5 regressor
param_regressor_type_cue={'SPE', 'RPE','uncertaintyM1','MAXinvFano12','Qsarsa','Qfwd','Qarb'};
param_regressor_type_cue={'SPE','RPE','uncertaintyM1','MAXinvFano12','SPE_BL','Qsarsa','Qfwd','Qarb'};
param_regressor_type_cue={'SPE','RPE','MAXinvFano12','uncertaintyM1','SPE_BL','Qsarsa','Qfwd','Qarb'};


% value test
% param_regressor_type_cue={'SPE', 'RPE', 'invFanoM1','Qfwd','Qsarsa','Qarb','dQbwdEnergy', 'weigtM1'};


reg_type_go_first=[1 1.5 2]; %[1 1.5] [CAUTION] The order should match with 'param_regressor_type_cue'.   [CAUTION] type"2" should go always last!!!
% ## the below method does not work - figure out why and fix it.
% reg_type_go_first=[1.5 1 2]; % FOR Q-test only!!! (19 NOV)
Do_create_regressors=1;
Is_save_files_local=1; % save optimization parameters and regressor filesjob
Is_save_files_cluster=1; % save optimization parameters and regressor files


%% OPTION - behaviroal analysis & display
Do_behavioral_analysis=[0]; % [NOTE] use list_sbj_included=[2:1:14 16:1:19 21 23:1:24] for better results
if(Do_behavioral_analysis(1)==1) % dont need to create regressors in behavioral analysis mode!
    Do_create_regressors=0;
    Is_save_files_cluster=0; % it is a analysis in local pc, so do not update files in neuroecon.
end



%% initialization
if(Is_save_files_local==0)    disp('### files will not be saved to your local PC.');     end
if(Is_save_files_cluster==0)    disp('### files will not be saved to the cluster PC.');     end
    
use_model_regressor_cue=0;  ind_regressor_total=[];   type_regressor=[];    ind_regressor_total_in_design_mat=[];
for ii=1:1:size(param_regressor_type_cue,2) % collect regressor information
    if(strcmp(param_regressor_type_cue{1,ii},'SPE')==1)    use_model_regressor_cue=1;  ind_chk=1;   end
    if(strcmp(param_regressor_type_cue{1,ii},'RPE')==1)    use_model_regressor_cue=1;  ind_chk=2;   end
    if(strcmp(param_regressor_type_cue{1,ii},'uncertaintyM1')==1)    use_model_regressor_cue=1;    ind_chk=3;   end
    if(strcmp(param_regressor_type_cue{1,ii},'uncertaintyM2')==1)    use_model_regressor_cue=1;    ind_chk=4;   end
    if(strcmp(param_regressor_type_cue{1,ii},'meanM1')==1)    use_model_regressor_cue=1;    ind_chk=5;   end
    if(strcmp(param_regressor_type_cue{1,ii},'meanM2')==1)    use_model_regressor_cue=1;    ind_chk=6;   end
    if(strcmp(param_regressor_type_cue{1,ii},'invFanoM1')==1)    use_model_regressor_cue=1;    ind_chk=7;   end
    if(strcmp(param_regressor_type_cue{1,ii},'invFanoM2')==1)    use_model_regressor_cue=1;    ind_chk=8;   end
    if(strcmp(param_regressor_type_cue{1,ii},'weigtM1')==1)    use_model_regressor_cue=1;    ind_chk=9;   end
    if(strcmp(param_regressor_type_cue{1,ii},'weigtM2')==1)    use_model_regressor_cue=1;    ind_chk=10;   end
    if(strcmp(param_regressor_type_cue{1,ii},'Qfwd')==1)    use_model_regressor_cue=1;    ind_chk=11;   end
    if(strcmp(param_regressor_type_cue{1,ii},'Qsarsa')==1)    use_model_regressor_cue=1;    ind_chk=12;   end
    if(strcmp(param_regressor_type_cue{1,ii},'Qarb')==1)    use_model_regressor_cue=1;    ind_chk=13;   end
    if(strcmp(param_regressor_type_cue{1,ii},'dQbwdEnergy')==1)    use_model_regressor_cue=1;    ind_chk=14;   end
    if(strcmp(param_regressor_type_cue{1,ii},'dQbwdMean')==1)    use_model_regressor_cue=1;    ind_chk=15;   end
    if(strcmp(param_regressor_type_cue{1,ii},'duncertaintyM1')==1)    use_model_regressor_cue=1;    ind_chk=16;   end
    if(strcmp(param_regressor_type_cue{1,ii},'dinvFanoM1')==1)    use_model_regressor_cue=1;    ind_chk=17;   end
    if(strcmp(param_regressor_type_cue{1,ii},'TR_alpha')==1)    use_model_regressor_cue=1;    ind_chk=18;   end
    if(strcmp(param_regressor_type_cue{1,ii},'TR_beta')==1)    use_model_regressor_cue=1;    ind_chk=19;   end
    if(strcmp(param_regressor_type_cue{1,ii},'dinvFano12')==1)    use_model_regressor_cue=1;    ind_chk=20;   end
    if(strcmp(param_regressor_type_cue{1,ii},'ABSdinvFano12')==1)    use_model_regressor_cue=1;    ind_chk=21;   end
    if(strcmp(param_regressor_type_cue{1,ii},'MAXinvFano12')==1)    use_model_regressor_cue=1;    ind_chk=22;   end
    if(strcmp(param_regressor_type_cue{1,ii},'CONFLICTinvFano12')==1)    use_model_regressor_cue=1;    ind_chk=23;   end    
    if(strcmp(param_regressor_type_cue{1,ii},'PMB')==1)    use_model_regressor_cue=1;    ind_chk=24;   end    
    if(strcmp(param_regressor_type_cue{1,ii},'invFanoM1_meancorrected')==1)    use_model_regressor_cue=1;    ind_chk=25;   end
    if(strcmp(param_regressor_type_cue{1,ii},'invFanoM2_meancorrected')==1)    use_model_regressor_cue=1;    ind_chk=26;   end
    if(strcmp(param_regressor_type_cue{1,ii},'SPE_BL')==1)    use_model_regressor_cue=1;  ind_chk=29;   end
    if(strcmp(param_regressor_type_cue{1,ii},'RPE_BL')==1)    use_model_regressor_cue=1;  ind_chk=30;   end
    
    % index of regressor in "SBJ" structure
    ind_regressor_total=[ind_regressor_total ind_chk];
    % regressor type
    type_regressor=[type_regressor TYPE_REGRESSOR(ind_chk)];
end
% make a regressor index matrix for 1st parametric modulations (normal)
% (1) make 'param_regressor_type_cue_abs_pos_in_design_mat'
reg_cnt=0; param_regressor_type_cue_abs_pos_in_design_mat=[];
ind_regressor_type_base{1,1}.ind_reg=ind_regressor_total(find(type_regressor==reg_type_go_first(1)));
reg_cnt=reg_cnt+1+length(ind_regressor_type_base{1,1}.ind_reg);   param_regressor_type_cue_abs_pos_in_design_mat=[param_regressor_type_cue_abs_pos_in_design_mat [2:1:reg_cnt]];
ind_regressor_type_base{1,2}.ind_reg=ind_regressor_total(find(type_regressor==reg_type_go_first(2)));
reg_cnt=reg_cnt+1+length(ind_regressor_type_base{1,2}.ind_reg);   param_regressor_type_cue_abs_pos_in_design_mat=[param_regressor_type_cue_abs_pos_in_design_mat [param_regressor_type_cue_abs_pos_in_design_mat(end)+2:1:reg_cnt]];
if(reg_type_go_first(3)==1.7)
    ind_regressor_type_base{1,3}.ind_reg=ind_regressor_total(find(type_regressor==reg_type_go_first(3)));
    reg_cnt=reg_cnt+1+length(ind_regressor_type_base{1,3}.ind_reg);   param_regressor_type_cue_abs_pos_in_design_mat=[param_regressor_type_cue_abs_pos_in_design_mat [param_regressor_type_cue_abs_pos_in_design_mat(end)+2:1:reg_cnt]];
end

% make a regressor index for 2nd parametric modulations (dummy)
ind_regressor_type_dummy.ind_reg=ind_regressor_total(find(type_regressor==2));
reg_cnt=reg_cnt+length(ind_regressor_type_dummy.ind_reg);   param_regressor_type_cue_abs_pos_in_design_mat=[param_regressor_type_cue_abs_pos_in_design_mat [param_regressor_type_cue_abs_pos_in_design_mat(end)+1:1:reg_cnt]];
for j=1:1:length(ind_regressor_type_dummy.ind_reg)
    ind_regressor_type_dummy.name{1,j}=LIST_REGRESSOR{1,ind_regressor_type_dummy.ind_reg(j)};
end
% (2)
ind_regressor_type_base{1,1}.abs_pos_in_design_mat=param_regressor_type_cue_abs_pos_in_design_mat(find(type_regressor==reg_type_go_first(1)));
ind_regressor_type_base{1,2}.abs_pos_in_design_mat=param_regressor_type_cue_abs_pos_in_design_mat(find(type_regressor==reg_type_go_first(2)));
if(reg_type_go_first(3)==1.7)
    ind_regressor_type_base{1,3}.abs_pos_in_design_mat=param_regressor_type_cue_abs_pos_in_design_mat(find(type_regressor==reg_type_go_first(3)));
end
ind_regressor_type_dummy.abs_pos_in_design_mat=param_regressor_type_cue_abs_pos_in_design_mat(find(type_regressor==2));


if(sum(abs(param_regressor_type_cue_abs_pos_in_design_mat-sort(param_regressor_type_cue_abs_pos_in_design_mat,'ascend')))~=0)
    error('- ERROR!!!!: the variable ''param_regressor_type_cue_abs_pos_in_design_mat'' should be in ascending order!!!');
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
    
    SBJ{1,i}.map_type=LIST_sbj_map_type(ind_sbj_included(i));
    
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
SBJ_event=SBJ;



%% model optimization
% param_in(1): myArbitrator.PE_tolerance_m1 (m1's threshold for zero PE)
% param_in(2): myArbitrator.PE_tolerance_m2 (m2's threshold for zero PE)
% param_in(3): myArbitrator.A_12
% param_in(x): myArbitrator.B_12 : based on A12
% param_in(4): myArbitrator.A_21
% param_in(x): myArbitrator.B_21 : based on A21
% param_in(5): myArbitrator.tau_softmax/param_sarsa.tau/param_fwd.tau : better to fix at 0,2. This should be determined in a way that maintains softmax values in a reasonable scale. Otherwise, this will drive the fitness value!
% param_in(6): % param_sarsa.alpha/param_fwd.alpha 0.01~0.2 to ensure a good "state_fwd.T" in phase 1

param_init=[0.4, 12, 2.0, 12, 0.2, 0.1]; 
param_BoundL=[0.1, 10, 0.5*param_init(3:1:4), 0.15, 0.01]; %[0.01, 3, 1.0, 2.0, 0.5, 5, 0.01, 0.01];
param_BoundU=[1.0/2, 40/2, 3*param_init(3:1:4), 0.25, 0.2]; %[1, 20, 6, 10, 4, 30, 0.5, 0.2];
mode.boundary_12=0.1;       mode.boundary_21=0.01; % boundary condition : gating fn(1)

mode.param_init=param_init; mode.param_BoundL=param_BoundL; mode.param_BoundU=param_BoundU;


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
        myFunc_bu = @(x) eval_ArbitrationRL3(x, SBJ_test, mode); % define a new anonymous function : eval_ArbitrationRL2(x, SBJ_test, mode) for full BayesArb
        [model_BayesArb.param, model_BayesArb.val]=fminsearchbnd(myFunc_bu, param_init, param_BoundL, param_BoundU, optimset('Display','iter','MaxIter',mode.max_iter));   % X0,LB,UB
        % [2-1] add regressor vector to SBJ
        mode.out=99;
        SBJ_test=eval_ArbitrationRL3(model_BayesArb.param,SBJ_test,mode); %  : eval_ArbitrationRL2(x, SBJ_test, mode) for full BayesArb
        % [3] Save
        model_BayesArb.mode=mode;
        SBJ_test{1,1}.model_BayesArb=model_BayesArb;
        SBJ{1,ind_sbj}=SBJ_test{1,1};
        save_file_name=['SBJ_structure.mat'];
        if(Is_save_files_local==1)
            eval(['save ' save_path_result save_file_name ' SBJ'])
        end
    end
    option_optimizing_model=2; % and then write regressors to SBJ structure based on this optimized parameter
end

if(option_optimizing_model==1)
    % ## (way2-batch) optimizing for *all* subjects and plug the result into each SBJ structure
    % [0] retrieve intial configuration for skipping pre-training
%     SBJ_keep=SBJ;
%     load_file_name=['SBJ_structure(backup,batch,Oct30_4).mat'];
%     eval(['load ' save_path_result load_file_name]);
%     for ff1=1:1:length(SBJ_keep)
%         SBJ_keep{1,ff1}.init_state_fwd=SBJ{1,ff1}.init_state_fwd;    SBJ_keep{1,ff1}.init_state_sarsa=SBJ{1,ff1}.init_state_sarsa;
%     end
%     SBJ=SBJ_keep;
    % [1] model optimization
    mode.out=1;
    myFunc_bu = @(x) eval_ArbitrationRL3(x, SBJ, mode); % define a new anonymous function : eval_ArbitrationRL2(x, SBJ_test, mode) for full BayesArb
    bunch_step=2; bunch=[2:bunch_step:mode.max_iter];
    for mm=[1:1:length(bunch)-1] % split optimization into bunch and save everytime
        max_iter_current=bunch(mm+1)-bunch(mm);
        [model_BayesArb.param, model_BayesArb.val]=fminsearchbnd(myFunc_bu, param_init, param_BoundL, param_BoundU, optimset('Display','iter','MaxIter',max_iter_current));   % X0,LB,UB
        param_init=model_BayesArb.param;
        % [2-1] add regressor vector to SBJ
        mode.out=99;
        SBJ=eval_ArbitrationRL3(model_BayesArb.param,SBJ,mode); %  : eval_ArbitrationRL2(x, SBJ_test, mode) for full BayesArb
        % [3] save
        model_BayesArb.mode=mode;
        for ind_sbj=1:1:size(SBJ,2) % plug in a identical parameter (because this is batch)
            SBJ{1,ind_sbj}.model_BayesArb=model_BayesArb;
        end
        disp('-saving intermediate optimization result to SBJ structure...');
        save_file_name=['SBJ_structure.mat'];
        if(Is_save_files_local==1)
            eval(['save ' save_path_result save_file_name ' SBJ'])
        end
        if(Is_save_files_cluster==1)
            eval(['save ' save_path_result save_file_name ' SBJ'])
        end
    end    
    option_optimizing_model=2; % and then write regressors to SBJ structure based on this optimized parameter
end


if(option_optimizing_model==2) % [NOTE] replace initial SBJ = just read SBJ from the "SBJ_structure.mat"
        
    % backup plan
    load_file_name=['SBJ_structure.mat']; % DO NOT CHANGE THIS. THIS WILL NOT USED.
    eval(['load ' save_path_result load_file_name])    
    SBJ_backpu_plan=SBJ;
    clear SBJ
    % 
    load_file_name=name_paramfile_to_test; % [THIS ONE APPLIES TO DATA!!!!]    
    eval(['load ' save_path_result load_file_name])
    % regressor part deleting and regenerating.
    for ff=1:1:length(list_sbj_included)
        disp(sprintf('- writing regressor to SBJ structure (SBJ%02d)...',list_sbj_included(ff)));
        % find my subject in "SBJ" strucure of the 'SBJ_structure.mat' file        
        did_find=0;
        if(size(findstr(load_file_name,'batch'),1)~=0) % if reading batch file
            did_find=1;
            SBJ0{1,1}=SBJ{1,1}; % get a template            
            SBJ0{1,1}.name=LIST_SBJ_included{1,ff}; % apply a name
        else % each param case
            for ss=1:1:size(SBJ,2)
                if(strcmp(SBJ{1,ss}.name,LIST_SBJ_included{1,ff})==1)
                    SBJ0{1,1}=SBJ{1,ss}; % SBJ : this includes SBJ structure for subjects to be included for this code
                    did_find=did_find+1;
                end
            end
        end
        % apply batch parameter to 'Thao' whose parameter is extremely unstable.
%         if(strcmp(LIST_SBJ_included{1,ff},'Thao')==1) % find 'Thao' in batch structure >> his parameter is highly unstable and produce a constant invF.
%             disp('- applying backup plan (using batch param for Thao...')
%             for ss=1:1:size(SBJ,2)
%                 if(strcmp(SBJ{1,ss}.name,LIST_SBJ_included{1,ff})==1)
%                     SBJ0{1,1}=SBJ_backup_plan{1,ss}; % SBJ : this includes SBJ structure for subjects to be included for this code
%                     did_find=1;
%                 end
%             end
%         end
        if(did_find~=1)            error('-ERROR:: no correponding subject found in the "SBJ_structure.mat" file!!!');   end
        
        if(isfield(SBJ0{1,1}, 'regressor')==1)
            SBJ0{1,1}=rmfield(SBJ0{1,1},'regressor'); %remove the regressor field
        end
        mode.out=99;
        model_BayesArb.param=SBJ0{1,1}.model_BayesArb.param;
%         if(isfield(SBJ0{1,1}.model_BayesArb.mode, 'boundary_12')==1) % parameter length=6 case
        if(length(SBJ0{1,1}.model_BayesArb.param)==6) % parameter length=6 case
            mode.boundary_12=SBJ0{1,1}.model_BayesArb.mode.boundary_12;
            mode.boundary_21=SBJ0{1,1}.model_BayesArb.mode.boundary_21;
            mode.param_length=6;
        end        
        if(length(SBJ0{1,1}.model_BayesArb.param)==4) % parameter length=8 case
            mode.boundary_12=SBJ0{1,1}.model_BayesArb.mode.boundary_12; %n/a indicator
            mode.boundary_21=SBJ0{1,1}.model_BayesArb.mode.boundary_21; %n/a indicator
            mode.param_length=4;
        end
        
        if(strcmp(name_paramfile_to_test,'SBJ_structure_new2.mat')==1)
            SBJ0=eval_ArbitrationRL_DPON2(model_BayesArb.param,SBJ0,mode); % refresh and add the regressor part  : eval_ArbitrationRL2(x, SBJ_test, mode) for full BayesArb
        else
            SBJ0=eval_ArbitrationRL3c(model_BayesArb.param,SBJ0,mode); % refresh and add the regressor part  : eval_ArbitrationRL2(x, SBJ_test, mode) for full BayesArb
        end
        SBJ1{1,ff}=SBJ0{1,1};
    end
    clear SBJ
    SBJ=SBJ1;
    if(update_SBJ_structure==1)        eval(['save ' save_path_result load_file_name ' SBJ']);  end
end


%test mode
if(option_optimizing_model==3)
    model_BayesArb.param=param_init;
    mode.out=99;
    SBJ=eval_ArbitrationRL3(model_BayesArb.param,SBJ,mode); %  : eval_ArbitrationRL2(x, SBJ_test, mode) for full BayesArb
end




%%
% if(Do_behavioral_analysis(1)==1)
% 
%     %% correlation between two signals
%     out_corr_tmp=[];
%     for jj=1:1:size(SBJ,2)
%         [c_val lag]=xcorr(SBJ{1,jj}.regressor{1,7}.value(8,:),SBJ{1,jj}.regressor{1,8}.value(8,:),'coeff');
%         out_corr_tmp=[out_corr_tmp c_val(find(lag==3))]; % zero-lag point = correlation
%     end
% 
% 
%     %% Behavioral analysis #0 (computational)
% 
%     d_step_size=1;
%         
%     % (1) for comparing alpha with beta in the dynamic model    
%     Arb_alpha=cell(1,size(SBJ,2)); %MF2MB-based on the performance of MF
%     Arb_beta=cell(1,size(SBJ,2)); %MB2MF-based on the performance of MB
%     Arb_PMB=cell(1,size(SBJ,2)); % for SBJ    
%     contribution_alpha=cell(1,size(SBJ,2)); contribution_beta=cell(1,size(SBJ,2));
%     contribution_alpha_total=[];    contribution_beta_total=[];
%     for jj2=1:1:size(SBJ,2)        % each subject
%         Arb_alpha{1,jj2}=SBJ{1,jj2}.regressor{1,18}.value(7,:);
%         d_Arb_alpha{1,jj2}=Arb_alpha{1,jj2}(1+d_step_size:end)-Arb_alpha{1,jj2}(1:end-d_step_size);
%         Arb_beta{1,jj2}=SBJ{1,jj2}.regressor{1,19}.value(7,:);
%         d_Arb_beta{1,jj2}=Arb_beta{1,jj2}(1+d_step_size:end)-Arb_beta{1,jj2}(1:end-d_step_size);
%         Arb_PMB{1,jj2}=SBJ{1,jj2}.regressor{1,9}.value(7,:);
%         d_Arb_PMB{1,jj2}=Arb_PMB{1,jj2}(1+d_step_size:end)-Arb_PMB{1,jj2}(1:end-d_step_size);        
%         
%         ind_incl=find(mod([1:1:length(d_Arb_PMB{1,jj2})],3)~=0); % exclude btw trials
%         
%         mat_test=d_Arb_PMB{1,jj2}(ind_incl)./d_Arb_alpha{1,jj2}(ind_incl);
%         ind_contribution_alpha=~isnan(mat_test);
%         mat_test2=mat_test(ind_contribution_alpha);
%         ind_contribution_alpha2=isfinite(mat_test2);
%         contribution_alpha{1,jj2}=mat_test2(ind_contribution_alpha2);
%         contribution_alpha_total=[contribution_alpha_total contribution_alpha{1,jj2}];
%         
%         mat_test=d_Arb_PMB{1,jj2}./d_Arb_beta{1,jj2};
%         ind_contribution_beta=~isnan(mat_test);
%         mat_test2=mat_test(ind_contribution_beta);
%         ind_contribution_beta2=isfinite(mat_test2);
%         contribution_beta{1,jj2}=(-1)*mat_test2(ind_contribution_beta2); % positive contribution if decreasing MB2MF causes increasing PMB
%         contribution_beta_total=[contribution_beta_total contribution_beta{1,jj2}];
%     end
%     contribution_alpha_total=contribution_alpha_total(abs(contribution_alpha_total)<100);
%     contribution_beta_total=contribution_beta_total(abs(contribution_beta_total)<100);    
% %     figure('Name','contribution index of alpha/beta to PMB')
% %     hb=boxplot([contribution_alpha_total contribution_beta_total],[ones(1,length(contribution_alpha_total)) 2*ones(1,length(contribution_beta_total))]);
% %     set(hb(7,:),'Visible','off')
% %     axis([0.5 2.5 -1.1 1.1])
% %     [h p] = ttest2(contribution_alpha_total,contribution_beta_total);
% %     title(sprintf('two-sample ttest p=%1.2e, contribution(beta-alpha):%1.2e',p,[mean(contribution_beta_total)-mean(contribution_alpha_total)]))
% 
%     
% %      save Arb_PMB_MBinvFfixed1 Arb_PMB
% %       save Arb_PMB_full Arb_PMB
% %       break;
%     
%     
%     
%     % (2) for comparing invF_M1 with invF_M2 in the dynamic model    
%     Arb_invF_M2=cell(1,size(SBJ,2));    Arb_invF_M2_meancorrected=cell(1,size(SBJ,2)); %reliability of MF
%     Arb_invF_M1=cell(1,size(SBJ,2));    Arb_invF_M1_meancorrected=cell(1,size(SBJ,2)); %reliability of MB
%     Arb_invF_maxRel=cell(1,size(SBJ,2));    Arb_invF_diffRel=cell(1,size(SBJ,2));
%     All_condition=cell(1,size(SBJ,2));
%     Arb_PMB=cell(1,size(SBJ,2));
%     Arb_PMB_backup=cell(1,size(SBJ,2)); % for SBJ_backup_plan    
%     contribution_invF_M2=cell(1,size(SBJ,2)); contribution_invF_M1=cell(1,size(SBJ,2));
%     contribution_invF_M2_total=[];    contribution_invF_M1_total=[];    
%     for jj2=1:1:size(SBJ,2)        % each subject
%         % condition        
%         tmp_mat_c=[];
%         for pp2=1:1:size(SBJ{1,jj2}.HIST_block_condition,2)            
%             tmp_mat_c=[tmp_mat_c reshape(ones(3,1)*SBJ{1,jj2}.HIST_block_condition{1,pp2}(2,:),[1,3*length(SBJ{1,jj2}.HIST_block_condition{1,pp2}(2,:))])];         
%         end
%         All_condition{1,jj2}=tmp_mat_c;
%         % arb signals
%         Arb_invF_M2{1,jj2}=SBJ{1,jj2}.regressor{1,8}.value(8,:);
%         Arb_invF_M2_meancorrected{1,jj2}=SBJ{1,jj2}.regressor{1,26}.value(8,:);
%         d_Arb_invF_M2{1,jj2}=Arb_invF_M2{1,jj2}(1+d_step_size:end)-Arb_invF_M2{1,jj2}(1:end-d_step_size);
%         Arb_invF_M1{1,jj2}=SBJ{1,jj2}.regressor{1,7}.value(8,:);
%         Arb_invF_M1_meancorrected{1,jj2}=SBJ{1,jj2}.regressor{1,25}.value(8,:);
%         Arb_invF_maxRel{1,jj2}=SBJ{1,jj2}.regressor{1,22}.value(8,:);
%         Arb_invF_diffRel{1,jj2}=SBJ{1,jj2}.regressor{1,20}.value(8,:);
%         d_Arb_invF_M1{1,jj2}=Arb_invF_M1{1,jj2}(1+d_step_size:end)-Arb_invF_M1{1,jj2}(1:end-d_step_size);
%         Arb_PMB{1,jj2}=SBJ{1,jj2}.regressor{1,9}.value(7,:);
%         d_Arb_PMB{1,jj2}=Arb_PMB{1,jj2}(1+d_step_size:end)-Arb_PMB{1,jj2}(1:end-d_step_size);
%         Arb_PMB_backup{1,jj2}=SBJ_backup_plan{1,jj2}.regressor{1,9}.value(7,:);
%         
%         ind_incl=find((mod([1:1:length(d_Arb_PMB{1,jj2})],3)~=0)&(d_Arb_PMB{1,jj2}~=0)); % exclude btw trials & measure only if PMB changes!        
%         
%         mat_test=d_Arb_PMB{1,jj2}(ind_incl)./d_Arb_invF_M2{1,jj2}(ind_incl);        
%         ind_contribution_invF_M2=~isnan(mat_test);
%         mat_test2=mat_test(ind_contribution_invF_M2);
%         ind_contribution_invF_M22=isfinite(mat_test2);
%         contribution_invF_M2{1,jj2}=(-1)*mat_test2(ind_contribution_invF_M22); % [NOTE] positive contribution if decreasing invF_M2 causes increasing PMB
%         contribution_invF_M2_total=[contribution_invF_M2_total contribution_invF_M2{1,jj2}];
%         
%         mat_test=d_Arb_PMB{1,jj2}(ind_incl)./d_Arb_invF_M1{1,jj2}(ind_incl);
%         ind_contribution_invF_M1=~isnan(mat_test);
%         mat_test2=mat_test(ind_contribution_invF_M1);
%         ind_contribution_invF_M12=isfinite(mat_test2);
%         contribution_invF_M1{1,jj2}=mat_test2(ind_contribution_invF_M12);
%         contribution_invF_M1_total=[contribution_invF_M1_total contribution_invF_M1{1,jj2}];
%         
%     end
%     
%    %% 
%     % example reliability trace
%     for lll2=[1]%[22:-1:1]
%         sbj_id_show=lll2; %1,17
%         hFig=figure('Name',['reliability trace' sprintf('-sbj#%d',sbj_id_show)]);
%         set(hFig, 'Position', [500 500 800 100])
%         tt=[1:1:length(All_condition{1,sbj_id_show})];
%         tt_show1=[2:3:length(All_condition{1,sbj_id_show})];    tt_show2=[3:3:length(All_condition{1,sbj_id_show})];
%         tt_show=sort([tt_show1 tt_show2]);
%         All_condition_show=All_condition{1,sbj_id_show};
%         [r_tmp c_goal]=find((All_condition_show==1)|(All_condition_show==2));
%         goal_condi_show=ones(1,length(All_condition_show))*(-0.5); goal_condi_show(c_goal)=(-0.4);
%         [r_tmp c_uncert]=find((All_condition_show==2)|(All_condition_show==3));
%         uncert_condi_show=ones(1,length(All_condition_show))*(-0.5); uncert_condi_show(c_uncert)=(-0.4);
% %         % (1) full version (w/ condition)
%         plot(tt(tt_show),goal_condi_show(tt_show),'g-',tt(tt_show),uncert_condi_show(tt_show),'m-',...
%             tt(tt_show),Arb_invF_M1_meancorrected{1,sbj_id_show}(tt_show),'r-',...
%             tt(tt_show),Arb_invF_M2_meancorrected{1,sbj_id_show}(tt_show),'b-');
%         axis([1 length(tt(tt_show)) -.6 .8]);
%         legend('flexible|goal','lowU|highU','MBreliability','MFreliability');
%         % (2) simple version (w/o condition)
% %         plot(tt(tt_show),Arb_invF_M1_meancorrected{1,sbj_id_show}(tt_show),'r-',...
% %             tt(tt_show),Arb_invF_M2_meancorrected{1,sbj_id_show}(tt_show),'b-');
% %         axis([1 length(tt(tt_show)) -.6 .8]);%         legend('MBreliability','MFreliability');
%     end
%     %%
%     
%     % plot MBrel vs MFrel according to the condition
%     a_relM1=cell(size(SBJ,2),4);  a_relM2=cell(size(SBJ,2),4); % display order
%     a_relM1_all=cell(1,4);  a_relM2_all=cell(1,4); % display order    
%     a_ttest_table_btwGoal=[];    a_ttest_table_btwUncertainty=[];
%     for sbj_id_show=1:1:size(SBJ,2)
%         tt=[1:1:length(All_condition{1,sbj_id_show})];
%         tt_show1=[2:3:length(All_condition{1,sbj_id_show})];    tt_show2=[3:3:length(All_condition{1,sbj_id_show})];
%         tt_show=sort([tt_show1 tt_show2]);
%         All_condition_show=All_condition{1,sbj_id_show};
%         [r_tmp c_goal]=find((All_condition_show==1)|(All_condition_show==2));
%         goal_condi_show=ones(1,length(All_condition_show))*2; goal_condi_show(c_goal)=0; % specific:0, flexible:2
%         [r_tmp c_uncert]=find((All_condition_show==2)|(All_condition_show==3));
%         uncert_condi_show=ones(1,length(All_condition_show))*1; uncert_condi_show(c_uncert)=2; %lowU:1, highU:2
%         
%         a_gc=goal_condi_show(tt_show); % specific:0, flexible:1
%         a_uc=uncert_condi_show(tt_show); %lowU:0, highU:1
%         a_relM1_tmp=Arb_invF_M1_meancorrected{1,sbj_id_show}(tt_show);
%         a_relM2_tmp=Arb_invF_M2_meancorrected{1,sbj_id_show}(tt_show);        
%         for a_ii=1:1:length(a_gc)
%             a_condi=goal_condi_show(a_ii)+uncert_condi_show(a_ii);
%             a_relM1{sbj_id_show,a_condi}=[a_relM1{sbj_id_show,a_condi}, a_relM1_tmp(a_ii)];
%             a_relM2{sbj_id_show,a_condi}=[a_relM2{sbj_id_show,a_condi}, a_relM2_tmp(a_ii)];
%         end  
%         for b_ii=1:1:4
%             a_relM1_all{1,b_ii}=[a_relM1_all{1,b_ii}, a_relM1{sbj_id_show,b_ii}];
%             a_relM2_all{1,b_ii}=[a_relM2_all{1,b_ii}, a_relM2{sbj_id_show,b_ii}];
%             a_mean1_withinSBJ(sbj_id_show,b_ii)=mean(a_relM1{sbj_id_show,b_ii});
%             a_mean2_withinSBJ(sbj_id_show,b_ii)=mean(a_relM2{sbj_id_show,b_ii});            
%             a_sem1_withinSBJ(sbj_id_show,b_ii)=std(a_relM1{sbj_id_show,b_ii})/sqrt(length(a_relM1{sbj_id_show,b_ii}));
%             a_sem2_withinSBJ(sbj_id_show,b_ii)=std(a_relM2{sbj_id_show,b_ii})/sqrt(length(a_relM2{sbj_id_show,b_ii}));                            
%         end
%         % ttest2
%         [a_ttest_table_btwUncertainty(sbj_id_show,1) p_tmp]=ttest2(a_relM1{sbj_id_show,1},a_relM1{sbj_id_show,2});
%         [a_ttest_table_btwUncertainty(sbj_id_show,2) p_tmp]=ttest2(a_relM1{sbj_id_show,3},a_relM1{sbj_id_show,4});
%     end    
%     a_mean1=[]; a_sem1=[];    a_mean2=[]; a_sem2=[];    
%     a_mean1_randomEff=[]; a_sem1_randomEff=[];    a_mean2_randomEff=[]; a_sem2_randomEff=[];    
%     a_mean1_withinSBJ_disp=[];  a_sem1_withinSBJ_disp=[];
%     a_mean2_withinSBJ_disp=[];  a_sem2_withinSBJ_disp=[];
%     a_ind_row=[1 1 2 2];    a_ind_col=[1 2 1 2];
%     for b_ii=1:1:4
%         a_mean1(a_ind_row(b_ii),a_ind_col(b_ii))=mean(a_relM1_all{1,b_ii});
%         a_sem1(a_ind_row(b_ii),a_ind_col(b_ii))=std(a_relM1_all{1,b_ii})/sqrt(length(a_relM1_all{1,b_ii}));
%         a_mean2(a_ind_row(b_ii),a_ind_col(b_ii))=mean(a_relM2_all{1,b_ii});
%         a_sem2(a_ind_row(b_ii),a_ind_col(b_ii))=std(a_relM2_all{1,b_ii})/sqrt(length(a_relM2_all{1,b_ii}));
%         a_mean1_randomEff(a_ind_row(b_ii),a_ind_col(b_ii))=mean(a_mean1_withinSBJ(:,b_ii));
%         a_sem1_randomEff(a_ind_row(b_ii),a_ind_col(b_ii))=std(a_mean1_withinSBJ(:,b_ii))./sqrt(length(a_mean1_withinSBJ(:,b_ii)));
%         a_mean2_randomEff(a_ind_row(b_ii),a_ind_col(b_ii))=mean(a_mean2_withinSBJ(:,b_ii));
%         a_sem2_randomEff(a_ind_row(b_ii),a_ind_col(b_ii))=std(a_mean2_withinSBJ(:,b_ii))./sqrt(length(a_mean2_withinSBJ(:,b_ii)));        
%     end
%     
%     % display a particular sbj
%     for lll=[1]
%     sbj_id_display=lll;   disp_margin_factor=0.4;
%     for b_ii=1:1:4
%         a_mean1_withinSBJ_disp(a_ind_row(b_ii),a_ind_col(b_ii))=a_mean1_withinSBJ(sbj_id_display,b_ii);
%         a_mean2_withinSBJ_disp(a_ind_row(b_ii),a_ind_col(b_ii))=a_mean2_withinSBJ(sbj_id_display,b_ii);
%         a_sem1_withinSBJ_disp(a_ind_row(b_ii),a_ind_col(b_ii))=a_sem1_withinSBJ(sbj_id_display,b_ii);
%         a_sem2_withinSBJ_disp(a_ind_row(b_ii),a_ind_col(b_ii))=a_sem2_withinSBJ(sbj_id_display,b_ii);
%     end
%     hFig = figure('Name',['MB,MFreliability' sprintf('-sbj#%d',sbj_id_display)]); %%%    
%     set(hFig, 'Position', [500 500 400 250])
%     colormap(summer)
%     yscale_min1=min(min([a_mean1_withinSBJ_disp])); yscale_min2=min(min([a_mean2_withinSBJ_disp]));
%     yscale_max1=max(max([a_mean1_withinSBJ_disp])); yscale_max2=max(max([a_mean2_withinSBJ_disp]));
%     subplot(1,2,1);    
%     barwitherr(a_sem1_withinSBJ_disp, a_mean1_withinSBJ_disp);    % Plot with errorbars
%     set(gca,'XTickLabel',{'Specific','Flexible'})
% %     legend('low uncertainty','high uncertainty')    
%     yscale0=yscale_max1-yscale_min1; axis([0.5 2.5 yscale_min1-yscale0*disp_margin_factor yscale_max1+yscale0*disp_margin_factor]);
%     title('MBreliability') 
%     subplot(1,2,2);    
%     barwitherr(a_sem2_withinSBJ_disp, a_mean2_withinSBJ_disp);    % Plot with errorbars
%     set(gca,'XTickLabel',{'Specific','Flexible'})
% %     legend('low uncertainty','high uncertainty')    
%     yscale0=yscale_max2-yscale_min2; axis([0.5 2.5 yscale_min2-yscale0*disp_margin_factor yscale_max2+yscale0*disp_margin_factor]);
%     title('MFreliability') 
%     end
%     
%     
%     figure('Name','signal trace2')
%     sbj_id_show=1; %1,17
%     tt=[1:1:length(All_condition{1,sbj_id_show})];
%     tt_show1=[2:3:length(All_condition{1,sbj_id_show})];    tt_show2=[3:3:length(All_condition{1,sbj_id_show})];
%     tt_show=sort([tt_show1 tt_show2]);
%     All_condition_show=All_condition{1,sbj_id_show};
%     [r_tmp c_goal]=find((All_condition_show==1)|(All_condition_show==2));  
%     goal_condi_show=ones(1,length(All_condition_show))*(-0.5); goal_condi_show(c_goal)=(-0.4);
%     [r_tmp c_uncert]=find((All_condition_show==2)|(All_condition_show==3));  
%     uncert_condi_show=ones(1,length(All_condition_show))*(-0.5); uncert_condi_show(c_uncert)=(-0.4);
%     plot(tt(tt_show),goal_condi_show(tt_show),'g-',tt(tt_show),uncert_condi_show(tt_show),'m-',...
%         tt(tt_show),Arb_invF_maxRel{1,sbj_id_show}(tt_show),'r-',...
%         tt(tt_show),Arb_invF_diffRel{1,sbj_id_show}(tt_show),'b-');
%     axis([1 length(tt(tt_show)) -.6 .8]);
%     legend('flexible|goal','lowU|highU','Max reliability','Diff reliability');
%     
%     figure('Name','signal trace3')
%     sbj_id_show=1; %1,17
%     tt=[1:1:length(All_condition{1,sbj_id_show})];
%     tt_show1=[2:3:length(All_condition{1,sbj_id_show})];    tt_show2=[3:3:length(All_condition{1,sbj_id_show})];
%     tt_show=sort([tt_show1 tt_show2]);
%     All_condition_show=All_condition{1,sbj_id_show};
%     [r_tmp c_goal]=find((All_condition_show==1)|(All_condition_show==2));  
%     goal_condi_show=ones(1,length(All_condition_show))*(-0.5); goal_condi_show(c_goal)=(-0.4);
%     [r_tmp c_uncert]=find((All_condition_show==2)|(All_condition_show==3));  
%     uncert_condi_show=ones(1,length(All_condition_show))*(-0.5); uncert_condi_show(c_uncert)=(-0.4);
%     plot(tt(tt_show),goal_condi_show(tt_show),'g-',tt(tt_show),uncert_condi_show(tt_show),'m-',...
%         tt(tt_show),Arb_PMB{1,sbj_id_show}(tt_show),'k-');
% %     axis([1 length(tt(tt_show)) -.6 .8]);
%     legend('flexible|goal','lowU|highU','PMB');
%     
%     %
%     contribution_invF_M2_total=contribution_invF_M2_total(abs(contribution_invF_M2_total)<100);
%     contribution_invF_M1_total=contribution_invF_M1_total(abs(contribution_invF_M1_total)<100);    
%     figure('Name','contribution index of reliability of MB/MF to PMB')
% %     hpb=boxplot([contribution_invF_M2_total contribution_invF_M1_total],[ones(1,length(contribution_invF_M2_total)) 2*ones(1,length(contribution_invF_M1_total))]); % [no whisker] 'whisker',0
% %     set(hpb(7,:),'Visible','off')
% %     axis([0.5 2.5 -3e-3 3e-3])
%     mat_disp=[mean(contribution_invF_M2_total) mean(contribution_invF_M1_total)]; % bin x state
%     err_y=[std(contribution_invF_M2_total)/sqrt(length(contribution_invF_M2_total)) std(contribution_invF_M1_total)/sqrt(length(contribution_invF_M1_total))];
%     barwitherr(err_y', mat_disp');    
% %     axis([0.5 2.5 -0.05 0.15])
%     [h p] = ttest2(contribution_invF_M2_total,contribution_invF_M1_total);
%     [ht pt] = ttest(contribution_invF_M1_total);
%     [ht2 pt2] = ttest(contribution_invF_M2_total);
%     title(sprintf('mean diff (p=%1.2e), contri. of invF1(p=%1.2e), of invF2(p=%1.2e)',p,pt,pt2))
%     % contribution(M1-M2): [mean(contribution_invF_M1_total)-mean(contribution_invF_M2_total)]
% 
%     % correlation btw PMB & PMB(backup) / correlation between invF1 & invF2
%     out_corr_PMB_new_old=zeros(1,size(SBJ,2));  out_R2_PMB_new_old=zeros(2,size(SBJ,2)); %[r2; rmse]
%     out_corr_invF_1_2=zeros(1,size(SBJ,2));  corr_btw_invF1_invF2=zeros(2,size(SBJ,2));        
%     for jj2=1:1:size(SBJ,2)        % each subject        
%         
%         % xcorr:Normalized cross-correlation so the autocorrelations at zero lag are identically 1.0.
%         [c_val lag]=xcorr(Arb_PMB{1,jj2},Arb_PMB_backup{1,jj2},'coeff');            out_corr_PMB_new_old(jj2)=c_val(find(lag==0)); % zero-lag point = correlation
%         [out_R2_PMB_new_old(1,jj2) out_R2_PMB_new_old(2,jj2)] = rsquare(Arb_PMB{1,jj2},Arb_PMB_backup{1,jj2});
%         
%         % xcorr:Normalized cross-correlation so the autocorrelations at zero lag are identically 1.0.
%         trial_effective=find(mod([1:1:length(Arb_invF_M1{1,jj2})],3)~=1);
%         [r_tmp p_tmp]=corrcoef(Arb_invF_M1{1,jj2}(trial_effective)',Arb_invF_M2{1,jj2}(trial_effective)');        
%         corr_btw_invF1_invF2(:,jj2)=[r_tmp(1,2); p_tmp(1,2)];
%         [c_val lag]=xcorr(Arb_invF_M1{1,jj2}(trial_effective),Arb_invF_M2{1,jj2}(trial_effective),'coeff');            out_corr_invF_1_2(jj2)=c_val(find(lag==0)); % zero-lag point = correlation
%         
%     end
%     
%     
%     
%     figure('Name','Correlation btw PMB(fullBayes) and PMB(new)')    
%     bar(out_corr_PMB_new_old)
% %     subplot(1,3,1)
% %     boxplot(out_corr_PMB_new_old)
% %     subplot(1,3,2)
% %     plot(Arb_PMB{1,3},Arb_PMB_backup{1,3},'b.');
% %     subplot(1,3,3)
% %     plot(Arb_PMB{1,4},Arb_PMB_backup{1,4},'b.');
%     
%     figure('Name','Correlation btw PMB(fullBayes) and PMB(new)')    
%     subplot(1,3,1)
%     boxplot(out_corr_PMB_new_old)
%     subplot(1,3,2)
%     plot(Arb_PMB{1,3},Arb_PMB_backup{1,3},'b.');
%     subplot(1,3,3)
%     plot(Arb_PMB{1,4},Arb_PMB_backup{1,4},'b.');
%     
%     figure('Name','Correlation btw invF1 and invF2')        
%     bar(corr_btw_invF1_invF2(1,:))
% %     bar(out_corr_invF_1_2);
%     
%     
%     %% Behavioral analysis #1 (computational)
%     
%     disp_opt.sbj_ind_for_goal_directed_trace=8; % pick only one subject to display as an example
%     
%     all_y_goal_amount_in_goalblock=[];  all_y_goal_amount_in_habitualblock=[];
%     all_out_corr=[];
%     all_y_mat_disp=cell(1,size(SBJ,2));
%     all_vec_gd_mb=[];       all_vec_gd_mf=[];   all_vec_gd_mn=[];
%     
%     
%     
%     %% 00. how well subjects are doing, relative to how well they could be doing in this task. the proportion of optimal choices (Figure S5) in each condition.
%     % [1:G(.9,.1),2:G'(.5,.5),4:H'(.9,.1),3:H(.5,.5)] 
%     all_opt_choice_rate=[];
%     all_opt_choice_rate_each_cond=[]; % [G,G',H',H]
%     all_RT_each_cond=[]; 
%     for jj2=1:1:size(SBJ,2)        % each subject
%         
%         % 00-1. proportiona of optimal choices
%         list_opt_action_block_identity=[[1 6]; [1 7]; [1 8];... % G
%             [2 6]; [2 7]; [2 8];... % G'
%             [3 -1]; [4 -1]; ]; % H and H'
%         list_opt_state_action{1,1}=[[1 2 3 4 5]' [2 -99 -99 2 2]']; % G40 [state action], 0: both are right, -99: both wrong
%         list_opt_state_action{1,2}=[[1 2 3 4 5]' [0 1 2 1 1]']; % G20 [state action]
%         list_opt_state_action{1,3}=[[1 2 3 4 5]' [1 2 1 -99 -99]']; % G10 [state action]
%         list_opt_state_action{1,4}=[[1 2 3 4 5]' [2 -99 -99 0 2]']; % G'40 [state action]
%         list_opt_state_action{1,5}=[[1 2 3 4 5]' [0 1 2 1 1]']; % G'20 [state action]
%         list_opt_state_action{1,6}=[[1 2 3 4 5]' [1 0 1 -99 -99]']; % G'10 [state action]
%         list_opt_state_action{1,7}=[[1 2 3 4 5]' [2 1 2 1 2]']; % H [state action]
%         list_opt_state_action{1,8}=[[1 2 3 4 5]' [2 1 2 2 1]']; % H' [state action]
%         val_score_total=0;  val_score=0;    
%         val_score_total_each_cond=zeros(1,4);  val_score_each_cond=zeros(1,4);        % G, G', H, H'
%         if(SBJ{1,jj2}.map_type==2) % new map only
%             nn_session=size(SBJ{1,jj2}.HIST_block_condition,2);
%             for ff2=1:1:nn_session
%                 nn_trial=size(SBJ{1,jj2}.HIST_block_condition{1,ff2},2);
%                 for qq2=1:1:nn_trial
%                     % identify the block condition
%                     id_blck=[SBJ{1,jj2}.HIST_block_condition{1,ff2}(2,qq2), SBJ{1,jj2}.HIST_behavior_info{1,ff2}(qq2,end)];
%                     tmp_m=sum(abs(list_opt_action_block_identity-ones(8,1)*id_blck),2);
%                     [r_idb c_idb]=find(tmp_m==0);
%                     % set the optimal state-action set
%                     opt_action_set=list_opt_state_action{1,r_idb};
%                     opt_action=[opt_action_set(1,2) opt_action_set(SBJ{1,jj2}.HIST_behavior_info{1,ff2}(qq2,5),2)]; % 1x2
%                     tmp_score=length(find((opt_action-SBJ{1,jj2}.HIST_behavior_info{1,ff2}(qq2,7:8))==0))+length(find(opt_action==0));                    
%                     tmp_score2=length(find((opt_action-SBJ{1,jj2}.HIST_behavior_info{1,ff2}(qq2,7:8))==0))+length(find(opt_action==0));                    
%                     val_score=val_score+tmp_score;                    
%                     val_score_total=val_score_total+2; % 2 actions each time
%                     if(r_idb<=3)    id_condi=1; end % G
%                     if((r_idb>3)&&(r_idb<=6))    id_condi=2; end % G'
%                     if(r_idb==7)    id_condi=4; end % H (highU)(change according to the order of the display)
%                     if(r_idb==8)    id_condi=3; end % H' (lowU)(change according to the order of the display)
%                     % (G,G',H',H) 
%                     val_score_each_cond(id_condi)=val_score_each_cond(id_condi)+tmp_score;   
%                     val_score_total_each_cond(id_condi)=val_score_total_each_cond(id_condi)+2;
%                 end
%             end
%             all_opt_choice_rate=[all_opt_choice_rate val_score/val_score_total];
%             all_opt_choice_rate_each_cond=[all_opt_choice_rate_each_cond; val_score_each_cond./val_score_total_each_cond];
%         end
%         
%         
%         % 00-2. RT
%         val_RT_totalcnt_each_cond=zeros(1,4);    val_RT_each_cond=zeros(1,4);        % G, G', H, H'
%         if(SBJ{1,jj2}.map_type==2) % new map only
%             nn_session=size(SBJ{1,jj2}.HIST_block_condition,2);
%             for ff2=1:1:nn_session
%                 nn_trial=size(SBJ{1,jj2}.HIST_block_condition{1,ff2},2);
%                 for qq2=1:1:nn_trial
%                     % identify the block condition
%                     id_blck=[SBJ{1,jj2}.HIST_block_condition{1,ff2}(2,qq2), SBJ{1,jj2}.HIST_behavior_info{1,ff2}(qq2,end)];
%                     tmp_m=sum(abs(list_opt_action_block_identity-ones(8,1)*id_blck),2);
%                     [r_idb c_idb]=find(tmp_m==0);
%                     % read reaction time
%                     if(r_idb<=3)    id_condi=1; end % G
%                     if((r_idb>3)&&(r_idb<=6))    id_condi=2; end % G'
%                     if(r_idb==7)    id_condi=4; end % H (change according to the order of the display)
%                     if(r_idb==8)    id_condi=3; end % H' (change according to the order of the display)
%                     % keep RT (G,G',H',H)
%                     val_RT_each_cond(id_condi)=val_RT_each_cond(id_condi)+sum(SBJ{1,jj2}.HIST_behavior_info{1,ff2}(qq2,9:10));
%                     val_RT_totalcnt_each_cond(id_condi)=val_RT_totalcnt_each_cond(id_condi)+2;
%                 end
%             end
%             all_RT_each_cond=[all_RT_each_cond; val_RT_each_cond./val_RT_totalcnt_each_cond];
%         end
%                 
%     end
%        
%     
%        
%     % bar plot
%     figure('Name','proportion of optimal choices');
%     colormap(summer)    
%     y_plot1=[mean(all_opt_choice_rate_each_cond(:,1:2)); mean(all_opt_choice_rate_each_cond(:,3:4))]; % (2 groups of 2 parameters)
% %     errY1 = zeros(2,2,2);    
%     errY1=[std(all_opt_choice_rate_each_cond(:,1:2))/sqrt(length(all_opt_choice_rate_each_cond(:,1)));...
%         std(all_opt_choice_rate_each_cond(:,3:4))/sqrt(length(all_opt_choice_rate_each_cond(:,1)))]; % SEM: error    
%     barwitherr(errY1, y_plot1);    % Plot with errorbars
%     set(gca,'XTickLabel',{'Specific Goal','Flexible goal'})
%     legend('low uncertainty','high uncertainty')
%     ylabel('proportion of optimal choices')
%     current_axis=axis;
%     axis([current_axis(1:2) 0 1]);
%     [h p]=ttest2(all_opt_choice_rate_each_cond(:,3),all_opt_choice_rate_each_cond(:,4));%ttest
%     
%     % bar plot
%     figure('Name','reaction time');
%     colormap(summer)    
%     y_plot1=[mean(all_RT_each_cond(:,1:2)); mean(all_RT_each_cond(:,3:4))]; % (2 groups of 2 parameters)
% %     errY1 = zeros(2,2,2);
%     errY1=[std(all_RT_each_cond(:,1:2))/sqrt(length(all_RT_each_cond(:,1)));...
%         std(all_RT_each_cond(:,3:4))/sqrt(length(all_RT_each_cond(:,1)))]; % SEM: error
%     barwitherr(errY1, y_plot1);    % Plot with errorbars
%     set(gca,'XTickLabel',{'Specific Goal','Flexible goal'})
%     legend('low uncertainty','high uncertainty')
%     ylabel('reaction time (sec)')
%     current_axis=axis;
% %     axis([current_axis(1:2) 0 1]);
%     [h p]=ttest2(all_RT_each_cond(:,3),all_RT_each_cond(:,4));%ttest
%     
%     
%     
%     %%
%     for jj2=1:1:size(SBJ,2)        % each subject
%         
%               
%         
%         %% 0. degree of goal-directed (Likelihood-ratio) "how much subjects are goal-directed?"
%         % we can use a likelihood ratio as an index, but we cannot do the
%         % likelihood ratio test because fwd and sarsa model are different
%         % model (i.e., one is not the special case of the other)        
%         % Test statistics of likelihood-ratio test
%         y_mat=2*[SBJ{1,jj2}.model_error{1,2}.value(7,:)-SBJ{1,jj2}.model_error{1,1}.value(7,:)]; % how much goal-directed (0:neutral) likelihood ratio
%         valid_col=find(~isnan(y_mat)); % y_mat is for state1/2 only; all state 3 has NaN, so eliminate them.
%         y_mat=y_mat(valid_col);
% %         y_mat=(y_mat-min(y_mat))/(max(y_mat)-min(y_mat));   y_mat_disp=(y_mat-0.5)*2;
%         y_mat=y_mat/max(abs(y_mat)); %[-1,1]  for display
%         y_mat_disp=y_mat;
%         all_y_mat_disp{1,jj2}=y_mat_disp;
%         
%         % 1: goal block, -1: habitual block
%         y_block_condi=(SBJ{1,jj2}.model_error{1,2}.value(5,valid_col)~=-1)-0.5; %[-.5,.5]
%         y_block_condi_goal=(y_block_condi>0); %[1,0]
%         y_block_condi_habit=(-1)*(y_block_condi<=0); %[0,-1]
%         % trace of degree of goal-directed (here, only for one subject as an example)
%         x_mat=[1:1:length(y_mat_disp)];        
%         if(jj2==disp_opt.sbj_ind_for_goal_directed_trace)
%             txt_disp=['deg. of goal-directed: ' '''' SBJ{1,jj2}.name ''''];
%             figure('Name',txt_disp);   plot(x_mat,y_mat_disp,'k-',x_mat,y_block_condi_goal,'r-.',x_mat,y_block_condi_habit,'b-.',x_mat,zeros(1,length(x_mat)),'-'); axis([min(x_mat) max(x_mat) -2 2])
%         end
%         % amount of goal-directed in goal-directed blocks and in habitual blocks
%         y_goal_amount_in_goalblock=sum(y_mat(find(y_block_condi_goal==1))); % sum of NegLogLik
%         y_goal_amount_in_habitualblock=sum(y_mat(find(y_block_condi_habit==-1))); % sum of NegLogLik
%         all_y_goal_amount_in_goalblock=[all_y_goal_amount_in_goalblock; y_goal_amount_in_goalblock];
%         all_y_goal_amount_in_habitualblock=[all_y_goal_amount_in_habitualblock; y_goal_amount_in_habitualblock];
%          
%         
%         
%         %% 1. correlation btw degree of goal-driected and each regressor
%         num_regressors=size(SBJ{1,jj2}.regressor,2);
%         list_testing_regressor=[1 2 3 4 7 8 9]; % index in SBJ.regressor : spe,rpe,U1,U2,invF1,invF2,W1
%         i_r0=0; vec_all=[];
%         vec_a=y_mat_disp;
%         for i_r=list_testing_regressor(1:2) % SPE, RPE regressor
%             i_r0=i_r0+1;
%             vec_b=SBJ{1,jj2}.regressor{1,i_r}.value(row_mat(i_r),valid_col);
%             vec_b=(vec_b-min(vec_b))/(max(vec_b)-min(vec_b))-0.5; % [-0.5,0.5]
%             list_legend{i_r0,1}=SBJ{1,jj2}.regressor{1,i_r}.name;
%             [c_val lag]=xcorr(vec_a,vec_b,'coeff');            out_corr(i_r0)=c_val(find(lag==0)); % zero-lag point = correlation
%         end
%         % U1 - U2
%         i_r0=i_r0+1;    list_legend{i_r0,1}='Uncertainty(M1-M2)';
%         vec_b=SBJ{1,jj2}.regressor{1,list_testing_regressor(3)}.value(row_mat(list_testing_regressor(3)),valid_col)-...
%             SBJ{1,jj2}.regressor{1,list_testing_regressor(4)}.value(row_mat(list_testing_regressor(4)),valid_col);
%         vec_b=vec_b/max(abs(vec_b)); %[-1,1]
%         [c_val lag]=xcorr(vec_a,vec_b,'coeff');        out_corr(i_r0)=c_val(find(lag==0));
%         % invF1 - invF2
%         i_r0=i_r0+1;    list_legend{i_r0,1}='invFano(M1-M2)';
%         vec_b=SBJ{1,jj2}.regressor{1,list_testing_regressor(5)}.value(row_mat(list_testing_regressor(5)),valid_col)-...
%             SBJ{1,jj2}.regressor{1,list_testing_regressor(6)}.value(row_mat(list_testing_regressor(6)),valid_col);
%         vec_b=vec_b/max(abs(vec_b)); %[-1,1]s 
%         [c_val lag]=xcorr(vec_a,vec_b,'coeff');        out_corr(i_r0)=c_val(find(lag==0));
%         % W1
%         i_r0=i_r0+1;
%         i_r=list_testing_regressor(7);
%         vec_b=SBJ{1,jj2}.regressor{1,i_r}.value(row_mat(i_r),valid_col); % [0,1] because this is probability
%         vec_b=2*vec_b-1; % [-1,1]
%         list_legend{i_r0,1}='Prob_{M1}';%SBJ{1,jj2}.regressor{1,i_r}.name;
%         [c_val lag]=xcorr(vec_a,vec_b,'coeff');            out_corr(i_r0)=c_val(find(lag==0)); % zero-lag point = correlation
%         % collecting values
%         all_out_corr=[all_out_corr; out_corr];
%         
%         
%         %% 1.1. P_MB vs  degree of goal-directed (Likelihood-ratio)
%         % have NaN if it did not find a match
%         vec_gd_mb=vec_a(find(vec_b>0.0)); % deg of goal-directed when PMB>0.5
%         all_vec_gd_mb=[all_vec_gd_mb mean(vec_gd_mb)];
%         vec_gd_mn=vec_a(find((vec_b>-0.4)&(vec_b<=0.4))); % deg of goal-directed when PMB middle
%         all_vec_gd_mn=[all_vec_gd_mn mean(vec_gd_mn)]; 
%         vec_gd_mf=vec_a(find(vec_b<=-0.0)); % deg of goal-directed when PMB<0.5
%         all_vec_gd_mf=[all_vec_gd_mf mean(vec_gd_mf)];
%         
%         
%     end
%     
%     % Removing outliers: all_vec_gd_mf
%     q1=prctile(all_vec_gd_mf,25);   q3=prctile(all_vec_gd_mf,75);   w=1.5; % default
%     thresh_outlier_h= q3 + w*(q3 - q1);     thresh_outlier_l=q1 -w*(q3-q1);
%     all_vec_gd_mf=all_vec_gd_mf(find((all_vec_gd_mf<=thresh_outlier_h)&(all_vec_gd_mf>=thresh_outlier_l)));
%     % Removing outliers: all_vec_gd_mb
%     q1=prctile(all_vec_gd_mb,25);   q3=prctile(all_vec_gd_mb,75);   w=1.5; % default
%     thresh_outlier_h= q3 + w*(q3 - q1);     thresh_outlier_l=q1 -w*(q3-q1);
%     all_vec_gd_mb=all_vec_gd_mb(find((all_vec_gd_mb<=thresh_outlier_h)&(all_vec_gd_mb>=thresh_outlier_l)));
%     
%     figure('Name','degree of goal-directed sbj to model''s prediction(Pmb) #1'); % barplot
%     hb=boxplot([all_vec_gd_mf all_vec_gd_mb],[ones(1,length(all_vec_gd_mf)) 2*ones(1,length(all_vec_gd_mb))]);
% %     hb=boxplot([all_vec_gd_mf all_vec_gd_mn all_vec_gd_mb],[ones(1,length(all_vec_gd_mf)) 2*ones(1,length(all_vec_gd_mn)) 3*ones(1,length(all_vec_gd_mb))]);
%     set(hb(7,:),'Visible','off')
%     %         axis([0.5 2.5 -1.1 1.1])
%     [h p] = ttest2(all_vec_gd_mf,all_vec_gd_mb);
%     [ht pt] = ttest(all_vec_gd_mf);
%     [ht pt2] = ttest(all_vec_gd_mb);
%     title(sprintf('mean diff p=1.2%e',p))
%     
%     figure('Name','degree of goal-directed sbj to model''s prediction(Pmb) #2'); % barplot
% %     mat_disp=[mean(all_vec_gd_mf) mean(all_vec_gd_mn) mean(all_vec_gd_mb)]; % bin x state
% %     err_y=[std(all_vec_gd_mf)/2 std(all_vec_gd_mn)/2 std(all_vec_gd_mb)/2];
%     mat_disp=[mean(all_vec_gd_mf) mean(all_vec_gd_mb)]; % bin x state
%     err_y=[std(all_vec_gd_mf)/sqrt(length(all_vec_gd_mf)) std(all_vec_gd_mb)/sqrt(length(all_vec_gd_mb))]; % SEM: Standard error of the mean
%     barwitherr(err_y', mat_disp');
%     title(sprintf('mean diff p=1.2%e | mf~=0 p=1.2%e | mb~=0 p=1.2%e',p,pt,pt2))
%     axis([0.5 2.5 -0.05 0.15])
%         
%     
%         
%         
%     figure('Name','degree of goal-directed in goal-directed/habit blocks');
%     % (1) error bar plot
% %     % mean of NegLogLik test
% %     y_plot1=[mean(all_y_goal_amount_in_goalblock) mean(all_y_goal_amount_in_habitualblock)];
% %     % 95% confidence interval
% %      [h,p_fwd,ci_fwd] = ttest(all_y_goal_amount_in_goalblock);  [h,p_sarsa,ci_sarsa] = ttest(all_y_goal_amount_in_habitualblock);
% %     errY1=[sum(abs(ci_fwd)) sum(abs(ci_sarsa))]/2;    
% %     barwitherr(errY1, y_plot1);
% %     set(gca,'XTickLabel',{'goal block','habit block'})
% %     ylabel('deg of goal-directed')
%     % (2) boxlot
%     y_plot1=[(all_y_goal_amount_in_goalblock) (all_y_goal_amount_in_habitualblock)];
%     [h_fwd,p_fwd,ci_fwd] = ttest(all_y_goal_amount_in_goalblock);  [h_sarsa,p_sarsa,ci_sarsa] = ttest(all_y_goal_amount_in_habitualblock);
%     pass_fwd='';    pass_sarsa='';
%     if(h_fwd==1)    pass_fwd='*'; end % add '*' if reject the null hypothesis at the 5% significance level. 
%     if(h_sarsa==1)    pass_sarsa='*'; end % add '*' if reject the null hypothesis at the 5% significance level. 
%     disp_txt{1,1}=['goal block' sprintf('\n') '(' pass_fwd sprintf('p=%1.1e)',p_fwd)]; disp_txt{1,2}=['habit block' sprintf('\n') '(' pass_sarsa sprintf('p=%1.1e)',p_sarsa)];
%     [h_diff,p_diff,ci_diff] = ttest(all_y_goal_amount_in_goalblock-all_y_goal_amount_in_habitualblock);
%     pass_diff='';
%     if(h_diff==1)    pass_diff='*'; end % add '*' if reject the null hypothesis at the 5% significance level.     
%     disp_title{1,1}=['goal block - control block' sprintf('\n') '(' pass_diff sprintf('p=%1.1e)',p_diff)];    
%     boxplot(y_plot1,disp_txt)    
%     ylabel('deg of goal-directed')
%     title(disp_title{1,1});
%     
%     figure('Name','Corr. btw (degree of goal-directed) and regressors'); % [note] error bar indicates std (all subjects)    
% %     % (1) error bar plot
% %     y_plot1=mean(all_out_corr);
% %     errY1=[];
% %     for j=1:1:size(all_out_corr,2)
% %         [h,p,ci_corr] = ttest(all_out_corr(:,j)-mean(all_out_corr(:,j))); % 95% confidence interval
% %         errY1=[errY1 sum(abs(ci_corr))/2];
% %     end
% %     barwitherr(errY1, y_plot1);
% %     set(gca,'XTickLabel',list_legend)
% %     ylabel('correlation')
%     % (2) boxlot    
%     clear p;
%     y_plot1=[all_out_corr];
%     for j=1:1:size(all_out_corr,2)
%         [h,p(j),ci_corr] = ttest(all_out_corr(:,j)); % 95% confidence interval        
%         pass_txt='';
%         if(h==1)    pass_txt='*'; end % add '*' if reject the null hypothesis at the 5% significance level. 
%         list_legend_disp{1,j}=[list_legend{j,1}, sprintf('\n') '(' pass_txt sprintf('p=%1.1e)',p(j))];
%     end
%     boxplot(y_plot1,list_legend_disp);    
%     ylabel('correlation')
%     
%     
%     %% Behavioral analysis #2 (intuitive version) - "how often sbj switches his choice in state4 in habitual block?"
%     if(1)
%         all_beh_perc=cell(1,size(SBJ,2));
%         all_beh_perc_determined=cell(1,size(SBJ,2));
%         cond_indi_list=[1, 6, 7, 8];
%         for jj2=1:1:size(SBJ,2)        % each subject
%             all_beh_perc0=zeros(4,5); % percentage of left choice in 4 conditions(h,g10,g20,g40) x 5states. thus '0' for always-left, '1' for always-right
%             all_beh_perc0_cnt=zeros(4,5);
%             num_sess=size(SBJ{1,jj2}.HIST_behavior_info,2);
%             for kk2=1:1:num_sess % each session
%                 mat_work0=SBJ{1,jj2}.HIST_behavior_info{1,kk2};
%                 for ii2=1:1:size(mat_work0,1) % each trial
%                     condi0=mat_work0(ii2,end); % condition
%                     for aa2=1:1:2 % 1st/2nd decision state
%                         i_c=find(cond_indi_list==condi0); % condition
%                         i_s=mat_work0(ii2,3+aa2); % state
%                         i_a=mat_work0(ii2,6+aa2); % action
%                         all_beh_perc0_cnt(i_c,i_s)=all_beh_perc0_cnt(i_c,i_s)+1;
%                         if(i_a==2) % +1 if right choice is made
%                             all_beh_perc0(i_c,i_s)=all_beh_perc0(i_c,i_s)+1;
%                         end
%                     end
%                 end
%             end
%             all_beh_perc0_cnt(find(all_beh_perc0_cnt==0))=1; % make the denominator 1 for states never visited.
%             all_beh_perc0=all_beh_perc0./all_beh_perc0_cnt;% normalized and percentized
%             all_beh_perc_determined0=abs(all_beh_perc0-0.5)*2; % 1: keep pressing the same button, 0: switching equally
%             
%             all_beh_perc{1,jj2}=all_beh_perc0;
%             all_beh_perc_determined{1,jj2}=all_beh_perc_determined0;
%         end
%         % collect habitual block's assessment
%         box0=[];
%         for jj2=1:1:size(SBJ,2)        % each subject
%             %         box0=[box0; all_beh_perc{1,jj2}(1,4)];
%             box0=[box0; all_beh_perc_determined{1,jj2}(1,4)];
%         end
%         figure('Name','percent of the consistent choice in state4 in habitual blocks');
%         [h,p(j),ci_corr] = ttest(box0); % 95% confidence interval
%         pass_txt='';
%         if(h==1)    pass_txt='*'; end % add '*' if reject the null hypothesis at the 5% significance level.
%         list_legend0{1,1}=[list_legend{j,1}, sprintf('\n') '(' pass_txt sprintf('p=%1.1e)',p(j))];
%         boxplot(box0,list_legend0)        
% %         barwitherr(var(box0), mean(box0));    % Plot with errorbars
% %         legend(list_legend0)
%         colormap(flipud(gray)) % flipped gray map
%         axis([0.5 1.5 0 1.1])
%     end
%     
%     
%     
%     %% Behavioral analysis #3 (how much arbitration signal reflects subject's actual choice switching?)
%     list_r_ind=[9];%[1:1:10]; % index of regressor to look (SPE,RPE,U1,U2,m1,m2,invF1,invF2,w1,w2)
%     list_r_test_row=[7 7 8 8 8 8 8 8 7 7]; % corresponding row in the regressor vector
%     list_r_test_trial_s={[2 3]',[2 3]',[2 3]',[2 3]',[2 3]',[2 3]',[2 3]',[2 3]',[1 2 3]',[1 2 3]'}; % trial_s set in which regressor signal to be extracted
%     list_block_test=[1 2 3 4]; % [1 2] - which block condition included for this analysis
%     list_block_test_next=[2 1 4 3]; % [2 1] - the block condition that should come right after the 'list_block_test'.
%     minimum_num_sample_for_analysis=10; %10
%     show_major_state_only=1; % show state 1&4 only
%     
%     % 1:percentile-based bin, 2:equal-sized bin (NOTE:percentile-based method is not ideal option because the corresponding threshold is driven by distribution bias!)
%     option_bin=2;
%     
%     list_percentile=[0 50 100]; % for option1
%     bin_size=3; % for option2, = # of threshold of bins (=# of bins +1)
%     list_bin_val=[0 .5 1]; % for option 3
%     
%     
%     % 1. obtain percentile-threshold from the regressor vector
%     all_r_val=cell(size(SBJ,2),length(list_r_ind));
%     all_r_threshold_percentile=cell(size(SBJ,2),length(list_r_ind));    all_r_threshold_equal_bin=cell(size(SBJ,2),length(list_r_ind)); all_r_threshold_specified_bin=cell(size(SBJ,2),length(list_r_ind));
%     for jj2=1:1:size(SBJ,2)        % each subject        
%         for rr2=1:1:length(list_r_ind) % each regressor
% %             disp(sprintf('[sbj-%02d/%02d] obtainning the threshold from the regressor vector (%02d/%02d)',jj2, size(SBJ,2),rr2,length(list_r_ind)));
%             % (1) threshold for percentile-based bind
%             line_ind=SBJ{1,jj2}.regressor{1,list_r_ind(rr2)}.value(4,:);
%             if(length(list_r_test_trial_s{1,list_r_ind(rr2)})==1)
%                 col_all0=find((line_ind==list_r_test_trial_s{1,list_r_ind(rr2)}(1))); % extract all corresponding trial_s's column
%             end
%             if(length(list_r_test_trial_s{1,list_r_ind(rr2)})==2)
%                 col_all0=find((line_ind==list_r_test_trial_s{1,list_r_ind(rr2)}(1))|(line_ind==list_r_test_trial_s{1,list_r_ind(rr2)}(2))); % extract all corresponding trial_s's column
%             end
%             if(length(list_r_test_trial_s{1,list_r_ind(rr2)})==3)
%                 col_all0=find((line_ind==list_r_test_trial_s{1,list_r_ind(rr2)}(1))|(line_ind==list_r_test_trial_s{1,list_r_ind(rr2)}(2))|(line_ind==list_r_test_trial_s{1,list_r_ind(rr2)}(3))); % extract all corresponding trial_s's column
%             end
%             all_r_val{jj2,rr2}=SBJ{1,jj2}.regressor{1,list_r_ind(rr2)}.value(list_r_test_row(list_r_ind(rr2)),col_all0);            
%             all_r_threshold_percentile{jj2,list_r_ind(rr2)}=prctile(all_r_val{jj2,rr2},list_percentile);
%             % (2) threshold for equal-sized bin
%             all_r_threshold_equal_bin{jj2,rr2}=[min(all_r_val{jj2,rr2}):(max(all_r_val{jj2,rr2})-min(all_r_val{jj2,rr2}))/(bin_size-1):max(all_r_val{jj2,rr2})];
%             % (3) specified bin
%             all_r_threshold_specified_bin{jj2,rr2}=list_bin_val;
%         end
%     end
%     
%     %2. measure the deg of choice switching 
%     if(option_bin==1)        all_r_threshold_use=all_r_threshold_percentile;    end
%     if(option_bin==2)        all_r_threshold_use=all_r_threshold_equal_bin;     end
%     if(option_bin==3)        all_r_threshold_use=all_r_threshold_specified_bin;     end
%     all_beh_perc_Arb=cell(length(list_r_ind),size(SBJ,2));
%     all_beh_perc_determined_Arb=cell(length(list_r_ind),size(SBJ,2));
%     all_beh_empty_cell_index=cell(length(list_r_ind),size(SBJ,2));
%     all_beh_choice_cnt=cell(length(list_r_ind),1);    all_beh_sample_cnt=cell(length(list_r_ind),1);    all_beh_perc_determined_Arb0=cell(length(list_r_ind),1);% just collect for all sbjs
%     % define the block set in which you measure choice consistency    
%     list_ref_blk=[list_block_test;  list_block_test_next];
%     for jj2=1:1:size(SBJ,2)        % each subject            
%         num_sess=size(SBJ{1,jj2}.HIST_behavior_info,2);
%         SBJ{1,jj2}.block_set=cell(1,num_sess);        
%         for kk2=1:1:num_sess % each session
%             tmp_mat_test=abs(SBJ{1,jj2}.HIST_block_condition{1,kk2}(2,2:end)-SBJ{1,jj2}.HIST_block_condition{1,kk2}(2,1:end-1));
%             [ind_chg0]=find(tmp_mat_test~=0); % means there is a change in col btw [ind_chg(x) ind_chg(x)+1].
%             ind_chg=[0 ind_chg0];
%             ind_valid=zeros(1,length(ind_chg));
%             for gg2=2:1:length(ind_chg)
%                 condi_chk=SBJ{1,jj2}.HIST_block_condition{1,kk2}(2,[ind_chg(gg2) ind_chg(gg2)+1]);
%                 oo=find(sum(abs(list_ref_blk-repmat(condi_chk',[1 size(list_ref_blk,2)])))==0);
%                 if(length(oo)>0) % valid change
%                     ind_valid(gg2)=1;
%                 end
%             end
%             ms_on=0;
%             block_set=[];
%             for gg2=[1:1:length(ind_chg)-1] % get start index and end index of block sets to be measured
%                 if((ind_valid(gg2)==0)&(ind_valid(gg2+1)==1)) % index start
%                     ind_st0=ind_chg(gg2)+1;                    ms_on=1;
%                 end
%                 if((ind_valid(gg2)==1)&(ind_valid(gg2+1)==0)) % index end
%                     ind_end0=ind_chg(gg2+1);                    ms_on=0;
%                     block_set=[block_set; [ind_st0 ind_end0]];
%                 end
%             end
%             SBJ{1,jj2}.block_set{1,kk2}=block_set; % # of testing blocks x 2. row: testing block index, col(1/2): starting/ending trial index of the testing block
%         end
%     end
%     
%     all_consistency_for_all_regressor_all_sbj=cell(size(SBJ,2),length(list_r_ind));
%     all_consistency_for_all_regressor=cell(length(list_r_ind),1);
%     for rr2=1:1:length(list_r_ind) % each regressor
%         all_beh_choice_cnt{rr2,1}=zeros(length(all_r_threshold_use{jj2,rr2})-1,5);        all_beh_sample_cnt{rr2,1}=zeros(length(all_r_threshold_use{jj2,rr2})-1,5);
%         all_beh_perc_Arb_name{1,rr2}=SBJ{1,1}.regressor{1,list_r_ind(rr2)}.name;                
%         all_consistency_r=cell(length(all_r_threshold_use{jj2,rr2})-1,5);
%         disp(sprintf('- measure the deg of choice switching based on percentile-threshold of regressor (%02d/%02d)...',rr2,length(list_r_ind)));
%         
%         for jj2=1:1:size(SBJ,2)        % each subject            
%             all_consistency=cell(length(all_r_threshold_use{jj2,rr2})-1,5);
%             all_beh_perc0=zeros(length(all_r_threshold_use{jj2,rr2})-1,5); % percentage of left choice in conditions(# of bins of regressor) x 5states. thus '0' for always-left, '1' for always-right
%             all_beh_perc0_cnt=zeros(length(all_r_threshold_use{jj2,rr2})-1,5);
%             all_beh_empty_ind=ones(length(all_r_threshold_use{jj2,rr2})-1,5); % zero if empty
%             num_sess=size(SBJ{1,jj2}.HIST_behavior_info,2);            
%             for kk2=1:1:num_sess % each session
%                 mat_work0=SBJ{1,jj2}.HIST_behavior_info{1,kk2};                
%                 for pp2=1:1:size(SBJ{1,jj2}.block_set{1,kk2},1) % each block set
%                     for ii2=[SBJ{1,jj2}.block_set{1,kk2}(pp2,1):1:SBJ{1,jj2}.block_set{1,kk2}(pp2,2)] % within each block in which consistency is computed
%                         for aa2=1:1:2 % 1st/2nd decision state
%                             % identify the condition (=bin# of the regressor signal)
%                             tmp_id=[kk2 mat_work0(ii2,1) mat_work0(ii2,2) list_r_test_trial_s{1,rr2}(aa2)]'; % [session block trial state]
%                             mycol=find(sum(abs(SBJ{1,jj2}.regressor{1,list_r_ind(rr2)}.value(1:4,:)-repmat(tmp_id,[1 length(SBJ{1,jj2}.regressor{1,list_r_ind(rr2)}.value(1,:))])))==0);
%                             val_r0=SBJ{1,jj2}.regressor{1,list_r_ind(rr2)}.value(list_r_test_row(list_r_ind(rr2)),mycol);
%                             mycol1=find((all_r_threshold_use{jj2,rr2}-val_r0)<=0);
%                             i_c=mycol1(end); % my condition = my bin#
%                             if(i_c==length(all_r_threshold_use{jj2,rr2}))                                i_c=i_c-1;                            end % the biggest value
%                             i_s=mat_work0(ii2,3+aa2); % state
%                             i_a=mat_work0(ii2,6+aa2); % action
%                             all_beh_perc0_cnt(i_c,i_s)=all_beh_perc0_cnt(i_c,i_s)+1;
%                             if(i_a==2) % +1 if right choice is made
%                                 all_beh_perc0(i_c,i_s)=all_beh_perc0(i_c,i_s)+1;
%                             end
%                         end
%                     end
%                     % compute choice consistency within each block set
%                     [r_valid c_valid]=find(all_beh_perc0_cnt>=minimum_num_sample_for_analysis); %  do not include the cell whose # of sample <20
%                     tmp1=all_beh_perc0./max(all_beh_perc0_cnt,1); % range [0,1]
%                     tmp2=abs(tmp1-0.5)*2; % range [0,1] : 0 means inconsistent and 1 means highly consistent.
%                     for tt2=1:1:length(r_valid) % collect choice consistency value
%                         all_consistency{r_valid(tt2),c_valid(tt2)}=[all_consistency{r_valid(tt2),c_valid(tt2)} tmp2(r_valid(tt2),c_valid(tt2))];
%                         all_consistency_r{r_valid(tt2),c_valid(tt2)}=[all_consistency_r{r_valid(tt2),c_valid(tt2)} tmp2(r_valid(tt2),c_valid(tt2))];
%                     end
%                 end                
%             end
%             all_consistency_for_all_regressor_all_sbj{jj2,rr2}=all_consistency; % all_consistency_for_all_regressor_all_sbj{sbj#,regressor#}{mybin#,state#}
%             
%         end
%         all_consistency_for_all_regressor{rr2,1}=all_consistency_r;
%         
% 
%         
%     end
%    
%      all_beh_perc_determined_Arb0=all_consistency_for_all_regressor;
%                 
%     % display
%      ind_r_disp=list_r_ind;%[1 2 3 5 7 9];%[1 2 3 4 5 6 7 8 9 10]; % index of regressor to look (SPE,RPE,U1,U2,m1,m2,invF1,invF2,w1,w2)
%      ind_state_disp=[1:1:5]; % state 1~5
%      label_bin=cell(1,(bin_size-1));
%      bin_ind=round([0:100/(bin_size-1):100]);     
%      for cc2=1:1:(bin_size-1)             label_bin{1,cc2}=sprintf('level%d',cc2);             end
% %      end
%      for rr2=1:1:length(ind_r_disp) % each regressor
%          
%          
%          figure('Name',['choice consistency -' all_beh_perc_Arb_name{1,rr2}]);
%          for ss2=1:1:length(ind_state_disp) % each state
%              label_disp{1,ss2}=sprintf('state%d',ind_state_disp(ss2));
%          end   
%          % ttest
%          ttest_str=cell(1,5); % 5 states
%          for aa2=1:1:(bin_size-1)
%              for bb2=aa2+1:1:(bin_size-1)
%                  for ss2=1:1:5      
%                      condi_valid=length(all_beh_perc_determined_Arb0{rr2,1}{aa2,ss2})>minimum_num_sample_for_analysis;
%                      if(condi_valid)
%                          [h,p] = ttest2(all_beh_perc_determined_Arb0{rr2,1}{aa2,ss2},all_beh_perc_determined_Arb0{rr2,1}{bb2,ss2}); % 95% confidence interval
%                          pass_txt='';
%                          if(h==1)    pass_txt='*'; end % add '*' if reject the null hypothesis at the 5% significance level.
%                          ttest_str{1,ss2}= [ttest_str{1,ss2} [sprintf('state%d(Lv%d-Lv%d) ',ss2,aa2,bb2) pass_txt sprintf('p=%1.1e',p)]];
%                      else
%                          ttest_str{1,ss2}='N/A (insufficient sample size)';
%                      end
%                  end
%              end
%          end
%          % compute mean consistency
%          mean_cons=zeros(size(all_beh_perc_determined_Arb0{rr2,1},1),size(all_beh_perc_determined_Arb0{rr2,1},2));
%          std_cons=zeros(size(all_beh_perc_determined_Arb0{rr2,1},1),size(all_beh_perc_determined_Arb0{rr2,1},2));
%          ind_valid_state=[];
%          for dd3=1:1:size(all_beh_perc_determined_Arb0{rr2,1},2) % state
%              condi_valid=1;
%              for dd2=1:1:size(all_beh_perc_determined_Arb0{rr2,1},1) % if any of bin in the state has small sample size, do not display                 
%                  condi_valid=condi_valid&(length(all_beh_perc_determined_Arb0{rr2,1}{dd2,dd3})>minimum_num_sample_for_analysis);
%              end             
%              if(condi_valid==1)             ind_valid_state=[ind_valid_state dd3];             end
%              for dd2=1:1:size(all_beh_perc_determined_Arb0{rr2,1},1) % bin                              
%                  if(condi_valid)                     
%                      mean_cons(dd2,dd3)=mean(all_beh_perc_determined_Arb0{rr2,1}{dd2,dd3});
%                      std_cons(dd2,dd3)=std(all_beh_perc_determined_Arb0{rr2,1}{dd2,dd3})/sqrt(length(all_beh_perc_determined_Arb0{rr2,1}{dd2,dd3})); % SEM
%                  else
%                       mean_cons(dd2,dd3)=0; std_cons(dd2,dd3)=0;
%                  end
%              end
%          end         
%          if(show_major_state_only==1)
%              ind_valid_state=[1 4 5]; % show only major state
%          end
%          mat_disp=mean_cons(:,ind_valid_state); % bin x state         
%          err_y=std_cons(:,ind_valid_state);         
%          barwitherr(err_y', mat_disp');          
%          for vv2=1:1:length(ind_valid_state)             label_disp2{1,vv2}=label_disp{1,ind_valid_state(vv2)};         end
%          set(gca,'XTickLabel',label_disp2)         
%          legend(label_bin);         
%          xlabel('state')
%          ylabel('choice consistency')
%          axis([0.5 length(ind_valid_state)+.5 0 1])
%          colormap(flipud(summer)) % flipped gray map
%          % add a text box
%          ttest_str1=[];
%          for ss2=1:1:size(ttest_str,2)             ttest_str1=[ttest_str1 sprintf('\n') ttest_str{1,ss2}];         end
%          annotation('textbox', [.6 .6 .2 .3], 'String', ttest_str1); % [x y w h] creates an editable text box annotation with its lower left corner at the point x,y, a width w, and a height h, specified in normalized figure units.
%          
%          
%          figure('Name',['Choice consistency (in major states 1,4,5) -' all_beh_perc_Arb_name{1,rr2}]);
%          state_interest=[1 4 5];
%          % compute mean consistency
%          data_collect=cell(1,(bin_size-1)); % consistency across bin         
%          mean_cons2=zeros((bin_size-1),1);   std_cons2=zeros((bin_size-1),1);
%          for dd2=1:1:(bin_size-1) % bin
%              for dd3=[state_interest] % state             
%                  data_collect{1,dd2}=[data_collect{1,dd2} all_beh_perc_determined_Arb0{rr2,1}{dd2,dd3}];
%              end
%              mean_cons2(dd2)=mean(data_collect{1,dd2});
%              std_cons2(dd2)=std(data_collect{1,dd2})/sqrt(length(data_collect{1,dd2})); % SEM
%          end
%          % ttest         
%          ss2_i=0;   clear ttest_str2
%          for aa2=1:1:(bin_size-1)
%              for bb2=(aa2+1):1:(bin_size-1)
%                  [h,p] = ttest2(data_collect{1,aa2},data_collect{1,bb2}); % 95% confidence interval
%                  pass_txt='';
%                  if(h==1)    pass_txt='*'; end % add '*' if reject the null hypothesis at the 5% significance level.
%                  ss2_i=ss2_i+1;
%                  ttest_str2{1,ss2_i}= [sprintf('(Lv%d-Lv%d) ',aa2,bb2) pass_txt sprintf('p=%1.1e',p)];
%              end
%          end
%          mat_disp=mean_cons2; % bin
%          err_y=std_cons2;         
%          barwitherr(err_y, mat_disp);          
%          set(gca,'XTickLabel',label_bin)                  
%          xlabel('signal level')
%          ylabel('choice consistency')
%          axis([0.5 bin_size-.5 0 1])
%          colormap(flipud(summer)) % flipped gray map
%          annotation('textbox', [.7 .7 .2 .2], 'String', ttest_str2); % [x y w h] creates an editable text box annotation with its lower left corner at the point x,y, a width w, and a height h, specified in normalized figure units.
%          
%      end
%      
%      
%      
%      %% Behavioral+model analysis #3-1 (how much arbitration signal reflects subject's actual choice behavior in terms of deg of goal-directed?)
%      % [] measure deg of goal-directed for each regressor bin
%     % extract state vector
%     all_s_mat=cell(1,size(SBJ,2));
%     for jj2=1:1:size(SBJ,2)    
%         s_mat=[];
%         for kk2=1:1:size(SBJ{1,jj2}.HIST_behavior_info,2)
%             for hh2=1:1:size(SBJ{1,jj2}.HIST_behavior_info{1,kk2},1)
%                 s_mat=[s_mat SBJ{1,jj2}.HIST_behavior_info{1,kk2}(hh2,[4 5])];
%             end
%         end
%         all_s_mat{1,jj2}=s_mat;
%     end    
%     all_deg_goal_Arb0=cell(length(list_r_ind),1);     all_deg_goal_Arb0_cnt=cell(length(list_r_ind),1);
%     for rr2=1:1:length(list_r_ind) % each regressor
%         adg=zeros(length(all_r_threshold_use{jj2,rr2})-1,5);    adg_cnt=zeros(length(all_r_threshold_use{jj2,rr2})-1,5); 
%         for jj2=1:1:size(SBJ,2)        % each subject   
%             i0_crt=0;
%             for kk2=1:1:size(SBJ{1,jj2}.HIST_behavior_info,2) % each session
%                 mat_work0=SBJ{1,jj2}.HIST_behavior_info{1,kk2};
%                 for ii2=1:1:size(SBJ{1,jj2}.HIST_behavior_info{1,kk2},1) % each trial
%                     for ss2=1:1:2 % each state (1,2)
%                         tmp_id=[kk2 mat_work0(ii2,1) mat_work0(ii2,2) ss2]'; % [session block trial state]
%                         mycol=find(sum(abs(SBJ{1,jj2}.regressor{1,list_r_ind(rr2)}.value(1:4,:)-repmat(tmp_id,[1 length(SBJ{1,jj2}.regressor{1,list_r_ind(rr2)}.value(1,:))])))==0);
%                         val_r0=SBJ{1,jj2}.regressor{1,list_r_ind(rr2)}.value(list_r_test_row(list_r_ind(rr2)),mycol); % regressor value                        
%                         mycol1=find((all_r_threshold_use{jj2,rr2}-val_r0)<=0); % binned
%                         i_c=mycol1(end); % my condition = my bin#
%                         if(i_c==length(all_r_threshold_use{jj2,rr2}))                                i_c=i_c-1;                            end % the biggest value
%                         i_s=mat_work0(ii2,3+ss2); % state                        
%                         i0_crt=i0_crt+1;
%                         adg(i_c,i_s)=adg(i_c,i_s)+all_y_mat_disp{1,jj2}(i0_crt);
%                         adg_cnt(i_c,i_s)=adg_cnt(i_c,i_s)+1;
%                     end
%                 end
%             end
%         end
%         all_deg_goal_Arb0{rr2,1}=adg./adg_cnt;
%         all_deg_goal_Arb0_cnt{rr2,1}=adg_cnt;
%     end
%     
%          
% end
%%




%% Create regressors
% state. 0.5: fixation mark on, 1: S1, 2: S2, 3: S3, 4: S4, 5: S5,
% 6(+/-)0.1: O1(with win/lost msg), 7(+/-)0.1: O2(with win/lost msg), 8(+/-)0.1: O3(with win/lost msg), 9: O4,
% 10:A1, 11:A2, 20: a short blank page display, -99:fail to choose in time limit, (-) when display off

if(Do_create_regressors==1)
    
    for jj2=1:1:size(SBJ,2)        % each subject
        disp(sprintf('##### creating regressor structures (sbj%02d/%02d) #######',jj2,size(SBJ,2)));
        SBJ{1,jj2}.HIST_behavior_info_Tag = SBJ_event{1,jj2}.HIST_behavior_info_Tag;
        SBJ{1,jj2}.HIST_event_info = SBJ_event{1,jj2}.HIST_event_info;
        SBJ{1,jj2}.HIST_event_info_Tag = SBJ_event{1,jj2}.HIST_event_info_Tag;
        SBJ{1,jj2}.HIST_block_condition_Tag = SBJ_event{1,jj2}.HIST_block_condition_Tag;
        for kk2=1:1:size(SBJ{1,jj2}.HIST_behavior_info,2)  % each main session
            
            
            
            mat_work=SBJ{1,jj2}.HIST_event_info{1,kk2};
            num_tot_events=size(SBJ{1,jj2}.HIST_behavior_info{1,kk2},1);
            
            
            
            %% 1. Regressor cue presentation - with parametric modulation (timing: stimulus onset)
            % [0-duration] SPE, RPE
            % [RT-duration] Q_fwd, Q_sarsa
            
            % regressor generation for each main session and save it to a single file that is compatible with SPM
            ind_reg=0;  % corresponds to the size of the structure
            ind_reg_abs=0; % actual number of regressors (including parametric)
            durations={};
            onsets={};
            names={};
            pmod=struct('name',{},'param',{},'poly',{});
            
            
            use_model_regressor_cue=1;
            
            % (1) durations, name, onset
            [tmp col_on]=find((mat_work(7,:)==1)|(mat_work(7,:)==2)|(mat_work(7,:)==3)|(mat_work(7,:)==4)|(mat_work(7,:)==5)...
                |(mat_work(7,:)==5.9)|(mat_work(7,:)==6.1)|(mat_work(7,:)==6.9)|(mat_work(7,:)==7.1)|...
                (mat_work(7,:)==7.9)|(mat_work(7,:)==8.1)|(mat_work(7,:)==9));
            if(length(col_on)~=num_tot_events*3)
                error('-ERROR: variable ''mat_work'' missed some event extraction. check!')
            end
            
            RT_mat=[];  onset_mat=[];
            prev_trial=0;  show_n_th_times_t=0;
            param_mat{1,1}=zeros(length(ind_regressor_type_base{1,1}.ind_reg),length(col_on));
            param_mat{1,2}=zeros(length(ind_regressor_type_base{1,2}.ind_reg),length(col_on));
            
            for ll2=1:1:length(col_on)
                
                if(ll2<length(col_on)) % usual case
                    pt_on=mat_work(4,col_on(ll2));
                    col_off=col_on(ll2)+1;
                    %                 col_off=col_on(ll2)-1+find(mat_work(7,[col_on(ll2):1:(col_on(ll2)+2)])==0.5); % find the next fixation mark presentation
                    pt_off=mat_work(4,col_off);
                    RT=pt_off-pt_on;
                else % last event in the session is the outcome presentation
                    RT=2.0;
                end
                RT_mat=[RT_mat RT];
                onset_t=mat_work(4,col_on(ll2));
                onset_mat=[onset_mat onset_t];
                
                % fill out regresssor values
                if(use_model_regressor_cue==1)
                    
                    % regressor type1
                    for nn=1:1:length(ind_regressor_type_base{1,1}.ind_reg)
                        mysession=kk2;
                        myblock=mat_work(1,col_on(ll2)); % block in session
                        mytrial=mat_work(2,col_on(ll2)); % trial in block
                        mytrial_s=mat_work(3,col_on(ll2)); % trial_s in trial (arbitration index: 1 at the second stage, 2 at the third stage)
                        mat_work_reg=SBJ{1,jj2}.regressor{1,ind_regressor_type_base{1,1}.ind_reg(nn)}.value(1:4,:);
                        identity_tmp=sum(abs(mat_work_reg-repmat([mysession myblock mytrial mytrial_s]',1,size(mat_work_reg,2))));
                        col_event=find(identity_tmp==0);
                        param_mat{1,1}(nn,ll2)=SBJ{1,jj2}.regressor{1,ind_regressor_type_base{1,1}.ind_reg(nn)}.value(row_mat(ind_regressor_type_base{1,1}.ind_reg(nn)),col_event);
                        
                    end
                    
                    
                    % regressor type2
                    for nn=1:1:length(ind_regressor_type_base{1,2}.ind_reg)
                        mysession=kk2;
                        myblock=mat_work(1,col_on(ll2)); % block in session
                        mytrial=mat_work(2,col_on(ll2)); % trial in block
                        mytrial_s=mat_work(3,col_on(ll2)); % trial_s in trial
                        mat_work_reg=SBJ{1,jj2}.regressor{1,ind_regressor_type_base{1,2}.ind_reg(nn)}.value(1:4,:);
                        identity_tmp=sum(abs(mat_work_reg-repmat([mysession myblock mytrial mytrial_s]',1,size(mat_work_reg,2))));
                        col_event=find(identity_tmp==0);
                        param_mat{1,2}(nn,ll2)=SBJ{1,jj2}.regressor{1,ind_regressor_type_base{1,2}.ind_reg(nn)}.value(row_mat(ind_regressor_type_base{1,2}.ind_reg(nn)),col_event);
                    end
                    
                end
                
            end
            
            
            % [[1st regressors]]
            ind_reg=ind_reg+1;
            onsets{1,ind_reg}=onset_mat;
            if(reg_type_go_first(1)==1)
                names{1,ind_reg}=['Cue_{0T}'];  durations{1,ind_reg}=0;
            end
            if(reg_type_go_first(1)==1.5)
                names{1,ind_reg}=['Cue_{1T}'];    durations{1,ind_reg}=RT_mat;
            end
            ind_reg_abs=1;      list_name_for_contrast{1,ind_reg_abs}=names{1,ind_reg}; % add to the global list of regresssors

            % (2) pmod: how many times each cue presented
            if(use_model_regressor_cue==1)
                ind_reg_param=0;
                for nn=1:1:length(ind_regressor_type_base{1,1}.ind_reg)
                    ind_reg_param=ind_reg_param+1;
                    pmod(1,ind_reg).name{1,ind_reg_param}=[SBJ{1,jj2}.regressor{1,ind_regressor_type_base{1,1}.ind_reg(ind_reg_param)}.name];
                    pmod(1,ind_reg).poly{1,ind_reg_param}=1;
                    pmod(1,ind_reg).param{1,ind_reg_param}=param_mat{1,1}(ind_reg_param,:);
                    ind_reg_abs=ind_regressor_type_base{1,1}.abs_pos_in_design_mat(ind_reg_param);      list_name_for_contrast{1,ind_reg_abs}=pmod(1,ind_reg).name{1,ind_reg_param};
                end
            end
            
            
            
            
            % [[2nd regressors]]
            ind_reg=ind_reg+1;
            onsets{1,ind_reg}=onset_mat;
            if(reg_type_go_first(2)==1)
                names{1,ind_reg}=['Cue_{0T}'];  durations{1,ind_reg}=0;
            end
            if(reg_type_go_first(2)==1.5)
                names{1,ind_reg}=['Cue_{1T}'];    durations{1,ind_reg}=RT_mat;
            end
            ind_reg_abs=ind_regressor_type_base{1,1}.abs_pos_in_design_mat(end)+1;      list_name_for_contrast{1,ind_reg_abs}=names{1,ind_reg}; % add to the global list of regresssors
            
            % (2) pmod:
            if(use_model_regressor_cue==1)
                ind_reg_param=0;
                for nn=1:1:length(ind_regressor_type_base{1,2}.ind_reg)
                    ind_reg_param=ind_reg_param+1;
                    pmod(1,ind_reg).name{1,ind_reg_param}=[SBJ{1,jj2}.regressor{1,ind_regressor_type_base{1,2}.ind_reg(ind_reg_param)}.name];
                    pmod(1,ind_reg).poly{1,ind_reg_param}=1;
                    pmod(1,ind_reg).param{1,ind_reg_param}=param_mat{1,2}(ind_reg_param,:);
                    ind_reg_abs=ind_regressor_type_base{1,2}.abs_pos_in_design_mat(ind_reg_param);      list_name_for_contrast{1,ind_reg_abs}=pmod(1,ind_reg).name{1,ind_reg_param};
                end
            end
            
            % Reserve spaces for extra regressors (will come from the dummy design matrix)
            if(reg_type_go_first(3)==2)
                for nn=1:1:length(ind_regressor_type_dummy.ind_reg)
                    ind_reg_param=ind_reg_param+1;
                    pmod(1,ind_reg).name{1,ind_reg_param}=ind_regressor_type_dummy.name{1,nn};
                    pmod(1,ind_reg).poly{1,ind_reg_param}=1;
                    pmod(1,ind_reg).param{1,ind_reg_param}=rand(1,length(param_mat{1,1}(1,:)));
                    ind_reg_abs=ind_regressor_type_dummy.abs_pos_in_design_mat(nn);      list_name_for_contrast{1,ind_reg_abs}=pmod(1,ind_reg).name{1,ind_reg_param};
                end
            end
            
            
            onset_event=onsets;
            
            
            %% 1-A. Regressor cue presentation - regressor type 3 only!!! with parametric modulation (timing: decision time)            
            % [0-duration] 
            % [RT-duration]
            
            if(reg_type_go_first(3)==1.7)
                
                use_model_regressor_cue=1;
                
                % (1) durations, name, onset
                [tmp col_on]=find((mat_work(7,:)==10)|(mat_work(7,:)==11));
                
                RT_mat=[];  onset_mat=[];
                prev_trial=0;  show_n_th_times_t=0;
                param_mat3{1,1}=zeros(length(ind_regressor_type_base{1,3}.ind_reg),length(col_on));
                
                
                for ll2=1:1:length(col_on)
                    
                    if(ll2<length(col_on)) % usual case
                        pt_on=mat_work(4,col_on(ll2));
                        col_off=col_on(ll2)+1;
                        %                 col_off=col_on(ll2)-1+find(mat_work(7,[col_on(ll2):1:(col_on(ll2)+2)])==0.5); % find the next fixation mark presentation
                        pt_off=mat_work(4,col_off);
                        RT=pt_off-pt_on;
                    else % last event in the session is the outcome presentation
                        RT=2.0;
                    end
                    RT_mat=[RT_mat RT];
                    onset_t=mat_work(4,col_on(ll2));
                    onset_mat=[onset_mat onset_t];
                    
                    % fill out regresssor values
                    if(use_model_regressor_cue==1)
                        
                        % regressor type3 only
                        for nn=1:1:length(ind_regressor_type_base{1,3}.ind_reg)
                            mysession=kk2;
                            myblock=mat_work(1,col_on(ll2)); % block in session
                            mytrial=mat_work(2,col_on(ll2)); % trial in block
                            mytrial_s=mat_work(3,col_on(ll2)); % trial_s in trial (arbitration index: 1 at the second stage, 2 at the third stage)
                            mat_work_reg=SBJ{1,jj2}.regressor{1,ind_regressor_type_base{1,3}.ind_reg(nn)}.value(1:4,:);
                            identity_tmp=sum(abs(mat_work_reg-repmat([mysession myblock mytrial mytrial_s]',1,size(mat_work_reg,2))));
                            col_event=find(identity_tmp==0);
                            param_mat3{1,1}(nn,ll2)=SBJ{1,jj2}.regressor{1,ind_regressor_type_base{1,3}.ind_reg(nn)}.value(row_mat(ind_regressor_type_base{1,3}.ind_reg(nn)),col_event);
                        end
                        
                    end
                    
                end
                
                
                % 3rd regressors
                ind_reg=ind_reg+1;
                onsets{1,ind_reg}=onset_mat;
                if(reg_type_go_first(3)==1.7)            names{1,ind_reg}=['Decision_{0T}'];  durations{1,ind_reg}=0;        end
                
                ind_reg_abs=ind_regressor_type_base{1,2}.abs_pos_in_design_mat(end)+1;      list_name_for_contrast{1,ind_reg_abs}=names{1,ind_reg}; % add to the global list of regresssors
                
                
                % (2) pmod:
                if(use_model_regressor_cue==1)
                    ind_reg_param=0;
                    for nn=1:1:length(ind_regressor_type_base{1,3}.ind_reg)
                        ind_reg_param=ind_reg_param+1;
                        pmod(1,ind_reg).name{1,ind_reg_param}=[SBJ{1,jj2}.regressor{1,ind_regressor_type_base{1,3}.ind_reg(ind_reg_param)}.name];
                        pmod(1,ind_reg).poly{1,ind_reg_param}=1;
                        pmod(1,ind_reg).param{1,ind_reg_param}=param_mat3{1,1}(ind_reg_param,:);
                        ind_reg_abs=ind_regressor_type_base{1,3}.abs_pos_in_design_mat(ind_reg_param);      list_name_for_contrast{1,ind_reg_abs}=pmod(1,ind_reg).name{1,ind_reg_param};
                    end
                    % Reserve spaces for extra regressors (will come from the dummy design matrix)
                    for nn=1:1:length(ind_regressor_type_dummy.ind_reg)
                        ind_reg_param=ind_reg_param+1;
                        pmod(1,ind_reg).name{1,ind_reg_param}=ind_regressor_type_dummy.name{1,nn};
                        pmod(1,ind_reg).poly{1,ind_reg_param}=1;
                        pmod(1,ind_reg).param{1,ind_reg_param}=rand(1,length(param_mat3{1,1}(1,:)));
%                         ind_reg_abs=ind_regressor_type_dummy.abs_pos_in_design_mat(nn);      list_name_for_contrast{1,ind_reg_abs}=pmod(1,ind_reg).name{1,ind_reg_param};
                    end
                end
                
            end
            
            %%
            
            
            % (3) Saving normal regressor file
            tot_num_myregressor=length(list_name_for_contrast);
            save_file_name=['Regressor--' SBJ{1,jj2}.name '_sess' sprintf('%02d.mat',kk2)];
            if(Is_save_files_local==1)
                eval(['save ' save_path_result save_file_name ' durations names onsets pmod'])
                eval(['save ' save_for_SPM save_file_name ' durations names onsets pmod'])
            end
            if(Is_save_files_cluster==1)
                eval(['save ' save_path_neuroecon save_file_name ' durations names onsets pmod'])
            end
            
            
            
            %% 2. Extra Regressors independent of cue presentation (timing: TR)
            % The regressors will be saved in separate .mat file.
            
            % regressor generation for each main session and save it to a single file that is compatible with SPM
            ind_reg=0;  % corresponds to the size of the structure
            ind_reg_abs=0; % actual number of regressors (including parametric)
            durations={};
            onsets={};
            names={};
            pmod=struct('name',{},'param',{},'poly',{});
            
            
            ind_reg=ind_reg+1;
            TR_CBIC=2.78; % (sec)
            % (1) duration,onset,name
            durations{1,ind_reg}=TR_CBIC;
            onsets{1,ind_reg}=[0:TR_CBIC:(mat_work(4,end)+20)];
            names{1,ind_reg}=['Dummy'];
            % the name will NOT be added to "list_name_for_contrast" because it is a dummy regressor.
            length_reg=length(onsets{1,1});
            % (2-1) determine the regressor values
            regressor_extra=zeros(length(ind_regressor_type_dummy.ind_reg),length_reg);
            for nn=1:1:length(ind_regressor_type_dummy.ind_reg) % for each regressor
                % find the previous onset point closest to the current time
                t_scan=0;
                for i_scan=1:1:length_reg
                    % find the recent onset (event) time on which the regressor value has been updated.
                    test_mat=onset_event{1,1}-t_scan;
                    test_ii=find(test_mat<0);
                    if(length(test_ii)==0) % first a few scans before the event
                        myscan_time=onset_event{1,1}(1);
                    else % use the value of previous event
                        myscan_time=onset_event{1,1}(min(length(onset_event{1,1}),test_ii(end)));
                    end
                    % find the corresponding {block#, trial#, and trial_s#} in event matrix
                    mycol=find(abs(mat_work(4,:)-myscan_time)<0.001);
                    mysession=kk2; % session
                    myblock=mat_work(1,mycol); % block in session
                    mytrial=mat_work(2,mycol); % trial in block
                    mytrial_s=mat_work(3,mycol); % trial_s in trial (mytrial_s=max(1,mat_work(3,mycol)-1);)
                    % find the corresponding regressor value in a regressor matrix
                    if(mytrial_s~=1)
                        mat_work_reg=SBJ{1,jj2}.regressor{1,ind_regressor_type_dummy.ind_reg(nn)}.value(1:4,:);
                        identity_tmp=sum(abs(mat_work_reg-repmat([mysession myblock mytrial mytrial_s]',1,size(mat_work_reg,2))));
                        col_event=find(identity_tmp==0);
                        regressor_extra(nn,i_scan)=SBJ{1,jj2}.regressor{1,ind_regressor_type_dummy.ind_reg(nn)}.value(row_mat(ind_regressor_type_dummy.ind_reg(nn)),col_event);
                    else % mytrial_s=1: the first state
                        if(i_scan==1) % for the very first scan, we simply read the first regressor value.
                            regressor_extra(nn,i_scan)=SBJ{1,jj2}.regressor{1,ind_regressor_type_dummy.ind_reg(nn)}.value(row_mat(ind_regressor_type_dummy.ind_reg(nn)),1);
                        else % for the first state in every trial, we simply take the regressor value at t-1 (because there is no update)
                            regressor_extra(nn,i_scan)=regressor_extra(nn,i_scan-1);
                        end
                    end
                    % compute the next scan time
                    t_scan=t_scan+TR_CBIC;
                end
            end
            % (2-2) pmod: parametric modulators
            for nn=1:1:length(ind_regressor_type_dummy.ind_reg)
                pmod(1,ind_reg).name{1,nn}=[SBJ{1,jj2}.regressor{1,ind_regressor_type_dummy.ind_reg(nn)}.name];
                pmod(1,ind_reg).poly{1,nn}=1;
                pmod(1,ind_reg).param{1,nn}=regressor_extra(nn,:);
                ind_reg_abs=ind_regressor_type_dummy.abs_pos_in_design_mat(nn);      
                list_name_for_contrast{1,ind_reg_abs}=pmod(1,ind_reg).name{1,nn}; % we do not add because it has been already added in the main design matrix
            end
            
            % (3) Saving dummy regressor file
            tot_num_myregressor=length(list_name_for_contrast);
            save_file_name=['Regressor_dummy--' SBJ{1,jj2}.name '_sess' sprintf('%02d.mat',kk2)];
            if(Is_save_files_local==1)
                eval(['save ' save_path_result save_file_name ' durations names onsets pmod ind_regressor_type_dummy'])
                eval(['save ' save_for_SPM save_file_name ' durations names onsets pmod ind_regressor_type_dummy'])
            end
            if(Is_save_files_cluster==1)
                eval(['save ' save_path_neuroecon save_file_name ' durations names onsets pmod ind_regressor_type_dummy'])
            end
            
        end
    end
    
    
    %% Saving Contrast file
    % [index of my regressors for contrast vector] : total main regressor=6, total regressors=7
    clear contrast_spm
    
    total_number_regressor=tot_num_myregressor+6; % # + 6 movements
    ind_contrast_vec=0;
    
    
    % individual : (ex) [0 1 0 0 0 0 0 0]
    for ii=1:1:tot_num_myregressor
        ind_contrast_vec=ind_contrast_vec+1;
        contrast=zeros(1,tot_num_myregressor);
        contrast(1,ii)=1;
        contrast_spm{1,ind_contrast_vec}.name=list_name_for_contrast{1,ii};
        contrast_spm{1,ind_contrast_vec}.vec=contrast;
    end
    
    % common Q : (ex) [0 1 0 1 0 0 0]
    
%     ind_lli3=[];    ind_lli2=[];
%     for ii=1:1:size(list_name_for_contrast,2) % collect regressor information
%         if(strcmp(list_name_for_contrast{1,ii},'Qsarsa')==1)
%             ind_lli3=[ind_lli3 ii];
%             ind_lli2=[ind_lli2 ii];
%         end
%         if(strcmp(list_name_for_contrast{1,ii},'Qfwd')==1)
%             ind_lli3=[ind_lli3 ii];
%             ind_lli2=[ind_lli2 ii];
%         end
%         if(strcmp(list_name_for_contrast{1,ii},'Qarb')==1)
%             ind_lli3=[ind_lli3 ii];
%         end
%     end
%     
%     ind_contrast_vec=ind_contrast_vec+1;
%     contrast=zeros(1,tot_num_myregressor);
%     contrast(1,ind_lli3)=1;
%     contrast_spm{1,ind_contrast_vec}.name='Qcommon3';
%     contrast_spm{1,ind_contrast_vec}.vec=contrast;
%     
%     ind_contrast_vec=ind_contrast_vec+1;
%     contrast=zeros(1,tot_num_myregressor);
%     contrast(1,ind_lli2)=1;
%     contrast_spm{1,ind_contrast_vec}.name='Qcommon2';
%     contrast_spm{1,ind_contrast_vec}.vec=contrast;
    
    % % difference : (ex) [0 0 0 0 1 -1 0 0]
    % combination_mat=combnk(param_regressor_type_cue_abs_pos_in_design_mat,2);
    % for kk=1:1:size(combination_mat,1)
    %     % A-B
    %     ind_contrast_vec=ind_contrast_vec+1;
    %     contrast=zeros(1,total_number_regressor);
    %     contrast(1,combination_mat(kk,1))=1;    contrast(1,combination_mat(kk,2))=-1;
    %     contrast_spm{1,ind_contrast_vec}.name=[list_name_for_contrast{1,combination_mat(kk,1)} '>' list_name_for_contrast{1,combination_mat(kk,2)}];
    %     contrast_spm{1,ind_contrast_vec}.vec=contrast;
    %     % B-A
    %     ind_contrast_vec=ind_contrast_vec+1;
    %     contrast=zeros(1,total_number_regressor);
    %     contrast(1,combination_mat(kk,2))=1;    contrast(1,combination_mat(kk,1))=-1;
    %     contrast_spm{1,ind_contrast_vec}.name=[list_name_for_contrast{1,combination_mat(kk,1)} '<' list_name_for_contrast{1,combination_mat(kk,2)}];
    %     contrast_spm{1,ind_contrast_vec}.vec=contrast;
    % end
    
    if(Is_save_files_local==1)
        eval(['save ' save_path_result 'contrast_spm.mat' ' contrast_spm ind_regressor_type_base ind_regressor_type_dummy param_regressor_type_cue_abs_pos_in_design_mat list_name_for_contrast'])
        eval(['save ' save_for_SPM 'contrast_spm.mat' ' contrast_spm ind_regressor_type_base ind_regressor_type_dummy param_regressor_type_cue_abs_pos_in_design_mat list_name_for_contrast'])
    end
    if(Is_save_files_cluster==1)
        eval(['save ' save_path_neuroecon 'contrast_spm.mat' ' contrast_spm ind_regressor_type_base ind_regressor_type_dummy param_regressor_type_cue_abs_pos_in_design_mat list_name_for_contrast'])
    end
    
end








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




















