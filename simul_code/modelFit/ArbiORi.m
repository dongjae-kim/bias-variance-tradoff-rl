% Arbitration

% subject list
LIST_SBJ={'subject001_ksy', 'subject03_keb', 'subject004_ksj'};

% reject abnormal(different session length)
list_sbj={LIST_SBJ{1:end}};
for i = 1 : length(list_sbj)
   sbj(i).name = list_sbj{i}; 
end
sbj(1).val = 2.588378999578160e+02;
sbj(1).param = [0.691069294004446	0.113920373959288	1.00186577280393	6.99698277480360	0.199125728205578	0.0899666159119893];
sbj(2).val = 2.536326261750106e+02;
sbj(2).param = [0.699795805527792	0.163006424890624	1.01680487798466	6.68556829786769	0.197619108764548	0.0934585886357233];
sbj(3).val = 2.391211004311237e+02;
sbj(3).param = [0.650309330934117	0.0595422552835219	1.00055028635363	1.01153200341886	0.199899421992754	0.0996512999753423];
% optimized parameters 예시.

param_ori=[0.3 + 0.04*randi(10) ,  0.02 + 0.018*randi(10) , randi(5)  , randi(7) ,  0.05 + 0.015 * randi(100) , 0.03 + 0.012*randi(10)];
param_BoundL = [0.3 0.02 1 1 0.05 0.03];
param_BoundU = [0.7 0.2 5 7 0.2 0.15];

%% mode
mode.param_length = size(param_ori,2);
mode.total_simul = 20;
mode.experience_sbj_events=[1 1]; % ones(1,2); % experience_sbj_events(2) = 1 일때 subject의 state transition을 이용.
mode.USE_FWDSARSA_ONLY=0; % 이게 1이면 fwd only고 이게 2이면 sarsa only. 둘다 아니게 하려고 이런 정보를 넣었지만 아예 그 if문을 둘다 빼버리는 방법도 있음.
mode.USE_BWDupdate_of_FWDmodel=1; % BWD update 사용하는 방식으로 학습하기로.
mode.DEBUG_Q_VALUE_CHG=0;
mode.simul_process_display=0;
mode.out=1;
mode.opt_ArbModel = 0; %
mode.boundary_12 = 0.1; % 위아래거 바뀔수도 있음 사실 뭐가 맞는지 모르겟음...
mode.boundary_21 = 0.01;
mode.max_iter=100;


%% valid subject

LIST_SBJ={'subject001_ksy', 'subject03_keb', 'subject004_ksj'};
% rejecting abnormal(different session length)
list_sbj={LIST_SBJ{1:end}};
%Surendra 4개, 10번째, Gjergji 3 16번째,  Domenick 4 20번째, Thao 4 21번째


%% data_in
% 지금 가지고 있는 가짜 데이터는 data_in 1~9명의 서브젝트에 대한 data_in
maxi=size(list_sbj,2);
batch_Save={};
tempval = {};
tempval2 = {};
tempval3={};
% 2016 11 07 full subject history

maxd=100; % seed parameter 몇개 사용해 볼 것인지.
exptime=0;
timedu=0;
exp1=clock;
exp0=clock;
outputval_oriparam=[];
disp('DATA loading...');

% Data preparation
for i = 1 : maxi
    Data_pre{i}= load([pwd '\result_save_sy\' list_sbj{i} '_pre_1.mat']);
    
    tt = dir([pwd '\result_save_sy']);
    tt = {tt.name};
    maxsess = sum(cell2mat(strfind(tt,[list_sbj{i} '_fmri_']))) - 1;
    
    for ii = 1 : maxsess
        Data_main{i,ii}=load([pwd '\result_save_sy\' list_sbj{i} '_fmri_' num2str(ii) '.mat']);
    end
end

outputval=cell(maxd,1);
outputparam=outputval;

for d= 1 : maxd
    clc;
    param_init=[0.3 + 0.04*randi(10) ,  0.02 + 0.018*randi(10) , randi(5)  , randi(7) ,  0.05 + 0.015 * randi(100) , 0.03 + 0.012*randi(10)];
    for i = 1 : 1 : maxi
        clc;
        timedu = etime (exp1,exp0);
        exp0=clock;
        % pre save
        data_pre=Data_pre{i};
        data_in{1,1}.HIST_behavior_info_pre{1,1}=data_pre.HIST_behavior_info{1,1};
        data_in{1,1}.HIST_block_condition_pre{1,1}=data_pre.HIST_block_condition{1,1};
        
        disp('############################################');
        disp('############################################');
        disp('############################################');
        fprintf('###   SUB_NUM: [%d / %d]\n',i,maxi);
        fprintf('### SIMUL_NUM: [%d / %d]\n',d,maxd);
        fprintf('### OPT_ITER : [%d]\n',mode.max_iter);
        fprintf('%d/%d/%d ###### %d시 %d분\n',t(1),t(2),t(3),t(4),t(5))
        disp('############################################');
        disp('############################################');
        disp('############################################');
        
        % main save
        
        % max sess eval
        
        temp = [];
        tt = dir([pwd '\result_save']);
        tt = {tt.name};
        maxsess = sum(cell2mat(strfind(tt,[list_sbj{i} '_fmri_']))) - 1;
        
        
        
        try
            for ii = 1 : maxsess
                data_main=Data_main{i,ii};
                data_in{1,1}.HIST_behavior_info{1,ii}=data_main.HIST_behavior_info0;
                data_in{1,1}.HIST_block_condition{1,ii}=data_main.HIST_block_condition{1,ii};
                disp(['     => current session num : ' num2str(ii)]);
            end
            data_in{1,1}.map_type=2;
        catch
            warning('MAX session error');
        end
        
        try
            myFunc_bu = @(x) eval_ArbitrationRL3c(x, data_in, mode);
            [model_BayesArb.param, model_BayesArb.val] = fminsearchbnd(myFunc_bu, param_init, param_BoundL, param_BoundU, optimset('Display','iter','MaxIter',mode.max_iter)); % X0,LB,UB
            outputval{d}=[outputval{d} model_BayesArb.val];
            outputparam{d}=[outputparam{d}; model_BayesArb.param];
        catch
            param_init=[0.3 + 0.04*randi(10) ,  0.02 + 0.018*randi(10) , randi(5)  , randi(7) ,  0.05 + 0.015 * randi(100) , 0.03 + 0.012*randi(10)]; % Initialized parameter에 문제가 있어서 그랬을 수 있으므로 이런식으로 다시 initialize
        end
        
        exp1=clock;
    end
    
    
end


disp('TIME_SPENT@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@');
disp(cputime-timeee);
value_param = struct('val', outputval, 'param', outputparam);
save([pwd '\RESULT.mat'], 'value_param');
