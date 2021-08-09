clear;

timeee=cputime;
%% param_in
%param_in 이건 stan에서 optimization 해줄수도 있을거같은데 조금 더 생각하고 알아봐야 할듯.
param_in_tag={};
param_in_tag{1,1} = 'a threshlod for defining zero SPE ; free parameter # 1 in sangwan lee’s 2014 neuron paper';
param_in_tag{1,2} = 'a threshlod for defining zero RPE ; free parameter # 1 in sangwan lee’s 2014 neuron paper';
param_in_tag{1,3} = 'inverse softmax temperature, tau of the softmax ; free parameter # 5 in sangwan lee’s 2014 neuron paper';
param_in_tag{1,4} = 'a leraning rate of the model-based.mode-free, alpha learning rate ; free parameter # 6 in sangwan lee’s 2014 neuron paper';
param_in=rand(1,4);

%% mode 
mode.param_length = size(param_in,2);
mode.total_simul = 1;
mode.experience_sbj_events=ones(1,2);
mode.USE_FWDSARSA_ONLY=99; % 이게 1이면 fwd only고 이게 2이면 sarsa only. 둘다 아니게 하려고 이런 정보를 넣었지만 아예 그 if문을 둘다 빼버리는 방법도 있음.
mode.USE_BWDupdate_of_FWDmodel=1; % BWD update 사용하는 방식으로 학습하기로.
mode.DEBUG_Q_VALUE_CHG=0;
mode.simul_process_display=1;
mode.out=1;
mode.opt_ArbModel = 1; % invFano model.


%% valid subject

LIST_SBJ={'Oliver', 'Hao', 'Breanna', 'Derek', 'Timothy', 'Teagan', 'Jeffrey', 'Seung', 'Carole', 'Tony', 'Surendra', 'Lark',...
    'Joaquin', 'DavidB', 'Christopher', 'Gjergji', 'Charles', 'Erin', 'Connor', 'Domenick', 'Thao', 'Arin', 'Pauline', 'Tho'};
% rejecting abnormal(different session length)
list_sbj={LIST_SBJ{2:10} LIST_SBJ{12:15} LIST_SBJ{17:19} LIST_SBJ{22:24}}; 


%% data_in
% 지금 가지고 있는 가짜 데이터는 data_in 1~9명의 서브젝트에 대한 data_in

maxi=size(list_sbj,2);


for i = 1 : 1 %maxi
    % pre save
   data_pre=load([pwd '\result_save\' list_sbj{i} '_pre_1.mat']); 
   data_in{1,i}.HIST_behavior_info_pre{1,1}=data_pre.HIST_behavior_info{1,1};
   data_in{1,i}.HIST_block_condition_pre{1,1}=data_pre.HIST_block_condition{1,1};   
   disp(['current sub num : ' num2str(i)]);
   % main save
   for ii = 1 : 5
    data_main=load([pwd '\result_save\' list_sbj{i} '_fmri_' num2str(ii) '.mat']); 
    data_in{1,i}.HIST_behavior_info{1,ii}=data_main.HIST_behavior_info{1,1};
    data_in{1,i}.HIST_block_condition{1,ii}=data_main.HIST_block_condition{1,1};
    disp(['     => current session num : ' num2str(ii)]);
   end
   data_in{1,i}.map_type=2;
end

data_in_tag{1,1}='data in의 열들은 subject 넘버. 그리고 그 내부의 HIST behavior info내부는 각각 세션을 의미';




% 
% for i = 1 : 1 : maxi
%     
%     
%     data_in{1,i}.HIST_behavior_info_pre{1,1}=PRE.HIST_behavior_info{1,1};
%     data_in{1,i}.HIST_block_condition_pre{1,1}=PRE.HIST_block_condition{1,1};
%     for ii= 1 : 4
%         data_in{1,i}.HIST_behavior_info{1,ii}=NPRE.HIST_behavior_info{1,1};
%         data_in{1,i}.HIST_block_condition{1,ii}=NPRE.HIST_block_condition{1,1};
%     end
%     data_in{1,i}.map_type=2; % this is for using myMap_new as map structure
% end




%% スタト
% Description : Arbitrator로 implement해줘야 하는 파라미터는 param_in(1:4) 각각의 의미는 앞에서
% 봤고.
outputval=eval_ArbitrationRL5(param_in,data_in,mode);
%outputval = eval_ArbitrationRL3_lessparam(param_in,data_in,mode);

map=outputval.map;
state_sarsa=outputval.sarsa;
state_fwd=outputval.fwd;
myArbitrator_top=outputval.Arbitrator;

disp('@@@@@@@@@@@@@@@@@@@@@@@@@TIME@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@');
disp(cputime-timeee);




save('DPplot', 'outputval');

DPploting;