function [output_info]=SIMUL_arbitration_fmri2_tmp(EXP_NAME_IN,index_num,session_opt)
% [NOTE] the goal is changing in each trial.
% (ex) SIMUL_arbitration_fmri('john',1,'pre')
% (ex) SIMUL_arbitration_fmri('john',1,'fmri')

% EXP_NAME_IN='testtest';
% index_num=1; % : session#
% session_opt='fmri'; %'pre','fmri'

% scanner key device : USB HHSC-2x4-C
% use "key_test.m" to see the actual number of keyboard input.


%% organizing all the cue presentation files before the practive session
% refer to ../result_save/List_subject_info.xlsx and update!
if(strcmp(session_opt,'pre')==1)
    if(index_num==1) % for the first session per each subject
        case_num=input('Enter the case number for this subject [1-6] (refer to ../result_save/List_subject_info.xlsx):');
        disp('- now starting to organize the cue presentation files...');
        SIMUL_arbitration_fmri2_init(case_num);
        disp('- organization of the cue presentation files completed.');
    end
end

%% function option check
okay_to_start=0;
if(strcmp(session_opt,'pre')==1)
    okay_to_start=1;
end
if(strcmp(session_opt,'fmri')==1)
    okay_to_start=1;
end
if(okay_to_start~=1)
    error('- ERROR:: check the "session_opt!"');
end

EXP_NAME=[EXP_NAME_IN '_' session_opt];

% seed_path0=['F:\0-Program\MATLABwork\work\One shot learning']; % desktop
seed_path0=pwd;%['D:\0-program\One shot learning']; % laptop
seed_path=[seed_path0 '\seed\'];
% save_path=[seed_path0 '\results_save\'];

%% check if session file exists, and other conditions
COND_NEW_FILE=1;
file_chk_name=[EXP_NAME sprintf('_%d.mat',index_num)];
file_name_check=[pwd '\result_save\' file_chk_name];
if(exist(file_name_check)==2)
    disp('$$$$$ ERROR: the corresponding file exists. try another name or session number. $$$$$$$$$$');
    COND_NEW_FILE=0;
    return;
end
if(index_num>1)
    file_chk_name2=[EXP_NAME sprintf('_%d.mat',index_num-1)];
    file_name_check2=[pwd '\result_save\' file_chk_name2];
    if(exist(file_name_check2)==0) % if the previous session file does not exist
        disp('$$$$$ ERROR: The previous session file does not exist. Check the previous session number. $$$$$$$$$$');
        COND_NEW_FILE=0;
        return;
    end
end
MAX_SESSION_NUM=5;
if(index_num>MAX_SESSION_NUM) % MAX session number =5 !
    disp(sprintf('$$$$$ ERROR: exceeds the max # of sessions =%d. $$$$$$$$$$',MAX_SESSION_NUM));
    COND_NEW_FILE=0;
    return;
end

warning('off')
close all
output_info=0;



%% options
KEEP_PICTURE_OPT=1; % always 1
DO_TAKE_SNAPSHOT=0;
SCREEN_RESOLUTION=[1024 768];
IMAGE_SIZE=[256 256]; %width,height - for cue presentation
OUTCOME_MSG_SIZE=[800 185]; %width,height - for reward message presentation
GOAL_IMG_SIZE=[256 245]; % width,height - for goal state presentation
disp_scale=1.0; % for cue presentation
disp_scale_goalimg=0.5; % for goal image presentation
disp_scale_outcome_msg=0.7; % for outcome msg presentation
disp_scale_scoring=0.7; % for mouse-clicking score submission display
Tot_session=1; % # of total sessions (fixed because this runs for each session)

if(strcmp(session_opt,'pre')==1) % pre-session
    Tot_block=(16+4); % # of total blocks (MUST be the multitude of 4(=#of conditions)
end
if(strcmp(session_opt,'fmri')==1) % fmri-session
    Tot_block=16; % 16; % # of total blocks (MUST be the multitude of 4(=#of conditions)
end

% session schedule
range_num_trials_G0_G1_H0_H1=[[3,5];[5,7];[5,7];[3,5]]; % each row: minmax # of trials of G,G',H,H'
time_estimation_trial_sec=13; %(sec)- used to estimate session time when scheduling
time_limit_session_min=17.8; %(min) - rescheduling until the time estimation meets thie criterion (limit)

if(strcmp(session_opt,'pre')==1) % pre-session
    ffw_speed=4;%4; % fast forward speed, 1: normal speed
end
if(strcmp(session_opt,'fmri')==1) % fmri-session
    ffw_speed=1;%1; % fast forward speed, 1: normal speed
end
% sec_stim_ready=.1; %(sec)
% sec_trial_ready=1.0/ffw_speed; %(sec)
sec_scanner_ready=5/ffw_speed; % sec for scanner stabilization
% sec_block_ready=0.5/ffw_speed; % sec for block ready signal
% sec_stim_display=0.0/ffw_speed;
sec_stim_interval=[1 4]/(ffw_speed);%1.5; %(sec)
sec_trial_interval=[1 4]/(ffw_speed);%1.5; %(sec)
sec_limit_decision=4;%/(ffw_speed); % time limit of decision (sec)
sec_jittered_blank_page=0.15/(ffw_speed); % (sec)
sec_reward_display=2/ffw_speed; % only for 'fmri' option. 1.5sec for 'pre' session

% text size
text_size_default=20; % font size (don't change)
text_size_reward=800; % height in pixel

% background color
BackgroundColor_block_intro=[130,130,130,150]; % gray
BackgroundColor_Cue_page=[210,210,210,150]; % light gray
BackgroundColor_Trial_ready_page=[210,210,210,150]; % light gray
BackgroundColor_Reward_page=[210,210,210,150]; % light gray
COLOR_FIXATION_MARK=[70,70,70,200]; % dark gray

%% key code

% fMRI scanner setting
KEY_L=49;
KEY_R=50;
KEY_Y=51;
KEY_N=52;
KEY_Q=81;
KEY_T=53;

% % laptop setting
% KEY_L=37;
% KEY_R=39;
% KEY_Y=89; %'y'
% KEY_N=78; %'n'
% KEY_Q=81; %'q'
% KEY_T=84; % 't', 5 in desktop, 84 in laptop



%% creating the mother of map and state
map_opt.transition_prob_seed=[0.9 0.1];
map_opt.reward_seed=[40 20 10 0];
[myMap N_state N_action N_transition]=Model_Map_Init2('sangwan2012b',map_opt);
% create my state
myState=Model_RL_Init(N_state,N_action,N_transition);


%% global variables (for each block)
HIST_event_info0=[];
HIST_behavior_info0=[];
if(index_num==1) % create the image usage matrix if this is the first session

    % event   : HIST_event_info{1,session#}
    HIST_event_info_Tag{1,1}='row1 - block#';    HIST_event_info_Tag{2,1}='row2 - trial# (in each block), 0 if outside of the trial';     HIST_event_info_Tag{3,1}='row3 - trial_s#, 0 if outside of the trial_s';
    HIST_event_info_Tag{4,1}='row4 - event time in session';   HIST_event_info_Tag{5,1}='row5 - event time in block';       HIST_event_info_Tag{6,1}='row6 - event time in trial';
    HIST_event_info_Tag{7,1}='row7 - state. 0.5: fixation mark on, 1: S1, 2: S2, 3: S3, 4: S4, 5: S5, 6(+/-)0.1: O1(with win/lost msg), 7(+/-)0.1: O2(with win/lost msg), 8(+/-)0.1: O3(with win/lost msg), 9: O4, 10:A1, 11:A2, 20: a short blank page display, -99:fail to choose in time limit, (-) when display off';
    HIST_event_info_Tag{8,1}='row8 - goal state=outcome state. 6: 40reward state, 7: 20reward state, 8: 10reward state, -1: universal state';

    % behavior : HIST_behavior_info{1,session#}    
    HIST_behavior_info_Tag{1,1}='col1 - block #';    HIST_behavior_info_Tag{1,2}='col2 - trial # (in each block)';
    HIST_behavior_info_Tag{1,3}='col3 - block condition - 1: G(with low uncertainty), 2: G''(with high uncertainty), 3:H(with high uncertainty), 4:H''(with low uncertainty)';
    HIST_behavior_info_Tag{1,4}='col4 - S1';
    HIST_behavior_info_Tag{1,5}='col5 - S2';        
    HIST_behavior_info_Tag{1,6}='col6 - S3';
    HIST_behavior_info_Tag{1,7}='col7 - A1 (action in S1)';     
    HIST_behavior_info_Tag{1,8}='col8 - A2 (action in S2)'; 
    HIST_behavior_info_Tag{1,9}='col9 - RT(A1)';        
    HIST_behavior_info_Tag{1,10}='col10 - RT(A2)';
    HIST_behavior_info_Tag{1,11}='col11 - onset (S1) from the trial start';      
    HIST_behavior_info_Tag{1,12}='col12 - onset (S2) from the trial start';      
    HIST_behavior_info_Tag{1,13}='col13 - onset (S3) from the trial start';
    HIST_behavior_info_Tag{1,14}='col14 - onset (A1) from the trial start';      
    HIST_behavior_info_Tag{1,15}='col15 - onset (A2) from the trial start';
    HIST_behavior_info_Tag{1,16}='col16 - reward amount (0/10/20/40) at S3';
    HIST_behavior_info_Tag{1,17}='col17 - total amount (total in the current session)';
    HIST_behavior_info_Tag{1,18}='col18 - goal state=outcome state. 6: 40reward state, 7: 20reward state, 8: 10reward state, -1: universal state';
    
else % load image usage matrix to update
    file_imgind_ld_name=[EXP_NAME '_info.mat'];
    file_name_ld=[pwd '\result_save\' file_imgind_ld_name];
    load(file_name_ld);
end

% seed image read
img_set=cell(1,N_state); % including othe output states
for img_ind=1:1:N_state
    file_full_path=[seed_path sprintf('s%03d.png',img_ind)];
    img_set{1,img_ind}=imresize(imread(file_full_path),[IMAGE_SIZE(2) IMAGE_SIZE(1)]); % get a frame
end
% reward message image read
img_set_reward_msg_win=cell(1,4); % msg for 40, 20, 10, 0
img_set_reward_msg_lost=cell(1,3); % msg for 40, 20, 10
for jj=1:1:4
    file_full_path=[seed_path sprintf('r%03d_win.png',jj)];
    img_set_reward_msg_win{1,jj}=imresize(imread(file_full_path),[text_size_reward NaN]);
end
for jj=1:1:3
    file_full_path=[seed_path sprintf('r%03d_lost.png',jj)];
    img_set_reward_msg_lost{1,jj}=imresize(imread(file_full_path),[text_size_reward NaN]);
end
% goal state image read
img_set_goal=cell(1,4); % 40, 20 ,10, universal
for jj=1:1:4
    file_full_path=[seed_path sprintf('g%03d.png',jj)];
    img_set_goal{1,jj}=imresize(imread(file_full_path),[GOAL_IMG_SIZE(2) GOAL_IMG_SIZE(1)]);
end


%% Scheduling (block sequence)
% G:1, G':2, H:3, H':4 
% # of trials : G,H'~[3,5], G',H~[5,7]
criterion=0;
HIST_block_condition_Tag{1,1}='row1: block #';
HIST_block_condition_Tag{2,1}='row2: block condition of each trial: 1:G, 2:G'', 3:H, 4:H''';
% HIST_block_condition_Tag{3,1}='row3: devaluation index: 1:normal trial, 0:devalued trial';
if(strcmp(session_opt,'pre')==1) % pre-session (5trials for each block)
    tmp_block_ind_mat=[1:1:Tot_block]'*ones(1,5);
    HIST_block_trial_index=[reshape(tmp_block_ind_mat',[1 5*Tot_block]); [3*ones(1,5*(Tot_block-2)) 2*ones(1,5*2)]];
%     HIST_block_trial_index=[HIST_block_trial_index; ones(1,size(HIST_block_trial_index,2))];
end
if(strcmp(session_opt,'fmri')==1) % fmri-session
    disp('- scheduling start...');
    while(criterion==0)
        HIST_block_trial_index=[]; % (2xtrial#) row1: block#, row2: the block condition of each trial
        block_index=0;
        for block_quad=1:4:Tot_block
            current_block_seq=randperm(size(range_num_trials_G0_G1_H0_H1,1)); % G:1, G':2, H:3, H':4
            if(block_quad>1)
                while(current_block_seq(1)==HIST_block_trial_index(end)) % do not allow the same subsequent condition.
                    current_block_seq=randperm(4); % G:1, G':2, H:3, H':4
                end
            end
            if((block_quad==1)&&(index_num==1)) % the first block of the first session should be G
                current_block_seq=[3 3 1 2];
            end

            for jj=1:1:length(current_block_seq)
                
                block_index=block_index+1;
                tmp=randperm(max(range_num_trials_G0_G1_H0_H1(jj,:))-min(range_num_trials_G0_G1_H0_H1(jj,:))+1); % pick integer
                num_trial=tmp(1)+min(range_num_trials_G0_G1_H0_H1(jj,:))-1;

                %                 condi_habit=(current_block_seq(jj)==3)||(current_block_seq(jj)==4);
                %                 if((rand>0.52)&&(condi_habit)) % add the devaluation index matrix (in the third row)
                %                     deval_mat0=[ones(1,num_trial-3) zeros(1,3)]; % devaluation for the last two trials
                %                 else
                %                     deval_mat0=[ones(1,num_trial)];
                %                 end
                %                 if((block_quad==1)&&(index_num==1)) % no devaluation for the first block of the first session
                %                     deval_mat0=[ones(1,num_trial)];
                %                 end
                %                 HIST_block_trial_index=[HIST_block_trial_index [ones(1,num_trial)*block_index; ones(1,num_trial)*current_block_seq(jj); deval_mat0]];
                HIST_block_trial_index=[HIST_block_trial_index [ones(1,num_trial)*block_index; ones(1,num_trial)*current_block_seq(jj)]];
                
            end
            
        end
        time_est=size(HIST_block_trial_index,2)*time_estimation_trial_sec/60; %min
        criterion=(time_est<time_limit_session_min);
        disp('- rescheduling to meet the session time limit criterion...')
    end
    str=sprintf('- scheduling done. # of trials = %d. Estimated session time = %02.1f. Proceed? (''n'' to quit, proceed otherwise)',size(HIST_block_trial_index,2),time_est); disp(str);
    [secs, keyCode] = KbPressWait;      [tmp tmp_key_code]=find(keyCode==1);
    if(tmp_key_code==KEY_N) % n pressed
        disp('### Experiment aborted as per user''s request. ###');    return;
    end
end
HIST_block_condition{1,index_num}=HIST_block_trial_index;
disp('- proceed to the experiment...');




%% Display initialization
whichScreen = 0;
wPtr  = Screen('OpenWindow',whichScreen);
[screenWidth, screenHeight] = Screen('WindowSize', wPtr);

white = WhiteIndex(wPtr); % pixel value for white
black = BlackIndex(wPtr); % pixel value for black
gray = (white+black)/2;
inc = white-gray;
inc_0=white-black;

imageArray={};

%% starting message
Screen('TextSize',wPtr, text_size_default);
% Screen('TextFont',wPtr, 'Times New Roman');
DrawFormattedText(wPtr, 'Are you ready for the experiment?\n(Press any key to wait for the trigger)', 'center', 'center');
Screen('Flip', wPtr);  
KbWait; % temporarily disabled for test APR 21
if(DO_TAKE_SNAPSHOT==1)
    snapshot=Screen(wPtr, 'GetImage', [1, 1, floor(1.0*SCREEN_RESOLUTION(1)), floor(1.0*SCREEN_RESOLUTION(2))]);
    imageArray=[imageArray; {snapshot}];
end


%% waiting for the trigger sign from the scanner
DrawFormattedText(wPtr, 'Now waiting for the trigger...', 'center', 'center');
Screen('Flip',wPtr);
% Look for trigger pulse
while 1 % temporarily disabled for test APR 21
    [ keyIsDown, timeSecs, keyCode ] = KbCheck;
    if keyIsDown
        [tmp tmp_key_code]=find(keyCode==1);

        %         if ((tmp_key_code==53)) % trigger
        if ((tmp_key_code==KEY_T)) % trigger
%         if (KbName(keyCode) == 5)
            break;
        end
        % If the user holds down a key, KbCheck will report multiple events.
        % To condense multiple 'keyDown' events into a single event, we wait until all
        % keys have been released.
        while KbCheck; end
    end
end


%% clock-start and then wait for another 5secs until the scanner stabilizes
session_clock_start = GetSecs;
WaitSecs(sec_scanner_ready);


%% block
for block=1:1:Tot_block % each block

    %% I. Stimuli presentation

    zzz=[];
    % block starts
    %     Screen('FillRect',wPtr,BackgroundColor_block_intro);
    %     str_block_intro=sprintf('*** Now block %d starts. ***',block);
    %     DrawFormattedText(wPtr, str_block_intro, 'center', 'center',[0, 0, 0, 255]);
    %     Screen(wPtr, 'Flip');
        block_clock_start = GetSecs;
    %     HIST_event_info0=[HIST_event_info0 [block; 0; 0; (GetSecs-session_clock_start); (GetSecs-block_clock_start); 0.5]]; % event save
    %     WaitSecs(sec_block_ready);

    trial_in_block=0;
    [tmp trial_set]=find(HIST_block_trial_index(1,:)==block);
    for trial=trial_set % each trial

        trial_in_block=trial_in_block+1;
        trial_clock_start = GetSecs;
                
        onset_state_from_trial_start=[];
        onset_action_from_trial_start=[];
        
        state_sbj=myState;
        map_sbj=myMap;

        
        %% [1] map preparation (G:1, G':2, H:3, H':4)
        block_condition=HIST_block_trial_index(2,trial);

        
        block_condition=3; 
        
        
        if(block_condition==1) % G
            % set T_prob
            prob_seed_mat=[0.9 0.1];
            % set Reward            
            tmp_mat=randperm(3);
            current_goal_state=map_sbj.goal_state_index(tmp_mat(1));
            map_sbj.reward=zeros(N_state,1);
            map_sbj.reward(current_goal_state)=map_sbj.reward_save(current_goal_state);
        end
        if(block_condition==2) % G'
            % set T_prob
            prob_seed_mat=[0.5 0.5];
            % set Reward            
            tmp_mat=randperm(3);
            current_goal_state=map_sbj.goal_state_index(tmp_mat(1));
            map_sbj.reward=zeros(N_state,1);
            map_sbj.reward(current_goal_state)=map_sbj.reward_save(current_goal_state);
        end
        if(block_condition==3) % H
            % set T_prob
            prob_seed_mat=[0.5 0.5];
            % set Reward
            map_sbj.reward=map_sbj.reward_save; % all the rwd state intact
            current_goal_state=-1;
        end
        if(block_condition==4) % H'
            % set T_prob
            prob_seed_mat=[0.9 0.1];
            % set Reward
            map_sbj.reward=map_sbj.reward_save; % all the rwd state intact
            current_goal_state=-1;
        end
        % T_prob encoding to the current map
        map_sbj.action(1,1).prob(1,[2 3])=prob_seed_mat;
        map_sbj.action(1,1).prob(2,[7 9])=prob_seed_mat;
        map_sbj.action(1,1).prob(3,[8 9])=prob_seed_mat;
        map_sbj.action(1,1).prob(4,[7 9])=prob_seed_mat;
        map_sbj.action(1,1).prob(5,[6 9])=prob_seed_mat;
        map_sbj.action(1,2).prob(1,[4 5])=prob_seed_mat;
        map_sbj.action(1,2).prob(2,[8 7])=prob_seed_mat;
        map_sbj.action(1,2).prob(3,[7 9])=prob_seed_mat;
        map_sbj.action(1,2).prob(4,[6 7])=prob_seed_mat;
        map_sbj.action(1,2).prob(5,[7 9])=prob_seed_mat;

        

        
        
        
        save map_tmp map_sbj
        qwe;
        
        
        
        
        
        
        
        while(map_sbj.IsTerminal(state_sbj.state_history(state_sbj.index))==0)

            
            current_state=state_sbj.state_history(state_sbj.index);


            %% [2] display
            
            % (1) display fixation mark and display off during the jittered interval
            Screen('FillRect',wPtr,BackgroundColor_Cue_page);
            DrawFormattedText(wPtr, '+', 'center', 'center', COLOR_FIXATION_MARK); % add 'o' mark at the click pt.
            Screen(wPtr, 'Flip');
            HIST_event_info0=[HIST_event_info0 [block; trial_in_block; state_sbj.index; ...
                (GetSecs-session_clock_start); (GetSecs-block_clock_start); (GetSecs-trial_clock_start); 0.5; current_goal_state]]; % event save
            sec_stim_interval0=rand*(max(sec_stim_interval)-min(sec_stim_interval))+min(sec_stim_interval);
            WaitSecs(sec_stim_interval0);
            
            % (2-1) add state image
            Screen('FillRect',wPtr,BackgroundColor_Cue_page);
            input_stim = Screen('MakeTexture', wPtr, img_set{1,current_state});
            xpos = round(screenWidth/2);    ypos = round(screenHeight/2);
            sx=floor(IMAGE_SIZE(1)*disp_scale);       sy=floor(IMAGE_SIZE(2)*disp_scale);
            destrect=[xpos-sx/2,ypos-sy/2,xpos+sx/2,ypos+sy/2];            
            Screen('DrawTexture', wPtr, input_stim,[],destrect);

            % (2-2) add goal image            
            if(current_goal_state==-1) % universal state
                current_goal_state_ind=4;
            else % other goal states
                current_goal_state_ind=current_goal_state-5;
            end
            input_stim2 = Screen('MakeTexture', wPtr, img_set_goal{1,current_goal_state_ind});
            sx2=floor(GOAL_IMG_SIZE(1)*disp_scale_goalimg);       sy2=floor(IMAGE_SIZE(2)*disp_scale_goalimg);
            xpos2=xpos; ypos2=ypos+sy/2+sy2/2+50;
            destrect2=[xpos2-sx2/2,ypos2-sy2/2,xpos2+sx2/2,ypos2+sy2/2];            
            Screen('DrawTexture', wPtr, input_stim2,[],destrect2);
            
            % for TEST only - APR 27, 2012
%             if(block_condition==1)                txt_test='G  (0.9, 0.1)';            end
%             if(block_condition==2)                txt_test='G'' (0.5, 0.5)';            end
%             if(block_condition==3)                txt_test='H  (0.5, 0.5)';            end
%             if(block_condition==4)                txt_test='H'' (0.9, 0.1)';            end
%             txt_test2=sprintf('current state=%d,',current_state);
%             DrawFormattedText(wPtr, [txt_test, ', ' txt_test2], 100, 200);
            
            % (2-3) display on
            Screen(wPtr, 'Flip'); % display on
            clock_time_limit_start=GetSecs;
            onset_state_from_trial_start=[onset_state_from_trial_start (GetSecs-trial_clock_start)];
            HIST_event_info0=[HIST_event_info0 [block; trial_in_block; state_sbj.index; ...
                (GetSecs-session_clock_start); (GetSecs-block_clock_start); onset_state_from_trial_start(end); current_state; current_goal_state]]; % event save

            %% [3] get chioce and update            
            decision_made=0;
            while(~decision_made)
                [secs, keyCode] = KbPressWait([], clock_time_limit_start+sec_limit_decision); % if no keyboard in time limit, then go ahead. if pressed earlier, then go ahead.
                onset_action_from_trial_start=[onset_action_from_trial_start (GetSecs-trial_clock_start)]; 
                state_sbj.RT(state_sbj.index)=GetSecs-clock_time_limit_start;                
                [tmp tmp_key_code]=find(keyCode==1);                
                if(tmp_key_code==KEY_L) % L pressed
                    state_sbj.action_history(state_sbj.index)=1;
                    decision_made=1;
                end
                if(tmp_key_code==KEY_R) % R pressed
                    state_sbj.action_history(state_sbj.index)=2;
                    decision_made=1;
                end
                if(tmp_key_code==KEY_Q) % 'q' pressed for aborting
                    state_sbj.action_history(state_sbj.index)=ceil(2*rand); % random select if fail to make a decision
                    HIST_event_info0=[HIST_event_info0 [block; trial_in_block; state_sbj.index; ...
                        (GetSecs-session_clock_start); (GetSecs-block_clock_start); onset_action_from_trial_start(end); -99; current_goal_state]]; % event save
                    qwe;
                    break;
                    decision_made=1;                    
                end                
                % check the time limit !@#$                
                if((state_sbj.RT(state_sbj.index)>sec_limit_decision)&&(decision_made==0)) % no decision made in time limit
                    state_sbj.action_history(state_sbj.index)=ceil(2*rand); % random select if failed to make a decision in time
                    decision_made=1;            Is_bet=1;
                    HIST_event_info0=[HIST_event_info0 [block; trial_in_block; state_sbj.index; ...
                        (GetSecs-session_clock_start); (GetSecs-block_clock_start); onset_action_from_trial_start(end); -99; current_goal_state]]; % event save
                else % in time, and decision made
                    if(decision_made==1)
                        HIST_event_info0=[HIST_event_info0 [block; trial_in_block; state_sbj.index; ...
                            (GetSecs-session_clock_start); (GetSecs-block_clock_start); onset_action_from_trial_start(end); 9+state_sbj.action_history(state_sbj.index); current_goal_state]]; % event save
                    end
                end
            end
            
            %% [3] moving to the next state
            [state_sbj map_sbj]=StateSpace(state_sbj,map_sbj);  % map&state index ++
            
            
            

        end % end of each choice
        
        
        %% [4] terminal state: display reward
        current_state=state_sbj.state_history(state_sbj.index);
        
        % (0) display fixation mark and display off during the jittered interval
        Screen('FillRect',wPtr,BackgroundColor_Cue_page);
        DrawFormattedText(wPtr, '+', 'center', 'center', COLOR_FIXATION_MARK); % add 'o' mark at the click pt.
        Screen(wPtr, 'Flip');
        HIST_event_info0=[HIST_event_info0 [block; trial_in_block; state_sbj.index; ...
            (GetSecs-session_clock_start); (GetSecs-block_clock_start); (GetSecs-trial_clock_start); 0.5; current_goal_state]]; % event save
        sec_stim_interval0=rand*(max(sec_stim_interval)-min(sec_stim_interval))+min(sec_stim_interval);
        WaitSecs(sec_stim_interval0);
        
        % (1) add state image
        Screen('FillRect',wPtr,BackgroundColor_Cue_page);
        input_stim = Screen('MakeTexture', wPtr, img_set{1,current_state});
        xpos = round(screenWidth/2);    ypos = round(screenHeight/2);
        sx=floor(IMAGE_SIZE(1)*disp_scale);       sy=floor(IMAGE_SIZE(2)*disp_scale);
        destrect=[xpos-sx/2,ypos-sy/2,xpos+sx/2,ypos+sy/2];
        Screen('DrawTexture', wPtr, input_stim,[],destrect);

        % (2) outcome message read
        outcome_state=state_sbj.state_history(3);
        original_rwd=map_sbj.reward_save(outcome_state);
        actual_rwd=state_sbj.reward_history(3);
        if(original_rwd==actual_rwd) % earn case = 40,20,10,0
            case_earn=1;
            % msg for 40, 20, 10, 0 added to the total
            input_stim2 = Screen('MakeTexture', wPtr, img_set_reward_msg_win{1,outcome_state-5});
        else % devalued case - lost case
            case_earn=-1;
            % msg for 40, 20, 10 did not added to the total.
            input_stim2 = Screen('MakeTexture', wPtr, img_set_reward_msg_lost{1,outcome_state-5});
        end
        % (3) add outcome message
        sx2=floor(OUTCOME_MSG_SIZE(1)*disp_scale_outcome_msg);       sy2=floor(OUTCOME_MSG_SIZE(2)*disp_scale_outcome_msg);
        xpos2=xpos; ypos2=ypos+sy/2+sy2/2+50;
        destrect2=[xpos2-sx2/2,ypos2-sy2/2,xpos2+sx2/2,ypos2+sy2/2];
        Screen('DrawTexture', wPtr, input_stim2,[],destrect2);

        % (4) display on
        Screen(wPtr, 'Flip'); % display on
        if(outcome_state~=9)
            value=outcome_state+0.1*case_earn;
        else
            value=outcome_state;
        end
        onset_state_from_trial_start=[onset_state_from_trial_start (GetSecs-trial_clock_start)];
        HIST_event_info0=[HIST_event_info0 [block; trial_in_block; state_sbj.index; ...
            (GetSecs-session_clock_start); (GetSecs-block_clock_start); onset_state_from_trial_start(end); value; current_goal_state]]; % event save

        if(strcmp(session_opt,'pre')==1) % pre-session
            WaitSecs(1.5);
        end
        if(strcmp(session_opt,'fmri')==1) % fmri-session
            WaitSecs(sec_reward_display);
        end

        %% Update HIST_event_info (overwrite at each block)
        HIST_event_info{1,index_num}=HIST_event_info0;
        
        %% update behavior matrix !@#$%
        if(isempty(HIST_behavior_info0))
            acc_rwd=0;
        else
            acc_rwd=HIST_behavior_info0(end,17);
        end
        mat_update=[block, trial_in_block, block_condition, ...
            state_sbj.state_history(1), state_sbj.state_history(2), state_sbj.state_history(3), ...
            state_sbj.action_history(1), state_sbj.action_history(2), ...
            state_sbj.RT(1), state_sbj.RT(2), ...
            onset_state_from_trial_start(1), onset_state_from_trial_start(2), onset_state_from_trial_start(3), ...
            onset_action_from_trial_start(1), onset_action_from_trial_start(2), ...
            state_sbj.reward_history(end), acc_rwd+state_sbj.reward_history(end), ...
            current_goal_state];
        HIST_behavior_info0=[HIST_behavior_info0; mat_update];
        
    end % end of each trial


    % behavior matrix update
    HIST_behavior_info{1,index_num}=HIST_behavior_info0;

    %% save the (updated) image usage matrix (overwriting)
    file_imgind_sv_name=[EXP_NAME '_info.mat'];
    file_name_sv=[pwd '\result_save\' file_imgind_sv_name];
    save(file_name_sv,'HIST_event_info','HIST_event_info_Tag','HIST_behavior_info','HIST_behavior_info_Tag','HIST_block_condition','HIST_block_condition_Tag');


end % end of each block

        
%% Ending message
str_end=sprintf('- Our experiments is over. Press any key to quit. -');
DrawFormattedText(wPtr, str_end, 'center', 'center');
Screen(wPtr, 'Flip');
KbWait; % temporarily disabled for test APR 21
% take a snapshot
if(DO_TAKE_SNAPSHOT==1)
    snapshot=Screen(wPtr, 'GetImage', [1, 1, floor(1.0*SCREEN_RESOLUTION(1)), floor(1.0*SCREEN_RESOLUTION(2))]);
    imageArray=[imageArray; {snapshot}];
end


%% save snapshots : CAUTION HEAVY PROCESS - might take a minute.
if(DO_TAKE_SNAPSHOT==1)
    for j=1:1:size(imageArray,1)
        str=sprintf('snapshot_dispay_exp_%03d.png',j);
        imwrite(imageArray{j},['snapshot\' str],'png');
    end
end


% behavior matrix update
HIST_behavior_info{1,index_num}=HIST_behavior_info0;


%% save the (updated) image usage matrix (overwriting)
file_imgind_sv_name=[EXP_NAME '_info.mat'];
file_name_sv=[pwd '\result_save\' file_imgind_sv_name];
save(file_name_sv,'HIST_event_info','HIST_event_info_Tag','HIST_behavior_info','HIST_behavior_info_Tag','HIST_block_condition','HIST_block_condition_Tag');

%% save all variables
file_sv_name=[EXP_NAME sprintf('_%d.mat',index_num)];
file_name_sv=[pwd '\result_save\' file_sv_name];
save(file_name_sv,'*');


%% session end sound
sec_dur_sound=1;
Beeper('med', 0.4, sec_dur_sound); WaitSecs(sec_dur_sound)
Beeper('high', 0.4, sec_dur_sound)



%% finish all
Screen('CloseAll');
clear mex
% clear Screen

disp('########################################################')
str_end1=sprintf('### session%d is done ############################',index_num);
disp(str_end1);
str_end2=sprintf('### next session = %d ############################',index_num+1);
disp(str_end2);
disp('########################################################')

% display the number of response failure
missed_count = length(find(HIST_event_info{1,index_num}(7,:)==-99));
disp(sprintf('- # of response failure in this session = %d. (will be penalized) ',missed_count));

output_info=1;
end


