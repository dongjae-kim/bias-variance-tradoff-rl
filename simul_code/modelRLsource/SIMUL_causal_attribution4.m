function [output_info]=SIMUL_causal_attribution4(EXP_NAME,index_num)
% EXP_NAME='test0'
% index_num=1; % : session#



% seed_path0=['F:\0-Program\MATLABwork\work\One shot learning']; % desktop
seed_path0=pwd;%['D:\0-program\One shot learning']; % laptop
seed_path=[seed_path0 '\seed\'];
save_path=[seed_path0 '\results_save\'];

%% check if session file exists, and other conditions
COND_NEW_FILE=1;
file_chk_name=[EXP_NAME sprintf('_%d.mat',index_num)];
file_name_check=[save_path file_chk_name];
if(exist(file_name_check)==2)
    disp('$$$$$ ERROR: the corresponding file exists. try another name or session number. $$$$$$$$$$');
    COND_NEW_FILE=0;
    return;
end
if(index_num>1)
    file_chk_name2=[EXP_NAME sprintf('_%d.mat',index_num-1)];
    file_name_check2=[save_path file_chk_name2];
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
KEEP_PICTURE_OPT=1; % keepting picture + bonus session
DO_TAKE_SNAPSHOT=0;
SCREEN_RESOLUTION=[1920 1200];
IMAGE_SIZE=[256 256]; %width,height - for cue presentation
disp_scale=1.5; % for cue presentation
disp_scale_scoring=0.8; % for mouse-clicking score submission display
Tot_session=1; % # of total sessions (fixed because this runs for each session)
Tot_block=8; % # of total blocks (must be EVEN number)
Tot_trial=5; % # of trials in each block
Tot_trial_s=5; % # of cue presentation in one trial [NOTE] Tot_trial*Tot_trial_s must be *even* number
tot_num_img_per_block=3;

ffw_speed=1; % fast forward speed, 1: normal speed
% sec_stim_ready=.1; %(sec)
% sec_trial_ready=1.0/ffw_speed; %(sec)
sec_scanner_ready=5/ffw_speed; % sec for scanner stabilization
sec_block_ready=0.5/ffw_speed; % sec for block ready signal
sec_stim_display=1.0/ffw_speed;
sec_stim_interval=[1 4]/(ffw_speed);%1.5; %(sec)
sec_trial_interval=[1 4]/(ffw_speed);%1.5; %(sec)
sec_limit_Q_rating=[4 4]/(ffw_speed); % time limit of [hedonic_rating causal_rating] (sec)
sec_limit_Bayesian_rating=8/(ffw_speed); % time limit of Bayesian triangular rating(sec)
sec_jittered_blank_page=0.15/(ffw_speed); % (sec)

earliest_novel_stimulus_show_in_block=0.6; % 0: any trial in block, 1: last?
earliest_reward_stimulus_show_in_block=0.6; % 0: any t3rial in block, 1: last?. must be in [reward_prob] range.
sec_reward_display=2/ffw_speed;
% how_early_novel_stimulus_show_in_trial=1; % 0: shown at 1-st trial, 1: any trial in block.
latest_novel_stimulus_show_in_trial_s=0.6; % 0: shown at 1-st trial_s, 1: any trial_s in trial.

% reward
reward_mag=[-50 10];
reward_prob=[0.25 0.75]; % must be 2nd col>1st colyy

% money
money_given=1; % point
% novelty-level for the non-novel stimuli. a novel stimulus will be shown only one time.
img_present_ratio=[2 1]; % novelty level :: size = tot_num_img_per_block-1.

% questions
Q{1,1}.text='*** How much do you like this picture? (hedonic rating) [L/R] to change, [y] to confirm. ***';
Q{1,1}.text2='(Dislike)';
Q{1,1}.text3=' (Like)';
Q{1,1}.text4='(Don''t know)';
Q{1,1}.color=[0, 0, 255, 255];
Q{1,1}.cursor_input=[-5:1:5];
Q{1,2}.text='*** How much do you think this causes a bad outcome?  [L/R] to change, [y] to confirm. ***';
Q{1,2}.text2='(Not at all)';
Q{1,2}.text3='(Very likely)';
Q{1,2}.text4='(Don''t know)';
Q{1,2}.color=[255, 0, 0, 255];
Q{1,2}.cursor_input=[0:1:10];

% text size
text_size_default=20; % font size (don't change)
text_size_reward=400; % height in pixel

% background color
BackgroundColor_block_intro=[130,130,130,150]; % gray
BackgroundColor_Cue_page=[210,210,210,150]; % light gray
BackgroundColor_Trial_ready_page=[210,210,210,150]; % light gray
BackgroundColor_Reward_page=[210,210,210,150]; % light gray
COLOR_FIXATION_MARK=[70,70,70,200]; % dark gray
% automatically determined
total_bet_trial=nchoosek(tot_num_img_per_block,3);
list_combination=combntns(1:tot_num_img_per_block,3);

% key code
KEY_L=37;
KEY_R=39;
KEY_Y=89; %'y'
KEY_N=78; %'n'
KEY_Q=81; %'q'
KEY_T=84; % 't', 5 in desktop, 84 in laptop


% recording variables (for each block)
HIST_bet_score_table_type1=cell(Tot_block,size(list_combination,1));
HIST_bet_money_table_type1=cell(Tot_block,size(list_combination,1));
HIST_Barycentric_coord_table_type1=cell(Tot_block,size(list_combination,1));
HIST_bet_score_table_type2=cell(Tot_block,size(Q,2));

HIST_bet_score_table_type1_bonus=cell(1,size(list_combination,1));
HIST_bet_money_table_type1_bonus=cell(1,size(list_combination,1));
HIST_Barycentric_coord_table_type1_bonus=cell(1,size(list_combination,1));
HIST_bet_score_table_type2_bonus=cell(1,size(Q,2));

HIST_event_info=[]; % row1=event time(in session), row2=event time(in block), row3=event type
HIST_event_type{1,1}='row1 - block#';    HIST_event_type{1,2}='row2 - trial#, 0 if outside of the trial';     HIST_event_type{1,3}='row3 - trial_s#, 0 if outside of the trial_s'; 
HIST_event_type{1,4}='row4 - event time in session';   HIST_event_type{1,5}='row5 - event time in block';
HIST_event_type{1,6}='row6 - 0.5: block start msg on, -0.5: block start msg off, 1: cue novelty level1, 2: cue novelty level2, 3: cue novelty level3, 4: reward delivery, 5: hedonic rating display, 6: hedonic rating answer, 7: causal rating display, 8: causal rating answer, 9: Bayesian rating display, 10: Bayesian rating answer, 20: a short blank page display, -99:fail to do ratings in time limit, (-) when display off';




%% seed image read
if(index_num==1) % create the image usage matrix if this is the first session
    Tot_num_img=120;%Tot_block*tot_num_img_per_block;
    Tot_img_info=[1:1:Tot_num_img];
    Tot_img_info=[Tot_img_info; zeros(8,Tot_num_img)]; % 2nd row: used index, 3rd row: session#, 4th row: block#, 5th row: novelty level (1:non-novel, 3:novel)
    Tot_img_info_Tag{1,1}='row1 - file #';    Tot_img_info_Tag{1,2}='row2 - is used';
    Tot_img_info_Tag{1,3}='row3 - session #';    Tot_img_info_Tag{1,4}='row4 - block #';
    Tot_img_info_Tag{1,5}='row5 - novelty level (1:non-novel, 3:novel)';        Tot_img_info_Tag{1,6}='row6 - ratings (hedonic)';
    Tot_img_info_Tag{1,7}='row7 - rating (bad outcome causality)';     Tot_img_info_Tag{1,8}='row8 - ratings (bayesian for, the bonus round)';
    Tot_img_info_Tag{1,9}='row9 - case 0: nov-pos, 1:nov-neg';
else % load image usage matrix to update
    file_imgind_ld_name=[EXP_NAME '_image_usage_info.mat'];
    file_name_ld=[save_path file_imgind_ld_name];
    load(file_name_ld);
end

img_set=cell(Tot_block,tot_num_img_per_block-1);
img0_set=cell(Tot_block,1);
img_set_bonus=cell(Tot_block,tot_num_img_per_block-1);
img0_set_bonus=cell(Tot_block,1);

HIST_current_pool=[];
for block=1:1:Tot_block
    
    % select all images to be used for each block (nonoverlap among blocks)
    [tmp current_img_pool_ind]=find(Tot_img_info(2,:)==0);
    current_pool=Tot_img_info(1,current_img_pool_ind);
    img_index_use_all=randsample(current_pool,tot_num_img_per_block,false);
    
    % update image usage matrix
    Tot_img_info(2,img_index_use_all)=1; % tagged as "used"
    Tot_img_info(3,img_index_use_all)=index_num; % record the current session#
    Tot_img_info(4,img_index_use_all)=block; % record the current block#
    Tot_img_info(5,img_index_use_all(1))=1; % record the novelty level
    Tot_img_info(5,img_index_use_all(2))=2; % record the novelty level
    Tot_img_info(5,img_index_use_all(3))=3; % record the novelty level
    
    % assign normal&novel stimulus
    img_index_use=img_index_use_all(1:end-1);
    img0_index_use=img_index_use_all(end);
    
    num_stim=length(img_index_use);
    for i=1:1:num_stim
        file_full_path=[seed_path sprintf('%03d.png',img_index_use(i))];
        img_set{block,i}=imresize(imread(file_full_path),[IMAGE_SIZE(2) IMAGE_SIZE(1)]); % get a frame
    end
    num_stim0=length(img0_index_use);
    for i=1:1:num_stim0
        file_full_path=[seed_path sprintf('%03d.png',img0_index_use(i))];
        img0_set{block,i}=imresize(imread(file_full_path),[IMAGE_SIZE(2) IMAGE_SIZE(1)]); % get a frame
    end
    HIST_current_pool=[HIST_current_pool; [img_index_use img0_index_use]];
    
end

% reward message image read
img_msg_set=cell(length(reward_mag),1);
for jj=1:1:length(reward_mag)
    if(reward_mag(jj)>=0)
        file_full_path=[seed_path sprintf('msg_p%02dpt.png',abs(reward_mag(jj)))];
    else
        file_full_path=[seed_path sprintf('msg_n%02dpt.png',abs(reward_mag(jj)))];
    end
    img_msg_set{jj,1}=imresize(imread(file_full_path),[text_size_reward NaN]);
end


%% Scheduling
img_present_freq=img_present_ratio/sum(img_present_ratio);

ind_mat=[];     num_acc=0;
tot_num=Tot_trial*Tot_trial_s;
HIST_schedule=cell(1,Tot_block);
HIST_reward=zeros(Tot_block,Tot_trial);
for j=1:1:tot_num_img_per_block-2
    numm=round((tot_num-1)*img_present_freq(j));
    num_acc=num_acc+numm;
    ind_mat=[ind_mat j*ones(1,numm)];
end
ind_mat=[ind_mat (tot_num_img_per_block-1)*ones(1,tot_num-1-num_acc)]; % second last column
for k=1:1:Tot_block
    rand_order=randperm(tot_num-1);
    schedule_index_mat=ind_mat(rand_order); % now we get index matrix without a novel stimulus
    % 1. plug in a single novel stimulus
    pos=min(max(1,ceil((tot_num*(1-earliest_novel_stimulus_show_in_block))*rand))+tot_num*earliest_novel_stimulus_show_in_block,tot_num);
    schedule_index_mat_final=[schedule_index_mat(1:pos-1) tot_num_img_per_block schedule_index_mat(pos:end)];
    HIST_schedule{1,k}=reshape(schedule_index_mat_final,Tot_trial_s,Tot_trial)';
    % 1.5. check if novel stimulus apears early in trial_s
    [b_trial b_trial_s]=find(HIST_schedule{1,k}==3);
    margin_s=max(1,floor(latest_novel_stimulus_show_in_trial_s*Tot_trial_s));
    if(b_trial_s>margin_s) % switch two cues
        b_new=max(1,ceil(margin_s*rand));
        tmp=HIST_schedule{1,k}(b_trial,b_trial_s);
        HIST_schedule{1,k}(b_trial,b_trial_s)=HIST_schedule{1,k}(b_trial,b_new);
        HIST_schedule{1,k}(b_trial,b_new)=tmp;
    end
    % 2. reward scheduling
    earliest_trial_r=round(earliest_reward_stimulus_show_in_block*Tot_trial);
    HIST_reward(k,1:earliest_trial_r)=reward_mag(2); % good reward
    % compute prob for the remaining trials
    remain_prob=[reward_prob(1) reward_prob(2)-earliest_reward_stimulus_show_in_block];
    remain_prob=remain_prob/sum(remain_prob);
    bad_reward_ind_0=round(remain_prob(1)*(Tot_trial-earliest_trial_r));
    reward_mat_0=[reward_mag(1)*ones(1,bad_reward_ind_0), reward_mag(2)*ones(1,Tot_trial-earliest_trial_r-bad_reward_ind_0)];
    rand_order_r=randperm(Tot_trial-earliest_trial_r);
    schedule_reward_mat_0=reward_mat_0(rand_order_r);
    HIST_reward(k,earliest_trial_r+1:Tot_trial)=schedule_reward_mat_0;
end

% balancing novel-positive_outcome VS novel-negative_outcome
nov_neg_block_index=[];     nov_pos_block_index=[];
for k=1:1:Tot_block
    [row1 col1]=find(HIST_schedule{1,k}==3); % identify the novel stimulus trial
    [row2 col2]=find(HIST_reward(k,:)==reward_mag(1)); % identify the negative outcome trial
    if(row1==col2) % if nov-neg
        nov_neg_block_index=[nov_neg_block_index k];
    else
        nov_pos_block_index=[nov_pos_block_index k];
    end
end
if(length(nov_neg_block_index)>length(nov_pos_block_index))
    block_picked=randsample(nov_neg_block_index,(length(nov_neg_block_index)-length(nov_pos_block_index))/2,false); % pick trials to be corrected
    for hh=1:1:length(block_picked) % correct
        [row1 col1]=find(HIST_schedule{1,block_picked(hh)}==3); % identify the novel stimulus trial                
        tmp_ind=HIST_reward(block_picked(hh),row1);
        if(row1==Tot_trial)
            HIST_reward(block_picked(hh),row1)=HIST_reward(block_picked(hh),row1-1);
            HIST_reward(block_picked(hh),row1-1)=tmp_ind;
        else
            HIST_reward(block_picked(hh),row1)=HIST_reward(block_picked(hh),row1+1);
            HIST_reward(block_picked(hh),row1+1)=tmp_ind;
        end
        nov_neg_block_index(find(nov_neg_block_index==block_picked(hh)))=[];
        nov_pos_block_index=[nov_pos_block_index block_picked(hh)];
    end
end
if(length(nov_pos_block_index)>length(nov_neg_block_index))
    block_picked=randsample(nov_pos_block_index,(length(nov_pos_block_index)-length(nov_neg_block_index))/2,false); % pick trials to be corrected
    for hh=1:1:length(block_picked) % correct
        [row1 col1]=find(HIST_schedule{1,block_picked(hh)}==3); % identify the novel stimulus trial
        HIST_reward(block_picked(hh),:)=reward_mag(2);
        HIST_reward(block_picked(hh),row1)=reward_mag(1);
        nov_pos_block_index(find(nov_pos_block_index==block_picked(hh)))=[];
        nov_neg_block_index=[nov_neg_block_index block_picked(hh)];
    end
end


%% Bonus block scheduling
if(KEEP_PICTURE_OPT==1)
    ind_mat=[];     num_acc=0;
    tot_num=Tot_trial*Tot_trial_s;
    HIST_schedule_bonus=cell(1,1);
    HIST_reward_bonus=zeros(1,Tot_trial);
    for j=1:1:tot_num_img_per_block-2
        numm=round((tot_num-1)*img_present_freq(j));
        num_acc=num_acc+numm;
        ind_mat=[ind_mat j*ones(1,numm)];
    end
    ind_mat=[ind_mat (tot_num_img_per_block-1)*ones(1,tot_num-1-num_acc)]; % second last column
    
    rand_order=randperm(tot_num-1);
    schedule_index_mat=ind_mat(rand_order); % now we get index matrix without a novel stimulus
    % 1. plug in a single novel stimulus
    pos=min(max(1,ceil((tot_num*(1-earliest_novel_stimulus_show_in_block))*rand))+tot_num*earliest_novel_stimulus_show_in_block,tot_num);
    schedule_index_mat_final=[schedule_index_mat(1:pos-1) tot_num_img_per_block schedule_index_mat(pos:end)];
    HIST_schedule_bonus{1,1}=reshape(schedule_index_mat_final,Tot_trial_s,Tot_trial)';
    % 1.5. check if novel stimulus apears early in trial_s
    [b_trial b_trial_s]=find(HIST_schedule_bonus{1,1}==3);
    margin_s=max(1,floor(latest_novel_stimulus_show_in_trial_s*Tot_trial_s));
    if(b_trial_s>=margin_s) % switch two cues
        b_new=max(1,floor(margin_s*rand));
        tmp=HIST_schedule_bonus{1,1}(b_trial,b_trial_s);
        HIST_schedule_bonus{1,1}(b_trial,b_trial_s)=HIST_schedule_bonus{1,1}(b_trial,b_new);
        HIST_schedule_bonus{1,1}(b_trial,b_new)=tmp;
    end
    % 2. reward scheduling
    earliest_trial_r=round(earliest_reward_stimulus_show_in_block*Tot_trial);
    HIST_reward_bonus(1,1:earliest_trial_r)=reward_mag(2); % good reward
    % compute prob for the remaining trials
    remain_prob=[reward_prob(1) reward_prob(2)-earliest_reward_stimulus_show_in_block];
    remain_prob=remain_prob/sum(remain_prob);
    bad_reward_ind_0=round(remain_prob(1)*(Tot_trial-earliest_trial_r));
    reward_mat_0=[reward_mag(1)*ones(1,bad_reward_ind_0), reward_mag(2)*ones(1,Tot_trial-earliest_trial_r-bad_reward_ind_0)];
    rand_order_r=randperm(Tot_trial-earliest_trial_r);
    schedule_reward_mat_0=reward_mat_0(rand_order_r);
    HIST_reward_bonus(1,earliest_trial_r+1:Tot_trial)=schedule_reward_mat_0;
end




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
HIST_key_press_button=ones(1,Tot_trial);
HIST_key_press_sec=ones(1,Tot_trial);

%% starting message
Screen('TextSize',wPtr, text_size_default);
% Screen('TextFont',wPtr, 'Times New Roman');
DrawFormattedText(wPtr, 'Are you ready for the experiment?\n(Press any key to wait for the trigger)', 'center', 'center');
Screen('Flip', wPtr);  
% KbWait; % temporarily disabled for test APR 21
if(DO_TAKE_SNAPSHOT==1)
    snapshot=Screen(wPtr, 'GetImage', [1, 1, floor(1.0*SCREEN_RESOLUTION(1)), floor(1.0*SCREEN_RESOLUTION(2))]);
    imageArray=[imageArray; {snapshot}];
end


%% waiting for the trigger sign from the scanner
DrawFormattedText(wPtr, 'Now waiting for the trigger...', 'center', 'center');
Screen('Flip',wPtr);
% Look for trigger pulse
while 0 % temporarily disabled for test APR 21
    [ keyIsDown, timeSecs, keyCode ] = KbCheck;
    if keyIsDown
%         [tmp tmp_key_code]=find(keyCode==1);
%         if ((KbName(keyCode) == 't')||((tmp_key_code==KEY_T)))
        if (KbName(keyCode) == 't')
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
img_set_all=cell(Tot_block,tot_num_img_per_block);
for block=1:1:Tot_block % each block
    
    %% I. Stimuli presentation
    
    zzz=[];    
    % block starts
    Screen('FillRect',wPtr,BackgroundColor_block_intro);
    str_block_intro=sprintf('*** Now block %d starts. ***',block);
    DrawFormattedText(wPtr, str_block_intro, 'center', 'center',[0, 0, 0, 255]);
    Screen(wPtr, 'Flip');
    block_clock_start = GetSecs;    
    HIST_event_info=[HIST_event_info [block; 0; 0; (GetSecs-session_clock_start); (GetSecs-block_clock_start); 0.5]]; % event save
    WaitSecs(sec_block_ready);
  
    % add fixation mark and display off during the jittered interval    
    DrawFormattedText(wPtr, '+', 'center', 'center', COLOR_FIXATION_MARK); % add 'o' mark at the click pt.
    Screen(wPtr, 'Flip');
    HIST_event_info=[HIST_event_info [block; 0; 0; (GetSecs-session_clock_start); (GetSecs-block_clock_start); -0.5]]; % event save
    sec_stim_interval0=rand*(max(sec_stim_interval)-min(sec_stim_interval))+min(sec_stim_interval);
    WaitSecs(sec_stim_interval0);

    % take a snapshot
    if(DO_TAKE_SNAPSHOT==1)
        snapshot=Screen(wPtr, 'GetImage', [1, 1, floor(1.0*SCREEN_RESOLUTION(1)), floor(1.0*SCREEN_RESOLUTION(2))]);
        imageArray=[imageArray; {snapshot}];
    end
    
    img0_index_showed=zeros(1,length(img0_index_use));
    
    
    for trial=1:1:Tot_trial % each trial
        
        %         novel_trial_s_ind=max(1,ceil((Tot_trial_s*how_early_novel_stimulus_show_in_trial_s)*rand));
        
        Screen('FillRect',wPtr,BackgroundColor_Trial_ready_page)
        
        % ready sign
%         str_ready=sprintf('Ready for the %d-th trial',trial);
%         DrawFormattedText(wPtr, str_ready, 'center', 'center',[0, 0, 0, 255]);
%         Screen(wPtr, 'Flip');
%         WaitSecs(sec_trial_ready);
        % take a snapshot
        if(DO_TAKE_SNAPSHOT==1)
            snapshot=Screen(wPtr, 'GetImage', [1, 1, floor(1.0*SCREEN_RESOLUTION(1)), floor(1.0*SCREEN_RESOLUTION(2))]);
            imageArray=[imageArray; {snapshot}];
        end
        
        
        for trial_s=1:1:Tot_trial_s % cue presentation
            
            %             % ready sign
            %             str_ready=sprintf('(%d/%d) cue ready',trial_s,Tot_trial_s);
            %             DrawFormattedText(wPtr, str_ready, 'center', 'center',[0, 0, 0, 255]);
            %             Screen(wPtr, 'Flip');
            %             WaitSecs(sec_stim_ready);
            
            % take a snapshot
            if(DO_TAKE_SNAPSHOT==1)
                snapshot=Screen(wPtr, 'GetImage', [1, 1, floor(1.0*SCREEN_RESOLUTION(1)), floor(1.0*SCREEN_RESOLUTION(2))]);
                imageArray=[imageArray; {snapshot}];
            end
            
            % Determine which stimlus will be presented.
            selected=HIST_schedule{1,block}(trial,trial_s);
            if(selected==tot_num_img_per_block) % the most novel cue
                input_stim = Screen('MakeTexture', wPtr, img0_set{block,1});
                img0_index_showed(1)=1;
            else
                input_stim = Screen('MakeTexture', wPtr, img_set{block,selected});
            end
            %             if((trial==novel_trial_ind)&&(trial_s==novel_trial_s_ind))
            %                 % novel trial
            %                 novel_trial_img_ind=ceil(length(img0_index_use)*rand);
            %                 selected=novel_trial_img_ind;
            %                 input_stim = Screen('MakeTexture', wPtr, img0_set{block,selected});
            %                 img0_index_showed(novel_trial_img_ind)=1;
            %             else
            %                 % normal trial
            %                 cumnum=cumsum(img_present_freq);
            %                 [temp1 selected]=histc(rand,[0 cumnum]);
            %                 input_stim = Screen('MakeTexture', wPtr, img_set{block,selected});
            %             end
            
            % image display
            xpos = round(screenWidth/2);    ypos = round(screenHeight/2);
            sx=floor(IMAGE_SIZE(1)*disp_scale);       sy=floor(IMAGE_SIZE(2)*disp_scale);
            destrect=[xpos-sx/2,ypos-sy/2,xpos+sx/2,ypos+sy/2];
            Screen('FillRect',wPtr,BackgroundColor_Cue_page);
            Screen('DrawTexture', wPtr, input_stim,[],destrect);
            %             % cue number display
            %             str_ready=sprintf('Cue (%d/%d)',trial_s,Tot_trial_s);
            %             DrawFormattedText(wPtr, str_ready, 'center', ypos+sy/2+20,[0, 0, 0, 255]);
            Screen(wPtr, 'Flip'); % display on 
            HIST_event_info=[HIST_event_info [block; trial; trial_s; (GetSecs-session_clock_start); (GetSecs-block_clock_start); selected]]; % event save
            WaitSecs(sec_stim_display);
            % take a snapshot
            if(DO_TAKE_SNAPSHOT==1)
                snapshot=Screen(wPtr, 'GetImage', [1, 1, floor(1.0*SCREEN_RESOLUTION(1)), floor(1.0*SCREEN_RESOLUTION(2))]);
                imageArray=[imageArray; {snapshot}];
            end
            % add fixation mark and display off during the jittered interval
            DrawFormattedText(wPtr, '+', 'center', 'center', COLOR_FIXATION_MARK); % add 'o' mark at the click pt.
            Screen(wPtr, 'Flip');
            HIST_event_info=[HIST_event_info [block; trial; trial_s; (GetSecs-session_clock_start); (GetSecs-block_clock_start); (-1)*selected]]; % event save
            sec_stim_interval0=rand*(max(sec_stim_interval)-min(sec_stim_interval))+min(sec_stim_interval);
            WaitSecs(sec_stim_interval0);
            % take a snapshot
            if(DO_TAKE_SNAPSHOT==1)
                snapshot=Screen(wPtr, 'GetImage', [1, 1, floor(1.0*SCREEN_RESOLUTION(1)), floor(1.0*SCREEN_RESOLUTION(2))]);
                imageArray=[imageArray; {snapshot}];
            end
            
        end
        
        % reward display
        %         Screen('TextSize',wPtr, text_size_reward);
        %         str_rwd=sprintf('### You have earned %d points. ###',HIST_reward(block,trial));
        %         DrawFormattedText(wPtr, str_rwd, 'center', 'center',[255, 0, 0, 255]);
        Screen('FillRect',wPtr,BackgroundColor_Reward_page);
        if(HIST_reward(block,trial)<0)
            input_stim_msg = Screen('MakeTexture', wPtr, img_msg_set{1,1});
            sx_msg=size(img_msg_set{1,1},2);    sy_msg=size(img_msg_set{1,1},1);
        else
            input_stim_msg = Screen('MakeTexture', wPtr, img_msg_set{2,1});
            sx_msg=size(img_msg_set{2,1},2);    sy_msg=size(img_msg_set{2,1},1);
        end
        destrect=[xpos-sx_msg/2,ypos-sy_msg/2,xpos+sx_msg/2,ypos+sy_msg/2];
        Screen('DrawTexture', wPtr, input_stim_msg,[],destrect);
%         str_rwd2=sprintf('### Press any button to continue. ###');
%         DrawFormattedText(wPtr, str_rwd2, 'center', round(screenHeight*4/5),[0, 0, 0, 255]);
        Screen(wPtr, 'Flip');
        HIST_event_info=[HIST_event_info [block; trial; 0; (GetSecs-session_clock_start); (GetSecs-block_clock_start); 4]]; % event save
        WaitSecs(sec_reward_display);
        
        % add fixation mark and display off during the jittered interval
        DrawFormattedText(wPtr, '+', 'center', 'center', COLOR_FIXATION_MARK); % add 'o' mark at the click pt.
        Screen(wPtr, 'Flip');
        HIST_event_info=[HIST_event_info [block; trial; 0; (GetSecs-session_clock_start); (GetSecs-block_clock_start); -4]]; % event save
        sec_trial_interval0=rand*(max(sec_trial_interval)-min(sec_trial_interval))+min(sec_trial_interval);
        WaitSecs(sec_trial_interval0);
        
        %         Screen('TextSize',wPtr, text_size_default);
        
        % take a snapshot
        if(DO_TAKE_SNAPSHOT==1)
            snapshot=Screen(wPtr, 'GetImage', [1, 1, floor(1.0*SCREEN_RESOLUTION(1)), floor(1.0*SCREEN_RESOLUTION(2))]);
            imageArray=[imageArray; {snapshot}];
        end

%         % button press - Decision-making
%         % 'q':81, esc: 27
%         decision_made=0;
%         exit_made=0;
%         while(~decision_made)
%             [secs, keyCode] = KbPressWait;
%             [tmp tmp_key_code]=find(keyCode==1);                zzz=[zzz tmp_key_code];
%             if(tmp_key_code==KEY_Q) % 'q' pressed for aborting
%                 clear mex
%                 decision_made=1;
%             end
%             decision_made=1;
%         end
%         WaitSecs(0.1);
%         % take a snapshot
%         if(DO_TAKE_SNAPSHOT==1)
%             snapshot=Screen(wPtr, 'GetImage', [1, 1, floor(1.0*SCREEN_RESOLUTION(1)), floor(1.0*SCREEN_RESOLUTION(2))]);
%             imageArray=[imageArray; {snapshot}];
%         end
        
    end
    
    
    % collect all the images used
    for h=1:1:length(img_index_use)
        img_set_all{block,h}=img_set{block,h};
    end
    [tmp ind_new]=find(img0_index_showed==1);
    for h=1:1:sum(img0_index_showed)
        img_set_all{block,h+length(img_index_use)}=img0_set{block,ind_new(h)};
    end
    
    
    
    
    
    %% III. Rating - type2
    
    img_toshow_index=randperm(tot_num_img_per_block);
    xpos = round(screenWidth/2);    ypos = round(screenHeight/2);
    length_bar=round(screenWidth*3/5);
    bar_pt_start=[round(screenWidth*1/5),ypos+250];
    bar_pt_end=[round(screenWidth*1/5)+length_bar,ypos+250];
    bet_score_table=[img_toshow_index; zeros(1,3)];
    bet_money_table=[img_toshow_index; zeros(1,3)];
    
    %     Screen('TextSize',wPtr, text_size_default);
    
    for q_ind=1:1:size(Q,2)
        
        ggg_ind=0;
        for m=1:1:tot_num_img_per_block
            ggg_ind=ggg_ind+1;
            
            
            Is_bet=0;
            bar_pt_cursor=(Q{1,q_ind}.cursor_input-min(Q{1,q_ind}.cursor_input))/(max(Q{1,q_ind}.cursor_input)-min(Q{1,q_ind}.cursor_input)); %[0,1]
            bar_pt_cursor=length_bar*bar_pt_cursor+bar_pt_start(1); %actual x-positions in display
            current_cursor_selection_ind=round(length(Q{1,q_ind}.cursor_input)/2);
            ind_ev=0;
            
            while(~Is_bet)
                
                ind_ev=ind_ev+1;
                
                % message
                str_block_intro=Q{1,q_ind}.text;
                DrawFormattedText(wPtr, str_block_intro, 'center', round(screenHeight/7) ,Q{1,q_ind}.color);
                % 2. Show Bars and images
                sx=floor(IMAGE_SIZE(1)*disp_scale);       sy=floor(IMAGE_SIZE(2)*disp_scale);
                % draw bar
                Screen('DrawLine', wPtr, Q{1,q_ind}.color, bar_pt_start(1), bar_pt_start(2), bar_pt_end(1), bar_pt_end(2),[3]);
                for jj=1:1:length(bar_pt_cursor)
                    Screen('DrawLine', wPtr, Q{1,q_ind}.color, bar_pt_cursor(jj), bar_pt_start(2)-5, bar_pt_cursor(jj), bar_pt_start(2)+5,[1]);
                    str_block_num=sprintf('%s',num2str(Q{1,q_ind}.cursor_input(jj)));
                    DrawFormattedText(wPtr, str_block_num, bar_pt_cursor(jj)-10, bar_pt_start(2)+25, Q{1,q_ind}.color);
                    if(jj==1) % like/dislike, notatall/verylikeli
                        DrawFormattedText(wPtr, Q{1,q_ind}.text2, bar_pt_cursor(jj)-50, bar_pt_start(2)+50, Q{1,q_ind}.color);
                        Screen('DrawLine', wPtr, Q{1,q_ind}.color, bar_pt_cursor(jj), bar_pt_start(2)-8, bar_pt_cursor(jj), bar_pt_start(2)+8,[3]);
                    end
                    if(jj==length(bar_pt_cursor)) % like/dislike, notatall/verylikeli
                        DrawFormattedText(wPtr, Q{1,q_ind}.text3, bar_pt_cursor(jj)-50, bar_pt_start(2)+50, Q{1,q_ind}.color);
                        Screen('DrawLine', wPtr, Q{1,q_ind}.color, bar_pt_cursor(jj), bar_pt_start(2)-8, bar_pt_cursor(jj), bar_pt_start(2)+8,[3]);
                    end
                    if(jj==round(length(Q{1,q_ind}.cursor_input)/2)) % don't know
                        DrawFormattedText(wPtr, Q{1,q_ind}.text4, bar_pt_cursor(jj)-80, bar_pt_start(2)+50, Q{1,q_ind}.color);
                        Screen('DrawLine', wPtr, Q{1,q_ind}.color, bar_pt_cursor(jj), bar_pt_start(2)-8, bar_pt_cursor(jj), bar_pt_start(2)+8,[3]);
                    end
                end
                % draw images
                input_stim = Screen('MakeTexture', wPtr, img_set_all{block,img_toshow_index(m)});
                img_pt_center=[xpos, bar_pt_start(2)-sy/2-100];
                destrect=[img_pt_center(1)-sx/2,img_pt_center(2)-sy/2,img_pt_center(1)+sx/2,img_pt_center(2)+sy/2];
                Screen('DrawTexture', wPtr, input_stim,[],destrect);
                % 4-2. add 'o' mark at the click pt.
                DrawFormattedText(wPtr, 'o', bar_pt_cursor(current_cursor_selection_ind)-8, bar_pt_start(2)-14, [0 0 0 255]); % add 'o' mark at the click pt.
                % display!
                Screen(wPtr, 'Flip');
                if(ind_ev==1)
                    HIST_event_info=[HIST_event_info [block; 0; 0; (GetSecs-session_clock_start); (GetSecs-block_clock_start); (5+2*(q_ind-1))]]; % event save
                    clock_time_limit_start=GetSecs;
                end
                % take a snapshot
                if(DO_TAKE_SNAPSHOT==1)
                    snapshot=Screen(wPtr, 'GetImage', [1, 1, floor(1.0*SCREEN_RESOLUTION(1)), floor(1.0*SCREEN_RESOLUTION(2))]);
                    imageArray=[imageArray; {snapshot}];
                end
                
                % 1. Get cursor (L:37 R:39)
                decision_made=0;
                while(~decision_made)                    
                    [secs, keyCode] = KbPressWait([], clock_time_limit_start+sec_limit_Q_rating(q_ind)); % if no keyboard in time limit, then go ahead. if pressed earlier, then go ahead.
                    [tmp tmp_key_code]=find(keyCode==1);                zzz=[zzz tmp_key_code];
                    if(tmp_key_code==KEY_L) % L pressed
                        current_cursor_selection_ind=max(1,current_cursor_selection_ind-1);
                        decision_made=1;
                    end
                    if(tmp_key_code==KEY_R) % R pressed
                        current_cursor_selection_ind=min(length(Q{1,q_ind}.cursor_input),current_cursor_selection_ind+1);
                        decision_made=1;
                    end
                    if(tmp_key_code==KEY_Y) % 'y' pressed
                        decision_made=1;            Is_bet=1;
                        bet_score_table(2,ggg_ind)=Q{1,q_ind}.cursor_input(current_cursor_selection_ind);                     
                        HIST_event_info=[HIST_event_info [block; 0; 0; (GetSecs-session_clock_start); (GetSecs-block_clock_start); (6+2*(q_ind-1))]]; % event save
                    end
                    if(tmp_key_code==KEY_Q) % 'q' pressed for aborting
                        clear mex
                        decision_made=1;
                    end
                    % check the time limit !@#$
                    if((GetSecs-clock_time_limit_start)>sec_limit_Q_rating(q_ind))
                        decision_made=1;            Is_bet=1;
                        dont_know_cursor_selection=(length(Q{1,q_ind}.cursor_input)+1)/2;
                        bet_score_table(2,ggg_ind)=Q{1,q_ind}.cursor_input(current_cursor_selection_ind);%Q{1,q_ind}.cursor_input(dont_know_cursor_selection);
                        bet_score_table(3,ggg_ind)=-99; % error code
                        HIST_event_info=[HIST_event_info [block; 0; 0; (GetSecs-session_clock_start); (GetSecs-block_clock_start); -99]]; % event save                        
                    end
                end
                
                
                % take a snapshot
                if(DO_TAKE_SNAPSHOT==1)
                    snapshot=Screen(wPtr, 'GetImage', [1, 1, floor(1.0*SCREEN_RESOLUTION(1)), floor(1.0*SCREEN_RESOLUTION(2))]);
                    imageArray=[imageArray; {snapshot}];
                end
                
            end
            
            % a short blank page : to give a sense of the page transition
            Screen(wPtr, 'Flip');
            HIST_event_info=[HIST_event_info [block; 0; 0; (GetSecs-session_clock_start); (GetSecs-block_clock_start); 20]]; % event save
            WaitSecs(sec_jittered_blank_page);
            
        end
        HIST_bet_score_table_type2{block,q_ind}=bet_score_table; % same as the above
    end
    
    
    
    
    
    if(KEEP_PICTURE_OPT==1)
        %% II. Rating - type1
        
        % compute points of a triangle
        
        tri_pt=zeros(3,2); % edges
        size_triangle=200;
        distance_limit=size_triangle+size_triangle*cos(pi/3);
        xpos = round(screenWidth/2);    ypos = round(screenHeight/2)+70;
        tri_pt(1,:)=[xpos,ypos-size_triangle];
        tri_pt(2,:)=[xpos-size_triangle*sin(pi/3),ypos+size_triangle*cos(pi/3)];
        tri_pt(3,:)=[xpos+size_triangle*sin(pi/3),ypos+size_triangle*cos(pi/3)];
        
        tri_pt_all=zeros(25,2); % all inner points including edges
        tri_pt_all(1,:)=tri_pt(1,:);    tri_pt_all(21,:)=tri_pt(2,:);    tri_pt_all(25,:)=tri_pt(3,:);  tri_pt_all(12,:)=[xpos ypos];
        tri_pt_all(4,:)=[tri_pt_all(1,1) tri_pt_all(1,2)+size_triangle/2];
        tri_pt_all(7,:)=[tri_pt_all(1,1) tri_pt_all(1,2)+3*size_triangle/4];
        tri_pt_all(17,:)=[tri_pt_all(12,1) tri_pt_all(12,2)+size_triangle/4];
        tri_pt_all(2,:)=[tri_pt_all(1,1) tri_pt_all(1,2)+size_triangle/4];
        th_ang=pi/3-atan(cos(pi/6));
        tri_pt_all(3,:)=[tri_pt_all(12,1)-sqrt(7)/4*size_triangle*sin(th_ang) tri_pt_all(12,2)-sqrt(7)/4*size_triangle*cos(th_ang)];
        tri_pt_all(5,:)=[tri_pt_all(12,1)+sqrt(7)/4*size_triangle*sin(th_ang) tri_pt_all(12,2)-sqrt(7)/4*size_triangle*cos(th_ang)];
        tri_pt_all(6,:)=[tri_pt_all(7,1)-size_triangle*tan(pi/3)/4 tri_pt_all(7,2)];
        tri_pt_all(8,:)=[tri_pt_all(7,1)+size_triangle*tan(pi/3)/4 tri_pt_all(7,2)];
        tri_pt_all(9,:)=[tri_pt_all(12,1)-size_triangle*cos(pi/6)/4 tri_pt_all(12,2)-size_triangle*sin(pi/6)/4];
        tri_pt_all(10,:)=[tri_pt_all(12,1)+size_triangle*cos(pi/6)/4 tri_pt_all(12,2)-size_triangle*sin(pi/6)/4];
        th_ang=atan(cos(pi/6))-pi/6;
        tri_pt_all(11,:)=[tri_pt_all(12,1)-sqrt(7)/4*size_triangle*cos(th_ang) tri_pt_all(12,2)+sqrt(7)/4*size_triangle*sin(th_ang)];
        tri_pt_all(13,:)=[tri_pt_all(12,1)+sqrt(7)/4*size_triangle*cos(th_ang) tri_pt_all(12,2)+sqrt(7)/4*size_triangle*sin(th_ang)];
        tri_pt_all(14,:)=[tri_pt_all(12,1)-size_triangle*cos(pi/6)/4 tri_pt_all(12,2)+size_triangle*sin(pi/6)/4];
        tri_pt_all(15,:)=[tri_pt_all(12,1)+size_triangle*cos(pi/6)/4 tri_pt_all(12,2)+size_triangle*sin(pi/6)/4];
        tri_pt_all(16,:)=[tri_pt_all(17,1)-size_triangle*cos(pi/6)/2 tri_pt_all(17,2)];
        tri_pt_all(18,:)=[tri_pt_all(17,1)+size_triangle*cos(pi/6)/2 tri_pt_all(17,2)];
        tri_pt_all(19,:)=[tri_pt_all(12,1)-size_triangle*3*cos(pi/6)/4 tri_pt_all(12,2)+size_triangle*3/8];
        tri_pt_all(20,:)=[tri_pt_all(12,1)+size_triangle*3*cos(pi/6)/4 tri_pt_all(12,2)+size_triangle*3/8];
        tri_pt_all(22,:)=[tri_pt_all(21,1)+size_triangle*cos(pi/6)/2 tri_pt_all(21,2)];
        tri_pt_all(23,:)=[tri_pt_all(12,1) tri_pt_all(12,2)+size_triangle/2];
        tri_pt_all(24,:)=[tri_pt_all(25,1)-size_triangle*cos(pi/6)/2 tri_pt_all(25,2)];
        tri_pt_all=round(tri_pt_all);
        tri_pt_arrange_tbl=[1 2 3 5 4 6 7 8 9 10 12 11 14 15 13 16 17 18 19 20 21:1:25];
        
        
        for kk=1:1:size(list_combination,1)
            
            img_toshow_index=list_combination(kk,randperm(3));
            
            
            Is_bet=0;
            bet_score_table=[img_toshow_index; zeros(1,3)];
            bet_money_table=[img_toshow_index; zeros(1,3)];
            Barycentric_coord_table=[img_toshow_index; zeros(1,3)];
            
            current_cursor_selection_ind_default=find(tri_pt_arrange_tbl==12);
            current_cursor_selection_ind=current_cursor_selection_ind_default;
            
            ind_ev=0;
            while(~Is_bet)
                
                ind_ev=ind_ev+1;
                % message
                str_block_intro=sprintf('*** Submit your rating for keeping pictures. ***');
                DrawFormattedText(wPtr, str_block_intro, 'center', round(screenHeight/8) ,[0, 0, 255, 255]);
                str_end=sprintf('[y]:submit, [n]:reset');
                DrawFormattedText(wPtr, str_end, 'center', ypos+size_triangle*cos(pi/3)+sy+50);
                
                % 2. Show triangle
                Screen('DrawLine', wPtr, [0, 0, 0, 255], tri_pt(1,1), tri_pt(1,2), tri_pt(2,1), tri_pt(2,2),[3]);
                Screen('DrawLine', wPtr, [0, 0, 0, 255], tri_pt(2,1), tri_pt(2,2), tri_pt(3,1), tri_pt(3,2),[3]);
                Screen('DrawLine', wPtr, [0, 0, 0, 255], tri_pt(3,1), tri_pt(3,2), tri_pt(1,1), tri_pt(1,2),[3]);
                Screen('DrawLine', wPtr, [0, 0, 0, 100], xpos, ypos, tri_pt(1,1), tri_pt(1,2),[1]);
                Screen('DrawLine', wPtr, [0, 0, 0, 100], xpos, ypos, tri_pt(2,1), tri_pt(2,2),[1]);
                Screen('DrawLine', wPtr, [0, 0, 0, 100], xpos, ypos, tri_pt(3,1), tri_pt(3,2),[1]);
                % show minor guideline in the triangle
                Screen('DrawLine', wPtr, [100, 100, 100, 100], tri_pt_all(12,1), tri_pt_all(12,2), tri_pt_all(6,1), tri_pt_all(6,2),[1]);
                Screen('DrawLine', wPtr, [100, 100, 100, 100], tri_pt_all(12,1), tri_pt_all(12,2), tri_pt_all(8,1), tri_pt_all(8,2),[1]);
                Screen('DrawLine', wPtr, [100, 100, 100, 100], tri_pt_all(12,1), tri_pt_all(12,2), tri_pt_all(23,1), tri_pt_all(23,2),[1]);
                % show all points of the triangle
                for ttt=1:1:size(tri_pt_all,1)
                    DrawFormattedText(wPtr, '+', tri_pt_all(ttt,1)-8, tri_pt_all(ttt,2)-14, [100 100 100 255]); % add '+' mark at the click pt.
                end
                % show images
                for j=1:1:3
                    selected=img_toshow_index(j);
                    input_stim = Screen('MakeTexture', wPtr, img_set_all{block,selected});
                    sx=floor(IMAGE_SIZE(1)*disp_scale_scoring);       sy=floor(IMAGE_SIZE(2)*disp_scale_scoring);
                    if(j==1)
                        xpos_img = tri_pt(j,1);    ypos_img = tri_pt(j,2)-sy/2;
                    end
                    if(j==2)
                        xpos_img = tri_pt(j,1)-sx/2;    ypos_img = tri_pt(j,2)+sy/2;
                    end
                    if(j==3)
                        xpos_img = tri_pt(j,1)+sx/2;    ypos_img = tri_pt(j,2)+sy/2;
                    end
                    destrect=[xpos_img-sx/2,ypos_img-sy/2,xpos_img+sx/2,ypos_img+sy/2];
                    Screen('DrawTexture', wPtr, input_stim,[],destrect);
                end
                
                x=tri_pt_all(tri_pt_arrange_tbl(current_cursor_selection_ind),1);   y=tri_pt_all(tri_pt_arrange_tbl(current_cursor_selection_ind),2);
                DrawFormattedText(wPtr, 'o', x-8,y-14, [0 0 255 255]); % add 'o' mark at the click pt.
                
                % 3-1. Barycentric_coordinate_system
                denorm=(tri_pt(2,2)-tri_pt(3,2))*(tri_pt(1,1)-tri_pt(3,1))+(tri_pt(3,1)-tri_pt(2,1))*(tri_pt(1,2)-tri_pt(3,2));
                Barycentric_coord_table(2,1)=((tri_pt(2,2)-tri_pt(3,2))*(x-tri_pt(3,1))+(tri_pt(3,1)-tri_pt(2,1))*(y-tri_pt(3,2)))/denorm;
                Barycentric_coord_table(2,1)=max(Barycentric_coord_table(2,1),0);   Barycentric_coord_table(2,1)=min(Barycentric_coord_table(2,1),1);
                Barycentric_coord_table(2,2)=((tri_pt(3,2)-tri_pt(1,2))*(x-tri_pt(3,1))+(tri_pt(1,1)-tri_pt(3,1))*(y-tri_pt(3,2)))/denorm;
                Barycentric_coord_table(2,2)=max(Barycentric_coord_table(2,2),0);   Barycentric_coord_table(2,2)=min(Barycentric_coord_table(2,2),1);
                Barycentric_coord_table(2,3)=1-Barycentric_coord_table(2,1)-Barycentric_coord_table(2,2);
                Barycentric_coord_table(2,:)=abs(Barycentric_coord_table(2,:));
                
                % 4-3. compute and add betting values from distance table
                bet_score_table(2,:)=Barycentric_coord_table(2,:);
                for j=1:1:3 % add the betting values to display
                    sx=floor(IMAGE_SIZE(1)*disp_scale_scoring);       sy=floor(IMAGE_SIZE(2)*disp_scale_scoring);
                    if(j==1)
                        xpos_score_disp = tri_pt(j,1)-40;    ypos_score_disp = tri_pt(j,2)-sy-35;
                    end
                    if(j==2)
                        xpos_score_disp = tri_pt(j,1)-sx;    ypos_score_disp = tri_pt(j,2)+sy+10;
                    end
                    if(j==3)
                        xpos_score_disp = tri_pt(j,1);    ypos_score_disp = tri_pt(j,2)+sy+10;
                    end
                    destrect=[xpos_score_disp-sx/2,ypos_score_disp-sy/2,xpos_score_disp+sx/2,ypos_score_disp+sy/2];
                    bet_money_table(2,j)=double(money_given*bet_score_table(2,j));
                    str_score=sprintf('%01.2f',bet_money_table(2,j));
                    DrawFormattedText(wPtr, str_score, xpos_score_disp, ypos_score_disp);
                end
                
                % display!
                Screen(wPtr, 'Flip');
                if(ind_ev==1)
                    HIST_event_info=[HIST_event_info [block; 0; 0; (GetSecs-session_clock_start); (GetSecs-block_clock_start); 9]]; % event save
                    clock_time_limit_start=GetSecs;
                end
                % take a snapshot
                if(DO_TAKE_SNAPSHOT==1)
                    snapshot=Screen(wPtr, 'GetImage', [1, 1, floor(1.0*SCREEN_RESOLUTION(1)), floor(1.0*SCREEN_RESOLUTION(2))]);
                    imageArray=[imageArray; {snapshot}];
                end
                
                % 1. Get cursor (L:37 R:39)
                decision_made=0;
                while(~decision_made)
                    [secs, keyCode] = KbPressWait([], clock_time_limit_start+sec_limit_Bayesian_rating); % if no keyboard in time limit, then go ahead. if pressed earlier, then go ahead.
                    [tmp tmp_key_code]=find(keyCode==1);                zzz=[zzz tmp_key_code];
                    if(tmp_key_code==KEY_L) % L pressed
                        current_cursor_selection_ind=max(1,current_cursor_selection_ind-1);
                        decision_made=1;
                    end
                    if(tmp_key_code==KEY_R) % R pressed
                        current_cursor_selection_ind=min(size(tri_pt_all,1),current_cursor_selection_ind+1);
                        decision_made=1;
                    end
                    if(tmp_key_code==KEY_Y) % 'y' pressed
                        decision_made=1;            Is_bet=1;
                        bet_score_table(2,:)=Barycentric_coord_table(2,:);
                        HIST_event_info=[HIST_event_info [block; 0; 0; (GetSecs-session_clock_start); (GetSecs-block_clock_start); 10]]; % event save
                    end
                    if(tmp_key_code==KEY_N) % 'n' pressed - reset
                        decision_made=1;            Is_bet=0;
                        current_cursor_selection_ind=current_cursor_selection_ind_default;                        
                    end
                    if(tmp_key_code==KEY_Q) % 'q' pressed for aborting
                        clear mex
                        decision_made=1;
                    end
                    % check the time limit !@#$
                    if((GetSecs-clock_time_limit_start)>sec_limit_Bayesian_rating)
                        decision_made=1;            Is_bet=1;
                        dont_know_cursor_selection=(length(Q{1,q_ind}.cursor_input)+1)/2;
                        bet_score_table(2,:)=Barycentric_coord_table(2,:);%[1/3 1/3 1/3];
                        bet_score_table(3,:)=-99; % error code
                        HIST_event_info=[HIST_event_info [block; 0; 0; (GetSecs-session_clock_start); (GetSecs-block_clock_start); -99]]; % event save                        
                    end
                    
                end
                
                
            end
            
            % a short blank page : to give a sense of the page transition
            Screen(wPtr, 'Flip');
            HIST_event_info=[HIST_event_info [block; 0; 0; (GetSecs-session_clock_start); (GetSecs-block_clock_start); 20]]; % event save
            WaitSecs(sec_jittered_blank_page);
            
            % save
            HIST_Barycentric_coord_table_type1{block,kk}=Barycentric_coord_table; % sum to one
            HIST_bet_score_table_type1{block,kk}=bet_score_table; % same as the above
            HIST_bet_money_table_type1{block,kk}=bet_money_table; % moeny bet amount
            
        end

    end

    % take a snapshot
    if(DO_TAKE_SNAPSHOT==1)
        snapshot=Screen(wPtr, 'GetImage', [1, 1, floor(1.0*SCREEN_RESOLUTION(1)), floor(1.0*SCREEN_RESOLUTION(2))]);
        imageArray=[imageArray; {snapshot}];
    end


    %     % block ends
    %     str_block_intro=sprintf('*** You have finished block %d. Press any key to continue. ***',block);
    %     DrawFormattedText(wPtr, str_block_intro, 'center', 'center',[0, 0, 0, 255]);
    %     Screen(wPtr, 'Flip');
    %     WaitSecs(0.3);
    %     KbWait;
    %     % take a snapshot
    %     if(DO_TAKE_SNAPSHOT==1)
    %         snapshot=Screen(wPtr, 'GetImage', [1, 1, floor(1.0*SCREEN_RESOLUTION(1)), floor(1.0*SCREEN_RESOLUTION(2))]);
    %         imageArray=[imageArray; {snapshot}];
    %     end


    %% Update Tot_img_info matrix (the ratings)
    if(length(find(nov_pos_block_index==block))==1)
        block_kind=0; % nov-pos
    else
        block_kind=1; % nov-neg
    end
    [r_tmp c_tmp]=find(Tot_img_info(3,:)==index_num);
    [r_tmp2 c_tmp2]=find(Tot_img_info(4,c_tmp)==block);
    absolute_index=c_tmp(c_tmp2);
    for hhh=1:1:length(absolute_index)
        novelty_lev=Tot_img_info(5,absolute_index(hhh)); % read out novelty level
        % update hedonic ratings
        [r_tmp3 c_tmp3]=find(HIST_bet_score_table_type2{block,1}(1,:)==novelty_lev);        
        Tot_img_info(6,absolute_index(hhh))=HIST_bet_score_table_type2{block,1}(2,c_tmp3);        
        % update causal ratings
        [r_tmp3 c_tmp3]=find(HIST_bet_score_table_type2{block,2}(1,:)==novelty_lev);        
        Tot_img_info(7,absolute_index(hhh))=HIST_bet_score_table_type2{block,2}(2,c_tmp3);                        
        % update bonus round ratings
        [r_tmp3 c_tmp3]=find(HIST_bet_score_table_type1{block,1}(1,:)==novelty_lev);        
        Tot_img_info(8,absolute_index(hhh))=HIST_bet_score_table_type1{block,1}(2,c_tmp3);       
        % update nov-pos=0, nov-neg=1
        Tot_img_info(9,absolute_index(hhh))=block_kind;
    end
    %% Update HIST_event_info_all (overwrite at each block)
    HIST_event_info_all{1,index_num}=HIST_event_info;

    
    %% save the (updated) image usage matrix (overwriting)
    file_imgind_sv_name=[EXP_NAME '_image_usage_info.mat'];
    file_name_sv=[save_path file_imgind_sv_name];
    % eval(['save ' file_name_sv ' Tot_img_info Tot_img_info_Tag']);
    save(file_name_sv,'Tot_img_info','Tot_img_info_Tag','HIST_event_info_all','HIST_event_type');
    
end













%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Bonus block %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if(0) % dont run bonus block here


    if(KEEP_PICTURE_OPT==1)
        block=1;
        img_set_all_bonus=cell(1,tot_num_img_per_block);
        %% image set read based on the choices in previous blocks

        % random-picking of 3 blocks
        block_pick=randperm(Tot_block);     block_pick=block_pick(1:3);
        current_pool_bonus=zeros(1,3);

        for hhh=1:1:3 % obtain a cue from each block based on prob. of keeping
            prob_tbl=HIST_bet_score_table_type1{block_pick(hhh),1}(2,:);
            prob_tbl=prob_tbl/sum(prob_tbl);
            cumnum=cumsum(prob_tbl);
            [temp1 selected]=histc(rand,[0 cumnum]);
            selected_file_index=HIST_bet_score_table_type1{block_pick(hhh),1}(1,selected);
            current_pool_bonus(1,hhh)=HIST_current_pool(block_pick(hhh),selected_file_index);
        end

        img_index_use_all=randsample(current_pool_bonus,tot_num_img_per_block,false);
        HIST_current_pool_bonus=img_index_use_all;
        Tot_img_info(2,img_index_use_all)=1; % tagged as "used"
        % assign normal&novel stimulus
        img_index_use=img_index_use_all(1:end-1);
        img0_index_use=img_index_use_all(end);

        num_stim=length(img_index_use);
        for i=1:1:num_stim
            file_full_path=[seed_path sprintf('%03d.png',img_index_use(i))];
            img_set_bonus{block,i}=imresize(imread(file_full_path),[IMAGE_SIZE(2) IMAGE_SIZE(1)]); % get a frame
        end
        num_stim0=length(img0_index_use);
        for i=1:1:num_stim0
            file_full_path=[seed_path sprintf('%03d.png',img0_index_use(i))];
            img0_set_bonus{block,i}=imresize(imread(file_full_path),[IMAGE_SIZE(2) IMAGE_SIZE(1)]); % get a frame
        end




        % I. Simuli presentation

        zzz=[];

        % block starts
        Screen('FillRect',wPtr,BackgroundColor_block_intro);
        str_block_intro=sprintf('*** Now Bonus block starts. Press any key to continue. ***');
        DrawFormattedText(wPtr, str_block_intro, 'center', 'center',[0, 0, 0, 255]);
        Screen(wPtr, 'Flip');
        WaitSecs(0.3);
        KbWait;

        % take a snapshot
        if(DO_TAKE_SNAPSHOT==1)
            snapshot=Screen(wPtr, 'GetImage', [1, 1, floor(1.0*SCREEN_RESOLUTION(1)), floor(1.0*SCREEN_RESOLUTION(2))]);
            imageArray=[imageArray; {snapshot}];
        end

        img0_index_showed=zeros(1,length(img0_index_use));


        for trial=1:1:Tot_trial % each trial

            %         novel_trial_s_ind=max(1,ceil((Tot_trial_s*how_early_novel_stimulus_show_in_trial_s)*rand));

            Screen('FillRect',wPtr,BackgroundColor_Trial_ready_page);

            % ready sign
            %         str_ready=sprintf('Ready for the %d-th trial',trial);
            %         DrawFormattedText(wPtr, str_ready, 'center', 'center',[0, 0, 0, 255]);
            %         Screen(wPtr, 'Flip');
            %         WaitSecs(sec_trial_ready);
            % take a snapshot
            if(DO_TAKE_SNAPSHOT==1)
                snapshot=Screen(wPtr, 'GetImage', [1, 1, floor(1.0*SCREEN_RESOLUTION(1)), floor(1.0*SCREEN_RESOLUTION(2))]);
                imageArray=[imageArray; {snapshot}];
            end


            for trial_s=1:1:Tot_trial_s % cue presentation

                %             % ready sign
                %             str_ready=sprintf('(%d/%d) cue ready',trial_s,Tot_trial_s);
                %             DrawFormattedText(wPtr, str_ready, 'center', 'center',[0, 0, 0, 255]);
                %             Screen(wPtr, 'Flip');
                %             WaitSecs(sec_stim_ready);

                % take a snapshot
                if(DO_TAKE_SNAPSHOT==1)
                    snapshot=Screen(wPtr, 'GetImage', [1, 1, floor(1.0*SCREEN_RESOLUTION(1)), floor(1.0*SCREEN_RESOLUTION(2))]);
                    imageArray=[imageArray; {snapshot}];
                end

                % Determine which stimlus will be presented.
                selected=HIST_schedule_bonus{1,block}(trial,trial_s);
                if(selected==tot_num_img_per_block) % the most novel cue
                    input_stim = Screen('MakeTexture', wPtr, img0_set_bonus{block,1});
                    img0_index_showed(1)=1;
                else
                    input_stim = Screen('MakeTexture', wPtr, img_set_bonus{block,selected});
                end
                %             if((trial==novel_trial_ind)&&(trial_s==novel_trial_s_ind))
                %                 % novel trial
                %                 novel_trial_img_ind=ceil(length(img0_index_use)*rand);
                %                 selected=novel_trial_img_ind;
                %                 input_stim = Screen('MakeTexture', wPtr, img0_set{block,selected});
                %                 img0_index_showed(novel_trial_img_ind)=1;
                %             else
                %                 % normal trial
                %                 cumnum=cumsum(img_present_freq);
                %                 [temp1 selected]=histc(rand,[0 cumnum]);
                %                 input_stim = Screen('MakeTexture', wPtr, img_set{block,selected});
                %             end

                % image display
                xpos = round(screenWidth/2);    ypos = round(screenHeight/2);
                sx=floor(IMAGE_SIZE(1)*disp_scale);       sy=floor(IMAGE_SIZE(2)*disp_scale);
                destrect=[xpos-sx/2,ypos-sy/2,xpos+sx/2,ypos+sy/2];
                Screen('FillRect',wPtr,BackgroundColor_Cue_page);
                Screen('DrawTexture', wPtr, input_stim,[],destrect);
                %             % cue number display
                %             str_ready=sprintf('Cue (%d/%d)',trial_s,Tot_trial_s);
                %             DrawFormattedText(wPtr, str_ready, 'center', ypos+sy/2+20,[0, 0, 0, 255]);
                Screen(wPtr, 'Flip');
                WaitSecs(sec_stim_display);
                % take a snapshot
                if(DO_TAKE_SNAPSHOT==1)
                    snapshot=Screen(wPtr, 'GetImage', [1, 1, floor(1.0*SCREEN_RESOLUTION(1)), floor(1.0*SCREEN_RESOLUTION(2))]);
                    imageArray=[imageArray; {snapshot}];
                end
                Screen(wPtr, 'Flip');
                sec_stim_interval0=rand*(max(sec_stim_interval)-min(sec_stim_interval))+min(sec_stim_interval);
                WaitSecs(sec_stim_interval0);
                % take a snapshot
                if(DO_TAKE_SNAPSHOT==1)
                    snapshot=Screen(wPtr, 'GetImage', [1, 1, floor(1.0*SCREEN_RESOLUTION(1)), floor(1.0*SCREEN_RESOLUTION(2))]);
                    imageArray=[imageArray; {snapshot}];
                end

            end

            % reward display
            %         Screen('TextSize',wPtr, text_size_reward);
            %         str_rwd=sprintf('You have earned %d points.',HIST_reward_bonus(1,trial));
            %         DrawFormattedText(wPtr, str_rwd, 'center', 'center',[255, 0, 0, 255]);
            Screen('FillRect',wPtr,BackgroundColor_Reward_page);
            if(HIST_reward_bonus(1,trial)<0)
                input_stim_msg = Screen('MakeTexture', wPtr, img_msg_set{1,1});
                sx_msg=size(img_msg_set{1,1},2);    sy_msg=size(img_msg_set{1,1},1);
            else
                input_stim_msg = Screen('MakeTexture', wPtr, img_msg_set{2,1});
                sx_msg=size(img_msg_set{2,1},2);    sy_msg=size(img_msg_set{2,1},1);
            end
            destrect=[xpos-sx_msg/2,ypos-sy_msg/2,xpos+sx_msg/2,ypos+sy_msg/2];
            Screen('DrawTexture', wPtr, input_stim_msg,[],destrect);
            str_rwd2=sprintf('### Press any button to continue. ###');
            DrawFormattedText(wPtr, str_rwd2, 'center', round(screenHeight*4/5),[0, 0, 0, 255]);
            Screen(wPtr, 'Flip');

            %         Screen('TextSize',wPtr, text_size_default);

            % take a snapshot
            if(DO_TAKE_SNAPSHOT==1)
                snapshot=Screen(wPtr, 'GetImage', [1, 1, floor(1.0*SCREEN_RESOLUTION(1)), floor(1.0*SCREEN_RESOLUTION(2))]);
                imageArray=[imageArray; {snapshot}];
            end
            % button press - Decision-making
            % 'q':81, esc: 27
            decision_made=0;
            exit_made=0;
            while(~decision_made)
                [secs, keyCode] = KbPressWait;
                [tmp tmp_key_code]=find(keyCode==1);                zzz=[zzz tmp_key_code];
                if(tmp_key_code==KEY_Q) % 'q' pressed for aborting
                    clear mex
                    decision_made=1;
                end
                decision_made=1;
            end
            WaitSecs(0.1);
            % take a snapshot
            if(DO_TAKE_SNAPSHOT==1)
                snapshot=Screen(wPtr, 'GetImage', [1, 1, floor(1.0*SCREEN_RESOLUTION(1)), floor(1.0*SCREEN_RESOLUTION(2))]);
                imageArray=[imageArray; {snapshot}];
            end

        end


        % collect all the images used
        for h=1:1:length(img_index_use)
            img_set_all_bonus{block,h}=img_set_bonus{block,h};
        end
        [tmp ind_new]=find(img0_index_showed==1);
        for h=1:1:sum(img0_index_showed)
            img_set_all_bonus{block,h+length(img_index_use)}=img0_set_bonus{block,ind_new(h)};
        end





        %% III. Rating - type2

        img_toshow_index=randperm(tot_num_img_per_block);
        xpos = round(screenWidth/2);    ypos = round(screenHeight/2);
        length_bar=round(screenWidth*3/5);
        bar_pt_start=[round(screenWidth*1/5),ypos+250];
        bar_pt_end=[round(screenWidth*1/5)+length_bar,ypos+250];
        bet_score_table=[img_toshow_index; zeros(1,3)];
        bet_money_table=[img_toshow_index; zeros(1,3)];

        %     Screen('TextSize',wPtr, text_size_default);

        for q_ind=1:1:size(Q,2)

            ggg_ind=0;
            for m=1:1:tot_num_img_per_block
                ggg_ind=ggg_ind+1;



                Is_bet=0;
                bar_pt_cursor=(Q{1,q_ind}.cursor_input-min(Q{1,q_ind}.cursor_input))/(max(Q{1,q_ind}.cursor_input)-min(Q{1,q_ind}.cursor_input)); %[0,1]
                bar_pt_cursor=length_bar*bar_pt_cursor+bar_pt_start(1); %actual x-positions in display
                current_cursor_selection_ind=round(length(Q{1,q_ind}.cursor_input)/2);

                while(~Is_bet)

                    % message
                    str_block_intro=Q{1,q_ind}.text;
                    DrawFormattedText(wPtr, str_block_intro, 'center', round(screenHeight/7) ,Q{1,q_ind}.color);
                    % 2. Show Bars and images
                    sx=floor(IMAGE_SIZE(1)*disp_scale);       sy=floor(IMAGE_SIZE(2)*disp_scale);
                    % draw bar
                    Screen('DrawLine', wPtr, Q{1,q_ind}.color, bar_pt_start(1), bar_pt_start(2), bar_pt_end(1), bar_pt_end(2),[3]);
                    for jj=1:1:length(bar_pt_cursor)
                        Screen('DrawLine', wPtr, Q{1,q_ind}.color, bar_pt_cursor(jj), bar_pt_start(2)-5, bar_pt_cursor(jj), bar_pt_start(2)+5,[1]);
                        str_block_num=sprintf('%s',num2str(Q{1,q_ind}.cursor_input(jj)));
                        DrawFormattedText(wPtr, str_block_num, bar_pt_cursor(jj)-10, bar_pt_start(2)+25, Q{1,q_ind}.color);
                        if(jj==1) % like/dislike, notatall/verylikeli
                            DrawFormattedText(wPtr, Q{1,q_ind}.text2, bar_pt_cursor(jj)-50, bar_pt_start(2)+50, Q{1,q_ind}.color);
                            Screen('DrawLine', wPtr, Q{1,q_ind}.color, bar_pt_cursor(jj), bar_pt_start(2)-8, bar_pt_cursor(jj), bar_pt_start(2)+8,[3]);
                        end
                        if(jj==length(bar_pt_cursor)) % like/dislike, notatall/verylikeli
                            DrawFormattedText(wPtr, Q{1,q_ind}.text3, bar_pt_cursor(jj)-50, bar_pt_start(2)+50, Q{1,q_ind}.color);
                            Screen('DrawLine', wPtr, Q{1,q_ind}.color, bar_pt_cursor(jj), bar_pt_start(2)-8, bar_pt_cursor(jj), bar_pt_start(2)+8,[3]);
                        end
                        if(jj==round(length(Q{1,q_ind}.cursor_input)/2)) % like/dislike, notatall/verylikeli
                            DrawFormattedText(wPtr, Q{1,q_ind}.text4, bar_pt_cursor(jj)-80, bar_pt_start(2)+50, Q{1,q_ind}.color);
                            Screen('DrawLine', wPtr, Q{1,q_ind}.color, bar_pt_cursor(jj), bar_pt_start(2)-8, bar_pt_cursor(jj), bar_pt_start(2)+8,[3]);
                        end
                    end
                    % draw images
                    input_stim = Screen('MakeTexture', wPtr, img_set_all_bonus{block,img_toshow_index(m)});
                    img_pt_center=[xpos, bar_pt_start(2)-sy/2-100];
                    destrect=[img_pt_center(1)-sx/2,img_pt_center(2)-sy/2,img_pt_center(1)+sx/2,img_pt_center(2)+sy/2];
                    Screen('DrawTexture', wPtr, input_stim,[],destrect);
                    % 4-2. add 'o' mark at the click pt.
                    DrawFormattedText(wPtr, 'o', bar_pt_cursor(current_cursor_selection_ind)-8, bar_pt_start(2)-14, [0 0 0 255]); % add 'o' mark at the click pt.
                    % display!
                    Screen(wPtr, 'Flip');
                    % take a snapshot
                    if(DO_TAKE_SNAPSHOT==1)
                        snapshot=Screen(wPtr, 'GetImage', [1, 1, floor(1.0*SCREEN_RESOLUTION(1)), floor(1.0*SCREEN_RESOLUTION(2))]);
                        imageArray=[imageArray; {snapshot}];
                    end

                    % 1. Get cursor (L:37 R:39)
                    decision_made=0;
                    while(~decision_made)
                        [secs, keyCode] = KbPressWait;
                        [tmp tmp_key_code]=find(keyCode==1);                zzz=[zzz tmp_key_code];
                        if(tmp_key_code==KEY_L) % L pressed
                            current_cursor_selection_ind=max(1,current_cursor_selection_ind-1);
                            decision_made=1;
                        end
                        if(tmp_key_code==KEY_R) % R pressed
                            current_cursor_selection_ind=min(length(Q{1,q_ind}.cursor_input),current_cursor_selection_ind+1);
                            decision_made=1;
                        end
                        if(tmp_key_code==KEY_Y) % 'y' pressed
                            decision_made=1;            Is_bet=1;
                            bet_score_table(2,ggg_ind)=Q{1,q_ind}.cursor_input(current_cursor_selection_ind);
                        end
                    end

                    % take a snapshot
                    if(DO_TAKE_SNAPSHOT==1)
                        snapshot=Screen(wPtr, 'GetImage', [1, 1, floor(1.0*SCREEN_RESOLUTION(1)), floor(1.0*SCREEN_RESOLUTION(2))]);
                        imageArray=[imageArray; {snapshot}];
                    end

                end

            end
            HIST_bet_score_table_type2_bonus{block,q_ind}=bet_score_table; % same as the above


        end





        if(KEEP_PICTURE_OPT==1)
            %% II. Rating - type1

            block=1;
            % compute points of a triangle

            tri_pt=zeros(3,2); % edges
            size_triangle=200;
            distance_limit=size_triangle+size_triangle*cos(pi/3);
            xpos = round(screenWidth/2);    ypos = round(screenHeight/2)+70;
            tri_pt(1,:)=[xpos,ypos-size_triangle];
            tri_pt(2,:)=[xpos-size_triangle*sin(pi/3),ypos+size_triangle*cos(pi/3)];
            tri_pt(3,:)=[xpos+size_triangle*sin(pi/3),ypos+size_triangle*cos(pi/3)];

            tri_pt_all=zeros(25,2); % all inner points including edges
            tri_pt_all(1,:)=tri_pt(1,:);    tri_pt_all(21,:)=tri_pt(2,:);    tri_pt_all(25,:)=tri_pt(3,:);  tri_pt_all(12,:)=[xpos ypos];
            tri_pt_all(4,:)=[tri_pt_all(1,1) tri_pt_all(1,2)+size_triangle/2];
            tri_pt_all(7,:)=[tri_pt_all(1,1) tri_pt_all(1,2)+3*size_triangle/4];
            tri_pt_all(17,:)=[tri_pt_all(12,1) tri_pt_all(12,2)+size_triangle/4];
            tri_pt_all(2,:)=[tri_pt_all(1,1) tri_pt_all(1,2)+size_triangle/4];
            th_ang=pi/3-atan(cos(pi/6));
            tri_pt_all(3,:)=[tri_pt_all(12,1)-sqrt(7)/4*size_triangle*sin(th_ang) tri_pt_all(12,2)-sqrt(7)/4*size_triangle*cos(th_ang)];
            tri_pt_all(5,:)=[tri_pt_all(12,1)+sqrt(7)/4*size_triangle*sin(th_ang) tri_pt_all(12,2)-sqrt(7)/4*size_triangle*cos(th_ang)];
            tri_pt_all(6,:)=[tri_pt_all(7,1)-size_triangle*tan(pi/3)/4 tri_pt_all(7,2)];
            tri_pt_all(8,:)=[tri_pt_all(7,1)+size_triangle*tan(pi/3)/4 tri_pt_all(7,2)];
            tri_pt_all(9,:)=[tri_pt_all(12,1)-size_triangle*cos(pi/6)/4 tri_pt_all(12,2)-size_triangle*sin(pi/6)/4];
            tri_pt_all(10,:)=[tri_pt_all(12,1)+size_triangle*cos(pi/6)/4 tri_pt_all(12,2)-size_triangle*sin(pi/6)/4];
            th_ang=atan(cos(pi/6))-pi/6;
            tri_pt_all(11,:)=[tri_pt_all(12,1)-sqrt(7)/4*size_triangle*cos(th_ang) tri_pt_all(12,2)+sqrt(7)/4*size_triangle*sin(th_ang)];
            tri_pt_all(13,:)=[tri_pt_all(12,1)+sqrt(7)/4*size_triangle*cos(th_ang) tri_pt_all(12,2)+sqrt(7)/4*size_triangle*sin(th_ang)];
            tri_pt_all(14,:)=[tri_pt_all(12,1)-size_triangle*cos(pi/6)/4 tri_pt_all(12,2)+size_triangle*sin(pi/6)/4];
            tri_pt_all(15,:)=[tri_pt_all(12,1)+size_triangle*cos(pi/6)/4 tri_pt_all(12,2)+size_triangle*sin(pi/6)/4];
            tri_pt_all(16,:)=[tri_pt_all(17,1)-size_triangle*cos(pi/6)/2 tri_pt_all(17,2)];
            tri_pt_all(18,:)=[tri_pt_all(17,1)+size_triangle*cos(pi/6)/2 tri_pt_all(17,2)];
            tri_pt_all(19,:)=[tri_pt_all(12,1)-size_triangle*3*cos(pi/6)/4 tri_pt_all(12,2)+size_triangle*3/8];
            tri_pt_all(20,:)=[tri_pt_all(12,1)+size_triangle*3*cos(pi/6)/4 tri_pt_all(12,2)+size_triangle*3/8];
            tri_pt_all(22,:)=[tri_pt_all(21,1)+size_triangle*cos(pi/6)/2 tri_pt_all(21,2)];
            tri_pt_all(23,:)=[tri_pt_all(12,1) tri_pt_all(12,2)+size_triangle/2];
            tri_pt_all(24,:)=[tri_pt_all(25,1)-size_triangle*cos(pi/6)/2 tri_pt_all(25,2)];
            tri_pt_all=round(tri_pt_all);
            tri_pt_arrange_tbl=[1 2 3 5 4 6 7 8 9 10 12 11 14 15 13 16 17 18 19 20 21:1:25];


            for kk=1:1:size(list_combination,1)

                img_toshow_index=list_combination(kk,randperm(3));


                Is_bet=0;
                bet_score_table=[img_toshow_index; zeros(1,3)];
                bet_money_table=[img_toshow_index; zeros(1,3)];
                Barycentric_coord_table=[img_toshow_index; zeros(1,3)];

                current_cursor_selection_ind_default=find(tri_pt_arrange_tbl==12);
                current_cursor_selection_ind=current_cursor_selection_ind_default;

                while(~Is_bet)

                    % message
                    str_block_intro=sprintf('*** Submit your rating for keeping pictures. ***');
                    DrawFormattedText(wPtr, str_block_intro, 'center', round(screenHeight/8) ,[0, 0, 255, 255]);
                    str_end=sprintf('[y]:submit, [n]:reset');
                    DrawFormattedText(wPtr, str_end, 'center', ypos+size_triangle*cos(pi/3)+sy+50);

                    % 2. Show triangle
                    Screen('DrawLine', wPtr, [0, 0, 0, 255], tri_pt(1,1), tri_pt(1,2), tri_pt(2,1), tri_pt(2,2),[3]);
                    Screen('DrawLine', wPtr, [0, 0, 0, 255], tri_pt(2,1), tri_pt(2,2), tri_pt(3,1), tri_pt(3,2),[3]);
                    Screen('DrawLine', wPtr, [0, 0, 0, 255], tri_pt(3,1), tri_pt(3,2), tri_pt(1,1), tri_pt(1,2),[3]);
                    Screen('DrawLine', wPtr, [0, 0, 0, 100], xpos, ypos, tri_pt(1,1), tri_pt(1,2),[1]);
                    Screen('DrawLine', wPtr, [0, 0, 0, 100], xpos, ypos, tri_pt(2,1), tri_pt(2,2),[1]);
                    Screen('DrawLine', wPtr, [0, 0, 0, 100], xpos, ypos, tri_pt(3,1), tri_pt(3,2),[1]);
                    % show minor guideline in the triangle
                    Screen('DrawLine', wPtr, [100, 100, 100, 100], tri_pt_all(12,1), tri_pt_all(12,2), tri_pt_all(6,1), tri_pt_all(6,2),[1]);
                    Screen('DrawLine', wPtr, [100, 100, 100, 100], tri_pt_all(12,1), tri_pt_all(12,2), tri_pt_all(8,1), tri_pt_all(8,2),[1]);
                    Screen('DrawLine', wPtr, [100, 100, 100, 100], tri_pt_all(12,1), tri_pt_all(12,2), tri_pt_all(23,1), tri_pt_all(23,2),[1]);
                    % show all points of the triangle
                    for ttt=1:1:size(tri_pt_all,1)
                        DrawFormattedText(wPtr, '+', tri_pt_all(ttt,1)-8, tri_pt_all(ttt,2)-14, [100 100 100 255]); % add '+' mark at the click pt.
                    end
                    % show images
                    for j=1:1:3
                        selected=img_toshow_index(j);
                        input_stim = Screen('MakeTexture', wPtr, img_set_all_bonus{block,selected});
                        sx=floor(IMAGE_SIZE(1)*disp_scale_scoring);       sy=floor(IMAGE_SIZE(2)*disp_scale_scoring);
                        if(j==1)
                            xpos_img = tri_pt(j,1);    ypos_img = tri_pt(j,2)-sy/2;
                        end
                        if(j==2)
                            xpos_img = tri_pt(j,1)-sx/2;    ypos_img = tri_pt(j,2)+sy/2;
                        end
                        if(j==3)
                            xpos_img = tri_pt(j,1)+sx/2;    ypos_img = tri_pt(j,2)+sy/2;
                        end
                        destrect=[xpos_img-sx/2,ypos_img-sy/2,xpos_img+sx/2,ypos_img+sy/2];
                        Screen('DrawTexture', wPtr, input_stim,[],destrect);
                    end

                    x=tri_pt_all(tri_pt_arrange_tbl(current_cursor_selection_ind),1);   y=tri_pt_all(tri_pt_arrange_tbl(current_cursor_selection_ind),2);
                    DrawFormattedText(wPtr, 'o', x-8,y-14, [0 0 255 255]); % add 'o' mark at the click pt.

                    % 3-1. Barycentric_coordinate_system
                    denorm=(tri_pt(2,2)-tri_pt(3,2))*(tri_pt(1,1)-tri_pt(3,1))+(tri_pt(3,1)-tri_pt(2,1))*(tri_pt(1,2)-tri_pt(3,2));
                    Barycentric_coord_table(2,1)=((tri_pt(2,2)-tri_pt(3,2))*(x-tri_pt(3,1))+(tri_pt(3,1)-tri_pt(2,1))*(y-tri_pt(3,2)))/denorm;
                    Barycentric_coord_table(2,1)=max(Barycentric_coord_table(2,1),0);   Barycentric_coord_table(2,1)=min(Barycentric_coord_table(2,1),1);
                    Barycentric_coord_table(2,2)=((tri_pt(3,2)-tri_pt(1,2))*(x-tri_pt(3,1))+(tri_pt(1,1)-tri_pt(3,1))*(y-tri_pt(3,2)))/denorm;
                    Barycentric_coord_table(2,2)=max(Barycentric_coord_table(2,2),0);   Barycentric_coord_table(2,2)=min(Barycentric_coord_table(2,2),1);
                    Barycentric_coord_table(2,3)=1-Barycentric_coord_table(2,1)-Barycentric_coord_table(2,2);
                    Barycentric_coord_table(2,:)=abs(Barycentric_coord_table(2,:));

                    % 4-3. compute and add betting values from distance table
                    bet_score_table(2,:)=Barycentric_coord_table(2,:);
                    for j=1:1:3 % add the betting values to display
                        sx=floor(IMAGE_SIZE(1)*disp_scale_scoring);       sy=floor(IMAGE_SIZE(2)*disp_scale_scoring);
                        if(j==1)
                            xpos_score_disp = tri_pt(j,1)-40;    ypos_score_disp = tri_pt(j,2)-sy-35;
                        end
                        if(j==2)
                            xpos_score_disp = tri_pt(j,1)-sx;    ypos_score_disp = tri_pt(j,2)+sy+10;
                        end
                        if(j==3)
                            xpos_score_disp = tri_pt(j,1);    ypos_score_disp = tri_pt(j,2)+sy+10;
                        end
                        destrect=[xpos_score_disp-sx/2,ypos_score_disp-sy/2,xpos_score_disp+sx/2,ypos_score_disp+sy/2];
                        bet_money_table(2,j)=double(money_given*bet_score_table(2,j));
                        str_score=sprintf('%01.2f',bet_money_table(2,j));
                        DrawFormattedText(wPtr, str_score, xpos_score_disp, ypos_score_disp);
                    end

                    % display!
                    Screen(wPtr, 'Flip');
                    % take a snapshot
                    if(DO_TAKE_SNAPSHOT==1)
                        snapshot=Screen(wPtr, 'GetImage', [1, 1, floor(1.0*SCREEN_RESOLUTION(1)), floor(1.0*SCREEN_RESOLUTION(2))]);
                        imageArray=[imageArray; {snapshot}];
                    end

                    % 1. Get cursor (L:37 R:39)
                    decision_made=0;
                    while(~decision_made)
                        [secs, keyCode] = KbPressWait;
                        [tmp tmp_key_code]=find(keyCode==1);                zzz=[zzz tmp_key_code];
                        if(tmp_key_code==KEY_L) % L pressed
                            current_cursor_selection_ind=max(1,current_cursor_selection_ind-1);
                            decision_made=1;
                        end
                        if(tmp_key_code==KEY_R) % R pressed
                            current_cursor_selection_ind=min(size(tri_pt_all,1),current_cursor_selection_ind+1);
                            decision_made=1;
                        end
                        if(tmp_key_code==KEY_Y) % 'y' pressed
                            decision_made=1;            Is_bet=1;
                            bet_score_table(2,:)=Barycentric_coord_table(2,:);
                        end
                        if(tmp_key_code==KEY_N) % 'n' pressed - reset
                            decision_made=1;            Is_bet=0;
                            current_cursor_selection_ind=current_cursor_selection_ind_default;

                        end

                    end


                end

                % save
                HIST_Barycentric_coord_table_type1_bonus{block,kk}=Barycentric_coord_table; % sum to one
                HIST_bet_score_table_type1_bonus{block,kk}=bet_score_table; % same as the above
                HIST_bet_money_table_type1_bonus{block,kk}=bet_money_table; % moeny bet amount

            end

        end


        % block ends
        str_block_intro=sprintf('*** You have finished Bonus block. Press any key to continue. ***');
        DrawFormattedText(wPtr, str_block_intro, 'center', 'center',[0, 0, 0, 255]);
        Screen(wPtr, 'Flip');
        WaitSecs(0.3);
        KbWait;
        % take a snapshot
        if(DO_TAKE_SNAPSHOT==1)
            snapshot=Screen(wPtr, 'GetImage', [1, 1, floor(1.0*SCREEN_RESOLUTION(1)), floor(1.0*SCREEN_RESOLUTION(2))]);
            imageArray=[imageArray; {snapshot}];
        end



    end %% Bonus block over




end % dont run bonus blo=1:0:0ck here




%% Ending message
str_end=sprintf('- Our experiments is over. Press any key to quit. -');
DrawFormattedText(wPtr, str_end, 'center', 'center');
Screen(wPtr, 'Flip');
% KbWait; % temporarily disabled for test APR 21
% take a snapshot
if(DO_TAKE_SNAPSHOT==1)
    snapshot=Screen(wPtr, 'GetImage', [1, 1, floor(1.0*SCREEN_RESOLUTION(1)), floor(1.0*SCREEN_RESOLUTION(2))]);
    imageArray=[imageArray; {snapshot}];
end













%% save all data

if(KEEP_PICTURE_OPT==1)
    data_save.HIST_bet_money_table_type1=cell(Tot_block,1);
    data_save.HIST_bet_score_table_type1=cell(Tot_block,1);
    data_save.HIST_Barycentric_coord_table_type1=cell(Tot_block,1);
    
end
data_save.HIST_bet_score_table_type2=cell(Tot_block,1);
for k=1:1:Tot_block
    if(KEEP_PICTURE_OPT==1)
        [tmp dd]=sort(HIST_bet_money_table_type1{k,1}(1,:)); % sort by index
        data_save.HIST_bet_money_table_type1{k,1}(:,:)=HIST_bet_money_table_type1{k,1}(:,dd);
        [tmp dd]=sort(HIST_bet_score_table_type1{k,1}(1,:)); % sort by index
        data_save.HIST_bet_score_table_type1{k,1}=HIST_bet_score_table_type1{k,1}(:,dd);
        [tmp dd]=sort(HIST_Barycentric_coord_table_type1{k,1}(1,:)); % sort by index
        data_save.HIST_Barycentric_coord_table_type1{k,1}=HIST_Barycentric_coord_table_type1{k,1}(:,dd);
    end
    [tmp dd]=sort(HIST_bet_score_table_type2{k,1}(1,:)); % sort by index
    data_save.HIST_bet_score_table_type2{k,1}=HIST_bet_score_table_type2{k,1}(:,dd);
end
data.HIST_mapping_index2filenum=HIST_current_pool;
if(KEEP_PICTURE_OPT==1)
%     data.data.HIST_mapping_index2filenum_bonus=HIST_current_pool_bonus;
%     [tmp dd]=sort(HIST_bet_money_table_type1_bonus{1,1}(1,:)); % sort by index
%     data_save.HIST_bet_money_table_type1_bonus{1,1}(:,:)=HIST_bet_money_table_type1_bonus{1,1}(:,dd);
%     [tmp dd]=sort(HIST_bet_score_table_type1_bonus{1,1}(1,:)); % sort by index
%     data_save.HIST_bet_score_table_type1_bonus{1,1}=HIST_bet_score_table_type1_bonus{1,1}(:,dd);
%     [tmp dd]=sort(HIST_Barycentric_coord_table_type1_bonus{1,1}(1,:)); % sort by index
%     data_save.HIST_Barycentric_coord_table_type1_bonus{1,1}=HIST_Barycentric_coord_table_type1_bonus{1,1}(:,dd);
%     [tmp dd]=sort(HIST_bet_score_table_type2_bonus{1,1}(1,:)); % sort by index
%     data_save.HIST_bet_score_table_type2_bonus{1,1}=HIST_bet_score_table_type2_bonus{1,1}(:,dd);
%     data_save.HIST_schedule_bonus=HIST_schedule_bonus;
%     data_save.HIST_reward_bonus=HIST_reward_bonus;
%     data_save.img_set_bonus=img_set_bonus;
%     data_save.img0_set_bonus=img0_set_bonus;
end
data_save.HIST_schedule=HIST_schedule;
data_save.HIST_reward=HIST_reward;
data_save.EXP_NAME=EXP_NAME;
data_save.Tot_block=Tot_block;
data_save.Tot_trial=Tot_trial;
data_save.Tot_trial_s=Tot_trial_s;
data_save.tot_num_img_per_block=tot_num_img_per_block;

data_save.img_set=img_set;
data_save.img0_set=img0_set;

% data_save.sec_stim_ready=sec_stim_ready; %(sec)
% data_save.sec_trial_ready=sec_trial_ready; %(sec)
data_save.sec_stim_interval=sec_stim_interval;%1.5; %(sec)
data_save.earliest_novel_stimulus_show_in_block=earliest_novel_stimulus_show_in_block; % 0: any trial in block, 1: last?

data_save.reward_mag=reward_mag;
data_save.reward_prob=reward_prob;

data_save.money_given=money_given; % point
data_save.img_present_ratio=img_present_ratio; % novelty level :: size = tot_num_img_per_block-1.



% eval(['save ' 'results\' EXP_NAME '.mat' ' data_save']);

%% save snapshots : CAUTION HEAVY PROCESS - might take a minute.
if(DO_TAKE_SNAPSHOT==1)
    for j=1:1:size(imageArray,1)
        str=sprintf('snapshot_dispay_exp_%03d.png',j);
        imwrite(imageArray{j},['snapshot\' str],'png');
    end
end






%% save the (updated) image usage matrix (overwriting)
file_imgind_sv_name=[EXP_NAME '_image_usage_info.mat'];
file_name_sv=[save_path file_imgind_sv_name];
% eval(['save ' file_name_sv ' Tot_img_info Tot_img_info_Tag']);
save(file_name_sv,'Tot_img_info','Tot_img_info_Tag','HIST_event_info_all','HIST_event_type');

%% save all variables
file_sv_name=[EXP_NAME sprintf('_%d.mat',index_num)];
file_name_sv=[save_path file_sv_name];
save(file_name_sv,'*');


%% session end sound
Beeper('low', 0.4, 0.4)
Beeper('medium', 0.4, 0.5)
Beeper('high', 0.4, 0.6)



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



output_info=1;
end


