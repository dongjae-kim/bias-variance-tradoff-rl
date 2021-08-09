% display subject's parameter in 'F:\0-Program\MATLABwork\work\ModelRL\result_simul\SBJ_structure_each/batch' and
% merge each SBJ structure and make a single one.
% display their NegLogLik value

clear all
close all

path_null='F:\0-Program\MATLABwork\work\ModelRL';
path00='F:\0-Program\MATLABwork\work\ModelRL\result_simul\';
path0='F:\0-Program\MATLABwork\work\ModelRL\result_simul\SBJ_structure_batch\';
path01='F:\0-Program\MATLABwork\work\ModelRL\result_simul\SBJ_structure_each\';

addpath('F:\0-Program\MATLABwork\work\spm8')


%% option

LIST_SBJ={'Oliver', 'Hao', 'Breanna', 'Derek', 'Timothy', 'Teagan', 'Jeffrey', 'Seung', 'Carole', 'Tony', 'Surendra', 'Lark',...
    'Joaquin', 'DavidB', 'Christopher', 'Gjergji', 'Charles', 'Erin', 'Connor', 'Domenick', 'Thao', 'Arin', 'Pauline', 'Tho'};
list_sbj_included=[2:1:19 21:1:24]; ind_sbj_included=list_sbj_included; sbj_included=list_sbj_included;

% 1. option for batch optimization
file_postfix_batch={'Fullbatch1'};
file_postfix_casename_batch={};{'Fullbatch1','Fullbatch2','Fullbatch3','Fullbatch4','Fullbatch5','Fullbatch6','Fullbatch7'}; %{batch case, each case1, 2, ,,,}
num_param_batch=[]; [4 4 4]; % option (for AIC,BIC)

% 2. option for each optimization
% [NOTE]
% '11_v29/27' best for individual opt
% '11_v53' best for visualizing the similarity btw the batch param and individual param set
% num_param: learning rate/invTemp are fixed for Arbitrators.
file_postfix={'11_v70','11_v77','11_v78','11_v79',...
    '11_v80','11_v81','11_v82','11_v86','11_v87','11_v92','11_v93','11_v94','11_v97','11_v98'};% simlified BayesArb - Bayesian_Arbitration_v4.m
%{'11_v15','11_v21','11_v22','11_v23',...
%     '11_v24','11_v25','11_v27','11_v28','11_v29','11_v30','11_v31','11_v32','11_v33','11_v34','11_v35','11_v36','11_v37','11_v38','11_v39'};
% {'11_v70','11_v77','11_v78','11_v79',...
%     '11_v80','11_v81','11_v82','11_v86','11_v87','11_v92','11_v93','11_v94','11_v97','11_v98'};% simlified BayesArb - Bayesian_Arbitration_v4.m
ind_file_postfix_for_parameter_plot=[-1];
file_postfix_casename_each=file_postfix; %{batch case, each case1, 2, ,,,}
num_param_each=4*ones(1,length(file_postfix));%[4 4 4 4 4 4 4]; % option (for AIC,BIC)

DO_EXCLUDE_OUTLIED_PERFORMANCE=0; 

% pick the best from all param for each sbj
%[NOTE] in this case all testing model should be of the same type!!!
DO_MAKE_BEST_SBJ_PARAM=1;

%% gating function display option
model_to_disp_gating_fn=[3]; % 0: use parameters from batch, -1: use mean of each(01 case), -2: use mean of each(11 case), 0>:index of sbj in 'sbj_included' to display
model_to_disp_gating_fn_opt2=[1]; % only if 'model_to_disp_gating_fn'>1. the index of case in 'file_postfix'
DRAW=[1 1 0];



% merge
num_param=[num_param_batch num_param_each]; % [[batch case] [each case]] number of parameter that correspond to each case(model)

%% load files (batch)
% extract each-SBJ parameter
mean_data_each=[];
all_data_each0=cell(1,size(file_postfix_batch,2));
val_AIC_each0=[];    val_AICc_each0=[];    val_BIC_each0=[];    val_NegLogLik_each0=[];

data_batch=[];
for i=1:1:size(file_postfix_batch,2) % for batch case
    load_file_name_batch=['SBJ_structure_batch_exp_' file_postfix_batch{1,i} '.mat'];
    eval(['load ' path0 load_file_name_batch])
    SBJ_batch=SBJ;
    
    %extract batch-SBJ parameter
    clear data_batch
    data_batch=SBJ_batch{1,1}.model_BayesArb.param;
    
    % create full SBJ structure from the batch param
    % regressor part deleting and regenerating.
    update_SBJ_structure=0; % 0: no update/just read and use, 1: update the changes to the saved SBJ file
    mode.opt_ArbModel=0; % 0: full arbitrator, 1: invF-based, 2: mean-based, 3: uncertainty-based arbitrator
    mode.USE_FWDSARSA_ONLY=0; % 0: arbitration, 1: use fwd only, 2: use sarsa only
    mode.USE_BWDupdate_of_FWDmodel=1; % 1: use the backward update for goal-directed model (fwd model), 0: do not use
    mode.DEBUG_Q_VALUE_CHG=0; % Debug option 1: show Q-value before/after whenever there is a goal change.
    mode.path_ext=path_null;
    mode.total_simul=1; % # of total simulation repetition per subject
    mode.simul_process_display=0; % 1: display model's process, 0: no diplay
    mode.experience_sbj_events=[1 1]; % [pre main]  +1: experience exactly the same events(decision,state) as subjects. 0: model's own experience -1: use saved setting
        
    num_sbj_included=length(ind_sbj_included);
    ind_included=ind_sbj_included;
    for k=1:1:num_sbj_included
        LIST_SBJ_included{1,k}=LIST_SBJ{1,ind_sbj_included(k)};
    end
    for ff=1:1:length(list_sbj_included)
        disp(sprintf('- computing and writing LogLik error to SBJ structure (SBJ%02d)...',list_sbj_included(ff)));
        % find my subject in "SBJ" strucure of the 'SBJ_structure.mat' file
        did_find=0;
        if(size(findstr(load_file_name_batch,'batch'),1)~=0) % if reading batch file
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
        model_BayesArb.param=data_batch;
        %         if(isfield(SBJ0{1,1}.model_BayesArb.mode, 'boundary_12')==1) % parameter length=6 case
        if(length(data_batch)==6) % parameter length=6 case
            mode.boundary_12=SBJ0{1,1}.model_BayesArb.mode.boundary_12;
            mode.boundary_21=SBJ0{1,1}.model_BayesArb.mode.boundary_21;
            mode.param_length=6;
        end
        if(length(data_batch)==8) % parameter length=8 case
            mode.boundary_12=NaN; %n/a indicator
            mode.boundary_21=NaN; %n/a indicator
            mode.param_length=8;
        end
        %     SBJ0=eval_ArbitrationRL3(model_BayesArb.param,SBJ0,mode); % refresh and add the regressor part  : eval_ArbitrationRL2(x, SBJ_test, mode) for full BayesArb
        
        % LogLik computing
        mode.out=1;
        SBJ0{1,1}.model_BayesArb.val=eval_ArbitrationRL3(model_BayesArb.param,SBJ0,mode); % refresh and add the regressor part  : eval_ArbitrationRL2(x, SBJ_test, mode) for full BayesArb
        
        % # of data        
        num_data=0;
        for kk=1:1:size(SBJ0{1,1}.HIST_behavior_info,2)
            num_data=num_data+size(SBJ0{1,1}.HIST_behavior_info{1,kk},1);
        end
        SBJ0{1,1}.num_data=num_data*2;
        
        SBJ1{1,ff}=SBJ0{1,1};
    end
    clear SBJ
    clear SBJ_batch
    SBJ_batch=SBJ1;
    
    % extracting info
    disp(sprintf('-processing batch case (%d/%d)...',i,size(file_postfix_batch,2)));
    
    data_each=[];   val_NegLogLik=[];
    for j=1:1:size(SBJ_batch,2)
        data_each=[data_each; SBJ_batch{1,j}.model_BayesArb.param];
    end
    mean_data_each=[mean_data_each; mean(data_each)];
    all_data_each0{1,i}=data_each;
    
    % extract each-SBJ LogLik value
    num_data_each=[];
    for j=1:1:size(SBJ_batch,2)
        val_NegLogLik=[val_NegLogLik; SBJ_batch{1,j}.model_BayesArb.val];
        num_data_each=[num_data_each; SBJ_batch{1,j}.num_data];
    end
    % compute AIC and BIC
    [v_aic v_bic]=aicbic( (-1)*val_NegLogLik , num_param_batch(i)*ones(size(val_NegLogLik,1),1) , num_data_each );
    val_NegLogLik_each0=[val_NegLogLik_each0 val_NegLogLik];
    val_AIC_each0=[val_AIC_each0 v_aic];%2*val_NegLogLik+2*num_param_batch(i)];
    val_AICc_each0=[val_AICc_each0 2*val_NegLogLik+2*num_param_batch(i)+2*num_param_batch(i)*(num_param_batch(i)+1)./(num_data_each-num_param_batch(i)-1)];
    val_BIC_each0=[val_BIC_each0 v_bic];%2*val_NegLogLik+num_param_batch(i)*log(num_data_each)];
end





%% load each

mean_data_each=[];
all_data_each=cell(1,size(file_postfix,2));
val_AIC_each=[];    val_AICc_each=[];    val_BIC_each=[];    val_NegLogLik_each=[];
for i=1:1:size(file_postfix,2) % for each case
    disp(sprintf('-processing each case (%d/%d)...',i,size(file_postfix,2)));
    % load files (each)
    for j=1:1:length(sbj_included)
        load_file_name_each=[sprintf('SBJ_structure_sbj%02d_each_exp_',sbj_included(j)) file_postfix{1,i} '.mat'];
        eval(['load ' path01 load_file_name_each])
        SBJ_each{1,j}=SBJ{1,1};
        num_data=0;
        for kk=1:1:size(SBJ_each{1,j}.HIST_behavior_info,2)
            num_data=num_data+size(SBJ_each{1,j}.HIST_behavior_info{1,kk},1);
        end
        SBJ_each{1,j}.num_data=num_data*2;
    end
    
    % extract each-SBJ parameter
    data_each=[];   val_NegLogLik=[];  
    for j=1:1:size(SBJ_each,2)
        data_each=[data_each; SBJ_each{1,j}.model_BayesArb.param];
    end
    mean_data_each=[mean_data_each; mean(data_each)];
    all_data_each{1,i}=data_each;
    
    % extract each-SBJ LogLik value
    num_data_each=[];
    for j=1:1:size(SBJ_each,2)
        val_NegLogLik=[val_NegLogLik; SBJ_each{1,j}.model_BayesArb.val];
        num_data_each=[num_data_each; SBJ_each{1,j}.num_data];
    end
    % compute AIC and BIC
    [v_aic v_bic]=aicbic((-1)*val_NegLogLik,num_param_each(i)*ones(size(val_NegLogLik,1),1),num_data_each);
    val_NegLogLik_each=[val_NegLogLik_each val_NegLogLik];
    val_AIC_each=[val_AIC_each v_aic];%2*val_NegLogLik+2*num_param_each(i)];
    val_AICc_each=[val_AICc_each 2*val_NegLogLik+2*num_param_each(i)+2*num_param_each(i)*(num_param_each(i)+1)./(num_data_each-num_param_each(i)-1)];
    val_BIC_each=[val_BIC_each v_bic];%2*val_NegLogLik+num_param_each(i)*log(num_data_each)];
        
    
    % save the collected structure
    clear SBJ;
    SBJ=SBJ_each;    
    eval(['save ' path00 'SBJ_structure_each_exp_' file_postfix{1,i} '.mat SBJ' ])
    
    if(i==ind_file_postfix_for_parameter_plot)
        figure('Name', ['parameter plot' file_postfix{1,i}])  
        ind_for_length6_param=[1 2 3 5 7 8];
        for a=1:1:size(data_each,2)
            subplot(1,size(data_each,2),a)
            hold on
            boxplot(data_each(:,a))
            if(length(data_batch)~=0)
                plot([1],[data_batch(1,ind_for_length6_param(a))],'-.b','linewidth',1,'marker','*')
                axis([0.5 1.5 [min([data_each(:,a); data_batch(1,ind_for_length6_param(a))])*0.8 max([data_each(:,a); data_batch(1,ind_for_length6_param(a))])*1.2]])
                % find percentile
                tmp=find((sort(data_each(:,a))-data_batch(1,ind_for_length6_param(a)))>0);
                percentile_val=(tmp(1)-0.5)/length(data_each(:,a))*100;
                disp_title=[sprintf('pcnt=%2.1f',percentile_val)];
                title(disp_title);
            end            
            hold off
        end
%         title(file_postfix_casename{1,i});        
    end

    
end




%% Analysis: interaction between parameters
% pick the best 10 and check the interaction between parameters
if(size(val_NegLogLik_each,2)>=10) % if # of models for the comparison >=10
    
    R_squared_each_sbj=cell(1,length(list_sbj_included));
    for lll=1:1:length(list_sbj_included)
        
        num_top_models_for_analysis=10; % max 10
        list_testing_sbj=[lll]; % [1:1:length(list_sbj_included)]
        
        ind_model_tops=[];      param_set_tops=[];
        param_set_tops_allsbj=cell(1,length(list_testing_sbj));
        for pp1=list_testing_sbj
            [tmp_v srt_ind]=sort(val_NegLogLik_each(pp1,:));
            ind_model_tops=[ind_model_tops; srt_ind(1:num_top_models_for_analysis)]; % [top1 ~ top10]
            param_set_tops=[]; % top10 params for sbj#-pp1
            for pp2=1:1:size(ind_model_tops,2) % for each top10 model
                ind_model_retrieve=ind_model_tops(end,pp2);
                param_set_tops=[param_set_tops; all_data_each{1,ind_model_retrieve}(pp1,:)];
            end
            param_set_tops_allsbj{1,pp1}=param_set_tops;
        end
        
        figure('Name','param interaction')
        if(size(param_set_tops,2)==8)    ind_param_show=[1 2 3 5 7 8];      num_params=6; end
        if(size(param_set_tops,2)==6)    ind_param_show=[1 2 3 4 5 6];      num_params=6; end
        R_squared=0*ones(num_params,num_params); %[0,1]: % 1 means perfect correlation, 0 means independent
        poly_fit_param=cell(num_params,num_params);
        R_squared_collect=[]; param_collect=[];
        for qq1=1:1:num_params
            zz=[];
            for pp1=list_testing_sbj % collect all sbj param
                zz=[zz; param_set_tops_allsbj{1,pp1}(:,ind_param_show(qq1))];
            end
            param_collect=[param_collect, zz];
            for qq2=[(qq1+1):1:num_params]
                xx=[];  yy=[];
                for pp1=list_testing_sbj % collect all sbj param
                    xx=[xx; param_set_tops_allsbj{1,pp1}(:,ind_param_show(qq1))];
                    yy=[yy; param_set_tops_allsbj{1,pp1}(:,ind_param_show(qq2))];
                end
                % plot
                subplot(num_params,num_params,qq2+num_params*(qq1-1));
                hold on
                scatter(xx,yy,'k');
                % [http://en.wikipedia.org/wiki/Coefficient_of_determination]
                [poly_f] = polyfit(xx,yy,1); poly_fit_param{qq1,qq2}=poly_f;
                val_f = polyval(poly_f,xx);
                [r2 rmse] = rsquare2(yy,val_f);  % http://www.mathworks.com/matlabcentral/fileexchange/34492-r-square-the-coefficient-of-determination
                R_squared(qq1,qq2)=r2;  R_squared_collect=[R_squared_collect r2];
                R_squared_each_sbj{1,lll}=R_squared;
                val_f_sorted = polyval(poly_f,sort(xx));
                plot(sort(xx),val_f_sorted,'b');    %axis([0 1 0 1]);
                xlabel(sprintf('parameter%01d',qq1));    ylabel(sprintf('parameter%01d',qq2));
                %             title(strcat(['R^{2}=' num2str(r2)]))
                title(sprintf('R^{2}=%0.3f',r2))
                hold off
            end
        end % OUTPUT:  "R_squared" & "poly_fit_param" (# param x # param matrix or cell)
        % correlation coefficient
        [r,p] = corrcoef(param_collect);
        r_test=r-eye(num_params);
        r_test=r_test.*(p<0.05);
        [max(max(r_test)) min(min(r_test))]; % max and min of significant correlation coeff (p<0.05)
        
        close all
    end
     
    % 3D bar plot
    R_squared_each_sbj_mean=mean(reshape(cell2mat(R_squared_each_sbj), [length(ind_param_show),length(ind_param_show),length(list_sbj_included)]), 3);
    R_squared_each_sbj_SEM=std(reshape(cell2mat(R_squared_each_sbj), [length(ind_param_show),length(ind_param_show),length(list_sbj_included)]), [], 3)/sqrt(length(list_sbj_included));
    R_squared_each_sbj_std=std(reshape(cell2mat(R_squared_each_sbj), [length(ind_param_show),length(ind_param_show),length(list_sbj_included)]), [], 3);        
    ThreeDBarWithErrorBars(R_squared_each_sbj_mean,R_squared_each_sbj_SEM);
end




%% pick the best from all param for each sbj and construct the best SBJ param
val_AIC_each_best=[];    val_AICc_each_best=[];    val_BIC_each_best=[];    val_NegLogLik_each_best=[];
val_NegLogLik_best=[];

num_data_each_best=[];  num_param_each_best=num_param_each(1);
data_each_best=[];

if(DO_MAKE_BEST_SBJ_PARAM==1)
    [a_val ind_best_model_each]=min(val_NegLogLik_each');
    
    clear SBJ;
    for jj=1:1:length(ind_best_model_each) % each sbj save the corresponding param file
        % load files (each)
        ind_sbj_collect=jj;
        ind_case_collect=ind_best_model_each(jj);
        % collect the best
        load_file_name_each=[sprintf('SBJ_structure_sbj%02d_each_exp_',list_sbj_included(ind_sbj_collect)) file_postfix{1,ind_case_collect} '.mat'];
        eval(['load ' path01 load_file_name_each])
        SBJ_eachbest{1,jj}=SBJ{1,1};
        num_data=0;
        for kk=1:1:size(SBJ_eachbest{1,jj}.HIST_behavior_info,2)
            num_data=num_data+size(SBJ_eachbest{1,jj}.HIST_behavior_info{1,kk},1);
        end
        SBJ_eachbest{1,jj}.num_data=num_data*2;
        data_each_best=[data_each_best; SBJ_eachbest{1,jj}.model_BayesArb.param];
        
        val_NegLogLik_best=[val_NegLogLik_best; SBJ_eachbest{1,jj}.model_BayesArb.val];
        num_data_each_best=[num_data_each_best; SBJ_eachbest{1,jj}.num_data];
                
    end
    
    % compute AIC and BIC
    [v_aic v_bic]=aicbic((-1)*val_NegLogLik_best,num_param_each_best*ones(size(val_NegLogLik_best,1),1),num_data_each_best);
    val_NegLogLik_each_best=[val_NegLogLik_each_best val_NegLogLik_best];
    val_AIC_each_best=[val_AIC_each_best v_aic];%2*val_NegLogLik+2*num_param_each_best];
    val_AICc_each_best=[val_AICc_each_best 2*val_NegLogLik+2*num_param_each_best+2*num_param_each_best*(num_param_each_best+1)./(num_data_each_best-num_param_each_best-1)];
    val_BIC_each_best=[val_BIC_each_best v_bic];%2*val_NegLogLik+num_param_each_best*log(num_data_each_best)];
    
    SBJ=SBJ_eachbest;
    eval(['save ' path00 'SBJ_structure_each_exp_BESTcollection.mat SBJ' ])
    
    
end

    
    
    
    

%% merge batch and each case
if(DO_MAKE_BEST_SBJ_PARAM==1)   best_add=1; else best_add=0;    end
all_data_each_tmp=cell(1,size(file_postfix_batch,2)+size(file_postfix,2)+best_add);
file_postfix_casename=cell(1,size(file_postfix_batch,2)+size(file_postfix,2)+best_add);
i_abs=0;
for i=1:1:size(file_postfix_batch,2) % batch
    i_abs=i_abs+1;
    all_data_each_tmp{1,i_abs}=all_data_each0{1,i};
    file_postfix_casename{1,i_abs}=file_postfix_casename_batch{1,i};
end
for i=1:1:size(file_postfix,2) % each
    i_abs=i_abs+1;
    all_data_each_tmp{1,i_abs}=all_data_each{1,i};
    file_postfix_casename{1,i_abs}=file_postfix_casename_each{1,i};
end
% best collection
i_abs=i_abs+1;
all_data_each_tmp{1,i_abs}=data_each_best;
file_postfix_casename{1,i_abs}='bestcollected';
    
val_AIC_each=[val_AIC_each0 val_AIC_each val_AIC_each_best];    
val_AICc_each=[val_AICc_each0 val_AICc_each val_AICc_each_best];    
val_BIC_each=[val_BIC_each0 val_BIC_each val_BIC_each_best];    
val_NegLogLik_each=[val_NegLogLik_each0 val_NegLogLik_each val_NegLogLik_each_best];



%% some post-processing

if(DO_EXCLUDE_OUTLIED_PERFORMANCE==1)
    NUM_outlier_exclude=2;    
    [aa ind_whole]=sort(val_NegLogLik_each(:,1),'descend');
    ind_nonoutliers=ind_whole(1+NUM_outlier_exclude:1:end);
    
    val_NegLogLik_each=val_NegLogLik_each(ind_nonoutliers,:);
    val_AIC_each=val_AIC_each(ind_nonoutliers,:);
    val_AICc_each=val_AICc_each(ind_nonoutliers,:);
    val_BIC_each=val_BIC_each(ind_nonoutliers,:);       
end
if(DO_EXCLUDE_OUTLIED_PERFORMANCE==2)    
    bar_percentile=[30 99];
    rnk=round(bar_percentile/100*size(val_NegLogLik_each,1)+1/2);
    [aa ind_whole]=sort(val_NegLogLik_each(:,1),'descend');    
    ind_nonoutliers=ind_whole(rnk(1):1:rnk(2));
    
    val_NegLogLik_each=val_NegLogLik_each(ind_nonoutliers,:);
    val_AIC_each=val_AIC_each(ind_nonoutliers,:);
    val_AICc_each=val_AICc_each(ind_nonoutliers,:);
    val_BIC_each=val_BIC_each(ind_nonoutliers,:);    
end


%% display the scores (NegLogLik, AIC, BIC)

% mean-based approach may not be an ideal because this is usually used for normal distributions. 
figure('Name','Performance (mean)')
subplot(1,4,1); bar(mean(val_NegLogLik_each));
title('NegLogLik')
subplot(1,4,2); bar(mean(val_AIC_each));
title('AIC')
subplot(1,4,3); bar(mean(val_AICc_each));
title('AICc')
subplot(1,4,4); bar(mean(val_BIC_each));
title('BIC')

% median-based approach  
figure('Name','Performance (median)')
subplot(1,4,1); bar(median(val_NegLogLik_each));
title('NegLogLik')
subplot(1,4,2); bar(median(val_AIC_each));
title('AIC')
subplot(1,4,3); bar(median(val_AICc_each));
title('AICc')
subplot(1,4,4); bar(median(val_BIC_each));
title('BIC')


% median-based approached might make more sense because it is used for
% skewed distribution and robust against outliers.
% These are *highly skewed* distribution; check to see "plot(ones(22,5),val_BIC_each(:,1),'*')"
figure('Name','Performance - box plot')
subplot(1,4,1); boxplot(val_NegLogLik_each);
title('NegLogLik')
subplot(1,4,2); boxplot(val_AIC_each);
title('AIC')
subplot(1,4,3); boxplot(val_AICc_each);
title('AICc')
subplot(1,4,4); boxplot(val_BIC_each);
title('BIC')


% figure('Name','Performance distribution')
% grid on
% 
% ind_model_compare=[1 2];
% 
% subplot(2,2,1); 
% hold on
% pl1=val_NegLogLik_each(:,ind_model_compare(1)); pl2=val_NegLogLik_each(:,ind_model_compare(2));
% plot(pl1,pl2,'o');   line([0 max([pl1; pl2])],[0 max([pl1; pl2])]);
% xlabel(file_postfix_casename{1,ind_model_compare(1)});  ylabel(file_postfix_casename{1,ind_model_compare(2)});
% title('NegLogLik')
% hold off
% subplot(2,2,2); 
% hold on
% pl1=val_AIC_each(:,ind_model_compare(1)); pl2=val_AIC_each(:,ind_model_compare(2));
% plot(pl1,pl2,'o');   line([0 max([pl1; pl2])],[0 max([pl1; pl2])]);
% xlabel(file_postfix_casename{1,ind_model_compare(1)});  ylabel(file_postfix_casename{1,ind_model_compare(2)});
% title('AIC')
% hold off
% subplot(2,2,3); 
% hold on
% pl1=val_AICc_each(:,ind_model_compare(1)); pl2=val_AICc_each(:,ind_model_compare(2));
% plot(pl1,pl2,'o');   line([0 max([pl1; pl2])],[0 max([pl1; pl2])]);
% xlabel(file_postfix_casename{1,ind_model_compare(1)});  ylabel(file_postfix_casename{1,ind_model_compare(2)});
% title('AICc')
% hold off
% subplot(2,2,4); 
% hold on
% pl1=val_BIC_each(:,ind_model_compare(1)); pl2=val_BIC_each(:,ind_model_compare(2));
% plot(pl1,pl2,'o');   line([0 max([pl1; pl2])],[0 max([pl1; pl2])]);
% xlabel(file_postfix_casename{1,ind_model_compare(1)});  ylabel(file_postfix_casename{1,ind_model_compare(2)});
% title('BIC')
% hold off


%% display the Bayesian model selection score
% [NOTE] .xp=exceedance probability indicates the probability that nonlinear hemodynamic models were better than linear models, regardless of any other aspect of model structure
% REFERENCE:
% Stephan KE, Penny WD, Daunizeau J, Moran RJ, Friston KJ (2009)
% Bayesian Model Selection for Group Studies. NeuroImage 46:1004-1017
% AIC/BIC used as an approximation of the log model evidence p(y|m) -
% but our AIC/BIC is based on "negative log-likelihood so we flip the sign
% for the first input to spm_BMS
BMSscore{1,1}.info='AIC approximation of log evidence';
[BMSscore{1,1}.alpha,BMSscore{1,1}.exp_r,BMSscore{1,1}.xp]=spm_BMS((-0.5)*val_AIC_each);
BMSscore{1,2}.info='AICc approximation of log evidence';
[BMSscore{1,2}.alpha,BMSscore{1,2}.exp_r,BMSscore{1,2}.xp]=spm_BMS((-0.5)*val_AICc_each);
BMSscore{1,3}.info='BIC approximation of log evidence';
[BMSscore{1,3}.alpha,BMSscore{1,3}.exp_r,BMSscore{1,3}.xp]=spm_BMS((-0.5)*val_BIC_each);
figure('Name','BMS score')
bar(BMSscore{1,3}.xp)

%% display gating function
if(model_to_disp_gating_fn==0)
    param_gating_display0=data_batch([3 4]);
end
if(model_to_disp_gating_fn<0)
    param_gating_display0=mean_data_each((-1)*model_to_disp_gating_fn,[3 4]);
end
if(model_to_disp_gating_fn>0)
    ind_s=find(sbj_included==model_to_disp_gating_fn);
    param_gating_display0=all_data_each{1,model_to_disp_gating_fn_opt2}(ind_s,[3 4])
end
if(size(all_data_each{1,1},2)==6) % constrainded version
    param_gating_display=[param_gating_display0(1) log(5*param_gating_display0(1)-1) param_gating_display0(2) log(10*param_gating_display0(2)-1)];
end
if(size(all_data_each{1,1},2)==8) % constrainded version
    param_gating_display=all_data_each{1,model_to_disp_gating_fn_opt2}(ind_s,[3 4 5 6]);
end


%% Draw dynamics of the arbitrator

% alpha : transition rate of SARSA->FWD by looking at inv_Fano(performance) of SARSA
% beta : transition rate of FWD->SARSA by looking at inv_Fano(performance) of FWD

% defines the rate of transition to the other "when both models behave badly"
% should be [alpha.A>beta.A] = "I would rather choose FWD when both go bad"

if(DRAW(1)==1)
    
    alpha.A=param_gating_display(1); % fwd 1->2
    beta.A=param_gating_display(3); % sarsa 2->1
    
    % defines "credibility" - how fastly do you give a trust when it behaves good.
    % should be [alpha.B>beta.B]
    alpha.B=param_gating_display(2); % fwd 1->2
    beta.B=param_gating_display(4); % sarsa 2->1
    
    cardD=20;
    
    inv_Fano_factor_equilibrium=log(alpha.A/beta.A)/(alpha.B-beta.B);
    
    out=Disp_gatingFn(alpha,beta, cardD);
    
end


%% Draw Beta distribution

if(DRAW(2)==1)
    
    x_in=[0:0.01:1];
    beta_param_mat=[1 1; 2 2; 3 10; 4 30];
    
    figure('Name','Beta','Position',[100 200 400 600]);
    for i=1:1:size(beta_param_mat,1)
        
        y_out = betapdf(x_in,beta_param_mat(i,1),beta_param_mat(i,2));                
        subplot(size(beta_param_mat,1),1,i)
        plot(x_in,y_out,'r')
        
        [mean var]=betastat(beta_param_mat(i,1),beta_param_mat(i,2));
        str_text=sprintf('[mean:%2.1f  variance:%1.3f  Fano^{-1}=%2.1f]',mean,var,mean/var);
        text('units','pixels','position',[10 20],'fontsize',10,'string',str_text)
    end
    
end

%% Draw the relationship between inverse Fano and posterior Beta

if(DRAW(3)==1)
    
    x_in=[0:0.01:1];
    mat_a=[1:1:15, 16:3:100];
    mat_b=[1:1:15, 16:3:100];
        
    inv_Fano=zeros(length(mat_a),length(mat_b));
    mean_mat=zeros(length(mat_a),length(mat_b));
    var_mat=zeros(length(mat_a),length(mat_b));
    for i=1:1:length(mat_a)
        for j=1:1:length(mat_b)
            [mean_mat(i,j) var_mat(i,j)]=betastat(mat_a(i),mat_b(j));
            inv_Fano(i,j)=mean_mat(i,j)/var_mat(i,j);            
        end
    end
    
    figure('Name','inverse Fano as a function of mean-var');
    surf(mean_mat,var_mat,log10(inv_Fano))
%     contour(mean_mat,var_mat,log10(inv_Fano),50)
%     grid on
    colorbar
    
    figure('Name','inverse Fano as a function of a,b');
    %     surf(mat_a,mat_b,log10(inv_Fano))    
    contour(mat_a,mat_b,log10(inv_Fano),100)
    grid on
    colorbar
    
    disp('done')
end

disp('done')


save score_BayesArb_mean val_NegLogLik_each_best val_AIC_each_best val_AICc_each_best val_BIC_each_best;