function SIMUL_Arbitration_v3_GA_Analysis_level2(maniX,maniY,Y_all,X_in_fwd,X_in_sarsa,X_in_fwd_mean,X_in_fwd_var,X_in_sarsa_mean,X_in_sarsa_var)

%% run this after SIMUL_Arbitration_v3_GA_Analysis.m
%% maniX and maniY should be placed in the workspace!!!

%% Further analysis following mani : SHOULD put in the ppt report!!!! (MAR 2)
%% [note] the case without outlier-type subjects shows also interesing!!!
%% [note] p=1 (expectation) hypothesis shows twice times the fitness of p=1e1 !!!
%% [note] do real pca and check this out again.

close all
figure_size=[600,500];
trial_array=[1:1:size(maniX,2)];

%% 0. EM clustering and display
num_cluster=3;
[label, model, llh] = emgm(maniX', num_cluster);
str_fig=sprintf('Projected beh. profile (EM clustered)');
f2=figure('Name',str_fig,'Position',[10+2*figure_size(1) 10 figure_size]);
PLOTCLR(maniY(:,1),maniY(:,2),label);


%% 1. averaged profile

% average profile for each subject
figure('Name','Average behavior profile for each subject')
hold on
avr_beh_profile=[];
each_sbj_ind0=0;
for each_sbj=0:1:max(Y_all(12,:))
    [sbj_line]=find(Y_all(12,:)==each_sbj);
    if(length(sbj_line)>0) %if correponding sub exists,
        each_sbj_ind0=each_sbj_ind0+1;
        avr_beh_profile=[avr_beh_profile; mean(maniX(sbj_line,:))];
    str_legend{1,each_sbj_ind0}=['sbj#' num2str(each_sbj)];
    end    
end
plot([1:1:size(maniX,2)],avr_beh_profile)
legend(str_legend)
hold off

%% 1.5. averaged invFano
figure('Name','Average inv Fano (Sbj0 use only)')
var_invFano_fwd=var(X_in_fwd);  var_invFano_sarsa=var(X_in_sarsa);
mean_invFano_fwd=mean(X_in_fwd);  mean_invFano_sarsa=mean(X_in_sarsa);
shade_mat_invFano_fwd=reshape(var_invFano_fwd',[size(var_invFano_fwd',1) 1 size(var_invFano_fwd',2)]);
shade_mat_invFano_sarsa=reshape(var_invFano_sarsa',[size(var_invFano_sarsa',1) 1 size(var_invFano_sarsa',2)]);
boundedline(trial_array,mean_invFano_fwd,shade_mat_invFano_fwd,'r',trial_array,mean_invFano_sarsa,shade_mat_invFano_sarsa,'b');
axis([min(trial_array) max(trial_array) 0 1]);


%% 2. projected distribution from mani
f1=figure('Name','Projected distribution (colorcode:sbj#)');
PLOTCLR(maniY(:,1),maniY(:,2),Y_all(12,:))

%% 3. click the intersting point # in the figure, and then plot the corresponding profile!

% parameters in Y_all
selection=8; % ###########[OPTION]################
string=cell(1,11);
string{1,1}='PE tolerance (fwd)';   string{1,11}='PE tolerance (sarsa)';
string{1,2}='Time step';
string{1,3}='A_12'; string{1,4}='B_12'; string{1,5}='A_21'; string{1,6}='B_21';
string{1,7}='learning rate';
string{1,8}='p';
string{1,9}='tau(softmax)';
string{1,10}='starting model'; string{1,12}='subject#';

% display color-coded by feature of interest
str_fig=sprintf('Projected beh. profile (colorcode: %s)',string{1,selection});
f2=figure('Name',str_fig,'Position',[50+2*figure_size(1) 10 figure_size]);
PLOTCLR(maniY(:,1),maniY(:,2),Y_all(selection,:));

% display the boxplots of p 
str_box=sprintf('distribution of %s',string{1,selection});
figure('Name',str_box)
boxplot(Y_all(selection,:)',Y_all(12,:)');
% and then display the corresponding function shape
if(selection==8) % display only for p
    str_box=sprintf('Q-value as a function of Q_(fwd) and Q_(sarsa) with p');
    f2_2=figure('Name',str_box);
    f_ax2_2 = axes('Parent',f2_2);
    range_minmax=[0 20];
    resolution=100;
    Q_range=[range_minmax(1):(range_minmax(2)-range_minmax(1))/(resolution-1):range_minmax(2)];
    [Q_x,Q_y] = meshgrid(Q_range, Q_range);
    pq=mean(Y_all(8,:));
    Q_z=(Q_x.^pq+Q_y.^pq).^(1/pq);
    surf(f_ax2_2 ,Q_x,Q_y,Q_z)
    view(2); colorbar;
end


% display for the control
str_fig=sprintf('Projected beh. profile (colorcode: data#)');
f3=figure('Name',str_fig,'Position',[10+2*figure_size(1) 10+figure_size(1) figure_size]);
PLOTCLR(maniY(:,1),maniY(:,2),[1:1:size(maniY,1)])
datacursormode on
dcm_obj=datacursormode(f3);
set(dcm_obj,'UpdateFcn',@plot_behavior_profile)

f_disp=figure('Name','Behavior profile display','Position',[10 10+figure_size(1) figure_size]);
f_ax = axes('Parent',f_disp);

f_disp2=figure('Name','P(fwd)','Position',[10+figure_size(1) 10+figure_size(1) figure_size]); %position-[left,bottom,w,h]
f_ax2 = axes('Parent',f_disp2);
cardinalityD=20;
range_minmax=[0 1]; % normalized
resolution=100;
inv_Fano=[range_minmax(1):(range_minmax(2)-range_minmax(1))/(resolution-1):range_minmax(2)];
[X_disp2,Y_disp2] = meshgrid(inv_Fano, inv_Fano);

f_disp3=figure('Name','transition rate: fwd->sarsa(red), sarsa->fwd(blue)','Position',[10+figure_size(1) 10 figure_size]);
f_ax3 = axes('Parent',f_disp3);

f_disp5=figure('Name','fwd:mean(red),var(red dot), sarsa:mean(blue),var(blue dot)','Position',[10+figure_size(1) 10 figure_size]);
f_ax5 = axes('Parent',f_disp5);

f_disp4=figure('Name','invFano: fwd(red), sarsa(blue)','Position',[10 10 figure_size]);
f_ax4 = axes('Parent',f_disp4);


    function [txt]=plot_behavior_profile(~,event_obj)
        
        %% get position
        pos = get(event_obj,'Position');
        pos_prev=pos;
        str_feature=sprintf('%s=',string{1,selection});
        txt = {['X=',num2str(pos(1))],['Y=',num2str(pos(2))],['dataID=',num2str(pos(3))],['sbjID=',num2str(Y_all(12,pos(3)))],[str_feature,num2str(Y_all(selection,pos(3)))]};
        
        %% plot behavior on display1
        plot(f_ax,trial_array,maniX(pos(3),:))
        
        
        %% plot n_inf = P(fwd) on display2
        param_alpha.A=Y_all(3,pos(3));
        param_alpha.B=Y_all(4,pos(3));
        param_beta.A=Y_all(5,pos(3));
        param_beta.B=Y_all(6,pos(3));
        Fn_alpha=param_alpha.A./(1.+exp(param_alpha.B*inv_Fano));
        Fn_beta=param_beta.A./(1.+exp(param_beta.B*inv_Fano));
        Fn_n_inf=param_beta.A./(1.+exp(param_beta.B*Y_disp2))./(param_alpha.A./(1.+exp(param_alpha.B*X_disp2))+param_beta.A./(1.+exp(param_beta.B*Y_disp2)));
        surf(f_ax2,X_disp2,Y_disp2,Fn_n_inf)
        %         view(2) % view point
        
        %% plot alpha/beta on display3
        plot(f_ax3,inv_Fano,Fn_alpha,'r',inv_Fano,Fn_beta,'b')
        
        %% plot mean/var of fwd(red) and sarsa(blue) on display5                
        cla(f_ax5,'reset')
        shade_mat_sarsa=reshape(X_in_sarsa_var(pos(3),:)',[size(X_in_sarsa_var(pos(3),:)',1) 1 size(X_in_sarsa_var(pos(3),:)',2)]);
        shade_mat_fwd=reshape(X_in_fwd_var(pos(3),:)',[size(X_in_fwd_var(pos(3),:)',1) 1 size(X_in_fwd_var(pos(3),:)',2)]);
%         boundedline(trial_array,X_in_sarsa_mean(pos(3),:),shade_mat_sarsa,'alpha','b',trial_array,X_in_fwd_mean(pos(3),:),shade_mat_fwd,'alpha','r',f_ax5);
        boundedline(trial_array,X_in_sarsa_mean(pos(3),:),shade_mat_sarsa,'b',trial_array,X_in_fwd_mean(pos(3),:),shade_mat_fwd,'r',f_ax5);        
%         plot(f_ax5,trial_array,X_in_fwd_mean(pos(3),:),'r-',X_in_fwd_var(pos(3),:),'r-')

        
        %% plot F_fwd (red), F_sarsa (blue) on display4
        plot(f_ax4,trial_array,X_in_fwd(pos(3),:),'r',trial_array,X_in_sarsa(pos(3),:),'b')
        

        
    end

%     function [out]=plot_behavior_profile(hObject,~)
%         pos=get(hObject,'CurrentPoint')
        % access to a child object
%         axesObjs = get(hObject, 'Children');
%         dataObjs = get(axesObjs(2), 'Children');
%         x1=get(dataObjs(1),'Xdata')
%         y1=get(dataObjs(1),'Ydata')
%         objTypes = get(dataObjs, 'Type')

end