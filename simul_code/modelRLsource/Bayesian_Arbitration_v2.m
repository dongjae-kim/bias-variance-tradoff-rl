function [myArbitrator myState1 myState2]=Bayesian_Arbitration_v2(myArbitrator, myState1, myState2, myMap)

%% let myState1:FWD, myState2:SARSA
% ind_active_model = 1 or 2: which model is currently active.

% myArbitrator.K=3; -> trichonomy of PE
% myArbitrator.T=20; -> half-life
% myArbitrator.T_current=0; -> # of accumulated events. cannot exceed T.

% myArbitrator.m1_threshold_PE=zeros(1,myArbitrator.K-1); myArbitrator.m2_threshold_PE=zeros(1,myArbitrator.K-1);
% myArbitrator.m1_thr_PE=zeros(1,myArbitrator.K-1); myArbitrator.m2_thr_PE=zeros(1,myArbitrator.K-1);

%% Bayesian part - PE history
% % row: first index = smallest PE -> last index = largest PE
% % column: first index = most recent -> last index = most past
% myArbitrator.m1_PE_history=zeros(myArbitrator.K,myArbitrator.T);  myArbitrator.m2_PE_history=zeros(myArbitrator.K,myArbitrator.T);
% myArbitrator.m1_PE_num=zeros(myArbitrator.K,1);   myArbitrator.m2_PE_num=zeros(myArbitrator.K,1);

%% Bayesian part - mean,var,inv_Fano
% myArbitrator.m1_mean=zeros(myArbitrator.K,1); myArbitrator.m2_mean=zeros(myArbitrator.K,1);
% myArbitrator.m1_var=zeros(myArbitrator.K,1);  myArbitrator.m2_var=zeros(myArbitrator.K,1);
% myArbitrator.m1_inv_Fano=zeros(myArbitrator.K,1); myArbitrator.m2_inv_Fano=zeros(myArbitrator.K,1);

%% Dynamic arbitration part
% myArbitrator.m1_wgt=[.1 .7 .2]; myArbitrator.m2_wgt=[.1 .7 .2];
% myArbitrator.A_12, myArbitrator.B_12
% myArbitrator.A_21, myArbitrator.B_21
% myArbitrator.transition_rate12, myArbitrator.transition_rate21
% myArbitrator.m1_prob_prev=0.5;

%% Q-integration part
% myArbitrator.p=1; % 1:expectation, inf:winner-take-all
% myArbitrator.tau_softmax=0.5; % use the same value as sarsa/fwd
% myArbitrator.action_history=zeros(N_transition,1);
% myArbitrator.action_prob=zeros(N_transition,1);
% myArbitrator.Q=zeros(N_state,N_action);
% myArbitrator.state_history=zeros(N_transition,1);    myArbitrator.state_history(1,1)=1;
% myArbitrator.action_history=zeros(N_transition,1);



%% Preparation
if(myArbitrator.ind_active_model==1)
    myState=myState1;
else
    myState=myState2;
end
%% index, state_history, action_history synchronization
% simply inherit because both should be in the same state
myArbitrator.index=myMap.index-1;
myArbitrator.state_history(myArbitrator.index)=myState.state_history(myArbitrator.index);
myArbitrator.action_history(myArbitrator.index)=myState.action_history(myArbitrator.index);


%% Hierarchical Bayesian Inference

if(myArbitrator.ind_active_model==1)
    
    myArbitrator.T_current1=min(myArbitrator.T_current1+1,myArbitrator.T); % update # of accumulated events
    
    % 1. model 1 (m1)
    % (1) find the corresponding row
    [tmp ind_neg]=find((myArbitrator.m1_thr_PE-myState1.RPE_history(myArbitrator.index))<0);
    ind_update=length(ind_neg)+1;
    % (2) update the current column(=1) in PE_history
    myArbitrator.m1_PE_history(:,2:end)=myArbitrator.m1_PE_history(:,1:end-1); % shift 1 column (toward past)
    myArbitrator.m1_PE_history(:,1)=zeros(myArbitrator.K,1); % empty the first column
    myArbitrator.m1_PE_history(ind_update,1)=1; % add the count 1 in the first column
    myArbitrator.m1_PE_num=myArbitrator.m1_PE_history*myArbitrator.discount_mat'; % compute discounted accumulated PE
    % (3) posterior mean & var
    sumK=sum(myArbitrator.m1_PE_num);
    sumK_excl=sumK-myArbitrator.m1_PE_num;
    myArbitrator.m1_mean=(1+myArbitrator.m1_PE_num)/(myArbitrator.K+sumK);
    myArbitrator.m1_var=((1+myArbitrator.m1_PE_num)/((myArbitrator.K+sumK)^2))/(myArbitrator.K+sumK+1).*(myArbitrator.K+sumK_excl-1);
    myArbitrator.m1_inv_Fano=myArbitrator.m1_mean./myArbitrator.m1_var;
    
else
    
    myArbitrator.T_current2=min(myArbitrator.T_current2+1,myArbitrator.T); % update # of accumulated events
    
    % 2. model 2 (m2)
    % (1) find the corresponding row
    [tmp ind_neg]=find((myArbitrator.m2_thr_PE-myState2.RPE_history(myState2.index))<0);
    ind_update=length(ind_neg)+1;
    % (2) update the current column(=1) in PE_history
    myArbitrator.m2_PE_history(:,2:end)=myArbitrator.m2_PE_history(:,1:end-1); % shift 1 column (toward past)
    myArbitrator.m2_PE_history(:,1)=zeros(myArbitrator.K,1); % empty the first column
    myArbitrator.m2_PE_history(ind_update,1)=1; % add the count 1 in the first column
    myArbitrator.m2_PE_num=myArbitrator.m2_PE_history*myArbitrator.discount_mat'; % compute discounted accumulated PE
    % (3) posterior mean & var
    sumK=sum(myArbitrator.m2_PE_num);
    sumK_excl=sumK-myArbitrator.m2_PE_num;
    myArbitrator.m2_mean=(1+myArbitrator.m2_PE_num)/(myArbitrator.K+sumK);
    myArbitrator.m2_var=((1+myArbitrator.m2_PE_num)/((myArbitrator.K+sumK)^2))/(myArbitrator.K+sumK+1).*(myArbitrator.K+sumK_excl-1);
    myArbitrator.m2_inv_Fano=myArbitrator.m2_mean./myArbitrator.m2_var;
    
end

%% minimum accumulation - prerequisite for arbitration
COND1=(myArbitrator.T_current1>=1);
COND2=(myArbitrator.T_current2>=1);

%% Dynamic Arbitration
if(myArbitrator.ind_active_model==1) % alpha, beta
    input0=myArbitrator.m1_inv_Fano;
    input=myArbitrator.m1_wgt'*myArbitrator.m1_inv_Fano;    
else
    input0=myArbitrator.m2_inv_Fano;
    input=myArbitrator.m2_wgt'*myArbitrator.m2_inv_Fano;
end

myArbitrator.temp=[myArbitrator.ind_active_model; input/sum(input0)];

myArbitrator.transition_rate12=myArbitrator.A_12/(1+exp(myArbitrator.B_12*input/sum(input0)));
% myArbitrator.transition_rate12=myArbitrator.A_12/(1+exp(myArbitrator.B_12*input));
myArbitrator.transition_rate12_prev=myArbitrator.transition_rate12;
% myArbitrator.transition_rate21=myArbitrator.transition_rate21_prev;

% myArbitrator.transition_rate12=myArbitrator.transition_rate12_prev;
myArbitrator.transition_rate21=myArbitrator.A_21/(1+exp(myArbitrator.B_21*input/sum(input0)));
% myArbitrator.transition_rate21=myArbitrator.A_21/(1+exp(myArbitrator.B_21*input));
myArbitrator.transition_rate21_prev=myArbitrator.transition_rate21;

myArbitrator.Tau=1/(myArbitrator.transition_rate12+myArbitrator.transition_rate21);
myArbitrator.m1_prob_inf=myArbitrator.transition_rate21*myArbitrator.Tau;
myArbitrator.m1_prob=myArbitrator.m1_prob_inf+(myArbitrator.m1_prob_prev-myArbitrator.m1_prob_inf)*exp((-1)*myArbitrator.Time_Step/myArbitrator.Tau);
myArbitrator.m1_prob_prev=myArbitrator.m1_prob;
myArbitrator.m2_prob=1-myArbitrator.m1_prob;

%% choice of the model
myArbitrator.ind_active_model_prev=myArbitrator.ind_active_model;
if(myArbitrator.m1_prob>0.5)
    myArbitrator.ind_active_model=1;
    myArbitrator.num_m1_chosen=myArbitrator.num_m1_chosen+1;
    % there is no Q-value hand-over because sarsa computes Q based on SPE.
else
    myArbitrator.ind_active_model=2;    
    myArbitrator.num_m2_chosen=myArbitrator.num_m2_chosen+1;
    % Q-value hand-over : sarsa uses RPE-based Q only, does not use SPE.
%     if((myArbitrator.ind_active_model_prev==1)&&(myArbitrator.num_m2_chosen==0))
%         myState2.Q=myState1.Q;
%     end    
end

%% Q-value computation and action choice
% (1) Q-value computing
state_current=myArbitrator.state_history(myArbitrator.index);
action_current=myArbitrator.action_history(myArbitrator.index);
% [myArbitrator.ind_active_model myArbitrator.ind_active_model_prev myArbitrator.index state_current action_current]
myArbitrator.Q(state_current,action_current)=...
    ((myArbitrator.m1_prob*myState1.Q(state_current,action_current))^myArbitrator.p+...
    (myArbitrator.m2_prob*myState2.Q(state_current,action_current))^myArbitrator.p)^(1/myArbitrator.p);
% (2) softmax decision making
check_cond0=(myMap.IsTerminal(state_current)~=1);
if(check_cond0) % if not in a terminal state
    var_exp=exp(myArbitrator.tau_softmax*myArbitrator.Q(state_current,:)); % (N_actionx1)
    Prob_action=var_exp/sum(var_exp); % (N_actionx1)
    myArbitrator.action_prob(myArbitrator.index)=Prob_action(1);
    if(rand<Prob_action(1))
        myArbitrator.action_history(myArbitrator.index)=1;
    else
        myArbitrator.action_history(myArbitrator.index)=2;
    end
end


%% if you want transition rate only, then look at
% myArbitrator.transition_rate12, myArbitrator.transition_rate21
%% if you want it to make a choice, the look at
% myArbitrator.action_history(myArbitrator.index)



end