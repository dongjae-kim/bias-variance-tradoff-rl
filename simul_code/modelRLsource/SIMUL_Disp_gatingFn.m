close all
clear all

DRAW=[1 0 0];

%% Draw dynamics of the arbitrator

% alpha : transition rate of SARSA->FWD by looking at inv_Fano(performance) of SARSA
% beta : transition rate of FWD->SARSA by looking at inv_Fano(performance) of FWD

% defines the rate of transition to the other "when both models behave badly"
% should be [alpha.A>beta.A] = "I would rather choose FWD when both go bad"

if(DRAW(1)==1)
    
    alpha.A=1.3*3; % fwd
    beta.A=11.7; % sarsa
    
    % defines "credibility" - how fastly do you give a trust when it behaves good.
    % should be [alpha.B>beta.B]
    alpha.B=2.67; % fwd
    beta.B=10; % sarsa
    
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








