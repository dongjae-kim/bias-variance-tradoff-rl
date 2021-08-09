% DP 써서 플랏하는 부분대해서.
clear; clc;
load('DPplot');
addpath([pwd '\DP']);
max_sub = size(outputval.full_arb.map,2);
SPE ={};
RPE ={};

for i = 1 : max_sub
    SPE{i} = outputval.full_arb.myAribtration(i).F_HIST_SPE(2,:); 
    RPE{i} = outputval.full_arb.myAribtration(i).F_HIST_RPE(2,:);
    
    
end
close all;

alphaSET= [0: 5 :100];
alphaSET(1)=1;
for a=1 : size(alphaSET,2)
    %% 2D case ㅜㅜ
    for i = 1 : max_sub
        einstein = 1;
        while(einstein == 1)
            try
                close all;
                X = [ SPE{i} ; RPE{i} ];
                label = ones(1,size(X,2));
                figure(1);
                title(['2-D case ' num2str(i) 'th subject']);
                plotClass(X,label);
                saveas(gcf,[pwd '\DP\DPrf\2D_CASE_BEFORE_' num2str(i) 'TH_SUBJECT']);
                saveas(gcf,[pwd '\DP\DPrf\2D_CASE_BEFORE_' num2str(i) 'TH_SUBJECT.png']);
                
                
                %param initialization
                [d,n]=size(X);
                mu = mean(X,2);
                Xo = bsxfun(@minus,X,mu);
                s = sum(Xo(:).^2)/(d*n);
                opt.kappa=1;
                opt.m=mean(X,2);
                opt.nu=d;
                opt.S = s*eye(d);
                
                
                opt.alpha= alphaSET(a);
                
                [y,model] = mixGaussGb(X, opt);
                save([pwd '\DP\rawmat\2Dcase_' num2str(i) 'th_SUB_' num2str(opt.alpha) 'A'], 'y','model');
                
                figure;
                plotClass(X,y);
                saveas(gcf,[pwd '\DP\DPrf\2D_CASE_A' num2str(opt.alpha) '_' num2str(i) 'TH_SUBJECT']);
                saveas(gcf,[pwd '\DP\DPrf\2D_CASE_A' num2str(opt.alpha) '_' num2str(i) 'TH_SUBJECT.png']);
                title(['# of cluster of subject ' num2str(i) ': ' num2str(max(y))]);
                
                
                disp('##########################');
                
                
                disp( ['SUBJECT #' num2str(i)  ' >>>> alpha value : ' num2str(a) ' >>>> TOTAL # of cluster : ' num2str(max(y))]);
                
                disp('##########################');
                
                
                for ii = 1 : max(y)
                    disp(['# of points in Cluster #' num2str(ii) ' : ' num2str(sum(y == ii))]);
                end
                figure;
                histogram2(X(1,:),X(2,:),100);
                title([num2str(i) 'th subject histogram2']);
                saveas(gcf,[pwd '\DP\DPrf\HIST2_' num2str(i) 'TH_SUBJECT']);
                saveas(gcf,[pwd '\DP\DPrf\HIST2_' num2str(i) 'TH_SUBJECT.png']);
                
                figure;
                histogram(X(1,:),100);
                title([num2str(i) 'th subject histogram of SPE']);
                saveas(gcf,[pwd '\DP\DPrf\HIST_' num2str(i) 'TH_SUBJECT_SPE']);
                saveas(gcf,[pwd '\DP\DPrf\HIST_' num2str(i) 'TH_SUBJECT_SPE.png']);
                
                figure;
                histogram(X(2,:),100);
                title([num2str(i) 'th subject histogram of RPE']);
                saveas(gcf,[pwd '\DP\DPrf\HIST_' num2str(i) 'TH_SUBJECT_RPE']);
                saveas(gcf,[pwd '\DP\DPrf\HIST_' num2str(i) 'TH_SUBJECT_RPE.png']);
                einstein = 0;
            catch
                einstein=1;
            end
        end
        


    end
end
%% 2D case with forgetting factor

%  あとでね



%% 1D case 
for a= 1 : size(alphaSET,2)
    %% 2D case ㅜㅜ
    for i = 1 : max_sub
        einstein = 1;
        while(einstein == 1)
            try
                close all;
                X = [ SPE{i} ; zeros(1,size(SPE{i},2)) ];
                label = ones(1,size(X,2));
                figure(1);
                title(['1-D SPE case ' num2str(i) 'th subject']);
                plotClass(X,label);
                saveas(gcf,[pwd '\DP\DPrf\1D_CASE_SPE_BEFORE_' num2str(i) 'TH_SUBJECT']);
                saveas(gcf,[pwd '\DP\DPrf\1D_CASE_SPE_BEFORE_' num2str(i) 'TH_SUBJECT.png']);
                
                
                %param initialization
                [d,n]=size(X);
                mu = mean(X,2);
                Xo = bsxfun(@minus,X,mu);
                s = sum(Xo(:).^2)/(d*n);
                opt.kappa=1;
                opt.m=mean(X,2);
                opt.nu=d;
                opt.S = s*eye(d);
                
                
                opt.alpha= alphaSET(a);
                
                [y,model] = mixGaussGb(X, opt);
                save([pwd '\DP\rawmat\1DcaseSPE_' num2str(i) 'th_SUB_' num2str(opt.alpha) 'A'], 'y','model');
                
                figure;
                plotClass(X,y);
                saveas(gcf,[pwd '\DP\DPrf\1D_CASE_SPE_A' num2str(opt.alpha) '_' num2str(i) 'TH_SUBJECT']);
                saveas(gcf,[pwd '\DP\DPrf\1D_CASE_SPE_A' num2str(opt.alpha) '_' num2str(i) 'TH_SUBJECT.png']);
                title(['# of cluster of subject, 1D_SPE ' num2str(i) ': ' num2str(max(y))]);
                
                
                disp('##########################');
                
                
                disp( ['1D_SPE_SUBJECT #' num2str(i)  ' >>>> alpha value : ' num2str(a) ' >>>> TOTAL # of cluster : ' num2str(max(y))]);
                
                disp('##########################');
                
                
                for ii = 1 : max(y)
                    disp(['1D_SPE # of points in Cluster #' num2str(ii) ' : ' num2str(sum(y == ii))]);
                end
                
                figure;
                histogram(X(1,:),100);
                title([num2str(i) 'th subject histogram of SPE']);
                saveas(gcf,[pwd '\DP\DPrf\1D_HIST_' num2str(i) 'TH_SUBJECT_SPE']);
                saveas(gcf,[pwd '\DP\DPrf\1D_HIST_' num2str(i) 'TH_SUBJECT_SPE.png']);
                
                einstein = 0;
            catch
                einstein=1;
            end
        end
        


    end
end


for a= 1 : size(alphaSET,2)
    %% 2D case ㅜㅜ
    for i = 1 : max_sub
        einstein = 1;
        while(einstein == 1)
            try
                close all;
                X = [ RPE{i} ; zeros(1,size(RPE{i},2)) ];
                label = ones(1,size(X,2));
                figure(1);
                title(['1-D RPE case ' num2str(i) 'th subject']);
                plotClass(X,label);
                saveas(gcf,[pwd '\DP\DPrf\1D_CASE_RPE_BEFORE_' num2str(i) 'TH_SUBJECT']);
                saveas(gcf,[pwd '\DP\DPrf\1D_CASE_RPE_BEFORE_' num2str(i) 'TH_SUBJECT.png']);
                
                
                %param initialization
                [d,n]=size(X);
                mu = mean(X,2);
                Xo = bsxfun(@minus,X,mu);
                s = sum(Xo(:).^2)/(d*n);
                opt.kappa=1;
                opt.m=mean(X,2);
                opt.nu=d;
                opt.S = s*eye(d);
                
                
                opt.alpha= alphaSET(a);
                
                [y,model] = mixGaussGb(X, opt);
                save([pwd '\DP\rawmat\1DcaseRPE_' num2str(i) 'th_SUB_' num2str(opt.alpha) 'A'], 'y','model');
                
                figure;
                plotClass(X,y);
                saveas(gcf,[pwd '\DP\DPrf\1D_CASE_RPE_A' num2str(opt.alpha) '_' num2str(i) 'TH_SUBJECT']);
                saveas(gcf,[pwd '\DP\DPrf\1D_CASE_RPE_A' num2str(opt.alpha) '_' num2str(i) 'TH_SUBJECT.png']);
                title(['# of cluster of subject, 1D_RPE ' num2str(i) ': ' num2str(max(y))]);
                
                
                disp('##########################');
                
                
                disp( ['1D_RPE_SUBJECT #' num2str(i)  ' >>>> alpha value : ' num2str(a) ' >>>> TOTAL # of cluster : ' num2str(max(y))]);
                
                disp('##########################');
                
                
                for ii = 1 : max(y)
                    disp(['1D_RPE, # of points in Cluster #' num2str(ii) ' : ' num2str(sum(y == ii))]);
                end
                
                figure;
                histogram(X(1,:),100);
                title([num2str(i) 'th subject histogram of SPE']);
                saveas(gcf,[pwd '\DP\DPrf\1D_HIST_' num2str(i) 'TH_SUBJECT_RPE']);
                saveas(gcf,[pwd '\DP\DPrf\1D_HIST_' num2str(i) 'TH_SUBJECT_RPE.png']);
                
                einstein = 0;
            catch
                einstein=1;
            end
        end
        


    end
end

% 1D with forgetting factor.