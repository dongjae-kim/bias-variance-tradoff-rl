function [out]=Disp_gatingFn(param_alpha,param_beta, cardinalityD)

% param.A, param.B
% cardinalityD: # of total trials for computing Fano factor = (# of rewarded trials + # of non rewarded trials)
min_v=(cardinalityD+2)*(cardinalityD+3)/cardinalityD;
max_v=(cardinalityD+2)*(cardinalityD+3);
range_minmax=[min_v max_v];
range_minmax=[0 1]; % normalized

resolution=100;

inv_Fano=[range_minmax(1):(range_minmax(2)-range_minmax(1))/(resolution-1):range_minmax(2)];

out.Fn_alpha=zeros(1,length(inv_Fano));
out.Fn_beta=zeros(1,length(inv_Fano));
out.Fn_Tau=zeros(1,length(inv_Fano));
out.Fn_n_inf=zeros(1,length(inv_Fano));

out.Fn_alpha=param_alpha.A./(1.+exp(param_alpha.B*inv_Fano));
out.Fn_beta=param_beta.A./(1.+exp(param_beta.B*inv_Fano));

[X,Y] = meshgrid(inv_Fano, inv_Fano);

out.Fn_Tau=1./(param_alpha.A./(1.+exp(param_alpha.B*X))+param_beta.A./(1.+exp(param_beta.B*Y)));


out.Fn_n_inf=param_beta.A./(1.+exp(param_beta.B*Y))./(param_alpha.A./(1.+exp(param_alpha.B*X))+param_beta.A./(1.+exp(param_beta.B*Y)));

n_half=0.5*ones(length(inv_Fano),length(inv_Fano));

figure(10)
subplot(1,2,1);
plot(inv_Fano,out.Fn_alpha)
title('Alpha');
subplot(1,2,2);
plot(inv_Fano,out.Fn_beta)
title('Beta');

figure(11)
plot(inv_Fano,out.Fn_alpha,'r',inv_Fano,out.Fn_beta,'b')
title('Alpha(red)/Beta(blue)');

figure(12)
surf(X,Y,out.Fn_Tau)
title('Tau');
view(-8,38) % view point

figure(13)
hold on
surf(X,Y,out.Fn_n_inf)
% surf(X,Y,n_half);
hold off
title('n_{inf}');
view(-11,31) % view point
grid on


end