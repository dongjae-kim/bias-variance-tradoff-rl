%% HIST_event_Info 
clear;
LIST_SBJ={'Oliver', 'Hao', 'Breanna', 'Derek', 'Timothy', 'Teagan', 'Jeffrey', 'Seung', 'Carole', 'Tony', 'Surendra', 'Lark',...
    'Joaquin', 'DavidB', 'Christopher', 'Gjergji', 'Charles', 'Erin', 'Connor', 'Domenick', 'Thao', 'Arin', 'Pauline', 'Tho'};
list_sbj_included=[2:1:15 17:1:19 21:1:24];
num_sbj_included = length(list_sbj_included);
for k=1:1:num_sbj_included
    LIST_SBJ_included{1,k}=LIST_SBJ{1,list_sbj_included(k)};
end

resi_path = [pwd '/result_simul/'];
yoto = load([resi_path 'SBJ_structure.mat']);

for sub_ind = 1 : size(list_sbj_included,2)
    dat_name = [LIST_SBJ_included{1, sub_ind} '_fmri_*.mat'];
    datlist = dir([pwd '/result_save/' dat_name]);
    maxsess = size(datlist,1)-1;
    et = load([pwd '/result_save/' datlist(maxsess).name]);
    for sei = 1 : maxsess
        yoto.SBJ2{1,sub_ind}.HIST_event_info{1,sei} = et.HIST_event_info{1, sei};
    end
end
SBJ = yoto.SBJ2;

save([resi_path 'SBJ_structure2.mat'],'SBJ');