addpath \\143.248.30.94\bmlsamba\kdj\repro_fmri\fmri_arbitration\modelRLsource

for i = 80:1:86
    switch i
        case 80 % 23 11 7
            param_regressor_type_cue={'SPE','RPE','MAXinvFano12','uncertaintyM1','MOTIV','Qsarsa','Qfwd','Qarb'};
            reg_type_go_first=[1 1.5 2];
        case 81 % 23 11 7
            param_regressor_type_cue={'SPE','RPE','MOTIV'};
            reg_type_go_first=[1 1.5 2];
        case 82 % 23 11 7
            param_regressor_type_cue={'SPE','RPE','MAXinvFano12','MOTIV','Qsarsa','Qfwd','Qarb'};
            reg_type_go_first=[1 1.5 2];
        case 83 % 23 11 7
            param_regressor_type_cue={'SPE','RPE','uncertaintyM1','MOTIV','MAXinvFano12','Qsarsa','Qfwd','Qarb'};
            reg_type_go_first=[1 1.5 2];
        case 84 % 23 11 7
            param_regressor_type_cue={'SPE','RPE','uncertaintyM1','MAXinvFano12','MOTIV','Qsarsa','Qfwd','Qarb'};
            reg_type_go_first=[1 1.5 2];
        case 85 % 23 11 7
            param_regressor_type_cue={'SPE','RPE','MAXinvFano12','MOTIV','uncertaintyM1','Qsarsa','Qfwd','Qarb'};
            reg_type_go_first=[1 1.5 2];
        case 86 % 23 11 7
            param_regressor_type_cue={'SPE','RPE','MAXinvFano12','uncertaintyM2','MOTIV','Qsarsa','Qfwd','Qarb'};
            reg_type_go_first=[1 1.5 2];
    end
    func_SIMUL_arbitration_regressor_generator_v10_kdj(i);
    func_job_dk_1thLevel()
    pause(60*10);
    job_dk_2ndLevel_passth()
    pause(60*1.5);
    copyfile([pwd '/fmri_arbitration/analysis_2nd_level'], [pwd '/fmri_arbitration/FINAL_BOSS_' list_dir '_tot'])
    
    
end
