% Fit model-1 RSA to single-trial multiclass probability profiles.
% The full class-probability matrix is treated as similarity. The reported
% cued models separate stimulus-response, rule, and dimension-general side
% information. Non-cued fits are computed but not used downstream; reported
% non-cued results come from model 3 (model 2 in the manuscript text).
% Author: Davide Gheza

function fits = fit_RSAmodels_multiclassRDM_sngltrl_m1__dim_gen_side__std(si,cls)

fitresiduals = 0; % use residuals from cued-model fitting to fit non-cued models
CLR = 1; % transform classification profile from [0 1] to centered log-ratio - to attenuate dependency in betas due to compositional data
local=0;

sn = [22,24:48];

if local==0
    load(['/home/gheza/EEG_data/RSA/distances_mltclass_sbj_' num2str(sn(si)) '.mat']);
else
    load(['/Users/davidegheza/Library/CloudStorage/OneDrive-Personal/Documenti/MATLAB/dg_workspace22/dimensions/analyses/preprocessed_data/rsa_data/RSA_pc/multiclass/singletrial/distances_mltclass_repCV_effHPT_bal_sbj_' num2str(sn(si)) '.mat']);
end
   keep.lda_prob_st(:,:,:,:) = cat(4,result.lda_prob_st(1).s,result.lda_prob_st(2).s,result.lda_prob_st(3).s,result.lda_prob_st(4).s);  
   keep.lda_acc_st(:,:,:) = cat(4,result.lda_acc_st(1).s,result.lda_acc_st(2).s,result.lda_acc_st(3).s,result.lda_acc_st(4).s);
   keep.y_test(:,:) = cat(2,squeeze(result.y_test(1).s),squeeze(result.y_test(2).s),squeeze(result.y_test(3).s),squeeze(result.y_test(4).s));
   keep.behav.b = result.behav;


% SR_model = eye(8);            
rule_model = [[1,1,zeros(1,6)]',[1,1,zeros(1,6)]', [0,0,1,1,zeros(1,4)]',[0,0,1,1,zeros(1,4)]', [zeros(1,4),ones(1,2),0,0]',[zeros(1,4),ones(1,2),0,0]', [zeros(1,6),1,1]',[zeros(1,6),1,1]'];
side_model = [repmat([repmat([1,-1],1,4); repmat([-1,1],1,4)]',1,4)];                       % dimension-general side
SR_model = rule_model.*side_model;
% formal way to standardize coefficients: (07/23/2024)
    % Calculate norms
    norm_rule_model = norm(rule_model, 'fro');
    norm_side_model = norm(side_model, 'fro');
    norm_SR_model = norm(SR_model, 'fro');

    % Define scaling factors
    scale_rule_model = norm_SR_model / norm_rule_model;  % normalize to SR_model 
    scale_side_model = norm_SR_model / norm_side_model;
    % Apply scaling
    rule_model = rule_model * scale_rule_model;
    side_model = side_model * scale_side_model;

column_map = [1,5,2,6,3,7,4,8];



    % LDA class probabilities (`cls` is no longer used; kept for existing calls)
    rdm = keep.lda_prob_st;
    rdm = permute(rdm, [4,3,2,1]); % ses, label_p, times, trials
    rdm = rdm(:,[1,5,2,6,3,7,4,8],:,:); % original order is ascending labels
    if CLR
        rdm = permute(rdm, [2,1,3,4]);
        sizes = size(rdm);
        rdm = reshape(rdm, 8, []);
        geo_mean = exp(mean(log(rdm),1, "omitnan"));
        rdm = log(rdm ./ geo_mean);
        rdm = reshape(rdm, sizes);
        rdm = permute(rdm, [2,1,3,4]);
    end
        
    % preallocate
    fits.slope_sr = nan(size(keep.lda_prob_st,[1,4,2]));
    fits.slope_ru = nan(size(keep.lda_prob_st,[1,4,2]));
    fits.slope_si = nan(size(keep.lda_prob_st,[1,4,2]));
    fits.intercept = nan(size(keep.lda_prob_st,[1,4,2]));
    fits.rsq = nan(size(keep.lda_prob_st,[1,4,2]));
    fits.slope_sr_t = nan(size(keep.lda_prob_st,[1,4,2]));      % added 01/18/2024
    fits.slope_ru_t = nan(size(keep.lda_prob_st,[1,4,2]));      % added 01/18/2024
    fits.slope_si_t = nan(size(keep.lda_prob_st,[1,4,2]));      % added 01/18/2024
    fits.slope_sr_dstr = nan(size(keep.lda_prob_st,[1,4,2]));
    fits.slope_ru_dstr = nan(size(keep.lda_prob_st,[1,4,2]));
    fits.slope_si_dstr = nan(size(keep.lda_prob_st,[1,4,2]));
    fits.slope_sr_dstr_t = nan(size(keep.lda_prob_st,[1,4,2]));      % added 01/18/2024
    fits.slope_ru_dstr_t = nan(size(keep.lda_prob_st,[1,4,2]));      % added 01/18/2024
    fits.slope_si_dstr_t = nan(size(keep.lda_prob_st,[1,4,2]));      % added 01/18/2024
    
    fits.slope_sr_dstr_con = nan(size(keep.lda_prob_st,[1,4,2]));
    fits.slope_ru_dstr_con = nan(size(keep.lda_prob_st,[1,4,2]));
    fits.slope_si_dstr_con = nan(size(keep.lda_prob_st,[1,4,2]));
    fits.slope_sr_dstr_inc = nan(size(keep.lda_prob_st,[1,4,2]));
    fits.slope_ru_dstr_inc = nan(size(keep.lda_prob_st,[1,4,2]));
    fits.slope_si_dstr_inc = nan(size(keep.lda_prob_st,[1,4,2]));
    fits.slope_sr_dstr_con_t = nan(size(keep.lda_prob_st,[1,4,2]));      % added 01/18/2024
    fits.slope_ru_dstr_con_t = nan(size(keep.lda_prob_st,[1,4,2]));      % added 01/18/2024
    fits.slope_si_dstr_con_t = nan(size(keep.lda_prob_st,[1,4,2]));      % added 01/18/2024
    fits.slope_sr_dstr_inc_t = nan(size(keep.lda_prob_st,[1,4,2]));      % added 01/18/2024
    fits.slope_ru_dstr_inc_t = nan(size(keep.lda_prob_st,[1,4,2]));      % added 01/18/2024
    fits.slope_si_dstr_inc_t = nan(size(keep.lda_prob_st,[1,4,2]));      % added 01/18/2024
     
    fits.slope_sr_singledstr = nan([size(keep.lda_prob_st,[1,4,2]),4]);  % added 01/19/2024:
    fits.slope_ru_singledstr = nan([size(keep.lda_prob_st,[1,4,2]),4]);
    fits.slope_si_singledstr = nan([size(keep.lda_prob_st,[1,4,2]),4]);
    fits.slope_sr_singledstr_t = nan([size(keep.lda_prob_st,[1,4,2]),4]);
    fits.slope_ru_singledstr_t = nan([size(keep.lda_prob_st,[1,4,2]),4]);
    fits.slope_si_singledstr_t = nan([size(keep.lda_prob_st,[1,4,2]),4]);

    for ses = 1:size(rdm,1)
        for st = 1:size(rdm,4)
            if ~isnan(keep.y_test(st,ses)) % go through all trials, for all sessions, although 3/4 are nan
                for t = 1:size(rdm,3)
                    tmp = squeeze(rdm(ses,:,t,st))'; % get 8 probabilities vector for a single trial, time, session, sbj
                    SR_model_st = SR_model(:,column_map==keep.y_test(st,ses));   % select model vector (column of model matrix) for cued-dim
                    rule_model_st = rule_model(:,column_map==keep.y_test(st,ses));
                    side_model_st = side_model(:,column_map==keep.y_test(st,ses));
                    lm = fitlm([SR_model_st, rule_model_st, side_model_st], tmp); %       ,'RobustOpts','on'
                    residuals = lm.Residuals.Raw;
                        slope_sr = lm.Coefficients.Estimate(2);
                        slope_ru = lm.Coefficients.Estimate(3);
                        slope_si = lm.Coefficients.Estimate(4);
                        intercept = lm.Coefficients.Estimate(1);
                        rsq = lm.Rsquared.Ordinary;
                        fits.slope_sr(st,ses,t) = slope_sr;
                        fits.slope_ru(st,ses,t) = slope_ru;
                        fits.slope_si(st,ses,t) = slope_si;
                        fits.intercept(st,ses,t) = intercept;
                        fits.rsq(st,ses,t) = rsq;
                        fits.slope_sr_t(st,ses,t) = lm.Coefficients.tStat(2);       % added 01/18/2024
                        fits.slope_ru_t(st,ses,t) = lm.Coefficients.tStat(3);
                        fits.slope_si_t(st,ses,t) = lm.Coefficients.tStat(4);       % added 01/18/2024
                    % non-cued dimensions representations (for each distractor, select model vector corresponding to left feature. NB: also cued-dimension label refers to left feature)
                    fns = fieldnames(keep.behav.b);
                    fns = fns(9:end); % dimensions' l r labels
                    slope_sr = nan(1,4);
                    slope_ru = nan(1,4);
                    slope_si = nan(1,4);
                    slope_sr_t = nan(1,4);
                    slope_ru_t = nan(1,4);
                    slope_si_t = nan(1,4);
                    for dim = 1:4
                        if keep.behav.b.congruency(st,dim)~=-1         % if not cued dim
                            dstr_lbl =  dim + (keep.behav.b.(fns{dim})(st,1)-1) * 4; % create label (dim and side, 1:8) for distractor - use left feature (as for creating labels for cued dimension in RSA_distances). NB: instructed feature is coded consistently across dims (forced side coding: either 1 or 2 for all dims)
                            SR_model_st_dstr = SR_model(:,column_map==dstr_lbl);
                            rule_model_st_dstr = rule_model(:,column_map==dstr_lbl);
                            side_model_st_dstr = side_model(:,column_map==dstr_lbl);
                            if fitresiduals
                                lm = fitlm([SR_model_st_dstr, rule_model_st_dstr, side_model_st_dstr], residuals); % nb: still fit model for 1 non-cued dim at a time - non-cued predictors across are also anti-correlated between each others
                            else
                                lm = fitlm([SR_model_st_dstr, rule_model_st_dstr, side_model_st_dstr], tmp);  % fit models for 1 non-cued dim only
                            end
%                             lm = fitlm([side_model_st, SR_model_st_dstr, rule_model_st_dstr, side_model_st_dstr], tmp); % fit models for 1 non-cued dim and also the side model for cued - but can't do: when congruent, the side models of cued and uncuded are identical, and one is not fitted. Get warning: Regression design matrix is rank deficient to within machine precision. 
                            slope_sr(dim) = lm.Coefficients.Estimate(2);
                            slope_ru(dim) = lm.Coefficients.Estimate(3);
                            slope_si(dim) = lm.Coefficients.Estimate(4);
                            slope_sr_t(dim) = lm.Coefficients.tStat(2);     % added 01/18/2024
                            slope_ru_t(dim) = lm.Coefficients.tStat(3);
                            slope_si_t(dim) = lm.Coefficients.tStat(4);
                        end
                    end
                    fits.slope_sr_dstr(st,ses,t) = mean(slope_sr,2,'omitnan');
                    fits.slope_ru_dstr(st,ses,t) = mean(slope_ru,2,'omitnan');
                    fits.slope_si_dstr(st,ses,t) = mean(slope_si,2,'omitnan');
                    fits.slope_sr_dstr_t(st,ses,t) = mean(slope_sr_t,2,'omitnan');      % added 01/18/2024
                    fits.slope_ru_dstr_t(st,ses,t) = mean(slope_ru_t,2,'omitnan');
                    fits.slope_si_dstr_t(st,ses,t) = mean(slope_si_t,2,'omitnan');
                    % different averages of slopes for congruent and incongruent non-cued dimensions
                    fits.slope_sr_dstr_con(st,ses,t) = mean(slope_sr(keep.behav.b.congruency(st,:)==1),2,'omitnan');
                    fits.slope_ru_dstr_con(st,ses,t) = mean(slope_ru(keep.behav.b.congruency(st,:)==1),2,'omitnan');
                    fits.slope_si_dstr_con(st,ses,t) = mean(slope_si(keep.behav.b.congruency(st,:)==1),2,'omitnan');
                    fits.slope_sr_dstr_inc(st,ses,t) = mean(slope_sr(keep.behav.b.congruency(st,:)==0),2,'omitnan');
                    fits.slope_ru_dstr_inc(st,ses,t) = mean(slope_ru(keep.behav.b.congruency(st,:)==0),2,'omitnan');
                    fits.slope_si_dstr_inc(st,ses,t) = mean(slope_si(keep.behav.b.congruency(st,:)==0),2,'omitnan');
                    fits.slope_sr_dstr_con_t(st,ses,t) = mean(slope_sr_t(keep.behav.b.congruency(st,:)==1),2,'omitnan');        % added 01/18/2024
                    fits.slope_ru_dstr_con_t(st,ses,t) = mean(slope_ru_t(keep.behav.b.congruency(st,:)==1),2,'omitnan');
                    fits.slope_si_dstr_con_t(st,ses,t) = mean(slope_si_t(keep.behav.b.congruency(st,:)==1),2,'omitnan');
                    fits.slope_sr_dstr_inc_t(st,ses,t) = mean(slope_sr_t(keep.behav.b.congruency(st,:)==0),2,'omitnan');
                    fits.slope_ru_dstr_inc_t(st,ses,t) = mean(slope_ru_t(keep.behav.b.congruency(st,:)==0),2,'omitnan');
                    fits.slope_si_dstr_inc_t(st,ses,t) = mean(slope_si_t(keep.behav.b.congruency(st,:)==0),2,'omitnan');
                    % save slopes of each non-cued, without averaging (01/19/2024)
                    fits.slope_sr_singledstr(st,ses,t,:) = slope_sr;
                    fits.slope_ru_singledstr(st,ses,t,:) = slope_ru;
                    fits.slope_si_singledstr(st,ses,t,:) = slope_si;
                    fits.slope_sr_singledstr_t(st,ses,t,:) = slope_sr_t;
                    fits.slope_ru_singledstr_t(st,ses,t,:) = slope_ru_t;
                    fits.slope_si_singledstr_t(st,ses,t,:) = slope_si_t;
                    
%                  fprintf('\tFitted st %g t %g\n', st, t)    
                end
            
            end
            
        end
         fprintf('\tFitted sbj %g ses %g\n', si, ses)
    end


if local==0
    fname = sprintf(('/scratch/gheza/RSA_fits/fits_singletrials_m1_sbj_%d.mat'), sn(si));
else
    if fitresiduals
        fname = sprintf(('/Users/davidegheza/Library/CloudStorage/OneDrive-Personal/Documenti/MATLAB/dg_workspace22/dimensions/analyses/preprocessed_data/rsa_data/fits_singletrials_m1_residualsNC_sbj_%d.mat'), sn(si));
    else
        fname = sprintf(('/Users/davidegheza/Library/CloudStorage/OneDrive-Personal/Documenti/MATLAB/dg_workspace22/dimensions/analyses/preprocessed_data/rsa_data/fits_singletrials_m1_sbj_%d.mat'), sn(si));
    end
end
save (fname, 'fits')

end
