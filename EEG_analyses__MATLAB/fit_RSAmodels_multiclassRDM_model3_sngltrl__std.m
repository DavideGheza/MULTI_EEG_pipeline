% Fit model-3 RSA to single-trial multiclass probability profiles.
% The full class-probability matrix is treated as similarity (unlike
% pairwise-classification accuracy, which is a dissimilarity). Model 3 fits
% feature-specific stimulus-response conjunctions for cued and non-cued
% dimensions. Side coding is omitted to avoid a singular design; `QRdec`
% optionally orthogonalizes the regression vectors.
% Author: Davide Gheza

function fits = fit_RSAmodels_multiclassRDM_model3_sngltrl__std(si,cls)

QRdec = 0; % apply QRdecomposition to orthogonalize the regression vectors
% stretch = 0; % transform classification profile from [0 1] to [-inf inf] -- work in progress with atanh   --> tmp = atanh(tmp-0.5)
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

         
rule_model = [[1,1,zeros(1,6)]',[1,1,zeros(1,6)]', [0,0,1,1,zeros(1,4)]',[0,0,1,1,zeros(1,4)]', [zeros(1,4),ones(1,2),0,0]',[zeros(1,4),ones(1,2),0,0]', [zeros(1,6),1,1]',[zeros(1,6),1,1]'];
side_model = [repmat([repmat([1,-1],1,4); repmat([-1,1],1,4)]',1,4)];
SR_model = rule_model.*side_model;
% standardize coefficient of SIDE based on SR (used only for QR dec). 
    % Calculate norms
    norm_side_model = norm(side_model, 'fro');
    norm_SR_model = norm(SR_model, 'fro');
    % Define scaling factors
    scale_side_model = norm_SR_model / norm_side_model;
    % Apply scaling
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
%     fits.slope_ru = nan(size(keep.lda_prob_st,[1,4,2]));
%     fits.slope_si = nan(size(keep.lda_prob_st,[1,4,2]));        % side model is dimension-general (cued and non-cued)
%     fits.intercept = nan(size(keep.lda_prob_st,[1,4,2]));
    fits.rsq = nan(size(keep.lda_prob_st,[1,4,2]));
    fits.slope_sr_t = nan(size(keep.lda_prob_st,[1,4,2]));      % added 01/18/2024
%     fits.slope_ru_t = nan(size(keep.lda_prob_st,[1,4,2]));      % added 01/18/2024
%     fits.slope_si_t = nan(size(keep.lda_prob_st,[1,4,2]));      % added 01/18/2024
    fits.slope_sr_dstr = nan(size(keep.lda_prob_st,[1,4,2]));
%     fits.slope_ru_dstr = nan(size(keep.lda_prob_st,[1,4,2]));
%     fits.slope_si_dstr = nan(size(keep.lda_prob_st,[1,4,2]));
    fits.slope_sr_dstr_t = nan(size(keep.lda_prob_st,[1,4,2]));      % added 01/18/2024
%     fits.slope_ru_dstr_t = nan(size(keep.lda_prob_st,[1,4,2]));      % added 01/18/2024
%     fits.slope_si_dstr_t = nan(size(keep.lda_prob_st,[1,4,2]));      % added 01/18/2024
    
    fits.slope_sr_dstr_con = nan(size(keep.lda_prob_st,[1,4,2]));
%     fits.slope_ru_dstr_con = nan(size(keep.lda_prob_st,[1,4,2]));
%     fits.slope_si_dstr_con = nan(size(keep.lda_prob_st,[1,4,2]));
    fits.slope_sr_dstr_inc = nan(size(keep.lda_prob_st,[1,4,2]));
%     fits.slope_ru_dstr_inc = nan(size(keep.lda_prob_st,[1,4,2]));
%     fits.slope_si_dstr_inc = nan(size(keep.lda_prob_st,[1,4,2]));
    fits.slope_sr_dstr_con_t = nan(size(keep.lda_prob_st,[1,4,2]));      % added 01/18/2024
%     fits.slope_ru_dstr_con_t = nan(size(keep.lda_prob_st,[1,4,2]));      % added 01/18/2024
%     fits.slope_si_dstr_con_t = nan(size(keep.lda_prob_st,[1,4,2]));      % added 01/18/2024
    fits.slope_sr_dstr_inc_t = nan(size(keep.lda_prob_st,[1,4,2]));      % added 01/18/2024
%     fits.slope_ru_dstr_inc_t = nan(size(keep.lda_prob_st,[1,4,2]));      % added 01/18/2024
%     fits.slope_si_dstr_inc_t = nan(size(keep.lda_prob_st,[1,4,2]));      % added 01/18/2024
     
    fits.slope_sr_singledstr = nan([size(keep.lda_prob_st,[1,4,2]),4]);  % added 01/19/2024:
%     fits.slope_ru_singledstr = nan([size(keep.lda_prob_st,[1,4,2]),4]);
%     fits.slope_si_singledstr = nan([size(keep.lda_prob_st,[1,4,2]),4]);
    fits.slope_sr_singledstr_t = nan([size(keep.lda_prob_st,[1,4,2]),4]);
%     fits.slope_ru_singledstr_t = nan([size(keep.lda_prob_st,[1,4,2]),4]);
%     fits.slope_si_singledstr_t = nan([size(keep.lda_prob_st,[1,4,2]),4]);

fits.R = nan([size(keep.lda_prob_st,[1,4]),5,5]);

    % non-cued dimensions representations (for each distractor, select model vector corresponding to left feature. NB: also cued-dimension label refers to left feature)
    fns = fieldnames(keep.behav.b);
    fns = fns(9:end); % dimensions' l r labels
    
    for ses = 1:size(rdm,1)
        for st = 1:size(rdm,4)
            if ~isnan(keep.y_test(st,ses)) % go through all trials, for all sessions, although 3/4 are nan
                SR_model_st = SR_model(:,column_map==keep.y_test(st,ses));   % select model vector (column of model matrix) for cued-dim
                rule_model_st = rule_model(:,column_map==keep.y_test(st,ses));
                side_model_st = side_model(:,column_map==keep.y_test(st,ses));
                SR_model_st_dstr = nan(8,4);
                for dim = 1:4   % prepare model vectors for the 3 non-cued dims
                    if keep.behav.b.congruency(st,dim)~=-1         % if not cued dim
                        dstr_lbl =  dim + (keep.behav.b.(fns{dim})(st,1)-1) * 4; % create label (dim and side, 1:8) for distractor - use left feature (as for creating labels for cued dimension in RSA_distances). NB: instructed feature is coded consistently across dims (forced side coding: either 1 or 2 for all dims)
                        SR_model_st_dstr(:,dim) = SR_model(:,column_map==dstr_lbl);
                    end
                end
                [~,cued_dim_idx] = find(isnan(SR_model_st_dstr));
                cued_dim_idx = cued_dim_idx(1);
                nc_idx = 1:4;
                nc_idx(cued_dim_idx) = [];
                if QRdec == 1
                    A = [side_model_st, SR_model_st, SR_model_st_dstr(:,nc_idx(1)), SR_model_st_dstr(:,nc_idx(2)), SR_model_st_dstr(:,nc_idx(3))];
                   [Q,R] = qr(A,'econ'); % order of vectors in A matrix matters!
                   for v=1:size(A,2)     % manually adjust sign to match A vectors (non-uniqueness of the QR decomposition)
                       first_nonzero_idx = find(Q(:, v) ~= 0, 1);
                       if sign(A(first_nonzero_idx,v)) ~= sign(Q(first_nonzero_idx,v))
                          Q(:,v) = -Q(:,v);
                          R(v,:) = -R(v,:);
                       end
                   end
                   fits.R(st,ses,:,:) = R; % save R trianguar matrix. Save once for trial and ses
                end
                for t = 1:size(rdm,3)
                    tmp = squeeze(rdm(ses,:,t,st))'; % get 8 probabilities vector for a single trial, time, session, sbj
                    if QRdec
                       tbl = table(Q(:,1), Q(:,2), Q(:,3), Q(:,4), Q(:,5), tmp, ... % matching order of QR decomposition!
                            'VariableNames', {'side_model_st', 'SR_model_st', 'SR_model_st_dstr1', 'SR_model_st_dstr2', 'SR_model_st_dstr3', 'tmp'});
                    else
                        tbl = table(side_model_st, SR_model_st, SR_model_st_dstr(:,nc_idx(1)), SR_model_st_dstr(:,nc_idx(2)), SR_model_st_dstr(:,nc_idx(3)), tmp,  ...
                            'VariableNames', {'side_model_st', 'SR_model_st', 'SR_model_st_dstr1', 'SR_model_st_dstr2', 'SR_model_st_dstr3', 'tmp'});
                    end
                    lm = fitlm(tbl, 'tmp ~ SR_model_st + SR_model_st_dstr1 + SR_model_st_dstr2 + SR_model_st_dstr3 -1'); % -1: specify model without intercept
%                         Coefficients = R \ lm.Coefficients.Estimate; % back-tranformed coefficients for unbiased effect sizes (use R to revert to original basis) 
%                         Coefficients_tStat = R \ lm.Coefficients.tStat;
                        Coefficients = lm.Coefficients.Estimate;
                        Coefficients_tStat = lm.Coefficients.tStat;
                        rsq = lm.Rsquared.Ordinary;
                        fits.slope_sr(st,ses,t) = Coefficients(1);
                        fits.rsq(st,ses,t) = rsq;
                        fits.slope_sr_t(st,ses,t) = Coefficients_tStat(1);

                        slope_sr = nan(1,4); % reset/allocate (to keep dims order, fill cued dim with nan)
                        slope_sr_t = nan(1,4); % reset/allocate

                        slope_sr(1,nc_idx(1)) = Coefficients(2);
                        slope_sr(1,nc_idx(2)) = Coefficients(3);
                        slope_sr(1,nc_idx(3)) = Coefficients(4);
                        slope_sr_t(1,nc_idx(1)) = Coefficients_tStat(2);
                        slope_sr_t(1,nc_idx(2)) = Coefficients_tStat(3);
                        slope_sr_t(1,nc_idx(3)) = Coefficients_tStat(4);



                    fits.slope_sr_dstr(st,ses,t) = mean(slope_sr,2,'omitnan');
                    fits.slope_sr_dstr_t(st,ses,t) = mean(slope_sr_t,2,'omitnan');
                    % different averages of slopes for congruent and incongruent non-cued dimensions
                    fits.slope_sr_dstr_con(st,ses,t) = mean(slope_sr(keep.behav.b.congruency(st,:)==1),2,'omitnan');
                    fits.slope_sr_dstr_inc(st,ses,t) = mean(slope_sr(keep.behav.b.congruency(st,:)==0),2,'omitnan');
                    fits.slope_sr_dstr_con_t(st,ses,t) = mean(slope_sr_t(keep.behav.b.congruency(st,:)==1),2,'omitnan');
                    fits.slope_sr_dstr_inc_t(st,ses,t) = mean(slope_sr_t(keep.behav.b.congruency(st,:)==0),2,'omitnan');
                    % save slopes of each non-cued, without averaging (01/19/2024)
                    fits.slope_sr_singledstr(st,ses,t,:) = slope_sr;
                    fits.slope_sr_singledstr_t(st,ses,t,:) = slope_sr_t;
                    
%                  fprintf('\tFitted st %g t %g\n', st, t)    
                end
            
            end
            
        end
%          fprintf('\tFitted sbj %g ses %g\n', si, ses)
    end


if local==0
    fname = sprintf(('/scratch/gheza/RSA_fits/fits_singletrials_m3_sbj_%d.mat'), sn(si));
else
    fname = sprintf(('/Users/davidegheza/Library/CloudStorage/OneDrive-Personal/Documenti/MATLAB/dg_workspace22/dimensions/analyses/preprocessed_data/rsa_data/fits_singletrials_m3_sbj_%d.mat'), sn(si));
end
save (fname, 'fits')

end
