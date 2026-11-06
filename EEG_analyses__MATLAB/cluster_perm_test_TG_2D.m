function [sig_mask, cluster_pvals, obs_cluster_stats, null_dist] = ...
    cluster_perm_test_TG_2D(data_3d, times, varargin)
% CLUSTER_PERM_TEST_TG_2D  Non-parametric cluster-based permutation test for
%   2-D temporal generalization (TG) matrices, following Maris & Oostenveld
%   (2007, J. Neurosci. Methods) and the multilevel-model variant applied to
%   continuous 2-D temporal spaces in EEG decoding studies.
%
%   Tests whether the group-mean TG map differs from zero using a one-sample
%   sign-permutation procedure across subjects in the 2-D (train × test)
%   time space.
%
%   USAGE
%     sig_mask = cluster_perm_test_TG_2D(data_3d, times)
%     [sig_mask, cluster_pvals, obs_cluster_stats, null_dist] = ...
%         cluster_perm_test_TG_2D(data_3d, times, 'n_perm', 5000, ...
%                                  'tg_xlim', [-300 1000])
%
%   INPUTS
%     data_3d          [n_times × n_times × n_subj]  TG maps per subject.
%                      One-sample test against zero is performed.
%     times            [1 × n_times]  Time vector in ms corresponding to
%                      both dimensions of data_3d.
%
%   OPTIONAL NAME-VALUE ARGUMENTS
%     'n_perm'         Number of Monte Carlo sign-permutation iterations.
%                      Default: 1000.
%     'alpha_cluster'  Cluster-forming threshold (two-tailed t-test p-value).
%                      Default: 0.01.
%     'alpha_sig'      Cluster-significance threshold (two-tailed).
%                      Default: 0.01.
%     'tail'           Test tail: 'two' (default), 'right' (positive only),
%                      or 'left' (negative only).
%     'tg_xlim'        [1×2] time window in ms to restrict the analysis.
%                      Default: full range of times.
%     'connectivity'   Pixel connectivity for bwconncomp: 4 or 8.
%                      Default: 8 (includes diagonal neighbours).
%     'report_top_k'   Print top K observed clusters and p-values.
%                      Default: 5; set 0 to disable.
%     'verbose'        Print progress to console. Default: true.
%
%   OUTPUTS
%     sig_mask          [n_times × n_times] logical array. True at (train,test)
%                       points belonging to a statistically significant cluster.
%     cluster_pvals     Struct with fields:
%                         .pos  [1 × n_pos_clusters]  two-tailed cluster p-values
%                               for each positive cluster (cluster mass > 0).
%                         .neg  [1 × n_neg_clusters]  two-tailed cluster p-values
%                               for each negative cluster (cluster mass < 0).
%     obs_cluster_stats Struct with fields:
%                         .pos  Observed positive cluster masses (sum of t).
%                         .neg  Observed negative cluster masses (sum of t).
%     null_dist         Struct with fields:
%                         .absmax [1 × n_perm] Null distribution of the maximum
%                                absolute cluster mass per permutation.
%                         .pos    [1 × n_perm] Null distribution of maximum
%                                positive cluster mass (for diagnostics).
%                         .neg    [1 × n_perm] Null distribution of minimum
%                                negative cluster mass (for diagnostics).
%
%   ALGORITHM
%     1. Compute one-sample t-statistic at every (train, test) time point.
%     2. Threshold at the cluster-forming level (p < alpha_cluster) with the
%        selected test tail, then identify candidate clusters using connected
%        components (bwconncomp from the Image Processing Toolbox).
%     3. Cluster mass = sum of t-values within each connected region.
%     4. Build null distribution: for n_perm iterations, randomly flip the
%        sign of each subject's entire map (equivalent to shuffling the
%        implicit two-class label in a one-sample test), recompute t-map,
%        find max positive and min negative cluster masses.
%     5. Cluster p-value = proportion of null max-absolute cluster masses
%        at least as extreme as the observed absolute cluster mass.
%     6. Clusters with p < alpha_sig are marked as significant.
%
%   REFERENCES
%     Maris, E. & Oostenveld, R. (2007). Nonparametric statistical testing of
%     EEG- and MEG-data. J. Neurosci. Methods, 164, 177-190.
%
%   NOTES
%     - Requires the MATLAB Image Processing Toolbox (bwconncomp).
%     - For large matrices / many permutations, consider parfor in the loop.
%     - tg_xlim restricts both the cluster analysis and the returned sig_mask
%       region; the remainder of sig_mask is set to false.
%

% ── Input parsing ──────────────────────────────────────────────────────────
p = inputParser;
addRequired(p,  'data_3d');
addRequired(p,  'times');
addParameter(p, 'n_perm',        1000,  @(x) isscalar(x) && x > 0);
addParameter(p, 'alpha_cluster', 0.01,  @(x) isscalar(x) && x > 0 && x < 1);
addParameter(p, 'alpha_sig',     0.01,  @(x) isscalar(x) && x > 0 && x < 1);
addParameter(p, 'tail',          'two', @(x) ischar(x) || isstring(x));
addParameter(p, 'tg_xlim',       [],    @(x) isempty(x) || (numel(x)==2 && x(1)<x(2)));
addParameter(p, 'connectivity',  8,     @(x) ismember(x,[4 8]));
addParameter(p, 'report_top_k',  5,     @(x) isscalar(x) && x >= 0);
addParameter(p, 'verbose',       true,  @islogical);
parse(p, data_3d, times, varargin{:});

n_perm        = p.Results.n_perm;
alpha_cluster = p.Results.alpha_cluster;
alpha_sig     = p.Results.alpha_sig;
tail          = char(lower(string(p.Results.tail)));
tg_xlim       = p.Results.tg_xlim;
connectivity  = p.Results.connectivity;
report_top_k  = p.Results.report_top_k;
verbose       = p.Results.verbose;

if ~ismember(tail, {'two','right','left'})
    error('cluster_perm_test_TG_2D: ''tail'' must be ''two'', ''right'', or ''left''.');
end

% ── Validate inputs ────────────────────────────────────────────────────────
[n_t, n_t2, n_subj] = size(data_3d);
assert(n_t == n_t2, ...
    'cluster_perm_test_TG_2D: data_3d must be square in the time dimensions.');
assert(n_t == numel(times), ...
    'cluster_perm_test_TG_2D: numel(times) must equal size(data_3d, 1).');
assert(n_subj >= 2, ...
    'cluster_perm_test_TG_2D: need at least 2 subjects for a t-test.');

% ── Time-window restriction ────────────────────────────────────────────────
if isempty(tg_xlim)
    t_idx = 1:n_t;
else
    t_idx = find(times >= tg_xlim(1) & times <= tg_xlim(2));
    if isempty(t_idx)
        error('cluster_perm_test_TG_2D: tg_xlim [%g %g] yields no time points.', ...
              tg_xlim(1), tg_xlim(2));
    end
end
n_t_win  = numel(t_idx);
data_win = data_3d(t_idx, t_idx, :);   % [n_t_win × n_t_win × n_subj]

% ── Cluster-forming t-threshold ───────────────────────────────────────────
df       = n_subj - 1;
if strcmp(tail, 'two')
    t_thresh = tinv(1 - alpha_cluster / 2, df);
else
    t_thresh = tinv(1 - alpha_cluster, df);
end

if verbose
    fprintf('Cluster perm test: n_subj=%d | df=%d | tail=%s | t_thresh=%.3f | n_perm=%d\n', ...
            n_subj, df, tail, t_thresh, n_perm);
    if ~isempty(tg_xlim)
        fprintf('  Time window: [%g %g] ms  (%d time points)\n', ...
                tg_xlim(1), tg_xlim(2), n_t_win);
    end
end

% ── Observed t-map and clusters ────────────────────────────────────────────
tmap_obs = tmap_from_data(data_win, n_subj);
[CC_pos_obs, CC_neg_obs, pos_masses_obs, neg_masses_obs] = ...
    find_clusters(tmap_obs, t_thresh, connectivity, tail);

obs_cluster_stats.pos = pos_masses_obs;
obs_cluster_stats.neg = neg_masses_obs;

if verbose
    fprintf('  Observed: %d positive cluster(s), %d negative cluster(s)\n', ...
            numel(pos_masses_obs), numel(neg_masses_obs));
end

% ── Null distribution (sign-flipping Monte Carlo) ──────────────────────────
% For each permutation: randomly assign ±1 to each subject's map, recompute
% t-map, and record the cluster-mass extremum required by the selected tail.
% When no clusters exceed threshold the neutral value 0 is stored.
null_max_pos = zeros(1, n_perm);
null_min_neg = zeros(1, n_perm);
null_max_abs = zeros(1, n_perm);

for perm_i = 1:n_perm
    % Random sign vector: ±1, broadcast over [n_t_win × n_t_win × n_subj]
    signs     = reshape(randi(2, 1, n_subj) * 2 - 3, 1, 1, n_subj);
    data_perm = data_win .* signs;
    tmap_perm = tmap_from_data(data_perm, n_subj);

    [~, ~, pm_pos, pm_neg] = find_clusters(tmap_perm, t_thresh, connectivity, tail);

    null_max_pos(perm_i) = max([0, pm_pos]);   % max positive mass (0 if none)
    null_min_neg(perm_i) = min([0, pm_neg]);   % min negative mass (0 if none)
    null_max_abs(perm_i) = max([0, pm_pos, abs(pm_neg)]); % strict two-tailed null
end

null_dist.absmax = null_max_abs;
null_dist.pos = null_max_pos;
null_dist.neg = null_min_neg;

% ── Cluster p-values and significance mask ─────────────────────────────────
sig_mask_win = false(n_t_win, n_t_win);
n_pos = numel(pos_masses_obs);
n_neg = numel(neg_masses_obs);

cluster_pvals.pos = ones(1, n_pos);
cluster_pvals.neg = ones(1, n_neg);

for k = 1:n_pos
    if strcmp(tail, 'two')
        p_val = mean(null_max_abs >= abs(pos_masses_obs(k)));
    else
        p_val = mean(null_max_pos >= pos_masses_obs(k));
    end
    cluster_pvals.pos(k) = p_val;
    if p_val < alpha_sig
        sig_mask_win(CC_pos_obs.PixelIdxList{k}) = true;
    end
end

for k = 1:n_neg
    if strcmp(tail, 'two')
        p_val = mean(null_max_abs >= abs(neg_masses_obs(k)));
    else
        p_val = mean(null_min_neg <= neg_masses_obs(k));
    end
    cluster_pvals.neg(k) = p_val;
    if p_val < alpha_sig
        sig_mask_win(CC_neg_obs.PixelIdxList{k}) = true;
    end
end

% Embed windowed mask back into full time grid
sig_mask = false(n_t, n_t);
sig_mask(t_idx, t_idx) = sig_mask_win;

if verbose
    fprintf('  Significant: %d pos cluster(s) | %d neg cluster(s)  (alpha_sig=%.2f)\n', ...
            sum(cluster_pvals.pos < alpha_sig), sum(cluster_pvals.neg < alpha_sig), alpha_sig);

    % Quick diagnostics: print strongest observed clusters and p-values.
    if report_top_k > 0
        report_cluster_table(pos_masses_obs, cluster_pvals.pos, 'pos', report_top_k);
        report_cluster_table(neg_masses_obs, cluster_pvals.neg, 'neg', report_top_k);
    end
end

end   % ── end main function ──────────────────────────────────────────────────


% ═══════════════════════════════════════════════════════════════════════════
%  Local helper functions
% ═══════════════════════════════════════════════════════════════════════════

function tmap = tmap_from_data(data_3d, n_subj)
% One-sample t-map of data_3d [... × n_subj] vs zero.
mn   = mean(data_3d, 3);
se   = std(data_3d, 0, 3) / sqrt(n_subj);
tmap = mn ./ (se + eps);   % eps avoids division by zero at flat time points
end


function [CC_pos, CC_neg, pos_masses, neg_masses] = find_clusters(tmap, t_thresh, connectivity, tail)
% Identify connected clusters above +t_thresh and below -t_thresh and
% compute their mass (sum of t-values).
pos_mask = false(size(tmap));
neg_mask = false(size(tmap));

if strcmp(tail, 'two') || strcmp(tail, 'right')
    pos_mask = tmap > t_thresh;
end
if strcmp(tail, 'two') || strcmp(tail, 'left')
    neg_mask = tmap < -t_thresh;
end

CC_pos = bwconncomp(pos_mask, connectivity);
CC_neg = bwconncomp(neg_mask, connectivity);

pos_masses = zeros(1, CC_pos.NumObjects);
for k = 1:CC_pos.NumObjects
    pos_masses(k) = sum(tmap(CC_pos.PixelIdxList{k}));
end

neg_masses = zeros(1, CC_neg.NumObjects);
for k = 1:CC_neg.NumObjects
    neg_masses(k) = sum(tmap(CC_neg.PixelIdxList{k}));   % will be negative
end
end


function report_cluster_table(masses, pvals, sign_label, top_k)
if isempty(masses)
    fprintf('  Top %s clusters: none above cluster-forming threshold.\n', sign_label);
    return;
end

[~, ord] = sort(abs(masses), 'descend');
n_show = min(top_k, numel(ord));
fprintf('  Top %s clusters (|mass|, p):\n', sign_label);
for i = 1:n_show
    idx = ord(i);
    fprintf('    #%d  mass=% .3f   p=%.4f\n', i, masses(idx), pvals(idx));
end
end
