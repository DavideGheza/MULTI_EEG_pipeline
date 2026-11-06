%% Conflict geometry — Encoding model figure panels
%
%  ── Similarity panels ────────────────────────────────────────────────────
%  Panel A:   12×12 axis-similarity heatmap (Cue / Location / Conflict domains;
%             averaged over the conflict window). Compact stylized layout with
%             per-cell dimension letters and block-domain labels.
%
%  Panel B:   Three 4×4 within-domain similarity blocks (Cue, Location, Conflict)
%             plotted as separate vertically stacked heatmaps, with shared colour
%             limits (p98 of off-diagonal values).
%
%  ── Orthogonality timecourse panels ──────────────────────────────────────
%  Panel C:   Within- vs across-dimension similarity over time for Location and
%             Conflict domains (group mean ± SEM; from the 16×16 RSA matrix).
%
%  Panel C.2: TG diagonal projections (t-train = t-test) of the fullmodel
%             temporal-generalisation matrices. Three stacked subplots:
%               (1) Location – cued:   feature processing (DSR) vs spatial
%                                      attention (DSI)  [cued trials]
%               (2) Location – uncued: feature processing (DSR) vs spatial
%                                      attention (DSI)  [noncued trials]
%               (3) Conflict:          dim-private (DSC) vs dim-shared (DGC)
%             Y-axis matched to the colour range of the corresponding TG heatmap.
%
%  ── Temporal-generalisation heatmap panels ───────────────────────────────
%  Panel D:   Four 2-D TG heatmaps (one per conflict distractor dimension) for
%             the dim-private conflict projection (noncued).
%
%  Panel D.2: Dimension-averaged TG for conflict: dim-private (DSC, left) and
%             dim-shared (DGC, right); shared symmetric colour limits.
%             Optional cluster-based permutation contours.
%
%  Panel D.3: Dimension-averaged TG for location (noncued): feature-processing
%             (DSR, left) and spatial attention (DSI, right); shared colour limits.
%             Optional cluster-based permutation contours.
%
%  Panel D.4: As D.3 but for cued trials (no distractor dimension).
%             Optional cluster-based permutation contours.
%
%  ── Conflict trajectory panels ───────────────────────────────────────────
%  Panel E:   3-D conflict trajectory per distractor dimension (2×2 subplots).
%             Axes: time (x), dim-shared projection (y), dim-private (z).
%             Trajectory coloured by time; SD ellipsoids across subjects.
%
%  Panel F:   3-D conflict trajectory averaged across distractor dimensions
%             (single plot). Same axis convention as Panel E.
%
%  Panel F.2: 2-D conflict trajectories (no SD ellipsoids, time-coded colour).
%             Figure 1: dimension-averaged (single plot).
%             Figure 2: per-distractor dimension (2×2 subplots).
%
%  ── Cluster permutation settings (D.2 / D.3 / D.4) ──────────────────────
%    run_cluster_perm     – toggle on/off
%    n_cluster_perm       – number of sign-permutation iterations
%    cluster_tail         – 'two' | 'right' | 'left'
%    cluster_alpha_form   – cluster-forming threshold (p-value)
%    cluster_alpha_sig    – cluster significance threshold (p-value)
%
%  Requires:
%    encoded_similarity  [16×16×n_times]   from encoding_fullmodel_repCV__similarity_HPC
%    cued/noncued structs                  from encoding_fullmodel_repCV__tempgeneralization_HPC
%    times, sbjIDX, sets already in workspace (re-uses your existing loading code)
%    cluster_perm_test_TG_2D.m             on the MATLAB path
%
%  Run after loading settings exactly as in plot_encoding_tempgeneralization.m
%
%  Author: Davide Gheza


%% 1. Load settings
% Keep preloaded workspace variables when available (sets, times, sbjIDX).
% This allows running the script as part of a larger analysis pipeline.

%Select onset or response events
onset = true;
if onset
    suffix = '';
else
    suffix = 'response_';   % response-locked
end

% Fallback subject list if not provided upstream.
if ~exist('sbjIDX','var') || isempty(sbjIDX)
    sbjIDX = [22,24:48];
end

% Fallback time vector if not provided upstream.
if ~exist('times','var') || isempty(times)
    times = -800:10:1498; % ms vector (230 points)
end

% Fallback settings load if not provided upstream.
if ~exist('sets','var') || ~isstruct(sets)
    if ispc
        load('C:\Users\gdavide\OneDrive\Documenti\MATLAB\dg_workspace22\dimensions\analyses\scripts_and_settings\preprocessing_settings.mat');
    else
        load('/Users/davidegheza/Library/CloudStorage/OneDrive-Personal/Documenti/MATLAB/dg_workspace22/dimensions/analyses/scripts_and_settings/preprocessing_settings_mac.mat');
    end
end

% Ensure RSA path is present.
sets.rsa_data_dir = [sets.main_saving_dir 'rsa_data/RSA_pc/multiclass/singletrial/'];



%% ── shared settings ──────────────────────────────────────────────────────

times      = -800:10:1498;          % ms vector (230 points)
n_times    = length(times);
nSubj      = length(sbjIDX);
tick_ms    = -800:300:1400;

% Domain index mapping (matches encoding_fullmodel_repCV__similarity_HPC)
cueIdx      = 1:4;
locIdx      = 5:8;
confIdx     = 9:12;
prevconfIdx = 13:16;

% Conflict window of interest (ms)
if onset
    conf_win_ms = [050 450];% Stimulus locked
else
    conf_win_ms = [-200 200];% Response locked
end
conf_win_t  = find(times >= conf_win_ms(1) & times <= conf_win_ms(2));

% Colours (consistent with plot_encoding_similarity.m)
cueCol      = [0.00 0.45 0.74];
locCol      = [0.85 0.33 0.10];
confCol     = [0.47 0.67 0.19];
prevconfCol = [0.49 0.18 0.56];
crossCol    = [0.50 0.50 0.50];
dimCols     = [0.12 0.47 0.71;   % colour dim
               0.20 0.63 0.17;   % shape dim
               0.89 0.10 0.11;   % edge dim
               1.00 0.50 0.00];  % pattern dim

fontsize = 14;


%% ── helper: diverging red-blue colormap ──────────────────────────────────


function cmap = redblue_colormap(n)
if nargin < 1, n = 256; end
half = floor(n/2);
% Diverging map: blue -> white -> red
neg = [linspace(0.00,1.00,half)', linspace(0.00,1.00,half)', ones(half,1)];
pos = [ones(n-half,1), linspace(1.00,0.00,n-half)', linspace(1.00,0.00,n-half)'];
cmap = [neg; pos];
end


%% ── helper: RdBu colormap (blue→white→red, ColorBrewer 11-stop) ─────────────
% Diverging: blue = negative projection, white = zero, red = positive.
% Center (white) maps to 0 when clim is symmetric around 0.

function cmap = rdBu(n)
if nargin < 1, n = 256; end
base = [ ...
      5  48  97; ...
     33 102 172; ...
     67 147 195; ...
    146 197 222; ...
    209 229 240; ...
    247 247 247; ...
    253 219 199; ...
    244 165 130; ...
    214  96  77; ...
    178  24  43; ...
    103   0  31] / 255;
x    = linspace(1, size(base,1), n);
cmap = interp1(1:size(base,1), base, x);
end


%% ── load similarity data (as in plot_encoding_similarity.m) ──────────────

groupS = nan(nSubj, 16, 16, n_times);

for sn = 1:nSubj
    load(fullfile(sets.rsa_data_dir, ...
        ['/encoding_fullmodel_similarity/encoding_fullmodel_similarity_' ...
         suffix 'sbj_' num2str(sbjIDX(sn)) '.mat']));
    groupS(sn,:,:,:) = encoded_similarity;
end

meanS = squeeze(mean(groupS, 1, 'omitnan'));   % [16×16×n_times]

%% ── load temporal-generalization data (as in plot_encoding_tempgeneralization.m) ──

% preallocate: 4 PCs × 230 train × 230 test × 26 subjects
noncued_DSC_all = nan(4, n_times, n_times, 4, nSubj);   % dim-private conflict
noncued_DGC_all = nan(4, n_times, n_times, 4, nSubj);   % dim-shared conflict
noncued_DSR_all = nan(4, n_times, n_times, 4, nSubj);   % dim-private location/side
noncued_DSI_all = nan(4, n_times, n_times, 4, nSubj);   % dim-shared location/side
% Cued SR/SI are single-dimension projections: [4 PCs x time x time]
cued_DSR_all    = nan(4, n_times, n_times, nSubj);      % dim-private location/side (cued)
cued_DSI_all    = nan(4, n_times, n_times, nSubj);      % dim-shared location/side (cued)

for sn = 1:nSubj
    load(fullfile(sets.rsa_data_dir, ...
        ['/encoding_fullmodel_proj_wtempgen/encoding_fullmodel_' suffix 'sbj_' num2str(sbjIDX(sn)) '.mat']));

    noncued_DSC_all(:,:,:,:,sn) = noncued.encoded_dim_spec_conf;
    noncued_DGC_all(:,:,:,:,sn) = noncued.encoded_dim_gen_conf;
    noncued_DSR_all(:,:,:,:,sn) = noncued.encoded_dim_spec_SR;
    noncued_DSI_all(:,:,:,:,sn) = noncued.encoded_dim_gen_SI;
    cued_DSR_all(:,:,:,sn)      = cued.encoded_dim_spec_SR;
    cued_DSI_all(:,:,:,sn)      = cued.encoded_dim_gen_SI;
end

% average over PCs (Neu1 + Neu2 = indices 1 and 4) and subjects
PC_idx = [1 4];
conf_DSC = squeeze(mean(mean(noncued_DSC_all(PC_idx,:,:,:,:), 1, 'omitnan'), 5, 'omitnan'));
%  → [n_times × n_times × 4 distractors]
conf_DGC = squeeze(mean(mean(noncued_DGC_all(PC_idx,:,:,:,:), 1, 'omitnan'), 5, 'omitnan'));
loc_DSR = squeeze(mean(mean(noncued_DSR_all(PC_idx,:,:,:,:), 1, 'omitnan'), 5, 'omitnan'));
loc_DSI = squeeze(mean(mean(noncued_DSI_all(PC_idx,:,:,:,:), 1, 'omitnan'), 5, 'omitnan'));
loc_cued_DSR = squeeze(mean(mean(cued_DSR_all(PC_idx,:,:,:), 1, 'omitnan'), 4, 'omitnan'));
loc_cued_DSI = squeeze(mean(mean(cued_DSI_all(PC_idx,:,:,:), 1, 'omitnan'), 4, 'omitnan'));


%% ══════════════════════════════════════════════════════════════════════════
%% PANEL A  –  12×12 similarity heatmap averaged over conflict window
%% ══════════════════════════════════════════════════════════════════════════
inches_out_a = [4 3];

S_win = squeeze(mean(meanS(1:12,1:12,conf_win_t), 3));   % [12×12], drop Prev.C domain

% domain boundary lines at 4 and 8
boundaries = [4.5, 8.5];

% Shared colour limit (p98 of off-diagonal values in 3 domain blocks) —
% computed once here so both Panel A and Panel B use the same range.
pA_blockIdx = {cueIdx, locIdx, confIdx};
panelD_vals = [];
for pA_b = 1:3
    pA_tmp = S_win(pA_blockIdx{pA_b}, pA_blockIdx{pA_b});
    pA_vals = pA_tmp(~eye(size(pA_tmp,1)));
    panelD_vals = [panelD_vals; pA_vals(:)];
end
panelD_vals = panelD_vals(~isnan(panelD_vals));
if isempty(panelD_vals)
    panelD_cmax = 0.2;
else
    panelD_cmax = prctile(abs(panelD_vals), 98);
    if panelD_cmax == 0, panelD_cmax = max(abs(panelD_vals)); end
    if panelD_cmax == 0, panelD_cmax = 0.2; end
end

figure('Name','Panel A – Similarity heatmap', 'Color','w', ...
       'Position',[100 100 340 320]);

imagesc(S_win);
axis square;
colormap(redblue_colormap());          % helper defined below; diverging -1->0->+1
cb = colorbar;
cb.Label.String = '';
caxis([-panelD_cmax panelD_cmax]);

% Stylized compact framework for tiny print: dimension letters per cell,
% block-domain labels in colour, strong domain separators.
dimLetters  = {'C','S','E','P'};
blockCols   = {cueCol, locCol, confCol};  % Cue / Location / Conflict
blockStarts = [1, 5, 9];                  % first index of each 4-cell block

hold on;
% Subtle inter-cell grid
for g = 0.5:1:12.5
    xline(g, '-', 'Color', [1 1 1]*0.82, 'LineWidth', 0.35);
    yline(g, '-', 'Color', [1 1 1]*0.82, 'LineWidth', 0.35);
end
% Strong domain separator lines
for b = boundaries
    xline(b, 'k-', 'LineWidth', 1.4);
    yline(b, 'k-', 'LineWidth', 1.4);
end

% Per-cell dimension letter on the x-axis (bottom) and y-axis (left)
ax = gca;
for bl = 1:3
    col = blockCols{bl};
    for d = 1:4
        ci = blockStarts(bl) + d - 1;    % column / row index (1..12)
        % X tick label (bottom)
        text(ax, ci, 13.0, dimLetters{d}, ...
            'HorizontalAlignment','center', 'VerticalAlignment','top', ...
            'FontSize', 6.6, 'Color', col, 'Clipping','off');
        % Y tick label (left)
        text(ax, 0.3, ci, dimLetters{d}, ...
            'HorizontalAlignment','right', 'VerticalAlignment','middle', ...
            'FontSize', 6.6, 'Color', col, 'Clipping','off');
    end
end

% Block-domain name labels above x-axis and to the left of y-axis
blockNames  = {'Cue', 'Loc', 'Conf'};
blockCentre = [2.5, 6.5, 10.5];
for bl = 1:3
    col = blockCols{bl};
    % above plot (x)
    text(ax, blockCentre(bl), 14.2, blockNames{bl}, ...
        'HorizontalAlignment','center', 'VerticalAlignment','top', ...
        'FontSize', 9.6, 'FontWeight','bold', 'Color', col, 'Clipping','off');
    % left of plot (y)
    text(ax, -1.1, blockCentre(bl), blockNames{bl}, ...
        'HorizontalAlignment','center', 'VerticalAlignment','middle', ...
        'FontSize', 9.6, 'FontWeight','bold', 'Color', col, 'Clipping','off', ...
        'Rotation', 90);
end

set(ax, 'XTick', [], 'YTick', [], 'FontSize', 7, 'TickLength', [0 0]);
title('Axis similarity', 'FontSize', fontsize-1);

hold off;
box on;

% Save out pdf.
set(gcf, 'Units', 'inches');
set(gcf, 'Position', [1 1 inches_out_a]);          % on screen
set(gcf, 'PaperUnits', 'inches');
set(gcf, 'PaperPosition', [0 0 inches_out_a]);     % for printed/exported output
set(gcf, 'PaperSize', inches_out_a);        

print(gcf, '-dpdf', 'cossim_full_matrix.pdf');

%% ══════════════════════════════════════════════════════════════════════════
%% PANEL B  –  Panel A split into domain blocks (4x4 each)
%% ══════════════════════════════════════════════════════════════════════════

blockIdx = {cueIdx, locIdx, confIdx};
blockTitles = {'Cue block', 'Location block', 'Conflict block'};
blockTitleCols = {cueCol, locCol, confCol};
blockAxisLabels = {'C','S','E','P'};

figure('Name','Panel D – Domain block similarity matrices', 'Color','w', ...
    'Position',[100 100 980 420]);
sgtitle(sprintf('Within-domain similarity blocks (avg %d–%d ms)', conf_win_ms(1), conf_win_ms(2)), ...
    'FontSize', fontsize+2);

% Typical diverging range for signed values: shared, symmetric around zero.
% Use robust limits from all off-diagonal entries across the 3 domain blocks.
% (panelD_cmax computed above, before Panel A, so both panels share the same range.)

for b = 1:3
    subplot(3,1,b);
    blockMat = S_win(blockIdx{b}, blockIdx{b});
    imagesc(blockMat);
    axis square;
    colormap(redblue_colormap());
    caxis([-panelD_cmax panelD_cmax]);
    colorbar;
    xticks(1:4);
    yticks(1:4);
    xticklabels(blockAxisLabels);
    yticklabels(blockAxisLabels);
    set(gca, 'FontSize', fontsize-2);
    title(blockTitles{b}, 'FontSize', fontsize, 'Color', blockTitleCols{b});
    hold on;
    for r = 1:4
        for c = 1:4
            v = blockMat(r,c);
            if isnan(v)
                txt = 'NaN';
                txt_col = [0 0 0];
            else
                txt = sprintf('%.3f', v);
                if abs(v) > 0.55*panelD_cmax
                    txt_col = [1 1 1];
                else
                    txt_col = [0 0 0];
                end
            end
            text(c, r, txt, 'HorizontalAlignment','center', 'FontSize', fontsize-5, ...
                 'FontWeight','bold', 'Color', txt_col);
        end
    end
    hold off;
end

%% ══════════════════════════════════════════════════════════════════════════
%% PANEL C  –  Within-vs-across similarity over time (Conflict and Location)
%% ══════════════════════════════════════════════════════════════════════════
inches_out_c = [4 6];

pairs = nchoosek(1:4, 2);                 % [6x2] all unique pairs
nPairs = size(pairs, 1);

domain_names = {'Location', 'Conflict'};
domain_indices = {locIdx, confIdx};
within_col = [0.85 0.15 0.15];
across_col = [0.25 0.25 0.25];

figure('Name','Panel C – Orthogonality summaries', 'Color','w', ...
       'Position',[100 100 980 420]);

axC = gobjects(1,2);
for dom = 1:2
    idx = domain_indices{dom};

    cross_sn = nan(nSubj, nPairs, n_times);
    self_sn  = nan(nSubj, 4, n_times);

    for sn = 1:nSubj
        S_sn = squeeze(groupS(sn,:,:,:));   % [16x16xn_times]
        for p = 1:nPairs
            i = idx(pairs(p,1));
            j = idx(pairs(p,2));
            cross_sn(sn,p,:) = squeeze(S_sn(i,j,:));
        end
        for d = 1:4
            self_sn(sn,d,:) = squeeze(S_sn(idx(d),idx(d),:));
        end
    end

    % Subject-level summaries before group averaging.
    cross_avg_sn = squeeze(mean(cross_sn, 2, 'omitnan')); % [nSubj x n_times]
    self_avg_sn  = squeeze(mean(self_sn,  2, 'omitnan')); % [nSubj x n_times]

    mean_cross = mean(cross_avg_sn, 1, 'omitnan');
    sem_cross  = std(cross_avg_sn, 0, 1, 'omitnan') / sqrt(nSubj);
    mean_self  = mean(self_avg_sn, 1, 'omitnan');
    sem_self   = std(self_avg_sn, 0, 1, 'omitnan') / sqrt(nSubj);

    axC(dom) = subplot(2,1,dom);
    hold on;
    title(sprintf('%s: within vs. across', domain_names{dom}), 'FontSize', fontsize);
    yline(0, 'k--', 'LineWidth', 1.2, 'HandleVisibility','off');

    hFill1 = fill([times fliplr(times)], ...
        [mean_self+sem_self, fliplr(mean_self-sem_self)], ...
        within_col, 'FaceAlpha', 0.20, 'EdgeColor','none');
    hFill1.HandleVisibility = 'off';
    hSelf = plot(times, mean_self, '-', 'Color', within_col, 'LineWidth', 2.2);

    hFill2 = fill([times fliplr(times)], ...
        [mean_cross+sem_cross, fliplr(mean_cross-sem_cross)], ...
        across_col, 'FaceAlpha', 0.20, 'EdgeColor','none');
    hFill2.HandleVisibility = 'off';
    hCross = plot(times, mean_cross, '-', 'Color', across_col, 'LineWidth', 2.2);

    yl = ylim;
    patch([conf_win_ms(1) conf_win_ms(2) conf_win_ms(2) conf_win_ms(1)], ...
          [yl(1) yl(1) yl(2) yl(2)], ...
          [0.8 0.8 0.8], 'FaceAlpha', 0.35, 'EdgeColor','none', 'HandleVisibility','off');

    xlabel('Time (ms)', 'FontSize', fontsize);
    ylabel('Cosine similarity', 'FontSize', fontsize);
    if dom==2
        legend([hSelf hCross], {'Within-dimension (on-diagonal)', 'Between-dimension (off-diagonal)'}, ...
            'Location','best', 'FontSize', fontsize-4);
    end
    xlim([-300 1000]);
    set(gca, 'FontSize', fontsize-2);
    box off;
    hold off;
end

% Enforce identical y-limits across conflict/location subplots.
yL1 = ylim(axC(1));
yL2 = ylim(axC(2));
yCommon = [min([yL1(1), yL2(1)]), max([yL1(2), yL2(2)])];
ylim(axC(1), yCommon);
ylim(axC(2), yCommon);

% sgtitle('Panel C: Within-vs-across similarity (orthogonality check)', 'FontSize', fontsize+2);

% Save out pdf.
set(gcf, 'Units', 'inches');
set(gcf, 'Position', [1 1 inches_out_c]);          % on screen
set(gcf, 'PaperUnits', 'inches');
set(gcf, 'PaperPosition', [0 0 inches_out_c]);     % for printed/exported output
set(gcf, 'PaperSize', inches_out_c);        

print(gcf, '-dpdf', 'cossim_traj.pdf');

%% ══════════════════════════════════════════════════════════════════════════
%% PANEL C.2  –  Within-vs-across from fullmodel TG projection (t_train=t_test)
%% ══════════════════════════════════════════════════════════════════════════
inches_out_c2 = [4 9]; % paper size of plots

diag_lin_idx = sub2ind([n_times n_times], 1:n_times, 1:n_times);

% Subject-level diagonal trajectories from fullmodel tempgen projections.
conf_DSC_bySubj = squeeze(mean(noncued_DSC_all(PC_idx,:,:,:,:), 1, 'omitnan')); % [time x time x dist x subj]
conf_DGC_bySubj = squeeze(mean(noncued_DGC_all(PC_idx,:,:,:,:), 1, 'omitnan')); % [time x time x dist x subj]
loc_DSR_bySubj  = squeeze(mean(noncued_DSR_all(PC_idx,:,:,:,:), 1, 'omitnan')); % [time x time x dist x subj]
loc_DSI_bySubj  = squeeze(mean(noncued_DSI_all(PC_idx,:,:,:,:), 1, 'omitnan')); % [time x time x dist x subj]
loc_cued_DSR_bySubj = squeeze(mean(cued_DSR_all(PC_idx,:,:,:), 1, 'omitnan'));  % [time x time x subj]
loc_cued_DSI_bySubj = squeeze(mean(cued_DSI_all(PC_idx,:,:,:), 1, 'omitnan'));  % [time x time x subj]

conf_within_sn = nan(nSubj, n_times);
conf_across_sn = nan(nSubj, n_times);
loc_within_sn  = nan(nSubj, n_times);
loc_across_sn  = nan(nSubj, n_times);
loc_cued_within_sn = nan(nSubj, n_times);
loc_cued_across_sn = nan(nSubj, n_times);

for sn = 1:nSubj
    conf_within_dist = nan(4, n_times);
    conf_across_dist = nan(4, n_times);
    loc_within_dist  = nan(4, n_times);
    loc_across_dist  = nan(4, n_times);

    for dist = 1:4
        tmp = squeeze(conf_DSC_bySubj(:,:,dist,sn));
        conf_within_dist(dist,:) = tmp(diag_lin_idx);
        tmp = squeeze(conf_DGC_bySubj(:,:,dist,sn));
        conf_across_dist(dist,:) = tmp(diag_lin_idx);

        tmp = squeeze(loc_DSR_bySubj(:,:,dist,sn));
        loc_within_dist(dist,:) = tmp(diag_lin_idx);
        tmp = squeeze(loc_DSI_bySubj(:,:,dist,sn));
        loc_across_dist(dist,:) = tmp(diag_lin_idx);
    end

    conf_within_sn(sn,:) = mean(conf_within_dist, 1, 'omitnan');
    conf_across_sn(sn,:) = mean(conf_across_dist, 1, 'omitnan');
    loc_within_sn(sn,:)  = mean(loc_within_dist, 1, 'omitnan');
    loc_across_sn(sn,:)  = mean(loc_across_dist, 1, 'omitnan');

    tmp = squeeze(loc_cued_DSR_bySubj(:,:,sn));
    loc_cued_within_sn(sn,:) = tmp(diag_lin_idx);
    tmp = squeeze(loc_cued_DSI_bySubj(:,:,sn));
    loc_cued_across_sn(sn,:) = tmp(diag_lin_idx);
end

mean_conf_within = mean(conf_within_sn, 1, 'omitnan');
sem_conf_within  = std(conf_within_sn, 0, 1, 'omitnan') / sqrt(nSubj);
mean_conf_across = mean(conf_across_sn, 1, 'omitnan');
sem_conf_across  = std(conf_across_sn, 0, 1, 'omitnan') / sqrt(nSubj);

mean_loc_within = mean(loc_within_sn, 1, 'omitnan');
sem_loc_within  = std(loc_within_sn, 0, 1, 'omitnan') / sqrt(nSubj);
mean_loc_across = mean(loc_across_sn, 1, 'omitnan');
sem_loc_across  = std(loc_across_sn, 0, 1, 'omitnan') / sqrt(nSubj);

mean_loc_cued_within = mean(loc_cued_within_sn, 1, 'omitnan');
sem_loc_cued_within  = std(loc_cued_within_sn, 0, 1, 'omitnan') / sqrt(nSubj);
mean_loc_cued_across = mean(loc_cued_across_sn, 1, 'omitnan');
sem_loc_cued_across  = std(loc_cued_across_sn, 0, 1, 'omitnan') / sqrt(nSubj);

% Match y-ranges in C.2 to the corresponding TG heatmap value ranges.
c2_conf_avg_map = mean(conf_DSC, 3, 'omitnan');
c2_conf_gen_map = mean(conf_DGC, 3, 'omitnan');
c2_conf_ylim = max(abs([c2_conf_avg_map(:); c2_conf_gen_map(:)]));
if ~isfinite(c2_conf_ylim) || c2_conf_ylim == 0, c2_conf_ylim = 0.1; end

c2_loc_uncued_spec_map = mean(loc_DSR, 3, 'omitnan');
c2_loc_uncued_gen_map  = mean(loc_DSI, 3, 'omitnan');
c2_loc_uncued_ylim = max(abs([c2_loc_uncued_spec_map(:); c2_loc_uncued_gen_map(:)]));
if ~isfinite(c2_loc_uncued_ylim) || c2_loc_uncued_ylim == 0, c2_loc_uncued_ylim = 0.1; end

c2_loc_cued_ylim = max(abs([loc_cued_DSR(:); loc_cued_DSI(:)]));
if ~isfinite(c2_loc_cued_ylim) || c2_loc_cued_ylim == 0, c2_loc_cued_ylim = 0.1; end

% Use the same displayed window as TG panels and ensure y-limits also cover
% the visible mean±SEM trajectories to avoid clipping at the top.
c2_xwin = [-300 1000];
c2_tmask = times >= c2_xwin(1) & times <= c2_xwin(2);

c2_loc_cued_curve_ylim = max(abs([
    mean_loc_cued_within(c2_tmask) + sem_loc_cued_within(c2_tmask), ...
    mean_loc_cued_within(c2_tmask) - sem_loc_cued_within(c2_tmask), ...
    mean_loc_cued_across(c2_tmask) + sem_loc_cued_across(c2_tmask), ...
    mean_loc_cued_across(c2_tmask) - sem_loc_cued_across(c2_tmask)
]));
if isfinite(c2_loc_cued_curve_ylim), c2_loc_cued_ylim = max(c2_loc_cued_ylim, c2_loc_cued_curve_ylim); end

c2_loc_uncued_curve_ylim = max(abs([
    mean_loc_within(c2_tmask) + sem_loc_within(c2_tmask), ...
    mean_loc_within(c2_tmask) - sem_loc_within(c2_tmask), ...
    mean_loc_across(c2_tmask) + sem_loc_across(c2_tmask), ...
    mean_loc_across(c2_tmask) - sem_loc_across(c2_tmask)
]));
if isfinite(c2_loc_uncued_curve_ylim), c2_loc_uncued_ylim = max(c2_loc_uncued_ylim, c2_loc_uncued_curve_ylim); end

c2_conf_curve_ylim = max(abs([
    mean_conf_within(c2_tmask) + sem_conf_within(c2_tmask), ...
    mean_conf_within(c2_tmask) - sem_conf_within(c2_tmask), ...
    mean_conf_across(c2_tmask) + sem_conf_across(c2_tmask), ...
    mean_conf_across(c2_tmask) - sem_conf_across(c2_tmask)
]));
if isfinite(c2_conf_curve_ylim), c2_conf_ylim = max(c2_conf_ylim, c2_conf_curve_ylim); end

figure('Name','Panel C.2 – Orthogonality from TG diagonal projections', 'Color','w', ...
       'Position',[100 420 980 620]);

axC2 = gobjects(1,3);

% Location (cued) on top
axC2(1) = subplot(3,1,1);
hold on;
title({
    'Location (cued):'
    'private vs. shared'
}, 'FontSize', fontsize);
yline(0, 'k--', 'LineWidth', 1.2, 'HandleVisibility','off');
hFill = fill([times fliplr(times)], [mean_loc_cued_within+sem_loc_cued_within, fliplr(mean_loc_cued_within-sem_loc_cued_within)], ...
    within_col, 'FaceAlpha', 0.20, 'EdgeColor','none');
hFill.HandleVisibility = 'off';
hLocCuedWithin = plot(times, mean_loc_cued_within, '-', 'Color', within_col, 'LineWidth', 2.2);
hFill = fill([times fliplr(times)], [mean_loc_cued_across+sem_loc_cued_across, fliplr(mean_loc_cued_across-sem_loc_cued_across)], ...
    across_col, 'FaceAlpha', 0.20, 'EdgeColor','none');
hFill.HandleVisibility = 'off';
hLocCuedAcross = plot(times, mean_loc_cued_across, '-', 'Color', across_col, 'LineWidth', 2.2);
yl = [-c2_loc_cued_ylim c2_loc_cued_ylim];
ylim(yl);
patch([conf_win_ms(1) conf_win_ms(2) conf_win_ms(2) conf_win_ms(1)], [yl(1) yl(1) yl(2) yl(2)], ...
    [0.8 0.8 0.8], 'FaceAlpha', 0.15, 'EdgeColor','none', 'HandleVisibility','off');
xlabel('Time (ms)', 'FontSize', fontsize);
ylabel('Projection', 'FontSize', fontsize);
legend([hLocCuedWithin hLocCuedAcross], {'Private', 'Shared'}, ...
    'Location','best', 'FontSize', fontsize-4);
xlim(c2_xwin);
set(gca, 'FontSize', fontsize-2);
box off;
hold off;

% Location (uncued) in the middle
axC2(2) = subplot(3,1,2);
hold on;
title({
    'Location (uncued):'
    'Private vs. Shared'
}, 'FontSize', fontsize);
yline(0, 'k--', 'LineWidth', 1.2, 'HandleVisibility','off');
hFill = fill([times fliplr(times)], [mean_loc_within+sem_loc_within, fliplr(mean_loc_within-sem_loc_within)], ...
    within_col, 'FaceAlpha', 0.20, 'EdgeColor','none');
hFill.HandleVisibility = 'off';
hLocWithin = plot(times, mean_loc_within, '-', 'Color', within_col, 'LineWidth', 2.2);
hFill = fill([times fliplr(times)], [mean_loc_across+sem_loc_across, fliplr(mean_loc_across-sem_loc_across)], ...
    across_col, 'FaceAlpha', 0.20, 'EdgeColor','none');
hFill.HandleVisibility = 'off';
hLocAcross = plot(times, mean_loc_across, '-', 'Color', across_col, 'LineWidth', 2.2);
yl = [-c2_loc_uncued_ylim c2_loc_uncued_ylim];
ylim(yl);
patch([conf_win_ms(1) conf_win_ms(2) conf_win_ms(2) conf_win_ms(1)], [yl(1) yl(1) yl(2) yl(2)], ...
    [0.8 0.8 0.8], 'FaceAlpha', 0.15, 'EdgeColor','none', 'HandleVisibility','off');
xlabel('Time (ms)', 'FontSize', fontsize);
ylabel('Projection', 'FontSize', fontsize);
legend([hLocWithin hLocAcross], {'Private', 'Shared'}, ...
    'Location','best', 'FontSize', fontsize-4);
xlim(c2_xwin);
set(gca, 'FontSize', fontsize-2);
box off;
hold off;

% Conflict on bottom
axC2(3) = subplot(3,1,3);
hold on;
title({'Conflict:'
    'Private vs. Shared'}, 'FontSize', fontsize);
yline(0, 'k--', 'LineWidth', 1.2, 'HandleVisibility','off');
hFill = fill([times fliplr(times)], [mean_conf_within+sem_conf_within, fliplr(mean_conf_within-sem_conf_within)], ...
    within_col, 'FaceAlpha', 0.20, 'EdgeColor','none');
hFill.HandleVisibility = 'off';
hConfWithin = plot(times, mean_conf_within, '-', 'Color', within_col, 'LineWidth', 2.2);
hFill = fill([times fliplr(times)], [mean_conf_across+sem_conf_across, fliplr(mean_conf_across-sem_conf_across)], ...
    across_col, 'FaceAlpha', 0.20, 'EdgeColor','none');
hFill.HandleVisibility = 'off';
hConfAcross = plot(times, mean_conf_across, '-', 'Color', across_col, 'LineWidth', 2.2);
yl = [-c2_conf_ylim c2_conf_ylim];
ylim(yl);
patch([conf_win_ms(1) conf_win_ms(2) conf_win_ms(2) conf_win_ms(1)], [yl(1) yl(1) yl(2) yl(2)], ...
    [0.8 0.8 0.8], 'FaceAlpha', 0.15, 'EdgeColor','none', 'HandleVisibility','off');
xlabel('Time (ms)', 'FontSize', fontsize);
ylabel('Projection', 'FontSize', fontsize);
legend([hConfWithin hConfAcross], {'Private', 'Shared'}, ...
    'Location','best', 'FontSize', fontsize-4);
xlim(c2_xwin);
set(gca, 'FontSize', fontsize-2);
box off;
hold off;

% sgtitle('Panel C.2: Within-vs-across from fullmodel TG diagonal (train=test)', 'FontSize', fontsize+2);

% Save out pdf.
set(gcf, 'Units', 'inches');
set(gcf, 'Position', [1 1 inches_out_c2]);          % on screen
set(gcf, 'PaperUnits', 'inches');
set(gcf, 'PaperPosition', [0 0 inches_out_c2]);     % for printed/exported output
set(gcf, 'PaperSize', inches_out_c2);        

print(gcf, '-dpdf', 'TGdiag_proj_attn_cnf_within_across.pdf');

%% ══════════════════════════════════════════════════════════════════════════
%% PANEL D  –  Within-dimension temporal generalization (4 distractor dims)
%% ══════════════════════════════════════════════════════════════════════════

dimNames = {'Color distractor', 'Shape distractor', ...
            'Edge distractor',  'Pattern distractor'};

figure('Name','Panel B – Conflict TG per distractor dimension', 'Color','w', ...
       'Position',[100 100 900 800]);

sgtitle('Dimension-Private conflict: temporal generalization', 'FontSize', fontsize+2);

for dist = 1:4

    subplot(2, 2, dist);
    hold on;

    map_dist = squeeze(conf_DSC(:,:,dist));
    % Per-subplot symmetric clim centered at 0 (tip 1: avoid asymmetric diverging map)
    delta_dist = max(abs(map_dist(:)));
    if ~isfinite(delta_dist) || delta_dist == 0, delta_dist = 0.1; end

    imagesc(times, times, map_dist);
    axis xy;
    axis square;
    xlim([times(1) times(end)]);
    ylim([times(1) times(end)]);
    colormap(rdBu());
    colorbar;
    clim([-delta_dist delta_dist]);

    % reference lines
    xline(0, 'k:', 'LineWidth', 1.5);
    yline(0, 'k:', 'LineWidth', 1.5);
    plot([times(1) times(end)], [times(1) times(end)], 'k--', 'LineWidth', 1);

    % conflict window box
    patch([conf_win_ms(1) conf_win_ms(2) conf_win_ms(2) conf_win_ms(1) conf_win_ms(1)], ...
          [conf_win_ms(1) conf_win_ms(1) conf_win_ms(2) conf_win_ms(2) conf_win_ms(1)], ...
          'k', 'FaceColor','none', 'EdgeColor','k', 'LineWidth', 1.5, 'LineStyle', '--');

    xticks(tick_ms);
    yticks(tick_ms);
    xtickangle(0);
    xlabel('Train time (ms)', 'FontSize', fontsize-2);
    ylabel('Test time (ms)', 'FontSize', fontsize-2);
    title(dimNames{dist}, 'FontSize', fontsize, 'Color', dimCols(dist,:));
    set(gca, 'FontSize', fontsize-2);

    hold off;
end


%% ══════════════════════════════════════════════════════════════════════════
%% PANEL D.2  –  Dimension-averaged temporal generalization (DSC vs DGC)
%% ══════════════════════════════════════════════════════════════════════════
inches_out = [9 6]; % paper size of plots

% ── Axis range and cluster-permutation settings ───────────────────────────
tg_plot_xlim    = [-300 1000];   % ms: display/analysis range for D.2 / D.3
run_cluster_perm = true;          % set false to skip (runtime increases with n_perm)
n_cluster_perm  = 5000;           % Monte Carlo sign-permutation iterations
cluster_tail    = 'right';          % 'two' (strict) | 'right' (directional +)
cluster_alpha_form = 0.01;        % cluster-forming threshold
cluster_alpha_sig  = 0.01;        % cluster significance threshold
cluster_report_topk = 5;          % print top-K observed clusters/p-values

% Average temporal-generalization matrices across distractor dimensions.
conf_DSC_avgDim = mean(conf_DSC, 3, 'omitnan');   % [n_times x n_times]
conf_DGC_avgDim = mean(conf_DGC, 3, 'omitnan');   % [n_times x n_times]

% Shared symmetric clim for left/right subplots in D.2.
panelD2_clim = max(abs([conf_DSC_avgDim(:); conf_DGC_avgDim(:)]));
if ~isfinite(panelD2_clim) || panelD2_clim == 0, panelD2_clim = 0.1; end

figure('Name','Panel D.2 – Avg TG across dimensions (DSC vs DGC)', 'Color','w', ...
       'Position',[120 120 980 430]);
sgtitle('Conflict orthogonalized', 'FontSize', fontsize+2);

subplot(1,2,1);
hold on;
imagesc(times, times, conf_DSC_avgDim);
axis xy;
axis square;
xlim(tg_plot_xlim);
ylim(tg_plot_xlim);
colormap(rdBu());
colorbar;
clim([-panelD2_clim panelD2_clim]);
xline(0, 'k:', 'LineWidth', 1.5);
yline(0, 'k:', 'LineWidth', 1.5);
plot([times(1) times(end)], [times(1) times(end)], 'k--', 'LineWidth', 1);
patch([conf_win_ms(1) conf_win_ms(2) conf_win_ms(2) conf_win_ms(1) conf_win_ms(1)], ...
      [conf_win_ms(1) conf_win_ms(1) conf_win_ms(2) conf_win_ms(2) conf_win_ms(1)], ...
      'k', 'FaceColor','none', 'EdgeColor','k', 'LineWidth', 1.5, 'LineStyle', '--');
xticks(tick_ms);
yticks(tick_ms);
xtickangle(0);   % force x‑axis labels to 0°
xlabel('Train time (ms)', 'FontSize', fontsize-2);
ylabel('Test time (ms)', 'FontSize', fontsize-2);
title('Dimension-Private conflict', 'FontSize', fontsize, 'Color', confCol);
set(gca, 'FontSize', fontsize-2);
box off;
hold off;

subplot(1,2,2);
hold on;
imagesc(times, times, conf_DGC_avgDim);
axis xy;
axis square;
xlim(tg_plot_xlim);
ylim(tg_plot_xlim);
colormap(rdBu());
colorbar;
clim([-panelD2_clim panelD2_clim]);
xline(0, 'k:', 'LineWidth', 1.5);
yline(0, 'k:', 'LineWidth', 1.5);
plot([times(1) times(end)], [times(1) times(end)], 'w--', 'LineWidth', 1);
patch([conf_win_ms(1) conf_win_ms(2) conf_win_ms(2) conf_win_ms(1) conf_win_ms(1)], ...
      [conf_win_ms(1) conf_win_ms(1) conf_win_ms(2) conf_win_ms(2) conf_win_ms(1)], ...
      'k', 'FaceColor','none', 'EdgeColor','k', 'LineWidth', 1.5, 'LineStyle', '--');
xticks(tick_ms);
yticks(tick_ms);
xtickangle(0);   % force x‑axis labels to 0°
xlabel('Train time (ms)', 'FontSize', fontsize-2);
ylabel('Test time (ms)', 'FontSize', fontsize-2);
title('Dimension-Shared conflict', 'FontSize', fontsize, 'Color', confCol);
set(gca, 'FontSize', fontsize-2);
box off;
hold off;

% ── Cluster-based permutation test: Panel D.2 ─────────────────────────────
if run_cluster_perm
    fprintf('\nPanel D.2 cluster permutation test (%d permutations)...\n', n_cluster_perm);
    % Subject-level dim-averaged TG maps: [n_times x n_times x nSubj]
    d_DSC_perm = squeeze(mean(conf_DSC_bySubj, 3, 'omitnan'));
    d_DGC_perm = squeeze(mean(conf_DGC_bySubj, 3, 'omitnan'));
    fprintf('  DSC...\n');
    [sig_DSC, cpval_DSC] = cluster_perm_test_TG_2D(d_DSC_perm, times, ...
        'n_perm', n_cluster_perm, 'tg_xlim', tg_plot_xlim, ...
        'tail', cluster_tail, 'alpha_cluster', cluster_alpha_form, ...
        'alpha_sig', cluster_alpha_sig, 'report_top_k', cluster_report_topk);
    fprintf('  DGC...\n');
    [sig_DGC, cpval_DGC] = cluster_perm_test_TG_2D(d_DGC_perm, times, ...
        'n_perm', n_cluster_perm, 'tg_xlim', tg_plot_xlim, ...
        'tail', cluster_tail, 'alpha_cluster', cluster_alpha_form, ...
        'alpha_sig', cluster_alpha_sig, 'report_top_k', cluster_report_topk);

    t_win_idx = find(times >= tg_plot_xlim(1) & times <= tg_plot_xlim(2));
    t_win = times(t_win_idx);

    subplot(1,2,1); hold on;
    if any(sig_DSC(:))
        contour(t_win, t_win, double(sig_DSC(t_win_idx, t_win_idx)), ...
            [0.5 0.5], 'k-', 'LineWidth', 2);
    end
    hold off;

    subplot(1,2,2); hold on;
    if any(sig_DGC(:))
        contour(t_win, t_win, double(sig_DGC(t_win_idx, t_win_idx)), ...
            [0.5 0.5], 'k-', 'LineWidth', 2);
    end
    hold off;
end

% Save out pdf.
set(gcf, 'Units', 'inches');
set(gcf, 'Position', [1 1 inches_out]);          % on screen
set(gcf, 'PaperUnits', 'inches');
set(gcf, 'PaperPosition', [0 0 inches_out]);     % for printed/exported output
set(gcf, 'PaperSize', inches_out);        

print(gcf, '-dpdf', 'TG_conflict.pdf');

%% ══════════════════════════════════════════════════════════════════════════
%% PANEL D.3  –  Dimension-averaged temporal generalization (SR vs SI)
%% ══════════════════════════════════════════════════════════════════════════

% Average temporal-generalization matrices across distractor dimensions.
loc_DSR_avgDim = mean(loc_DSR, 3, 'omitnan');   % [n_times x n_times]
loc_DSI_avgDim = mean(loc_DSI, 3, 'omitnan');   % [n_times x n_times]

% Shared symmetric clim for left/right subplots in D.3.
panelD3_clim = max(abs([loc_DSR_avgDim(:); loc_DSI_avgDim(:)]));
if ~isfinite(panelD3_clim) || panelD3_clim == 0, panelD3_clim = 0.1; end

figure('Name','Panel D.3 – Avg TG across dimensions (SR vs SI)', 'Color','w', ...
       'Position',[120 570 980 430]);
sgtitle('Location orthogonalized', 'FontSize', fontsize+2);

subplot(1,2,1);
hold on;
imagesc(times, times, loc_DSR_avgDim);
axis xy;
axis square;
xlim(tg_plot_xlim);
ylim(tg_plot_xlim);
colormap(rdBu());
colorbar;
clim([-panelD3_clim panelD3_clim]);
xline(0, 'k:', 'LineWidth', 1.5);
yline(0, 'k:', 'LineWidth', 1.5);
plot([times(1) times(end)], [times(1) times(end)], 'k--', 'LineWidth', 1);
patch([conf_win_ms(1) conf_win_ms(2) conf_win_ms(2) conf_win_ms(1) conf_win_ms(1)], ...
      [conf_win_ms(1) conf_win_ms(1) conf_win_ms(2) conf_win_ms(2) conf_win_ms(1)], ...
      'k', 'FaceColor','none', 'EdgeColor','k', 'LineWidth', 1.5, 'LineStyle', '--');
xticks(tick_ms);
yticks(tick_ms);
xtickangle(0);   % force x‑axis labels to 0°
xlabel('Train time (ms)', 'FontSize', fontsize-2);
ylabel('Test time (ms)', 'FontSize', fontsize-2);
title('Feature processing - Uncued', 'FontSize', fontsize, 'Color', locCol);
set(gca, 'FontSize', fontsize-2);
box off;
hold off;

subplot(1,2,2);
hold on;
imagesc(times, times, loc_DSI_avgDim);
axis xy;
axis square;
xlim(tg_plot_xlim);
ylim(tg_plot_xlim);
colormap(rdBu());
colorbar;
clim([-panelD3_clim panelD3_clim]);
xline(0, 'k:', 'LineWidth', 1.5);
yline(0, 'k:', 'LineWidth', 1.5);
plot([times(1) times(end)], [times(1) times(end)], 'w--', 'LineWidth', 1);
patch([conf_win_ms(1) conf_win_ms(2) conf_win_ms(2) conf_win_ms(1) conf_win_ms(1)], ...
      [conf_win_ms(1) conf_win_ms(1) conf_win_ms(2) conf_win_ms(2) conf_win_ms(1)], ...
      'k', 'FaceColor','none', 'EdgeColor','k', 'LineWidth', 1.5, 'LineStyle', '--');
xticks(tick_ms);
yticks(tick_ms);
xtickangle(0);   % force x‑axis labels to 0°
xlabel('Train time (ms)', 'FontSize', fontsize-2);
ylabel('Test time (ms)', 'FontSize', fontsize-2);
title('Spatial attention - Uncued', 'FontSize', fontsize, 'Color', locCol);
set(gca, 'FontSize', fontsize-2);
box off;
hold off;

% ── Cluster-based permutation test: Panel D.3 ─────────────────────────────
if run_cluster_perm
    fprintf('\nPanel D.3 cluster permutation test (%d permutations)...\n', n_cluster_perm);
    d_DSR_perm = squeeze(mean(loc_DSR_bySubj, 3, 'omitnan'));
    d_DSI_perm = squeeze(mean(loc_DSI_bySubj, 3, 'omitnan'));
    fprintf('  DSR...\n');
    [sig_DSR, cpval_DSR] = cluster_perm_test_TG_2D(d_DSR_perm, times, ...
        'n_perm', n_cluster_perm, 'tg_xlim', tg_plot_xlim, ...
        'tail', cluster_tail, 'alpha_cluster', cluster_alpha_form, ...
        'alpha_sig', cluster_alpha_sig, 'report_top_k', cluster_report_topk);
    fprintf('  DSI...\n');
    [sig_DSI, cpval_DSI] = cluster_perm_test_TG_2D(d_DSI_perm, times, ...
        'n_perm', n_cluster_perm, 'tg_xlim', tg_plot_xlim, ...
        'tail', cluster_tail, 'alpha_cluster', cluster_alpha_form, ...
        'alpha_sig', cluster_alpha_sig, 'report_top_k', cluster_report_topk);

    t_win_idx = find(times >= tg_plot_xlim(1) & times <= tg_plot_xlim(2));
    t_win = times(t_win_idx);

    subplot(1,2,1); hold on;
    if any(sig_DSR(:))
        contour(t_win, t_win, double(sig_DSR(t_win_idx, t_win_idx)), ...
            [0.5 0.5], 'k-', 'LineWidth', 2);
    end
    hold off;

    subplot(1,2,2); hold on;
    if any(sig_DSI(:))
        contour(t_win, t_win, double(sig_DSI(t_win_idx, t_win_idx)), ...
            [0.5 0.5], 'k-', 'LineWidth', 2);
    end
    hold off;
end

% Save out pdf.
set(gcf, 'Units', 'inches');
set(gcf, 'Position', [1 1 inches_out]);          % on screen
set(gcf, 'PaperUnits', 'inches');
set(gcf, 'PaperPosition', [0 0 inches_out]);     % for printed/exported output
set(gcf, 'PaperSize', inches_out);        

print(gcf, '-dpdf', 'TG_attn_noncued.pdf');

%% ══════════════════════════════════════════════════════════════════════════
%% PANEL D.4  –  Dimension-averaged temporal generalization (SR vs SI, cued)
%% ══════════════════════════════════════════════════════════════════════════

% Cued maps are already single-dimension TG projections after PC/subj averaging.
loc_cued_DSR_map = loc_cued_DSR;   % [n_times x n_times]
loc_cued_DSI_map = loc_cued_DSI;   % [n_times x n_times]

% Shared symmetric clim for left/right subplots in D.4.
panelD4_clim = max(abs([loc_cued_DSR_map(:); loc_cued_DSI_map(:)]));
if ~isfinite(panelD4_clim) || panelD4_clim == 0, panelD4_clim = 0.1; end

figure('Name','Panel D.4 – Avg TG across dimensions (SR vs SI, cued)', 'Color','w', ...
       'Position',[1120 570 980 430]);
sgtitle('Location orthogonalized (cued)', 'FontSize', fontsize+2);

subplot(1,2,1);
hold on;
imagesc(times, times, loc_cued_DSR_map);
axis xy;
axis square;
xlim(tg_plot_xlim);
ylim(tg_plot_xlim);
colormap(rdBu());
colorbar;
clim([-panelD4_clim panelD4_clim]);
xline(0, 'k:', 'LineWidth', 1.5);
yline(0, 'k:', 'LineWidth', 1.5);
plot([times(1) times(end)], [times(1) times(end)], 'k--', 'LineWidth', 1);
patch([conf_win_ms(1) conf_win_ms(2) conf_win_ms(2) conf_win_ms(1) conf_win_ms(1)], ...
      [conf_win_ms(1) conf_win_ms(1) conf_win_ms(2) conf_win_ms(2) conf_win_ms(1)], ...
      'k', 'FaceColor','none', 'EdgeColor','k', 'LineWidth', 1.5, 'LineStyle', '--');
xticks(tick_ms);
yticks(tick_ms);
xtickangle(0);   % force x‑axis labels to 0°
xlabel('Train time (ms)', 'FontSize', fontsize-2);
ylabel('Test time (ms)', 'FontSize', fontsize-2);
title('Feature processing - Cued', 'FontSize', fontsize, 'Color', locCol);
set(gca, 'FontSize', fontsize-2);
box off;
hold off;

subplot(1,2,2);
hold on;
imagesc(times, times, loc_cued_DSI_map);
axis xy;
axis square;
xlim(tg_plot_xlim);
ylim(tg_plot_xlim);
colormap(rdBu());
colorbar;
clim([-panelD4_clim panelD4_clim]);
xline(0, 'k:', 'LineWidth', 1.5);
yline(0, 'k:', 'LineWidth', 1.5);
plot([times(1) times(end)], [times(1) times(end)], 'k--', 'LineWidth', 1);
patch([conf_win_ms(1) conf_win_ms(2) conf_win_ms(2) conf_win_ms(1) conf_win_ms(1)], ...
      [conf_win_ms(1) conf_win_ms(1) conf_win_ms(2) conf_win_ms(2) conf_win_ms(1)], ...
      'k', 'FaceColor','none', 'EdgeColor','k', 'LineWidth', 1.5, 'LineStyle', '--');
xticks(tick_ms);
yticks(tick_ms);
xtickangle(0);   % force x‑axis labels to 0°
xlabel('Train time (ms)', 'FontSize', fontsize-2);
ylabel('Test time (ms)', 'FontSize', fontsize-2);
title('Spatial attention - Cued', 'FontSize', fontsize, 'Color', locCol);
set(gca, 'FontSize', fontsize-2);
box off;
hold off;

% ── Cluster-based permutation test: Panel D.4 ─────────────────────────────
if run_cluster_perm
    fprintf('\nPanel D.4 cluster permutation test (%d permutations)...\n', n_cluster_perm);
    % Subject-level cued maps after PC averaging: [time x time x subj]
    d_cued_DSR_perm = squeeze(mean(cued_DSR_all(PC_idx,:,:,:), 1, 'omitnan'));
    d_cued_DSI_perm = squeeze(mean(cued_DSI_all(PC_idx,:,:,:), 1, 'omitnan'));

    fprintf('  Cued DSR...\n');
    [sig_cued_DSR, cpval_cued_DSR] = cluster_perm_test_TG_2D(d_cued_DSR_perm, times, ...
        'n_perm', n_cluster_perm, 'tg_xlim', tg_plot_xlim, ...
        'tail', cluster_tail, 'alpha_cluster', cluster_alpha_form, ...
        'alpha_sig', cluster_alpha_sig, 'report_top_k', cluster_report_topk);
    fprintf('  Cued DSI...\n');
    [sig_cued_DSI, cpval_cued_DSI] = cluster_perm_test_TG_2D(d_cued_DSI_perm, times, ...
        'n_perm', n_cluster_perm, 'tg_xlim', tg_plot_xlim, ...
        'tail', cluster_tail, 'alpha_cluster', cluster_alpha_form, ...
        'alpha_sig', cluster_alpha_sig, 'report_top_k', cluster_report_topk);

    t_win_idx = find(times >= tg_plot_xlim(1) & times <= tg_plot_xlim(2));
    t_win = times(t_win_idx);

    subplot(1,2,1); hold on;
    if any(sig_cued_DSR(:))
        contour(t_win, t_win, double(sig_cued_DSR(t_win_idx, t_win_idx)), ...
            [0.5 0.5], 'k-', 'LineWidth', 2);
    end
    hold off;

    subplot(1,2,2); hold on;
    if any(sig_cued_DSI(:))
        contour(t_win, t_win, double(sig_cued_DSI(t_win_idx, t_win_idx)), ...
            [0.5 0.5], 'k-', 'LineWidth', 2);
    end
    hold off;
end

% Save out pdf.
set(gcf, 'Units', 'inches');
set(gcf, 'Position', [1 1 inches_out]);          % on screen
set(gcf, 'PaperUnits', 'inches');
set(gcf, 'PaperPosition', [0 0 inches_out]);     % for printed/exported output
set(gcf, 'PaperSize', inches_out);        

print(gcf, '-dpdf', 'TG_attn_cued.pdf');
