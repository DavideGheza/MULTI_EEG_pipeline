%% Plot congruency-sequence effects in time-frequency power
% Loads single-trial time-frequency data, computes stimulus- or
% response-locked CSE summaries, and produces time-frequency maps,
% topographies, and cluster plots.
% Author: Davide Gheza

clc
clear

%% 1. Load settings

freqshown = [3 8];
controlfortrialnumber = 0;
singletrial_baseline = 0;

%Select onset or response events
onset = false;
if onset
    suffix_eeg = 'onset_';
    suffix = '';
    timeshown = [200 500];
else
    timeshown = [-75 75]; % peri-response
    suffix = 'response_';   % response-locked requires precomputing stimulus-locked to get a fixed pre-stimulus baseline
    suffix_eeg = 'response_';
end

%Struct 'sets' is loaded from the 'scripts_and_settings' folder:
if ispc
    load('C:\Users\gdavide\OneDrive\Documenti\MATLAB\dg_workspace22\dimensions\analyses\scripts_and_settings\preprocessing_settings.mat');
else
    load('/Users/davidegheza/Library/CloudStorage/OneDrive-Personal/Documenti/MATLAB/dg_workspace22/dimensions/analyses/scripts_and_settings/preprocessing_settings_mac.mat');
end

% load an eeg dataset to fetch info
sampleEEG = pop_loadset('filename',['22' suffix_eeg 'trialrej.set'],...
                'filepath',[sets.clean_data_dir]...
                );
                
chanlocs = sampleEEG.chanlocs;
labels = {sampleEEG.chanlocs.labels};

%% Load normalized power 
sn=[22 24:48]; % v3


% sbj, chans, freqs, times, cond
power = zeros(length(sn),61,30,230,4);
power_si = zeros(4,61,30,230);

for si = 1:length(sn)
    
    load([sets.sessions_data_tf_dir...
            'data_' suffix 'tf_singletrial_NB_sbj_'...
            sets.sub_idx{sn(si)}...
            '.mat']);
    
    % load behavioral data
    behav = load([sets.sessions_data_dir '/V3_behav_seqinfo_trimmed_' suffix 'sbj_' num2str(sn(si))]);
    behav = behav.behav_trim;
    
    % define vectors for indexing trials by PC
    Neu_1_2 = nan(length(behav.blkwise_PCs),2);     %preallocate - two neutral dimensions per trial
    MI = nan(length(behav.blkwise_PCs),1);   
    MC = nan(length(behav.blkwise_PCs),1);
    % find task based on PC
    for t = 1:length(behav.blkwise_PCs)      
        Neu_1_2(t,:) = find(behav.blkwise_PCs(t,:) == 1/2);
        MI(t,:) = find(behav.blkwise_PCs(t,:) == 1/3);
        MC(t,:) = find(behav.blkwise_PCs(t,:) == 2/3);
        if Neu_1_2(t,1) == 1 && Neu_1_2(t,2) == 4               % Neu_1 and Neu_2 are artificially segregated as neutral dimension on left vs right column (of behav.blkwise_PCs). Hence, if Neu1 happens to be contiguous to Neu2, Neu1 never contains cond 4, and Neu2 never contains cond 1.
            Neu_1_2(t,1) = 4;                                   % To restore original balance, swap neu left and neu right, when they are respectively dimension 1 or 4 (i.e., follow right-shift logic used for associating taskID-PCs over blocks)
            Neu_1_2(t,2) = 1;
        elseif Neu_1_2(t,1) == 1 && Neu_1_2(t,2) == 3 || Neu_1_2(t,1) == 2 && Neu_1_2(t,2) == 4 % Alternatively, if Neu1 and Neu2 are not contiguous, Neu1 will always contain cond 1 and 2, and Neu2 3 and 4. 
            if Neu_1_2(t,1) == 2 && MI(t) == 1
                Neu_1_2(t,1) = 4;                                   % To restore original balance, swap neu left and neu right, based on relative position of N1 and MI (i.e., follow right-shift logic used for associating taskID-PCs over blocks)
                Neu_1_2(t,2) = 2;
            elseif Neu_1_2(t,1) == 1 && MI(t) == 4
                Neu_1_2(t,1) = 3;                                   % To restore original balance, swap neu left and neu right, based on relative position of N1 and MI  (i.e., follow right-shift logic used for associating taskID-PCs over blocks)
                Neu_1_2(t,2) = 1;
            end
        end
    end
    Neu_1 = Neu_1_2(:,1);
    Neu_2 = Neu_1_2(:,2);
    PC_matrix = [MI,MC,Neu_1,Neu_2];

    % define vectors for indexing trial: logical index for congruency
    CNG0 = sum(behav.congruency,2)==-1;
    CNG1 = sum(behav.congruency,2)==0;
    CNG2 = sum(behav.congruency,2)==1;
    CNG3 = sum(behav.congruency,2)==2;
    CNG_matrix = [CNG0,CNG1,CNG2,CNG3];


    % % define vectors for indexing trial: CSE. Combine pairs CNG levels
    % II = (behav.prev_sum_con==0 | behav.prev_sum_con==1) & (behav.sum_con==0 | behav.sum_con==1)  & behav.taskcount>1 & behav.post_correct==1;
    % IC = (behav.prev_sum_con==0 | behav.prev_sum_con==1) & (behav.sum_con==2 | behav.sum_con==3)  & behav.taskcount>1 & behav.post_correct==1;
    % CI = (behav.prev_sum_con==2 | behav.prev_sum_con==3) & (behav.sum_con==0 | behav.sum_con==1)  & behav.taskcount>1 & behav.post_correct==1;
    % CC = (behav.prev_sum_con==2 | behav.prev_sum_con==3) & (behav.sum_con==2 | behav.sum_con==3)  & behav.taskcount>1 & behav.post_correct==1;
    % define vectors for indexing trial: CSE. SELECT EXTREME CNG levels
    II = behav.prev_sum_con==0 & behav.sum_con==0  & behav.taskcount>1 & behav.post_correct==1;
    IC = behav.prev_sum_con==0 & behav.sum_con==3  & behav.taskcount>1 & behav.post_correct==1;
    CI = behav.prev_sum_con==3 & behav.sum_con==0  & behav.taskcount>1 & behav.post_correct==1;
    CC = behav.prev_sum_con==3 & behav.sum_con==3  & behav.taskcount>1 & behav.post_correct==1;
    CSE_matrix = [II,IC,CI,CC];
    
    % define PC filter 
    PC_filter = sum(data.labels_ses1==PC_matrix(:,3:4), 2); % neutral PC only

    if onset==false
        data.baseline = permute(data.baseline,[2,3,1]);      %  trials, chans, freq (permuted one too many in decoding_tf_sessions_singletrial)
    end
    trial_avg_baseline = mean(data.baseline,1); % trial average
    
    if singletrial_baseline
            
        

        % testing different single-trial baseline-normalization approaches
        
        % 1) single trial baseline, but baseline obtained as average across trials and times
        % (approach used for decoding data - should produce negative skewed results)
        % This is because power is right skewed, and baseline freq values
        % (mean over trials) tend to be larger than single-trial data. The
        % log of values<0 is negative.
        % % expanded_baseline = reshape(data.baseline, [1, 61, 30, 1]); % 09/22/2025 This is needed for obsolete single trial t/f. Now baseline is saved trial-wise and averaged above
        % power_normalized = 10*log10(bsxfun(@rdivide,data.ses1,trial_avg_baseline)); % baseline correction on single trials (baseline obtained as trial average)
        % for CSE = 1:size(CSE_matrix,2)
        %     power_si(CSE,:,:,:) = squeeze(mean( power_normalized(CSE_matrix(:,CSE) & PC_filter, :,:,:) ,1));    % mean over trials
        % end
        % permuted_power = permute(power_si,[2,3,4,1]);
        % power(si,:,:,:,:) = permuted_power; % catenate sbj data along 1st dim (sbj, chan, freq, time, cond)


        % 2) single trial baseline normalization (gain model)
        % positively skewed ("large positive post-baseline values are dominating")
        % base_time_idxs = (dsearchn(data.times',-800):dsearchn(data.times',-700)); % for stimulus locked: select the times -800 -700
        % st_baseline = mean(data.ses1(:,:,:,base_time_idxs),4); % mean over baseline times
        % power_normalized_percent = (bsxfun(@rdivide,data.ses1,st_baseline)); % baseline correction on single trials (w/ single trial baseline)
        % for CSE = 1:size(CSE_matrix,2)
        %     power_si(CSE,:,:,:) = squeeze(mean(power_normalized_percent(CSE_matrix(:,CSE) & PC_filter, :,:,:) ,1));    % mean over trials
        % end
        % power_normalized = 10*log10(power_si);
        % permuted_power = permute(power_normalized,[2,3,4,1]);
        % power(si,:,:,:,:) = permuted_power; % catenate sbj data along 1st dim (sbj, chan, freq, time, cond)


        % 3) dB full-trial single trial baseline normalization (gain-model - recommended by Grandchamp)
        % the 2nd step (post trial average) allows to normalize across
        % conditions. However, using single-trial in decoding the 2nd step does not apply. 
        % Does the full epoch, within trial normalization remove condition differences? OPEN Q --> yes, see below
        st_baseline = mean(data.ses1(:,:,:,:),4); % mean over baseline times (full epoch)
        power_normalized_st = (bsxfun(@rdivide,data.ses1,st_baseline)); % baseline correction on single trials (w/ single trial 'baseline' corersponding to the entire epoch mean)
        for CSE = 1:size(CSE_matrix,2)
            power_si(CSE,:,:,:) = squeeze(mean(power_normalized_st(CSE_matrix(:,CSE) & PC_filter, :,:,:) ,1));    % mean over trials
        end 
        % now perform 'classic' baseline normalization over the trial-averaged normalized power estimates - note: assumes across-conditions baseline
        base_time_idxs = (dsearchn(data.times',-600):dsearchn(data.times',-400)); % for stimulus locked: select the times -600 -400
        st_baseline_2 = mean(power_si(:,:,:,base_time_idxs),[1,4]);      % mean across times and conditions
        power_normalized_percent = bsxfun(@rdivide,power_si,st_baseline_2);

        power_normalized = 10*log10(power_normalized_percent);
        permuted_power = permute(power_normalized,[2,3,4,1]);
        power(si,:,:,:,:) = permuted_power; % catenate sbj data along 1st dim (sbj, chan, freq, time, cond)
        
        
        % 4) z-score full-trial single trial baseline normalization (additive-model - recommended by Grandchamp)
        % mu_st_baseline = mean(data.ses1(:,:,:,:),4); % mean over baseline times (full epoch)
        % sd_st_baseline = std(data.ses1(:,:,:,:),0,4); % std over baseline times (full epoch); 0 -> normalizes by N-1, on 4th dimension (times)
        % power_normalized_st = ((data.ses1-mu_st_baseline)./sd_st_baseline); % baseline NORMALIZATION on single trials (ADDITIVE model, w/ single trial 'baseline' corresponding to the entire epoch mean)
        % for CSE = 1:size(CSE_matrix,2)
        %     power_si(CSE,:,:,:) = squeeze(mean(power_normalized_st(CSE_matrix(:,CSE) & PC_filter, :,:,:) ,1));    % mean over trials
        % end 
        % % now perform 'classic' baseline normalization over the trial-averaged normalized power estimates - note: assumes across-conditions baseline
        % base_time_idxs = (dsearchn(data.times',-800):dsearchn(data.times',-700)); % for stimulus locked: select the times -800 -700
        % mu_st_baseline_2 = mean(power_si(:,:,:,base_time_idxs),[1,4]);     % mean across times and conditions
        % sd_st_baseline_2 = std(power_si(:,:,:,base_time_idxs),0,[1,4]);    % std over baseline times (full epoch); 0 -> normalizes by N-1, on 4th dimension (times) and conditions
        % power_normalized_z = ((power_si-mu_st_baseline_2)./sd_st_baseline_2);
        % 
        % permuted_power = permute(power_normalized_z,[2,3,4,1]);
        % power(si,:,:,:,:) = permuted_power; % catenate sbj data along 1st dim (sbj, chan, freq, time, cond)

        
        % 5) z-score normalization (additive model) - but NOT full epoch normalization. Therefore, no 2nd step. (keep the values positively skewed)
        % issue: 1) positive skew post baseline, 2) no way of eliminating condition-level biases in baseline
        % base_time_idxs = (dsearchn(data.times',-800):dsearchn(data.times',-700)); % for stimulus locked: select the times -800 -700
        % mu_st_baseline = mean(data.ses1(:,:,:,base_time_idxs),4); % mean over baseline times (full epoch)
        % sd_st_baseline = std(data.ses1(:,:,:,base_time_idxs),0,4); % std over baseline times (full epoch); 0 -> normalizes by N-1, on 4th dimension (times)
        % power_normalized_st = ((data.ses1-mu_st_baseline)./sd_st_baseline); % baseline NORMALIZATION on single trials (ADDITIVE model, w/ single trial 'baseline' corresponding to the entire epoch mean)
        % for CSE = 1:size(CSE_matrix,2)
        %     power_si(CSE,:,:,:) = squeeze(mean(power_normalized_st(CSE_matrix(:,CSE) & PC_filter, :,:,:) ,1));    % mean over trials
        % end 
        % 
        % permuted_power = permute(power_si,[2,3,4,1]);
        % power(si,:,:,:,:) = permuted_power; % catenate sbj data along 1st dim (sbj, chan, freq, time, cond)

        
        % Final note: 
        % Grandchamp & Delorme (2011) recommend a 2-step procedure where first, the single trial normalization is applied using the entire trial epoch, and then a 'classic' across trial baseline normalziation is applied
        % In the 2nd step, it is important to compute a baseline across conditions (accounts for noise condition-differences at baseline, and helps stabilize the baseline). This is possible with 3) and 4), with the 2nd 'classic' normalization that follows trial average.
        % However, for single-trial decoding or regression analyses, the second step of across-trial or across-condition baseline normalization DOES NOT APPLY
        % So, using full-epoch normalization for single-trial analyese is NOT recommended: systematic difference in power between conditions across the trial gets partially "absorbed" into mu at 1st step normalization (full epoch includes the post-stimulus activity that differs by condition). 
        % Therefore, condition contrasts are attenuated and effects can be masked.
        % For single-trial analysis, prefer normalization using only baseline time window: 2) and 5) (single trial baseline, giving positive skew. Probably better than across trials baseline, giving negative skew and losing trial-specific variability)



        
    else

        for CSE = 1:size(CSE_matrix,2)
            if controlfortrialnumber == 0
            power_si(CSE,:,:,:) = squeeze(mean( data.ses1(CSE_matrix(:,CSE) & PC_filter, :,:,:) ,1));    % mean over trials
            else
            power_temp =  data.ses1(CSE_matrix(:,CSE), :,:,:) ;   
            power_temp(20:end,:,:,:) =  []; % cap trial number
            power_si(CSE,:,:,:) = squeeze(mean( power_temp,1));    % mean over trials
            end
        end
        power(si,:,:,:,:) = 10*log10(bsxfun(@rdivide,permute(power_si,[2,3,4,1]),squeeze(trial_avg_baseline))); % baseline normalization after trial averaging, catenate sbj data along 1st dim (sbj, chan, freq, time, cond)

    end
    




    fprintf('\tLoaded sbj %g / %g\n', si, length(sn))
end
 
% Time-frequency data use the 1.5--40 Hz range from the decoding settings.
% min_freq = 1.5;
% max_freq = 40;
% num_frex = 30


%% plotting t/f stimulus

% plotting sets
chantoplot = {'Cz' 'FCz' 'FC1' 'FC2' 'C1' 'C2' 'CPz' 'CP1' 'CP2'};
chan = ismember(labels(:), chantoplot);
timein = -700;   % avoid edge artifacts at -800
timeout = 700; % for resp-locked. Max: 1498
freqin = 1.5;           
freqout = 40;  % 100         
clim = [-1.8 1.8];   
ymin_idx = dsearchn(data.freqs', freqin);
ymax_idx = dsearchn(data.freqs', freqout);
ylim([ymin_idx ymax_idx]);
tick_idx = round(linspace(ymin_idx,ymax_idx,8));
valid_idx = tick_idx(tick_idx>=ymin_idx & tick_idx<=ymax_idx);
%

figure

    
% 
% 
% splabels = {strjoin(['II - ' (string(chantoplot))]), ...
%             strjoin(['IC - ' (string(chantoplot))]), ...
%             strjoin(['CI - ' (string(chantoplot))]), ...
%             strjoin(['CC - ' (string(chantoplot))])};
splabels = {('II'), ...
            ('IC'), ...
            ('CI'), ...
            ('CC')};

spdata = {};            
for cnd = 1:size(power,5)
    spdata = [spdata, mean(power(:,chan,:,:,cnd),2)];
end


for sp = (1:4)

subplot(4,1,sp);
imagesc(data.times,[],squeeze(mean(spdata{1,sp}(:,:,:,:,:),1)));     
cbh = colorbar;
set(get(cbh,'xlabel'),'string','dB','FontSize',10);
set(gca,'clim',clim,'xlim',[timein timeout],'YDir','normal'); 
if sp>3
    xlabel('Time (ms)','FontSize',17,'FontName','Arial')
end
ylabel('Frequency (hz)','FontSize',10,'FontName','Arial')
set(gca, 'YTick', valid_idx, 'YTickLabel', round(data.freqs(valid_idx),1),'FontSize',8);
title(splabels(sp),'FontSize',12,'FontName','Arial')





% For horizontal 3 Hz line -> convert to index
idx3 = dsearchn(data.freqs', 3);
line([timein timeout],[idx3 idx3], 'LineWidth',1);


% For a vertical 0 ms line:
line_x = [0 0]; 
line_y = [ymin_idx ymax_idx];    % indices
line(line_x,line_y,'Color',[0.5 0.5 0.5],'LineStyle',':','LineWidth',2)
% set freq band line
idx1 = dsearchn(data.freqs', 3);
line([timein timeout],[idx1 idx1],'Color',[0.4 0.4 0.4],'LineStyle','-','LineWidth',1)
idx2 = dsearchn(data.freqs', 7);
line([timein timeout],[idx2 idx2],'Color',[0.4 0.4 0.4],'LineStyle','-','LineWidth',1)



end 

% save out figure
aspectRatio = 0.5;
set(gcf,'Units','inches');
screenposition = get(gcf,'Position');
newWidth = screenposition(4) * aspectRatio;  % width and height based on aspect ratio
newHeight = screenposition(4);  % Keep the height the same

set(gcf, ...
    'PaperPosition', [0 0 newWidth newHeight], ...
    'PaperSize', [newWidth newHeight]);
if onset == 1
    filename = 'CSE_stim_locked';
else
    filename = 'CSE_resp_locked';
end
print('-dpdf', '-vector', filename);



%% Scalp topography - cue

cmap = colormap(parula);
colormap parula
% prepare topo (sub plot) data

timeindexes = dsearchn(data.times',timeshown(1)):dsearchn(data.times',timeshown(2));
freqindexes = dsearchn(data.freqs',freqshown(1)):dsearchn(data.freqs',freqshown(2));

toposplabels = {['II'], ...
                ['IC'], ...
                ['CI'], ...
                ['CC']};
topospdata = {};            
for cnd = 1:size(power,5)
    topospdata = [topospdata, power(:,:,freqindexes,timeindexes,cnd)];
end



figure;
% preparing topo plot
EEG.chanlocs = chanlocs;

for toposp = 1:4
    
subplot(4,1,toposp);
% average on the 1st (sbj), 4th (freq) and 5th (time) dimensions
topoplot(mean(mean(mean(topospdata{1,toposp},1),4),5),...
            EEG.chanlocs,...
            'plotchans',[1:61],...
            'maplimits',clim, ...
            'colormap', cmap);    
% if rem(toposp,4)==0   # commented out since already included in tf plot (same clim)
%     cbh = colorbar; 
%     set(get(cbh,'xlabel'),'string','dB','FontSize',16);
%     set(gca,'FontSize', 16)
%     cbh.Position = cbh.Position + [0.10 -0.05 0 0.1]; % Adjust position of colorbar [x, y, width, height] 
% end
title(toposplabels(toposp),'Fontsize', 16,'FontName','Arial')

end

% sgtitle([ num2str(timeshown) 'ms - ' num2str(freqshown) 'Hz'])
% sgtitle('Average over 0-500 ms & 3-7 Hz','FontSize',18)


% save out figure
set(gcf,'Units','inches');
screenposition = get(gcf,'Position');
newWidth = screenposition(4) * aspectRatio;  % width and height based on aspect ratio
newHeight = screenposition(4);  % Keep the height the same

set(gcf, ...
    'PaperPosition', [0 0 newWidth newHeight], ...
    'PaperSize', [newWidth newHeight]);

if onset == 1
    filename = 'topo_CSE_stim_locked';
else
    filename = 'topo_CSE_resp_locked';
end
print('-dpdf', '-vector', filename);



%% contrasts

% prepare topo (sub plot) data
% cond order: [II, IC, CI, CC]

timeindexes = dsearchn(data.times',timeshown(1)):dsearchn(data.times',timeshown(2));
freqindexes = dsearchn(data.freqs',freqshown(1)):dsearchn(data.freqs',freqshown(2));

toposplabels = {['Ix-Cx  ' num2str(timeshown) ' - ' num2str(freqshown) ], ...
            ['xI - xC  ' num2str(timeshown) ' - ' num2str(freqshown) ]...
            ['CI-II  ' num2str(timeshown) ' - ' num2str(freqshown) ]...
            ['CC-IC  ' num2str(timeshown) ' - ' num2str(freqshown) ]};
power_selection = squeeze(power(:,:,freqindexes,timeindexes,:));

topospdata = {(mean(power_selection(:,:,[1,2]),3) - mean(power_selection(:,:,[3,4]),3)), ...
              (mean(power_selection(:,:,[1,3]),3) - mean(power_selection(:,:,[2,4]),3)), ...
              (power_selection(:,:,3) - power_selection(:,:,1)), ...
              (power_selection(:,:,4) - power_selection(:,:,2))};      


figure;
% preparing topo plot
EEG.chanlocs = chanlocs;


for toposp = 1:4
    
subplot(1,4,toposp);
% average on the 1st (sbj), 4th (freq) and 5th (time) dimensions
topoplot(mean(mean(mean(topospdata{1,toposp},1),4),5),EEG.chanlocs,'plotchans',[1:61],'maplimits',[-0.5 0.5], 'colormap', cmap);   
cbh = colorbar;
set(get(cbh,'xlabel'),'string','dB','FontSize',12);
set(gca,'FontSize', 12)
title(toposplabels(toposp),'Fontsize', 14,'FontName','Arial')

end





%% average for single plot 
cmap = colormap(parula);
colormap parula

figure

chantoplot = {'Cz' 'FCz' 'FC1' 'FC2' 'C1' 'C2' 'CPz' 'CP1' 'CP2'};
chan = ismember(labels(:), chantoplot);
timein = -700;   % avoid edge artifacts at -800
timeout = 700; % for resp-locked. Max: 1498
freqin = 1.5;           
freqout = 40;  % 100         
clim = [-1.8 1.8];   
ymin_idx = dsearchn(data.freqs', freqin);
ymax_idx = dsearchn(data.freqs', freqout);
ylim([ymin_idx ymax_idx]);
tick_idx = round(linspace(ymin_idx,ymax_idx,8));
valid_idx = tick_idx(tick_idx>=ymin_idx & tick_idx<=ymax_idx);



          
spdata = mean(power(:,chan,:,:,:),[2,5]);



imagesc(data.times,[],squeeze(mean(spdata,1)));   
cbh = colorbar;
set(get(cbh,'xlabel'),'string','dB','FontSize',16);
set(gca,'clim',clim,'xlim',[timein timeout],'YDir','normal'); 
xlabel('Time (ms)','FontSize',17,'FontName','Arial')
ylabel('Frequency (hz)','FontSize',17,'FontName','Arial')
set(gca, 'YTick', valid_idx, 'YTickLabel', round(data.freqs(valid_idx),1),'FontSize', 16);
title(('Grandaverage'),'FontSize',20,'FontName','Arial')



% For horizontal 3 Hz line -> convert to index
idx3 = dsearchn(data.freqs', 3);
line([timein timeout],[idx3 idx3], 'LineWidth',1);


% For a vertical 0 ms line:
line_x = [0 0]; 
line_y = [ymin_idx ymax_idx];    % indices
line(line_x,line_y,'Color',[0.5 0.5 0.5],'LineStyle',':','LineWidth',2)
% set freq band line
idx1 = dsearchn(data.freqs', 3);
line([timein timeout],[idx1 idx1],'Color',[0.4 0.4 0.4],'LineStyle','-','LineWidth',1)
idx2 = dsearchn(data.freqs', 7);
line([timein timeout],[idx2 idx2],'Color',[0.4 0.4 0.4],'LineStyle','-','LineWidth',1)












%% stats - 11/13/2025

%% Cluster-based permutation test: Previous Incongruent vs Previous Congruent

% assumes power(subj,chan,freq,time,cond)
% cond order: [II, IC, CI, CC]

opt1=1;

if opt1==1
    % Opt 1) running stats on the channel-averaged data for CENTRO-FRONTAL (test for current CNG)
    chantoplot = {'Cz' 'FCz' 'FC1' 'FC2' 'C1' 'C2' 'CPz' 'CP1' 'CP2'};
else
    % Opt 2) running stats on the channel-averaged data for PARIETO-OCCIPITAL (test for previou CNG)
    chantoplot = {'Pz','POz','CPz','CP1','CP2', ...
     'P1','P3','P5','PO3','PO7','O1','P9', ...
     'P2','P4','P6','PO4','PO8','O2','P10'};
end
chans = ismember(labels(:), chantoplot);

disp('Running cluster-based permutation test ...');

% average conditions and channels
if opt1==1
    cond1 = squeeze(mean(power(:,chans,:,:, [2 4]), [2,5])); % IC + CC = Con
    cond2 = squeeze(mean(power(:,chans,:,:, [1 3]), [2,5])); % II + CI = Inc
else
    cond1 = squeeze(mean(power(:,chans,:,:, [3 4]), [2,5])); % CI + CC %  = previous Cong
    cond2 = squeeze(mean(power(:,chans,:,:, [1 2]), [2,5])); % II + IC %  = previous Inc
end
diffdata = cond2 - cond1;  % subj x chan x freq x time
cond1 = permute(cond1,[2,3,1]); % the last dimension is the one permuted with statcond - make it subj (freq, times, sbj)
cond2 = permute(cond2,[2,3,1]);

% contrast

nSubj = size(power,1);
nFreq = size(power,3);
nTime = size(power,4);

% initialize mask
sigMask_tf = false(nFreq, nTime);
pval_tf = nan(nFreq, nTime);

% run permutation cluster-based test per channel
naccu = 5000; % number of permutations
alpha = 0.01;

%  test a freq x time matrix

[stats, df, pvals, surrog] = statcond({cond1, cond2}, ...
    'mode','perm','naccu',naccu,'paired','on');
% store
pval_tf(:,:) = pvals;
sigMask_tf(:,:) = pvals < 0.01;


disp('Permutation tests complete.');

% save('sigMask_tf_clustermass.mat','sigMask_tf','pval_tf','-v7.3');





%% plotting contrast and statistical mask
figure

clim = [-1.2 1.2];       


tfmap = squeeze(mean(diffdata,1)); % mean over sbjs
sig = sigMask_tf;


timein = -700;   % avoid edge artifacts at -800
timeout = 700; % for resp-locked. Max: 1498
freqin = 1.5;           
freqout = 40;  % 100         
ymin_idx = dsearchn(data.freqs', freqin);
ymax_idx = dsearchn(data.freqs', freqout);
ylim([ymin_idx ymax_idx]);
tick_idx = round(linspace(ymin_idx,ymax_idx,8));
valid_idx = tick_idx(tick_idx>=ymin_idx & tick_idx<=ymax_idx);



imagesc(data.times, [], tfmap);     % do not specify freqs!
cbh = colorbar;
set(get(cbh,'xlabel'),'string','dB','FontSize',22);
set(gca,'clim',clim,'xlim',[timein timeout],'YDir','normal');    
ylabel('Frequency (hz)','FontSize',22,'FontName','Arial')
xlabel('Time (ms)','FontSize',22,'FontName','Arial')
set(gca, 'YTick', valid_idx, 'YTickLabel', round(data.freqs(valid_idx),1),'FontSize', 22);



hold on
freq_idx = 1:length(data.freqs);   % convert freq axis → indices
contour(data.times, freq_idx, sig, 1, 'k', 'LineWidth', 4);
hold off




% For horizontal 3 Hz line -> convert to index
idx3 = dsearchn(data.freqs', 3);
line([timein timeout],[idx3 idx3], 'LineWidth',1);


% For a vertical 0 ms line:
line_x = [0 0]; 
line_y = [ymin_idx ymax_idx];    % indices
line(line_x,line_y,'Color',[0.5 0.5 0.5],'LineStyle',':','LineWidth',2)
% set freq band line
idx1 = dsearchn(data.freqs', 3);
line([timein timeout],[idx1 idx1],'Color',[0.4 0.4 0.4],'LineStyle','-','LineWidth',1)
idx2 = dsearchn(data.freqs', 7);
line([timein timeout],[idx2 idx2],'Color',[0.4 0.4 0.4],'LineStyle','-','LineWidth',1)


% save out figure
set(gcf,'Units','inches');
screenposition = get(gcf,'Position');
set(gcf,...
    'PaperPosition',[0 0 screenposition(3:4)],...
    'PaperSize',[screenposition(3:4)]);
if opt1 == 1
    filename = 'CSE_curr_CNG_cluster';
else
    filename = 'CSE_prev_CNG_cluster';
end
print('-dpdf', '-vector', filename);






%% using fieldtrip clustering method over chan space

% -------------------------
% Spatio-temporal cluster test (channels x freq x time)
% -------------------------
% ensure FieldTrip is on path

% cond order: [II, IC, CI, CC]

% build FieldTrip-style neighbours from EEG chanlocs
cfg_neighb = [];
cfg_neighb.method = 'triangulation'; % or 'distance' if triangulation fails
elec = [];
for ic = 1:length(chanlocs)
    elec.label{ic} = chanlocs(ic).labels;
    if isfield(chanlocs(ic),'X') && isfield(chanlocs(ic),'Y')
        elec.chanpos(ic,1:2) = [chanlocs(ic).X chanlocs(ic).Y];
        elec.elecpos(ic,1:2) = [chanlocs(ic).X chanlocs(ic).Y];
    else
        elec.chanpos(ic,1:2) = [0 0]; elec.elecpos(ic,1:2) = [0 0];
    end
end
cfg_neighb.elec = elec;
neighbours = ft_prepare_neighbours(cfg_neighb, elec);

% prepare data cells for statcondfieldtrip:
% statcondfieldtrip expects each cell to be a data matrix with dimensions:
% (channels x freq x time x subjects) OR other variations. Here we will
% pass channels as explicit dim by setting 'chandim'.
% We'll construct: data{1} = subj x chan x freq x time -> then makefieldtripdata handles it.
cond1   = mean(power(:,:,:,:, [3 4]), 5);   % subj x chan x freq x time === Prev Con
cond2 = mean(power(:,:,:,:, [1 2]), 5);   % subj x chan x freq x time === Prev Inc

% permute to chan x freq x time x subj
dataA = permute(cond2, [2 3 4 1]);
dataB = permute(cond1,   [2 3 4 1]);

% call with chandim=1 (channels are in dimension 1)
stats = statcondfieldtrip({dataA, dataB}, ...
    'method','permutation', ...
    'naccu', 50, ...
    'neighbours', neighbours, ...
    'alpha', 0.05, ...
    'chandim',1, ...
    'freqdim',2, ...
    'timdim',3, ...
    'structoutput','on', ...
    'mcorrect','cluster', ...
    'clusteralpha',0.05, ...
    'clusterstatistic','maxsum'); % clusters only over channels

% outputs:
% stats  -> t (channels x freq x time)
% pvals  -> p (channels x freq x time)
mask = pvals < 0.05; % logical mask across chan x freq x time

% Example: plot grand-mean TF (ROI or global) and mask for a channel set
% Average contrast over some freq/time window, then mark significant channels
tidx = data.times >= 200 & data.times <= 400;
fidx = data.freqs >= 8 & data.freqs <= 12;

% mean effect per channel (mean over freq/time)
effect_chan = squeeze(mean(mean( (cond2 - cond1), 3), 4 )); % subj x chan -> but we want grandmean
grandmean_chan = squeeze(mean(effect_chan,1)); % chan x 1

% which channels have any significant cluster in window?
sig_chan = squeeze(any(any(mask(:, fidx, tidx), 2), 3)); % chan x 1 logical

figure;
topoplot(mean(cond2 - cond1,1), chanlocs, 'maplimits','absmax'); % or use your average
hold on;
topoplot(double(sig_chan).*grandmean_chan, chanlocs, 'style','blank', 'electrodes','pts', ...
         'emarker2', {find(sig_chan),'o','k',1,1});
title('Channels significant (8-12 Hz, 200-400 ms)');




