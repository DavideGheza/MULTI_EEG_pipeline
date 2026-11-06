%% Plot group-level time-frequency/RT regression results
% Loads the saved 2025 regression slopes and produces time-frequency maps
% and topographies for stimulus- or response-locked data.
% Author: Davide Gheza

clc
clear

%% 1. Load settings


%Select onset or response events
onset = false;
if onset
    suffix_eeg = 'onset_';
    suffix = '';
else
    suffix_eeg = 'response_';
    suffix = 'response_';  
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

%% Load RT~tf regression betas ('power')

load([sets.sessions_data_tf_dir...
        'rt_power_regression_resplocked2025_slope.mat']);
% sbj, chans, freqs, times - 26    61    30    76
if onset
    times_regression = -50:10:700; % for stim-locked data
else
    times_regression = -700:10:50; % for response-locked data. 
end
data.freqs = logspace(log10(1.5),log10(40),30); % WARNING - hardcoded as per tf 2025 sets




%% average for single plot - 2025

chantoplot = {'Cz' 'FCz' 'FC1' 'FC2' 'C1' 'C2' 'CPz' 'CP1' 'CP2'};  % fronto-central channels show positive betas
% chantoplot = {'CPz' 'CP1' 'CP2' 'Pz' 'P1' 'P2' 'POz' 'P3' 'P4'}; % parieto-occipital channels show negative betas
chan = ismember(labels(:), chantoplot);

% tf plot
cmap = colormap(parula);
colormap parula

figure

timein = -700;   
timeout = 50;  
freqin = 1.5;           
freqout = 10;         
clim = [-0.3 0.3];
% clim = [-1 1];
ymin_idx = dsearchn(data.freqs', freqin);
ymax_idx = dsearchn(data.freqs', freqout);
tick_idx = round(linspace(ymin_idx,ymax_idx,8));
valid_idx = tick_idx(tick_idx>=ymin_idx & tick_idx<=ymax_idx);

%%%
spdata = mean(power(:,chan,:,:,:),[2,5]); % subset of chans
% spdata = mean(power(:,:,:,:,:),[2,5]); % (all channels)

imagesc(times_regression,[],squeeze(mean(spdata,1)));   
cbh = colorbar;
set(get(cbh,'xlabel'),'string','Regression Coefficient ({\it\beta})','FontSize',16);
set(gca,'clim',clim,'xlim',[timein timeout],'YDir','normal'); 
xlabel('Time (ms)','FontSize',17,'FontName','Arial')
ylabel('Frequency (hz)','FontSize',17,'FontName','Arial')
set(gca, 'YTick', valid_idx, 'YTickLabel', round(data.freqs(valid_idx),1),'FontSize', 16);
ylim([ymin_idx ymax_idx]);
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
%%%

set(gcf,'Units','inches');
screenposition = get(gcf,'Position');
set(gcf,...
    'PaperPosition',[0 0 screenposition(3:4)],...
    'PaperSize',[screenposition(3:4)]);
print -dpdf -vector RT_reg_resplocked



         
%%




% scalp topography
% Theta
freqshown = [3 8];
timeshown = [-200 -0];

timeindexes = dsearchn(times_regression',timeshown(1)):dsearchn(times_regression',timeshown(2));
freqindexes = dsearchn(data.freqs',freqshown(1)):dsearchn(data.freqs',freqshown(2));
        
topospdata = mean(power(:,:,freqindexes,timeindexes,:),5);

figure;

% average on the 1st (sbj), 4th (freq) and 5th (time) dimensions
topoplot(mean(mean(topospdata,1),4),...
            chanlocs,...
            'plotchans',[1:61],...
            'maplimits',[-0.03 0.03], ...
            'colormap', cmap);    %,'maplimits',[-4 4] / 'maxmin'  ,'electrodes','ptslabels'
cbh = colorbar;
set(get(cbh,'xlabel'),'string','Regression Coefficient ({\it\beta})','FontSize',16);
set(gca,'FontSize', 16)
% title(['Grandaverage within ' num2str(timeshown) 'ms & ' num2str(freqshown) 'Hz'],'Fontsize', 16,'FontName','Arial')
title(['Grandaverage within 3-8 Hz'],'Fontsize', 20,'FontName','Arial')

set(gcf,'Units','inches');
screenposition = get(gcf,'Position');
set(gcf,...
    'PaperPosition',[0 0 screenposition(3:4)],...
    'PaperSize',[screenposition(3:4)]);
print -dpdf -vector RT_reg_resplocked_topo_theta_min200_0



% Alpha
freqshown = [8 13];

timeindexes = dsearchn(times_regression',timeshown(1)):dsearchn(times_regression',timeshown(2));
freqindexes = dsearchn(data.freqs',freqshown(1)):dsearchn(data.freqs',freqshown(2));
        
topospdata = mean(power(:,:,freqindexes,timeindexes,:),5);

figure;

% average on the 1st (sbj), 4th (freq) and 5th (time) dimensions
topoplot(mean(mean(topospdata,1),4),...
            chanlocs,...
            'plotchans',[1:61],...
            'maplimits',[-0.05 0.05], ...
            'colormap', cmap);    %,'maplimits',[-4 4] / 'maxmin'  ,'electrodes','ptslabels'
cbh = colorbar;
set(get(cbh,'xlabel'),'string','Regression Coefficient ({\it\beta})','FontSize',16);
set(gca,'FontSize', 16)
% title(['Grandaverage within ' num2str(timeshown) 'ms & ' num2str(freqshown) 'Hz'],'Fontsize', 16,'FontName','Arial')
title(['Grandaverage within 8-13 Hz'],'Fontsize', 17,'FontName','Arial')





