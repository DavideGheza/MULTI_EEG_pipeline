 %% cluster based permutation testing 

 function hplot = RSA_permutest(betas,times,stratify)

 %  % if downsample
 % betas = movmean(betas,[3,3],2);
 % betas=downsample(betas',3)'; % downsample by a factor (columnwise)
 % times = downsample(times,3);
 
 if nargin>2
     if strcmp(stratify,'PC')
        PC1 = 3;
        PC2 = 2;
        PC_labels = [{'Neu1'}, {'MI'}, {'MC'}, {'Neu2'}];
     end
 else
    PC1 = 2;
    PC2 = 3;
    PC_labels = [{'Neu1'}, {'MI'}, {'MC'}, {'Neu2'}];
 end
    r = [];
    times2save = times; % from decoding_tf_on_sessions.m
    times2test = [-400 1200];  % subset to smaller time window (can increase sensitivity)
    timestested = times2save(dsearchn(times2save',times2test(1)):dsearchn(times2save',times2test(2)));
    cond1 = squeeze(betas(:,PC1,(dsearchn(times2save',times2test(1)):dsearchn(times2save',times2test(2)))))';
%         cond2 = squeeze(mean([betas(:,1,(dsearchn(times2save',times2test(1)):dsearchn(times2save',times2test(2)))),betas(:,4,(dsearchn(times2save',times2test(1)):dsearchn(times2save',times2test(2))))],2))'; %mean over Neu_1 and Neu_2
    cond2 = squeeze(betas(:,PC2,(dsearchn(times2save',times2test(1)):dsearchn(times2save',times2test(2)))))';
%     cond2 = ones(size(cond1,1),size(cond1,2))*50; % a flat line at 50, for all participant. Use to test cond1 against chance
%     cond2 = zeros(size(cond1,1),size(cond1,2)); % a flat line at 0, for all participant.
    [r.clusters, r.p_values, r.t_sums, r.perm_dstrb] = permutest(cond1,cond2,true,0.05,[],true); %2nd logical: false=one sided (positive)
    for i = 1:length(r.clusters)
        r.clusters_timing(i) = {timestested(r.clusters{i})};
    end

    %plot
    figure;
    title(['RSA - ' PC_labels{PC1} '-' PC_labels{PC2}], 'FontSize', 18)
    hold on
    set(gca,'XTick',timestested(1):200:timestested(end),'FontSize', 14');
    plot(timestested,mean(cond1,2)'-mean(cond2,2)','LineWidth',2,'LineStyle','-','color', 'blue');
%         plot([timestested(1) timestested(end)], [50 50], 'k-')
    ylim([-0.1 0.4])
    xlabel('Time [ms]', 'FontSize', 18)
    ylabel('\Deltabeta [RSA]', 'FontSize', 18)
    hold off
    dasz = [0 0];
    dasw = ylim; % set as the current y-axis limits
    line(dasz,dasw,'Color',[0.7 0.7 0.7],'LineStyle','--','LineWidth',2)    % set dashed y-axis vertical lines (onset cue -200ms, onset stim 0ms)
    dasz = [-200 -200];
    line(dasz,dasw,'Color',[0.7 0.7 0.7],'LineStyle','--','LineWidth',2)
    dasx = [timestested(1) timestested(end)];                % set dashed x-axis horizontal line
    dasy = [0 0];
    line(dasx,dasy,'Color',[0.7 0.7 0.7],'LineStyle','--','LineWidth',2)
    if r.p_values(1)<0.05
        sig_line = [r.clusters_timing{1,1}(1) r.clusters_timing{1,1}(end)];
        line(sig_line,dasy,'Color',[0.1 0.1 0.1],'LineStyle','-','LineWidth',3)
    end
        


 end