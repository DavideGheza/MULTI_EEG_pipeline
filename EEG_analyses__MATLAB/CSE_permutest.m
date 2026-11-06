 %% cluster based permutation testing 
 % based off RSA_permutest - 04/22/2025
 
 function hplot = CSE_permutest(betas,times,labels)

 betas = permute(betas,[3,1,2]); % sbj, cond, times
 % Note: for CSE we test for suppression, i.e., one-sided HP. Flip sign of
 % timeseries since permutest test for positive deflection

 %  % if downsample
 % betas = movmean(betas,[3,3],2);
 % betas=downsample(betas',3)'; % downsample by a factor (columnwise)
 % times = downsample(times,3);
 
    cnd1 = 1;
    cnd2 = nan;

    if nargin<3
        labels = [{'On'}, {'Off'}, {'zero'}];
        if isnan(cnd2)
            labels = labels([1,3,2]);
        end
    end


    r = [];
    times2save = times; % from decoding_tf_on_sessions.m
    times2test = [-400 1200];  % subset to smaller time window (can increase sensitivity)
    timestested = times2save(dsearchn(times2save',times2test(1)):dsearchn(times2save',times2test(2)));
    cond1 = squeeze(betas(:,cnd1,(dsearchn(times2save',times2test(1)):dsearchn(times2save',times2test(2)))))';
    if isnan(cnd2)
        cond2 = zeros(size(betas,3),size(betas,1)); % a flat line at 0, for all participants
        cond2 = squeeze(cond2(dsearchn(times2save',times2test(1)):dsearchn(times2save',times2test(2)),:));
    else
        cond2 = squeeze(betas(:,cnd2,(dsearchn(times2save',times2test(1)):dsearchn(times2save',times2test(2)))))';
    end

    [r.clusters, r.p_values, r.t_sums, r.perm_dstrb] = permutest(cond1*-1,cond2*-1,true,0.05,[],false); %1st logical: true=paired. 2nd logical: false=one sided (positive)
    for i = 1:length(r.clusters)
        r.clusters_timing(i) = {timestested(r.clusters{i})};
    end

    %plot
    endtitle = "";
    if(~isnan(cnd2)) 
        endtitle = labels{cnd2};
    end
    figure;
    title(['projection - ' labels{cnd1} '-' endtitle], 'FontSize', 18)
    hold on
    set(gca,'XTick',timestested(1):200:timestested(end),'FontSize', 14');
    plot(timestested,mean(cond1,2)'-mean(cond2,2)','LineWidth',2,'LineStyle','-','color', 'blue');
%         plot([timestested(1) timestested(end)], [50 50], 'k-')
    ylim([-0.02 0.01])
    xlabel('Time [ms]', 'FontSize', 18)
    ylabel('\Delta', 'FontSize', 18)
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