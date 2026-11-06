% Plot fixed effects from model-3 RSA or encoding/RT regressions.
% Accepts a fitted-regression structure and a t-value threshold. The
% model-3 RSA-by-congruency configuration exports the Figure 5 panel; other
% selections visualize exploratory full-model encoding regressions.
% Author: Davide Gheza


function plt = plot_RT_regress_RSAm3orEncoding(reg,tvalues)



if tvalues
    tolegend = ["", "intercept", "cued", "task_2", "task_3", "task_4", "side_2", "dstr-avg", "cng", "task2Xside_2", "task3Xside_2", "task4Xside_2", "cuedXcng", "dstrXcng"];

else
    tolegend = ["intercept", "cued", "task_2", "task_3", "task_4", "side_2", "dstr-avg", "cng", "task2Xside_2", "task3Xside_2", "task4Xside_2", "cuedXcng", "dstrXcng"]; 
end

toplot = [];
fn = fieldnames(reg);
for i = 1:2:length(fn)
        toplot = cat(3,toplot,reg.(fn{i+tvalues})); % 2nd dimension is prev_cng (0==post_inc, 1==post_con, 2=all)
end




    cmap = colororder();

    rng(10);


    times = [-800:10:1498]; % tf
 
    % critical t values for a 95% confidence interval centered on zero, given the sample n = 26 
    % i.e., (100(1 – α/2) percentile of the t-distribution with n degrees of freedom)
    t_95 = tinv(0.975, 26 - 1);
    
    
    





       
       




   %% cued and dstr - main RSA effect and interactions - REPORTED 'all' trials (ignore stratification by previous cng)
        
fontsize = 18;


if tvalues
    tolegend = ["", "cued", "... x congruency", "non-cued (average)", "... x congruency"];

else
    tolegend = ["cued", "... x congruency", "non-cued (average)", "... x congruency"]; 
end

toplot = [];
fn = fieldnames(reg);
for i = [3,23,13,25]
        toplot = cat(3,toplot,reg.(fn{i+tvalues})); % 2nd dimension is prev_cng (0==post_inc, 1==post_con, 2=all)
end

% match Fig 5D colors
custom_colors = [0 0 0;        % cued - black
                 0 0 0;        % cuedXcng - black
                 0 0.45 0.74;  % dstr - blue
                 0 0.45 0.74]; % dstrXcng - blue
        

titles = {'previous low congruent (0,1)', 'previous high congruent (2,3)', ''};
for prev_cng = 1:3

        figure
        sgtitle(['RT ~ target and distractors strength ' titles{prev_cng}], fontsize=20)
        if tvalues==1
            t_95_series = repmat(t_95,1,length(times));
            fill([times fliplr(times)], [t_95_series -t_95_series], [0.9 0.9 0.9], 'EdgeColor', 'none');
        end
        hold on
        for sbp = 1:size(toplot,3)
            if ismember(sbp,[2,4])
                plot(times, toplot(:,prev_cng,sbp), 'linewidth', 2, 'LineStyle',':', ...% dstr effects have dashed line
                    'Color', custom_colors(sbp,:));
            else
                plot(times, toplot(:,prev_cng,sbp), 'linewidth', 2, ...
                    'Color', custom_colors(sbp,:));
            end
        end
        set(gca,'FontSize', fontsize)
        plot(times,zeros(length(times),1),'color','black','LineStyle',':','linewidth',2)
        xlim([times(1) times(end)])
        if tvalues
            ylim([-12 10])
        else
            ylim([-0.05 0.05])
        end
        xlabel('Time (ms)')
        if tvalues  
            ylabel('RT~projection (t-value)')
        else
            ylabel('RT~projection (beta)')
        end
        legend(tolegend, 'Location', 'southeast', fontsize=12)
        legend('boxoff')
        legend('AutoUpdate','off')
        plot([times(1) times(end)], [0 0], ':k', 'LineWidth', 2);
        xline(0, ':k', 'LineWidth', 2, 'HandleVisibility', 'off');
        hold off
        

end
set(gca, 'Box', 'off')
        
    
% save pdf
fig = gcf;
fig.Units = 'inches';
fig.Position = [1 1 7 5.5];     
fig.PaperUnits = 'inches';
fig.PaperPosition = [0 0 7 5.5];
fig.PaperSize = [7 5.5];
if tvalues
    pdf_name = 'RT_regress_RSAm3__Fig5.pdf';
else
    pdf_name = 'RT_regress_RSAm3__betas.pdf';
end
exportgraphics(fig, pdf_name, 'ContentType', 'vector');
  


       



  %% cued and dstr - difference post-high - post-low cng
        
        fontsize = 18;
        fontsize2 = 16;
        % for poster
        fontsize = 22;
        fontsize2 = 19;
        cs = [];


if tvalues
    tolegend = ["", "cued", "dstr-avg", "cuedXcng", "dstrXcng"];

else
    tolegend = ["cued", "dstr-avg", "cuedXcng", "dstrXcng"]; 
end

toplot = [];
fn = fieldnames(reg);
for i = [3,13,23,25]
        toplot = cat(3,toplot,reg.(fn{i+tvalues})); % 2nd dimension is prev_cng (0==post_inc, 1==post_con, 2=all)
end

toplot = squeeze(toplot(:,2,:) - toplot(:,1,:)); 


figure
sgtitle('RT ~ dim-specific axis -- previous congruency high(2,3) - low(0,1)', fontsize=20)
if tvalues==1
    t_95_series = repmat(t_95,1,length(times));
    fill([times fliplr(times)], [t_95_series -t_95_series], [0.9 0.9 0.9], 'EdgeColor', 'none');
end
hold on
for sbp = 1:size(toplot,2)
    if ismember(sbp,[2,4])
        plot(times, toplot(:,sbp), 'linewidth', 2, 'LineStyle',':'); % dstr effects have dashed line
    else
        plot(times, toplot(:,sbp), 'linewidth', 2);
    end
end
set(gca,'FontSize', fontsize)
plot(times,zeros(length(times),1),'color','black','LineStyle',':','linewidth',2)
xlim([times(1) times(end)])
if tvalues
    ylim([-12 9])
else
    ylim([-0.05 0.05])
end
xlabel('Time (ms)')
ylabel('RT~projection')
legend(tolegend)
legend('AutoUpdate','off')
plot([times(1) times(end)], [0 0], ':k', 'LineWidth', 2);
hold off

        


end

