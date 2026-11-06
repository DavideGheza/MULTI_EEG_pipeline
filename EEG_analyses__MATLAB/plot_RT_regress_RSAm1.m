% Plot fixed effects from model-1 RSA/RT multilevel regressions.
% Accepts a fitted-regression structure and a t-value threshold and exports
% the Figure 4 panel when the corresponding option is selected below.
% Author: Davide Gheza

function plt = plot_RT_regress_RSAm1(reg,tvalues)

    
    cmap = colororder();

    rng(10);


    times = [-800:10:1498]; % tf
 
    % critical t values for a 95% confidence interval centered on zero, given the sample n = 26 
    % i.e., (100(1 – α/2) percentile of the t-distribution with n degrees of freedom)
    t_95 = tinv(0.975, 26 - 1);
    
    
    
    % beta values
    if tvalues==0
        intercept = reg.intercept;
        sr = reg.sr;         
        sr_cng = reg.slope_srXcng;
        ru = reg.ru;         
        ru_cng = reg.slope_ruXcng;
        si = reg.si;         
        si_cng = reg.slope_siXcng;
        
        cng = reg.cng;
    % t values
    else
        intercept = reg.intercept_t;
        sr = reg.sr_t;         
        sr_cng = reg.slope_srXcng_t;
        ru = reg.ru_t;         
        ru_cng = reg.slope_ruXcng_t;
        si = reg.si_t;         
        si_cng = reg.slope_siXcng_t;
        
        cng = reg.cng_t;

    end




        %% cued - main RSA effect and interactions

        tolegend = ["Task", "Spatial attention", "Feature", "... x congruency"];
        
        fontsize = 18;
        fontsize2 = 16;
        % for poster
        fontsize = 22;
        fontsize2 = 19;
        cs = [];

        custom_colors = [                   % red, yellow, blue
             1, 0, 0;        % RU
             1, 0.753, 00;   % SI
             0, 0.439, 0.753;% SR
             1, 0, 0;        % RU x congruency
             1, 0.753, 0;    % SI x congruency
             0, 0.439, 0.753]; % SR x congruency

        
        figure
        sgtitle('RT ~ RSA - m1 cued', fontsize=20)
        if tvalues==1
            t_95_series = repmat(t_95,1,length(times));
            fill([times fliplr(times)], [t_95_series -t_95_series], [0.9 0.9 0.9], 'EdgeColor', 'none', 'HandleVisibility', 'off');
        end
        hold on
        toplot = [ru, si, sr, ru_cng, si_cng, sr_cng];
        h_main = gobjects(1,3);
        h_int = gobjects(1,1);
        for sbp = 1:size(toplot,2)
            if ismember(sbp,[4,5,6])
                htmp = plot(times, toplot(:,sbp), 'linewidth', 2, 'LineStyle', ':', 'Color', custom_colors(sbp,:));
                if sbp == 4
                    h_int = htmp;
                end
            else
                h_main(sbp) = plot(times, toplot(:,sbp), 'linewidth', 2, 'Color', custom_colors(sbp,:));
            end
        end
        h_int = plot(nan, nan, 'LineWidth', 2, 'LineStyle', ':', 'Color', [0.5 0.5 0.5]);
        set(gca,'FontSize', fontsize)
        plot(times,zeros(length(times),1),'color','black','LineStyle',':','linewidth',2,'HandleVisibility','off')
        xlim([times(1) times(end)])
        if tvalues
            ylim([-14 4])
        else
            ylim([-0.025 0.01])
        end
        xlabel('Time (ms)')
        if tvalues
            ylabel('RT~RSA (t-value)')
        else
            ylabel('RT~RSA (beta)')
        end
        legend([h_main h_int], tolegend, 'Location', 'southwest')
        legend('boxoff')
        legend('AutoUpdate','off')
        plot([times(1) times(end)], [0 0], ':k', 'LineWidth', 2, 'HandleVisibility', 'off');
        xline(0, ':k', 'LineWidth', 2, 'HandleVisibility', 'off');
        hold off
        set(gca, 'Box', 'off')
        % save pdf
        fig = gcf;
        fig.Units = 'inches';
        fig.Position = [1 1 9 6];      % width=9, height=6
        fig.PaperUnits = 'inches';
        fig.PaperPosition = [0 0 9 6];
        fig.PaperSize = [9 6];
        if tvalues
            pdf_name = 'RT_regress_RSAm1__Fig4.pdf';
        else
            pdf_name = 'RT_regress_RSAm1__betas.pdf';
        end
        exportgraphics(fig, pdf_name, 'ContentType', 'vector');

        

       


      


       






end

