function d = getISCofRespresentation(cfg, d)

% evaluate input
if ~isfield(cfg, 'plot_sem'); cfg.plot_sem = true;end

cfg.timepoints = d.pairRep.included_time;
cfg.all_timepoints = d.pairRep.all_time;

% get number of time points
ntimepoints=numel(cfg.timepoints);

% loop through categories
for cate_num = 1:numel(cfg.categories)
    category = char(cfg.categories{cate_num});

    % loop through time points
    for tp=1:ntimepoints

        % preallocate RDM matrix
        RDMmat = nan(nchoosek(cfg.nTrials/2, 2), cfg.n);

        % loop through subjects
        for iSub = 1:length(cfg.subNums)
            subID = sprintf('sub-%0.3d', cfg.subNums(iSub));
            subID2 = strrep(subID, '-', '');

            % make a matrix with vectorized RDMs
            rdm = squeeze(d.pairRep.(subID2).rdm(:,:,tp));
            rdm(eye(size(rdm)) == 1) = 0;

            % filter for category
            if strcmp(category, 'bathroom')
                rdm = rdm(1:cfg.nTrials/2, 1:cfg.nTrials/2);
            else
                rdm = rdm(cfg.nTrials/2 + 1:end, cfg.nTrials/2 + 1:end);
            end

            RDMmat(:, iSub) = squareform(rdm);

        end

        % make IS-RDM
        [~, mat_out, ~] = make_RDM(RDMmat, cfg);
        if cfg.dissimilarity
            median_mat_out = 1 - mat_out;
        else
            median_mat_out = mat_out;
        end

        % take median and standard error
        median_mat_out(eye(size(median_mat_out)) == 1) = 0;
        medianISC = median(squareform(median_mat_out), 'omitnan');
        seISC = std(squareform(median_mat_out), 'omitnan') / sqrt(length(squareform(median_mat_out)));

        % store in structure
        d.ISC.([category,'_RDM']).pairRep(tp).RDM = mat_out;
        d.ISC.([category,'_RDM']).pairRep(tp).color = [0, 0, 0];
        d.ISC.([category,'_RDM']).pairRep(tp).name = (num2str(tp));
        d.ISC.medianISC.(category).pairRep(tp) = medianISC;
        d.ISC.medianISC_SE.(category).pairRep(tp) = seISC;

    end % time point
end % category

%% plotting
if cfg.plotting

    %% median ISC plot
    set(0, 'DefaultTextFontSize', cfg.FontSize);
    set(0, 'DefaultAxesFontSize', cfg.FontSize);
    set(0, 'DefaultTextFontName', cfg.FontName);
    set(0, 'DefaultAxesFontName', cfg.FontName);
    
    c = lines(numel(cfg.categories));
    x = d.pairRep.all_time(d.pairRep.included_time);

    fig=figure;
    hold on;

    for i=1:length(cfg.categories)
        y = d.ISC.medianISC.(cfg.categories{i}).pairRep;
        se = d.ISC.medianISC_SE.(cfg.categories{i}).pairRep;
        if cfg.plot_sem
            x2 = [x, fliplr(x)];
            inBetween_bath = [y+ se, fliplr(y - se)];
            fill(x2, inBetween_bath, c(i, :), 'FaceAlpha', 0.2, 'EdgeColor', 'none');
            hold on;
        end

        % plot median ISC
        h(i)=plot(x, y, 'color', c(i, :), 'LineWidth', 2);
        hold on;
    end

    set(gca, 'box', 'off');
    yline(0, '--', 'LineWidth', 1.5);
    xlabel('Time (s)');
    ylabel('Median ISC');
    xlim([min(x)-min(x)-0.2, max(x)+0.01])
    xline(0, '--k');  % Mark stimulus onset
    title('ISC of representation');
    legend(h, cfg.categories);

    if cfg.saving
        save_plot(fig, 'ISCofRepresentation', cfg.figPath);
    end

end
end