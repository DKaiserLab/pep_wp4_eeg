function d = getISCofRespresentation(cfg, d)

% evaluate input
if ~isfield(cfg, 'plot_sem'); cfg.plot_sem = true;end

cfg.timepoints = d.pairRep.included_time;
cfg.all_timepoints = d.pairRep.all_time;

% get number of ROI masks
ntimepoints=numel(cfg.timepoints);

% loop through categories
for cate_num = 1:numel(cfg.categories)
    category = char(cfg.categories{cate_num});

    % preallocate
    meanAcc = nan(length(d.pairRep.included_time), cfg.n);

    % get RDM of ISC type and category
    % loop through timepoints
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
            if strcmp(category, cfg.categories{1})
                rdm = rdm(1:cfg.nTrials/2, 1:cfg.nTrials/2);
            else
                rdm = rdm(cfg.nTrials/2 + 1:end, cfg.nTrials/2 + 1:end);
            end

            RDMmat(:, iSub) = squareform(rdm);
            meanAcc(tp, iSub) = mean(squareform(rdm), 'omitnan');

        end

        % make IS-RDM
        [~, mat_out, ~] = make_RDM(RDMmat, cfg);
        if cfg.dissimilarity
            median_mat_out = 1 - mat_out;
        else
            median_mat_out = mat_out;
        end

        median_mat_out(eye(size(median_mat_out)) == 1) = 0;
        medianISC = median(squareform(median_mat_out), 'omitnan');
        seISC = std(squareform(median_mat_out), 'omitnan') / sqrt(length(squareform(median_mat_out)));

        % store in structure
        d.ISC.([category,'_RDM']).pairRep(tp).RDM = mat_out;
        d.ISC.([category,'_RDM']).pairRep(tp).color = [0, 0, 0];
        d.ISC.([category,'_RDM']).pairRep(tp).name = (num2str(tp));
        d.ISC.medianISC.(category).pairRep(tp) = medianISC;
        d.ISC.medianISC_SE.(category).pairRep(tp) = seISC;
        d.meanAcc.(category).pairRep = meanAcc;

    end % timepoint
end % category

%% plotting
if cfg.plotting

    %% median ISC plot
    set(0, 'DefaultTextFontSize', cfg.FontSize);
    set(0, 'DefaultAxesFontSize', cfg.FontSize);
    set(0, 'DefaultTextFontName', cfg.FontName);
    set(0, 'DefaultAxesFontName', cfg.FontName);
    c = lines(numel(cfg.categories));
    x = 1:ntimepoints;

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
    xlabel('time');
    ylabel('median ISC');
    xlim([min(x)-min(x)-0.01, max(x)+0.01])
    xticks([1, 10:10:60]);
    xticklabels(cfg.all_timepoints(cfg.timepoints([1, 10:10:60])));
    xline(0, '--k');  % Mark stimulus onset
    title('ISC of representation');
    %     xticks(1:length(cfg.timepoints));
    %     xticklabels(cfg.all_timepoints(cfg.timepoints));
    legend(h, cfg.categories);

    if cfg.saving
        save_plot(fig, 'ISC of Representation', cfg.figPath);
    end


    %% mean accuracy plot
    fig = figure;

    % get group stats
    x = d.pairRep.all_time(d.pairRep.included_time);
    y1 = mean(d.meanAcc.bathroom.pairRep, 2);
    y2 = mean(d.meanAcc.kitchen.pairRep, 2);
    e1 = std(d.meanAcc.bathroom.pairRep, [], 2)/sqrt(size(d.meanAcc.bathroom.pairRep, 2));
    e2 = std(d.meanAcc.kitchen.pairRep, [], 2)/sqrt(size(d.meanAcc.kitchen.pairRep, 2));

    % plot
    b(1) = boundedline(x, y1, e1, 'color', c(1, :), 'alpha', 'LineWidth', 2);
    b(2) = boundedline(x, y2, e2,'color', c(2, :), 'alpha', 'LineWidth', 2);

    xlim([min(x)-min(x)-0.01, max(x)+0.01])
    xline(0, '--k');  % Mark stimulus onset
    yline(0.5, '--k');  % Mark cahnce level
    xlabel('Time (s)');
    ylabel('Accuracy');
    legend(b, {'Bathroom', 'Kitchen'})
    title('Mean pairwise decoding');

    if cfg.saving
        save_plot(fig, 'Mean-pairwise-decoding', cfg.figPath);
    end

end
end