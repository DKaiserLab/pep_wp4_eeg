function d = getISCofRespresentation(cfg, d)
if ~isfield(cfg, 'plot_sem'); cfg.plot_sem = true;end
if ~isfield(cfg, 'correlation_type'); cfg.correlation_type = 'Spearman';end
if ~isfield(cfg, 'save_name'); cfg.save_name = 'compare_tp_RDMs_to_stimuli_RDMs';end
if ~isfield(cfg, 'alpha'); cfg.alpha = 0.05;end
if ~isfield(cfg, 'n_permutations'); cfg.n_permutations = 5000;end
if ~isfield(cfg, 'permutationtest'); cfg.permutationtest = false;end
if ~isfield(cfg, 'smoothing_window'); cfg.smoothing_window = 1;end


%first_pos_tp = find(d.pairRep.all_time(d.pairRep.included_time)>=0, 1, 'first');
cfg.plot_rdm = false;
ntimepoints= length(d.pairRep.included_time);
n = cfg.n;
rng(1);
signs_all = randi([0 1], n, cfg.n_permutations)*2 - 1;
null_distrib = struct;
d.ISC.stats = struct;
d.ISC.stats.cluster = struct;

if isempty(gcp('nocreate'))
    parpool(8);
end

% loop through categories
for cate_num = 1:numel(cfg.categories)
    category = char(cfg.categories{cate_num});
    null_distrib.(category) = nan(ntimepoints, cfg.n_permutations);
    d.ISC.stats.cluster.(category) = struct;
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
        medianISC = median(squareform(median_mat_out), 'omitnan'); % median VZ-permutieren, ganze zeile+spalte
        seISC = std(squareform(median_mat_out), 'omitnan') / sqrt(length(squareform(median_mat_out)));

        if cfg.permutationtest
            perm_median = nan(1, cfg.n_permutations);
            parfor perm = 1:cfg.n_permutations
                signs = signs_all(:,perm);
                signs_mat = diag(signs)
                perm_median_mat_out = signs_mat * median_mat_out * signs_mat;
                perm_medianISC = median(squareform(perm_median_mat_out), 'omitnan');
                perm_median(perm) = perm_medianISC;
            end
            null_distrib.(category)(tp,:) = perm_median;
        end
        

        % store in structure
        d.ISC.([category,'_RDM']).pairRep(tp).RDM = mat_out;
        d.ISC.([category,'_RDM']).pairRep(tp).color = [0, 0, 0];
        d.ISC.([category,'_RDM']).pairRep(tp).name = (num2str(tp));
        d.ISC.medianISC.(category).pairRep(tp) = medianISC;
        d.ISC.medianISC_SE.(category).pairRep(tp) = seISC;

    end % time point
    if cfg.permutationtest
        d.ISC.stats.cluster.(category) = stats_cluster(d.ISC.medianISC.(category).pairRep, null_distrib.(category), cfg, d);
    end
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
    h = gobjects(1, length(cfg.categories));
    max_ISC = round(max([d.ISC.medianISC.bathroom.pairRep, d.ISC.medianISC.kitchen.pairRep]), 2);
    min_ISC = round(min([d.ISC.medianISC.bathroom.pairRep, d.ISC.medianISC.kitchen.pairRep]), 2);
    fig=figure;
    hold on;

    for i=1:length(cfg.categories)
        category = cfg.categories{i};
        %y = d.ISC.medianISC.(category).pairRep;
        y = smoothdata(d.ISC.medianISC.(category).pairRep,'movmean', cfg.smoothing_window);
        se = d.ISC.medianISC_SE.(category).pairRep;

        % plot median ISC
        h(i)=plot(x, y, 'color', c(i, :), 'LineWidth', 2);
        if cfg.plot_sem
            x2 = [x, fliplr(x)];
            inBetween_bath = [y+ se, fliplr(y - se)];
            fill(x2, inBetween_bath, c(i, :), 'FaceAlpha', 0.2, 'EdgeColor', 'none');
            %hold on;
        end
        
        sigMat = false(1, ntimepoints);
        if ~isempty(d.ISC.stats.cluster.(category).obs_cluster_stats)
            for p=1:length(d.ISC.stats.cluster.(category).pvals_cluster)
                if d.ISC.stats.cluster.(category).pvals_cluster(p) <= cfg.alpha
                    idx = d.ISC.stats.cluster.(category).obs_clusters.PixelIdxList{p};
                    sigMat(idx)=1;
                end
            end
        end

        if cfg.permutationtest
            pos = max_ISC + i*0.002;
            plot(x(sigMat), repmat(pos, 1, sum(sigMat)), ...
                'color', c(i, :), 'marker' ,'O', 'MarkerFaceColor', c(i, :) ,'MarkerSize', 5, 'LineStyle','none');
        end
    end

    set(gca, 'box', 'off');
    yline(0, '--', 'LineWidth', 1.5);
    xlabel('Time (s)');
    ylabel('Median ISC');
    xlim([min(x)-min(x)-0.2, max(x)+0.01])
    ylim([min_ISC, max_ISC+0.006]);
    xline(0, '--k');  % Mark stimulus onset
    title('ISC of representation');
    legend(h, cfg.categories);

    set(gca, 'LineWidth', 1, 'FontName', cfg.FontName, 'FontSize', cfg.FontSize, 'FontWeight', 'bold');
    set(gca, 'box', 'off');

    if cfg.saving
        save_plot(fig, 'ISCofRepresentation', cfg.figPath);
    end

end

if ~isempty(gcp('nocreate'))
    delete(gcp('nocreate'));
end

end