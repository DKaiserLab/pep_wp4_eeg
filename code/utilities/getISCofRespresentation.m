function d = getISCofRespresentation(cfg, d)

if ~isfield(cfg, 'plot_sem'); cfg.plot_sem = true;end

cfg.timepoints = d.pairRep.included_time;
cfg.all_timepoints = d.pairRep.all_time;

% get number of ROI masks
ntimepoints=numel(cfg.timepoints);

% loop through categories
for cate_num = 1:numel(cfg.categories)
    category = char(cfg.categories{cate_num});

    % loop through ISC types
    for ISC_type = cfg.ISC_types
        cfg.ISC_type = char(ISC_type);

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
                rdm = squeeze(d.pairRep.(subID2).rdm(:,:,tp)); %rdm = squeeze(d.(cfg.ISC_type).(subID2).(tp).rdm);
                rdm(eye(size(rdm)) == 1) = 0;

                % filter for category
                if strcmp(category, cfg.categories{1})
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
                median_mat_out(eye(size(median_mat_out)) == 1) = 0;
                medianISC = median(squareform(median_mat_out), 'omitnan');
                seISC = std(squareform(median_mat_out), 'omitnan') / sqrt(length(squareform(median_mat_out)));
            else
                median_mat_out = mat_out;
                median_mat_out(eye(size(median_mat_out)) == 1) = 0;
                medianISC = median(squareform(median_mat_out), 'omitnan');
                seISC = std(squareform(median_mat_out), 'omitnan') / sqrt(length(squareform(median_mat_out)));
            end

            % store in structure 
            d.ISC.([category,'_RDM']).(cfg.ISC_type)(tp).RDM = mat_out;
            d.ISC.([category,'_RDM']).(cfg.ISC_type)(tp).color = [0, 0, 0];
            d.ISC.([category,'_RDM']).(cfg.ISC_type)(tp).name = (num2str(tp));
            d.ISC.medianISC.(category).(cfg.ISC_type)(tp) = medianISC;
            d.ISC.medianISC_SE.(category).(cfg.ISC_type)(tp) = seISC;

        end % timepoint
    end % isc type
end % category

if cfg.plotting

    set(0, 'DefaultTextFontSize', cfg.FontSize);
    set(0, 'DefaultAxesFontSize', cfg.FontSize);
    set(0, 'DefaultTextFontName', cfg.FontName);
    set(0, 'DefaultAxesFontName', cfg.FontName);
    c = {'red', 'blue'};
    x = 1:ntimepoints;
    figure;
    hold on;

    for i=1:length(cfg.categories)


        y = d.ISC.medianISC.(cfg.categories{i}).pairRep;
        se = d.ISC.medianISC_SE.(cfg.categories{i}).pairRep;
        if cfg.plot_sem
            x2 = [x, fliplr(x)];
            inBetween_bath = [y+ se, fliplr(y - se)];
            fill(x2, inBetween_bath, c{i}, 'FaceAlpha', 0.2, 'EdgeColor', 'none');
            hold on;
        end
        % plot median ISC
        h(i)=plot(x, y, 'color', c{i}, 'LineWidth', 2);
        hold on;
    end

    set(gca, 'box', 'off');
    yline(0, '--', 'LineWidth', 1.5);
    xlabel('time');
    ylabel('median ISC');
    xticks([1, 10:10:60]);
    xticklabels(cfg.all_timepoints(cfg.timepoints([1, 10:10:60])));
%     xticks(1:length(cfg.timepoints));
%     xticklabels(cfg.all_timepoints(cfg.timepoints));
    legend(h, cfg.categories)
end
end