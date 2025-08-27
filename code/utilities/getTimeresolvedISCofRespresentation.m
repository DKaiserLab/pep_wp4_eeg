function d = getTimeresolvedISCofRespresentation(cfg, d)

cfg.plotting = true;
cfg.ISC_types = {'pairRep'};

% loop through ISC types
for ISC_type = cfg.ISC_types
    cfg.ISC_type = char(ISC_type);

    % loop through categories
    for cate_num = 1:numel(cfg.categories)
        category = char(cfg.categories{cate_num});

        % clear struct
        d.([category,'_RDM'])= struct;

        % preallocate
        medianISC = nan(1, length(d.pairRep.included_time));
        tpISC_RDM = nan(cfg.n, cfg.n, length(d.pairRep.included_time));
        meanAcc = nan(length(d.pairRep.included_time), cfg.n);

        for iTp = 1:length(d.pairRep.included_time)

            % preallocate RDM matrix
            RDMmat = nan(nchoosek(cfg.nTrials/2, 2), cfg.n);

            % loop through subjects
            for iSub = 1:length(cfg.subNums)
                subID = sprintf('sub-%0.3d', cfg.subNums(iSub));
                subID2 = strrep(subID, '-', '');

                % make a matrix with vectorized RDMs
                rdm = squeeze(d.(cfg.ISC_type).(subID2).rdm(:, :, iTp));
                rdm(eye(size(rdm)) == 1) = 0;

                % filter for category
                if strcmp(category, cfg.categories{1})
                    rdm = rdm(1:cfg.nTrials/2, 1:cfg.nTrials/2);
                else
                    rdm = rdm(cfg.nTrials/2 + 1:end, cfg.nTrials/2 + 1:end);
                end
                RDMmat(:, iSub) = squareform(rdm);
                meanAcc(iTp, iSub) = mean(squareform(rdm), 'omitnan');

            end

            % make IS-RDM
            [~, matOut, ~] = make_RDM(RDMmat, cfg);
            if cfg.dissimilarity
                median_mat_out = 1 - matOut;
                median_mat_out(eye(size(median_mat_out)) == 1) = 0;
                medianISC(iTp) = median(squareform(median_mat_out), 'omitnan');
            else
                median_mat_out = matOut;
                median_mat_out(eye(size(median_mat_out)) == 1) = 0;
                medianISC(iTp) = median(squareform(median_mat_out), 'omitnan');
            end

            % store in timeresloved ISC matrix
            tpISC_RDM(:, :, iTp) = matOut;

        end

        % store in structure
        d.([category,'_RDM']).(cfg.ISC_type).RDM = tpISC_RDM;
        d.([category,'_RDM']).(cfg.ISC_type).color = [0, 0, 0];
        d.([category,'_RDM']).(cfg.ISC_type).name = 'timeresolvedISC';
        d.medianISC.(category).(cfg.ISC_type) = medianISC;
        d.timePoints.(category).(cfg.ISC_type) = d.pairRep.included_time;
        d.meanAcc.(category).(cfg.ISC_type) = meanAcc;


    end % category
    if cfg.plotting
        % plot pairwise decoding
        if strcmp(cfg.ISC_type, 'pairRep')
            figure
            % get group stats
            x = d.pairRep.all_time(d.pairRep.included_time);
            y1 = mean(d.meanAcc.bathroom.(cfg.ISC_type), 2);
            y2 = mean(d.meanAcc.kitchen.(cfg.ISC_type), 2);
            e1 = std(d.meanAcc.bathroom.(cfg.ISC_type), [], 2)/sqrt(size(d.meanAcc.bathroom.(cfg.ISC_type), 2));
            e2 = std(d.meanAcc.kitchen.(cfg.ISC_type), [], 2)/sqrt(size(d.meanAcc.kitchen.(cfg.ISC_type), 2));

            % plot
            b(1) = boundedline(x, y1, e1, 'cmap', [1, 0, 0], 'alpha');
            b(2) = boundedline(x, y2, e2,'cmap', [0, 0, 1], 'alpha');

            xlim([min(x)-min(x)-0.01, max(x)+0.01])
            xline(0, '--k');  % Mark stimulus onset
            yline(0.5, '--k');  % Mark cahnce level
            xlabel('Time (s)');
            ylabel('Accuracy');
            legend(b, {'Bathroom', 'Kitchen'})
            title('Mean pairwise decoding');
        end

        % plot median ISC
        figure
        hold on
        x = d.pairRep.all_time(d.pairRep.included_time);
        y1 = d.medianISC.bathroom.(cfg.ISC_type);
        y2 = d.medianISC.kitchen.(cfg.ISC_type);

        p(1) = plot(x, y1, 'r');
        p(2) = plot(x, y2, 'b');
        xlim([min(x)-min(x)-0.01, max(x)+0.01])
        xline(0, '--k');  % Mark stimulus onset
        yline(0, '--k');  % Mark zero
        xlabel('Time (s)');
        ylabel('Median ISC');
        legend(p, {'Bathroom', 'Kitchen'})
        title('Mean pairwise decoding');
    end
end % isc type

end