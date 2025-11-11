% ============================================
% PHASE III — VISUALISATION FOR GNU OCTAVE
% No GUI popups, all PNG only
% ============================================

% Load dataset
raw = csvread("horse_No_scaling.csv");

cols = {"surgery", "age", "rectal_temp", "pulse", "respiratory_rate", ...
        "temp_of_extremities", "peripheral_pulse", "mucous_membrane", ...
        "capillary_refill_time", "pain", "peristalsis", ...
        "abdominal_distention", "nasogastric_tube", ...
        "nasogastric_reflux", "nasogastric_reflux_ph", ...
        "rectal_exam_feces", "abdomen", "packed_cell_volume", ...
        "total_protein", "abdomo_appearance", "abdomo_protein", ...
        "outcome", "surgical_lesion", "cp_data", ...
        "lesion_1_site", "lesion_1_type", "lesion_1_subtype", ...
        "lesion_1_specific", "num_lesions"};

for i = 1:length(cols)
  data.(cols{i}) = raw(:, i);
end

% --------------------------------------------
% FOLDER CREATION
% --------------------------------------------

mkdir("central_tendency_png");
mkdir("histograms_png");
mkdir("kde_png");
mkdir("heatmap_png");

% --------------------------------------------
% CENTRAL TENDENCY + DISPERSION
% --------------------------------------------

numeric_cols = {"rectal_temp", "pulse", "respiratory_rate", ...
                "packed_cell_volume", "total_protein", "abdomo_protein"};

fid = fopen("central_tendency_png/stats.txt", "w");

fprintf(fid, "Central Tendency and Dispersion\n\n");

for i = 1:length(numeric_cols)
    col = data.(numeric_cols{i});
    fprintf(fid, "%s: mean=%.3f, median=%.3f, var=%.3f, std=%.3f, range=%.3f\n", ...
            numeric_cols{i}, mean(col), median(col), var(col), std(col), max(col)-min(col));
end

fclose(fid);

% --------------------------------------------
% HISTOGRAMS (NO POPUPS)
% --------------------------------------------

for i = 1:length(numeric_cols)
    figure("visible", "off");
    hist(data.(numeric_cols{i}), 20);
    title(["Histogram of ", numeric_cols{i}]);
    xlabel(numeric_cols{i});
    ylabel("Frequency");
    print(strcat("histograms_png/", numeric_cols{i}, "_hist.png"), "-dpng");
end

% --------------------------------------------
% MANUAL KDE FUNCTION
% --------------------------------------------

function [xgrid, pdf_vals] = manual_kde(x)
    x = x(:);
    n = length(x);
    sigma = std(x) * (4/(3*n))^(1/5);             
    xgrid = linspace(min(x), max(x), 300);
    pdf_vals = zeros(size(xgrid));

    for j = 1:length(xgrid)
        pdf_vals(j) = mean(normpdf(xgrid(j), x, sigma));
    end
end

% --------------------------------------------
% KDE PLOT FOR TOP 5 FEATURES
% --------------------------------------------

top5 = {"packed_cell_volume", "num_lesions", "total_protein", "pulse", "abdominal_distention"};
colors = lines(length(top5));

figure("visible", "off"); hold on;

for i = 1:length(top5)
    col = data.(top5{i});
    [xg, pdfv] = manual_kde(col);
    plot(xg, pdfv, "linewidth", 2, "color", colors(i,:));
end

legend(top5, "location", "northeast");
title("KDE Distribution of Top 5 Discriminant Features");
xlabel("Value");
ylabel("Density");
hold off;

print("kde_png/kde_top5.png", "-dpng");

% --------------------------------------------
% CORRELATION HEATMAP (NO POPUPS)
% --------------------------------------------

C = corrcoef(raw);

figure("visible", "off");
imagesc(C);
colorbar;
colormap("jet");
title("Correlation Heatmap (All Encoded Numeric Features)");

set(gca, "XTick", 1:length(cols), "XTickLabel", cols);
set(gca, "YTick", 1:length(cols), "YTickLabel", cols);
xtickangle(45);

print("heatmap_png/correlation_heatmap.png", "-dpng");



% ============================================
%         VIOLIN STYLE PLOTS (OUTCOME GROUPS)
% ============================================

mkdir("violin_plots");

% manual KDE function
function [xgrid, pdf_vals] = kde_manual(x)
    x = x(:);
    n = length(x);
    sigma = std(x) * (4/(3*n))^(1/5);  
    xgrid = linspace(min(x), max(x), 200);
    pdf_vals = zeros(size(xgrid));
    for j = 1:length(xgrid)
        pdf_vals(j) = mean(normpdf(xgrid(j), x, sigma));
    end
end

% outcome column index
label_idx = find(strcmp(cols, "outcome"));

% Loop over all features EXCEPT "outcome"
for f = 1:length(cols)
    
    if strcmp(cols{f}, "outcome")
        continue;
    end

    feat_name = cols{f};
    x = raw(:, f);
    y = raw(:, label_idx);

    % Split by class
    g0 = x(y == 0);   % died
    g1 = x(y == 1);   % euthanized
    g2 = x(y == 2);   % lived

    % KDE for each group
    [x0, p0] = kde_manual(g0);
    [x1, p1] = kde_manual(g1);
    [x2, p2] = kde_manual(g2);

    % Normalize width so shapes are comparable
    p0 = p0 / max(p0);
    p1 = p1 / max(p1);
    p2 = p2 / max(p2);

    % Create violin-like mirroring
    figure("visible", "off"); hold on;

    fill([p0 -p0(end:-1:1)], [x0 x0(end:-1:1)], [1 0.5 0.5], "facealpha", 0.7); % red-ish
    fill([p1 -p1(end:-1:1)] + 3, [x1 x1(end:-1:1)], [1 0.8 0.3], "facealpha", 0.7); % yellow-ish
    fill([p2 -p2(end:-1:1)] + 6, [x2 x2(end:-1:1)], [0.4 0.8 1], "facealpha", 0.7); % blue-ish

    title(sprintf("Violin Style Distribution – %s", feat_name));
    ylabel(feat_name);
    set(gca, "XTick", [0 3 6], "XTickLabel", {"Died", "Euthanized", "Lived"});

    hold off;

    % Save file
    fname = strcat("violin_plots/", feat_name, "_violin.png");
    print(fname, "-dpng");

end

disp("✅ Violin-style plots saved in folder: violin_plots/");
