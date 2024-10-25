% Multivariate Signal Feature Extraction and Analysis

%% Initialization and Path Setup
clc; clear; close all;
basePath = 'D:\arman\';
dataPath = fullfile(basePath, 'arman3', '8');
outputPath = fullfile(basePath, 'output');
mkdir(outputPath); % Create output directory if it doesn't exist

%% Load Data from Excel Files
cd(dataPath);
fileList = dir(fullfile(dataPath, '*.xlsx'));
numFiles = numel(fileList);
dataColumns = cell(1, numFiles);

for i = 1:numFiles
    M = xlsread(fileList(i).name);
    dataColumns{i} = M(:, 2); % Extract the second column of data
end

%% Select First Data Column and Segment for Analysis
M2_eye = dataColumns{1}; % Example column data for eye signal
numSegments = 10; % Define number of segments for analysis
segmentLength = 1740; % Segment length

% Initialize feature matrix
features_eye = zeros(numSegments, 11);

% Calculate features for each segment
for i = 1:numSegments
    segmentData = M2_eye((i-1)*segmentLength + 1:segmentLength*i);
    features_eye(i, :) = ex_features(segmentData);
end

%% Visualization of Features
close all;
for i = 1:11
    figure;
    plot(features_eye(:, i));
    title(['Feature ', num2str(i), ' Across Segments']);
    
    % Save the figure as a JPG file
    saveas(gcf, fullfile(outputPath, ['Feature_', num2str(i), '.jpg']));
    close(gcf); % Close the figure to prevent overlap
end

%% Extract Features for Two Sets of Data
numPairs = floor(numFiles / 2);
featuresStart = zeros(11, numPairs);
featuresEnd = zeros(11, numPairs);

for i = 1:numPairs
    featuresStart(:, i) = ex_features(dataColumns{2 * i - 1});
    featuresEnd(:, i) = ex_features(dataColumns{2 * i});
end

%% Visualize Data from First Two Columns (Scatter Plot)
figure;
hold on;
for idx = 1:2
    data = dataColumns{idx};
    plot(data(1:end-1), data(2:end), 'DisplayName', ['Data ', num2str(idx)]);
end
title('Scatter Plot of First Two Data Columns');
legend;
hold off;

%% Visualize Sample Data for Comparison
figure;
plot(dataColumns{3}, 'DisplayName', 'Data Column 3');
hold on;
plot(dataColumns{4}(20:end), 'DisplayName', 'Data Column 4');
title('Comparison of Data Columns 3 and 4');
legend;
hold off;

%% Classification Using Feature Data
FF = [featuresStart; featuresEnd];
L = [zeros(11, 1); ones(11, 1)];
classificationLearner(FF, L);

%% Feature Comparison Using Bar Graphs
for i = 1:11
    figure;
    bar([featuresStart(i, :); featuresEnd(i, :)]);
    title(['Comparison of Feature ', num2str(i)]);
    legend('Start', 'End');
end

%% Mean and Variance Analysis for Different Feature Sets
% Assume hip_ans1 through hip_ans9 are stored as variables
AAA(:,:,1) = [features_eye(1,:)];
AAA(:,:,2) = [features_eye(4,:)];
AAA(:,:,3) = [features_eye(5,:)];
AAA(:,:,4) = [features_eye(6,:)];
AAA(:,:,5) = [features_eye(8,:)];
AAA(:,:,6) = [features_eye(9,:)];

% Calculate means and variances along the third dimension
meanValues = mean(AAA, 3);
varianceValues = var(AAA, 0, 3); % Variance calculation

% X-axis values (assuming a simple range)
x = 1:size(AAA, 1);

% Plot mean and variance lines
figure;
hold on;
for i = 1:size(meanValues, 2)
    plot(meanValues(:, i), '--', 'LineWidth', 2, 'DisplayName', ['Mean ', num2str(i)]);
end

% Add variance lines with scaling factor
scaleFactor = 0.5;
for i = 1:size(varianceValues, 2)
    plot(x, meanValues(:, i) + sqrt(varianceValues(:, i)) * scaleFactor, '-');
    plot(x, meanValues(:, i) - sqrt(varianceValues(:, i)) * scaleFactor, '-');
end

xlabel('X-axis');
ylabel('Y-axis');
title('Mean and Variance Analysis');
legend;
grid on;
hold off;

%% Feature Extraction Function (ex_features)
function [features] = ex_features(data)
    % Calculate PSD and statistical features for each column of data
    length_data_col = size(data, 2);
    features = zeros(11, length_data_col);

    for i = 1:length_data_col
        % Define shifted segments for each column of data
        an_data_col = data(:, i);
        an_1_data_col = an_data_col(2:end);
        an_2_data_col = an_1_data_col(2:end);
        an_3_data_col = an_2_data_col(2:end);
        an_4_data_col = an_3_data_col(2:end);
        
        len_an_1 = length(an_1_data_col);
        len_an_2 = length(an_2_data_col);
        len_an_3 = length(an_3_data_col);
        len_an_4 = length(an_4_data_col);

        % Calculate features SCCA, SSVL, TACR, and others
        scca = sum(pi * ((an_1_data_col - an_data_col(1:end-1)).^2 + (an_2_data_col - an_1_data_col).^2));
        ssvl = sum(sqrt((an_1_data_col - an_data_col(1:end-1)).^2 + (an_2_data_col - an_1_data_col).^2));
        tacr = sum(abs(det([0 an_data_col(1:end-2)'; an_1_data_col'; an_2_data_col'])));
        
        % SDHC Calculation
        sdhc = sum(sqrt(((an_1_data_col + an_2_data_col + an_3_data_col) / 3 - ...
                         (an_data_col(1:end-3) + an_1_data_col + an_2_data_col) / 3).^2 + ...
                         ((an_2_data_col + an_3_data_col + an_4_data_col) / 3 - ...
                         (an_1_data_col + an_2_data_col + an_3_data_col) / 3).^2));

        % SH45 and SH135 Calculation
        sh45 = sum(abs(an_1_data_col - an_data_col(1:end-1)) / sqrt(2));
        sh135 = sum(abs(an_1_data_col + an_data_col(1:end-1)) / sqrt(2));

        d_det_SHCA = 0;
        d_det_SCTA_s = 0;
        for j = 1:len_an_3
            d_det_SCTA = det([an_data_col(j), an_1_data_col(j), an_2_data_col(j); ...
                              an_1_data_col(j), an_2_data_col(j), an_3_data_col(j); ...
                              1, 1, 1]);
            d_det_SCTA_s = 0.5 * (d_det_SCTA_s + d_det_SCTA);
            sqss = sqrt((an_1_data_col(j) - an_data_col(j)).^2 + (an_2_data_col(j) - an_1_data_col(j)).^2) + ...
                   sqrt((an_2_data_col(j) - an_data_col(j)).^2 + (an_3_data_col(j) - an_1_data_col(j)).^2) + ...
                   sqrt((an_2_data_col(j) - an_1_data_col(j)).^2 + (an_3_data_col(j) - an_2_data_col(j)).^2);
            d_det_SHCA = pi * (d_det_SCTA / sqss).^2 + d_det_SHCA;
        end

        % Store calculated features
        features(1, i) = scca;
        features(2, i) = ssvl;
        features(3, i) = 0.5 * tacr;
        features(4, i) = sdhc;
        features(5, i) = sh45;
        features(6, i) = sh135;
        features(7, i) = d_det_SCTA_s;
        features(8, i) = d_det_SHCA;

        % Additional Features: SDTC, SABP, SCRA
        SDTC = sum(sqrt(an_1_data_col.^2 + an_data_col(1:end-1).^2));
        SCRA = 0.5 * sum(abs(an_1_data_col + an_data_col(1:end-1)) .* abs(an_1_data_col - an_data_col(1:end-1)));

        SABP = 0;
        for j2 = 1:len_an_3
            sss1 = sqrt((an_1_data_col(j2) - an_data_col(j2)).^2 + (an_2_data_col(j2) - an_1_data_col(j2)).^2);
            sss2 = sqrt((an_2_data_col(j2) - an_1_data_col(j2)).^2 + (an_3_data_col(j2) - an_2_data_col(j2)).^2);
            SABP = SABP + ((an_1_data_col(j2) - an_data_col(j2)) * (an_2_data_col(j2) - an_1_data_col(j2)) + ...
                          (an_2_data_col(j2) - an_1_data_col(j2)) * (an_3_data_col(j2) - an_2_data_col(j2))) / (sss1 + sss2);
        end

        features(9, i) = SDTC;
        features(10, i) = SABP;
        features(11, i) = SCRA;
    end
end