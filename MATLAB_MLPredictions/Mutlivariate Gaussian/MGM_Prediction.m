%% Load Data (High run time)
filename = 'ProcessedData_AMLBD_small.csv'
%Import manually,IMPORT TYPE: numerical matrix
importdata(filename);
%% Normalize Data
% file should be a 2D numerical array, row by case, column by feature
data = FeatureNormalize(ProcessedAMLBDsmall);
%% Data Prep
% Select 2 Random Genes for Comparison
col1 = randi(size(data,2)); col2 = randi(size(data,2));
plot(data(:,col1), data(:, col2), 'bx');
%axis([-3 5 -3 5]);   % use for loglog plots
xlabel(['Gene No. :',num2str(col1)]);
ylabel(['Gene No. :', num2str(col2)]);
title('Exploratory: expression of 2 genes')
%% Visualize MGM
X = [data(:,col1), data(:,col2)];
[mu, sigma2] = estimateGaussian(X);

%  Returns the density of the multivariate normal at each data point (row) of X
p = multivariateGaussian(X, mu, sigma2);

%  Visualize the fit
visualizeFit(X,  mu, sigma2);
xlabel(['Gene No. :', num2str(col1)]);5
ylabel(['Gene No. :', num2str(col2)]);
%%
Xval = [data(:,col1) data(:,col2)]; 
yval = ones(size(Xval,1),1);
%% Optimize Threshold
pval = multivariateGaussian(Xval, mu, sigma2);
[epsilon, F1] = selectThreshold(yval, pval);
fprintf('Best epsilon found using cross-validation: %e\n', epsilon);
fprintf('Best F1 on Cross Validation Set:  %f\n', F1);

%  Find the outliers in the training set and plot the
outliers = find(p < epsilon);

%  Visualize the fit
visualizeFit(X,  mu, sigma2);
xlabel('Latency (ms)');
ylabel('Throughput (mb/s)');
%  Draw a red circle around those outliers
hold on
plot(X(outliers, 1), X(outliers, 2), 'ro', 'LineWidth', 2, 'MarkerSize', 10);
hold off
%% ================== Multidimensional Outliers ===================
%% High Dimension Dataset
%  Apply the same steps to the larger dataset
%% Select sets
X = data;    % All data
Xval = data(round(size(X,1)/3):end,:);    %subset of data(cross-validation)
yval = ones(size(Xval,1),1)
[mu, sigma2] = estimateGaussian(X);
%%
%  Training set 
p = multivariateGaussian(X, mu, sigma2);

%  Cross-validation set
pval = multivariateGaussian(Xval, mu, sigma2);

%  Find the best threshold
[epsilon, F1] = selectThreshold(yval, pval);
fprintf('Best epsilon found using cross-validation: %e\n', epsilon);

fprintf('Best F1 on Cross Validation Set:  %f\n', F1);
fprintf('# Outliers found: %d\n', sum(p < epsilon));