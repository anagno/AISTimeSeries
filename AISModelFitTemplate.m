% This script file demonstrates the fitting of a AIS 
% model on the monthly CO2 dataset.

% Perform clearing operations
clc
clear all
% Load montly CO2 data from corresponding .mat file.
% This operation results in loading variable co2pppm.
load ('MaunaLoaMonthlyCO2.mat');
% Transform matrix data into a corresponding row vector.
co2ppm = reshape(co2ppm',1,numel(co2ppm));
% Remove zero values.
co2ppm = co2ppm(co2ppm~=0);
% Set the time interval for each observation to be exactly one month.
dt = 1/12;
% Set the corresponding time limits;
MinYear = 1958;
MaxYear = 2014;
% Set the corresponding time interval.
t_min = MinYear + 3*dt;
t_max = MaxYear + 2*dt;
t = [t_min:dt:t_max];
% Scale time.
t_scaled = (t-mean(t))/std(t);
% Set the proportion of training data.
training_percentage = 0.90;
% Get the corresponding training and testing data points.
[x_train,x_test,y_train,y_test] = SplitTrainingTestingData(t_scaled,co2ppm,training_percentage);






% Set up the linear fit model parameters.
model_name = 'Linear Model ';
coefs = {'A','B'};
parameters_number = length(coefs);
% Set up the actual functional form of the fitting model.
functional_form = 'A*x+B';
ft = fittype(functional_form,'coef',coefs);
% Set up the fitting options.
opts = fitoptions(ft);
opts.Algorithm = 'Levenberg-Marquard';
opts.StartPoint = rand(1,parameters_number);
% Fit model to training data.
[fitresult,gof] = fit(x_train,y_train,ft,opts);
% Report goodness of fit measures.
fprintf('Goodness of fit measures on training data:\n');
fprintf('------------------------------------------\n');
fprintf('RMSE(train): %d\n',gof.rmse);
fprintf('RSQUARE(train): %d\n',gof.rsquare);

% Evaluate the linear model fit on the test data.
y_test_est = feval(fitresult,x_test);
% Calculate goodness of fit measures on the testing data.
rmse_test = sqrt(mean((y_test - y_test_est).^2));
y_test_mean = mean(y_test);
SStot = sum((y_test-y_test_mean).^2);
SSres = sum((y_test-y_test_est).^2);
rsquare_test = 1 - (SSres/SStot);
% Report goodness of fit measures on testing data.
fprintf('Goodness of fit measures on testing data:\n');
fprintf('------------------------------------------\n');
fprintf('RMSE(test): %d\n',rmse_test);
fprintf('RSQUARE(test): %d\n',rsquare_test);

% Visualize the fitted model.
PlotFitModel(t_scaled,co2ppm,training_percentage,model_name,fitresult);
% Visualize the fitted model residuals.
PlotFitModelResiduals(t_scaled,co2ppm,training_percentage,model_name,fitresult);




