% This script file demonstrates the fitting of a constructive non-linear 
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
%--------------------------------------------------------------------------
% 1st REGRESSION COMPONENT
%--------------------------------------------------------------------------
% Set up the single sinusoidal fit model parameters.
model_name1 = 'Single Sinusoidal Model ';
coefs1 = {'A1','P1','T'};
parameters_number1 = length(coefs1);
% Set up the actual functional form of the fitting model.
functional_form1 = 'A1*sin((2*pi/T)*(x+P1))';
ft1 = fittype(functional_form1,'coef',coefs1);
% Set up the fitting options.
opts1 = fitoptions(ft1);
opts1.Algorithm = 'Levenberg-Marquard';
opts1.StartPoint = rand(1,parameters_number1);
% Fit model to training data.
[fitresult1,gof1] = fit(x_train,y_train,ft1,opts1);
% Report goodness of fit measures.
fprintf('RMSE on training data: %d\n',gof1.rmse);
fprintf('------------------------------------------\n');

% Evaluate the single sinusoidal model fit on the test data.
y_test_est1 = feval(fitresult1,x_test);
% Calculate goodness of fit measures on the testing data.
rmse_test1 = sqrt(mean((y_test - y_test_est1).^2));
% Report goodness of fit measures on testing data.
fprintf('RMSE on testing data: %d\n',rmse_test1);
fprintf('------------------------------------------\n');

% Visualize the fitted model.
PlotFitModel(t_scaled,co2ppm,training_percentage,model_name1,fitresult1);
% Visualize the fitted model residuals.
PlotFitModelResiduals(t_scaled,co2ppm,training_percentage,model_name1,fitresult1);

%--------------------------------------------------------------------------
% 2nd REGRESSION COMPONENT
%--------------------------------------------------------------------------
% Set up the double sinusoidal fit model parameters.
model_name2 = 'Double Sinusoidal Model ';
coefs2 = {'A1','P1','A2','P2','T'};
parameters_number2 = length(coefs2);
% Set up the actual functional form of the fitting model.
functional_form2 = 'A1*sin((2*pi/T)*(x+P1))+A2*sin((2*pi*2/T)*(x+P2))';
ft2 = fittype(functional_form2,'coef',coefs2);
% Set up the fitting options.
opts2 = fitoptions(ft2);
opts2.Algorithm = 'Levenberg-Marquard';
opts2.StartPoint = rand(1,parameters_number2);
% Fit model to training data.
[fitresult2,gof2] = fit(x_train,y_train,ft2,opts2);
% Report goodness of fit measures.
fprintf('RMSE on training data: %d\n',gof2.rmse);
fprintf('------------------------------------------\n');

% Evaluate the double sinusoidal model fit on the test data.
y_test_est2 = feval(fitresult2,x_test);
% Calculate goodness of fit measures on the testing data.
rmse_test2 = sqrt(mean((y_test - y_test_est2).^2));
% Report goodness of fit measures on testing data.
fprintf('RMSE on testing data: %d\n',rmse_test2);
fprintf('------------------------------------------\n');

% Visualize the fitted model.
PlotFitModel(t_scaled,co2ppm,training_percentage,model_name2,fitresult2);
% Visualize the fitted model residuals.
PlotFitModelResiduals(t_scaled,co2ppm,training_percentage,model_name2,fitresult2);

%--------------------------------------------------------------------------
% 3rd REGRESSION COMPONENT
%--------------------------------------------------------------------------
% Set up the triple sinusoidal fit model parameters.
model_name3 = 'Triple Sinusoidal Model ';
coefs3 = {'A1','P1','A2','P2','A3','P3','T'};
parameters_number3 = length(coefs3);
% Set up the actual functional form of the fitting model.
functional_form3 = 'A1*sin((2*pi/T)*(x+P1))+A2*sin((2*pi*2/T)*(x+P2))+A3*sin((2*pi*3/T)*(x+P3))';
ft3 = fittype(functional_form3,'coef',coefs3);
% Set up the fitting options.
opts3 = fitoptions(ft3);
opts3.Algorithm = 'Levenberg-Marquard';
opts3.StartPoint = rand(1,parameters_number3);
% Fit model to training data.
[fitresult3,gof3] = fit(x_train,y_train,ft3,opts3);
% Report goodness of fit measures.
fprintf('RMSE on training data: %d\n',gof3.rmse);
fprintf('------------------------------------------\n');

% Evaluate the triple sinusoidal model fit on the test data.
y_test_est3 = feval(fitresult3,x_test);
% Calculate goodness of fit measures on the testing data.
rmse_test3 = sqrt(mean((y_test - y_test_est3).^2));
% Report goodness of fit measures on testing data.
fprintf('RMSE on testing data: %d\n',rmse_test3);
fprintf('------------------------------------------\n');

% Visualize the fitted model.
PlotFitModel(t_scaled,co2ppm,training_percentage,model_name3,fitresult3);
% Visualize the fitted model residuals.
PlotFitModelResiduals(t_scaled,co2ppm,training_percentage,model_name3,fitresult3);

%--------------------------------------------------------------------------
% 4th REGRESSION COMPONENT
%--------------------------------------------------------------------------
% Set up the triple sinusoidal exponential fit model parameters.
model_name4 = 'Triple Sinusoidal Exponential Model ';
coefs4 = {'A1','P1','A2','P2','A3','P3','T','A','B','C'};
parameters_number4 = length(coefs4);
% Set up the actual functional form of the fitting model.
functional_form4 = 'A1*sin((2*pi/T)*(x+P1))+A2*sin((2*pi*2/T)*(x+P2))+A3*sin((2*pi*3/T)*(x+P3))+A*exp(-B*x)+C';
ft4 = fittype(functional_form4,'coef',coefs4);
% Set up the fitting options.
opts4 = fitoptions(ft4);
opts4.Algorithm = 'Levenberg-Marquard';
%opts4.Lower = -Inf * ones(1,parameters_number4);
%opts.Upper = Inf * ones(1,parameters_number4);
opts4.StartPoint = rand(1,parameters_number4);
% Fit model to training data.
[fitresult4,gof4] = fit(x_train,y_train,ft4,opts4);
% Report goodness of fit measures.
fprintf('RMSE on training data: %d\n',gof4.rmse);
fprintf('RSQUARE(train): %d\n',gof4.rsquare);
fprintf('------------------------------------------\n');

% Evaluate the triple sinusoidal exponential model fit on the test data.
y_test_est4 = feval(fitresult4,x_test);
% Calculate goodness of fit measures on the testing data.
rmse_test4 = sqrt(mean((y_test - y_test_est4).^2));
y_test_mean = mean(y_test);
SStot = sum((y_test-y_test_mean).^2);
SSres = sum((y_test-y_test_est4).^2);
rsquare_test = 1 - (SSres/SStot);
% Report goodness of fit measures on testing data.
fprintf('RMSE on testing data: %d\n',rmse_test4);
fprintf('RSQUARE(test): %d\n',rsquare_test);
fprintf('------------------------------------------\n');

% Visualize the fitted model.
PlotFitModel(t_scaled,co2ppm,training_percentage,model_name4,fitresult4);
% Visualize the fitted model residuals.
PlotFitModelResiduals(t_scaled,co2ppm,training_percentage,model_name4,fitresult4);