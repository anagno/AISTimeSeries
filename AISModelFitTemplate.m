% This script file demonstrates the fitting of a AIS 
% model on the monthly CO2 dataset.

% Perform clearing operations
clc
clear all
% Load montly CO2 data from corresponding .mat file.
% This operation results in loading variable co2pppm.
load ('MaunaLoaMonthlyCO2.mat');

% Remove incoplete data
co2ppm = co2ppm(2:end-1,:);

MinYear = 1959;
MaxYear = MinYear + size(co2ppm,1);
dt = 1/size(co2ppm,2);
t_min = MinYear;
t_max = MaxYear;
t = [t_min:dt:t_max];
% Due to the fact that the 2014 is not included.
t = t(1:end-1);

training_percentage = 0.8;
cutoff_index = round(size(co2ppm,1) * training_percentage);

training_data = co2ppm(1:cutoff_index-1,:);
testing_data = co2ppm(cutoff_index-1:end-1,:);

% !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!%
% Set the parameters.
threshold = 0.00355;
relax_threshold = 0.1;
max_iterations = 5000 ;
beta = 0.04 ;

% Set up the AIS Long model forecasting
[forecast_long, confidence_long, antibodies_long,iter_long, ... 
            total_time_long, train_rmse_long, train_errors_long]= ...
            AISLongForecasting(training_data, size(testing_data,1), ...
            threshold, relax_threshold, max_iterations, beta,[],true);
       
% Set up the AIS Short model forecasting
[forecast, confidence, antibodies,iter, total_time, antigens, ...
            train_rmse, train_errors]= ...
            AISShortForecasting(training_data, threshold, relax_threshold, ...
            max_iterations, beta,testing_data,antibodies_long,false,true);
                            
% !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!%                            

% To show the differences
training_data(2:end,:) = training_data(2:end,:) + train_errors;
training_data_long = training_data;
training_data_long(2:end,:) = training_data(2:end,:) + train_errors_long;

train_data_reshape = reshape(training_data',numel(training_data),1);
train_data_long_reshape = reshape(training_data_long', ...
                            numel(training_data_long),1);
forecast_reshape = reshape(forecast',numel(forecast),1);
forecast_long_reshape = reshape(forecast_long',numel(forecast),1);
co2ppm_reshape = reshape(co2ppm',numel(co2ppm),1);

% Visualize the fitted model.
figure('Name','Fit Visualization');
hold on
plot(t(1:size(co2ppm_reshape))',co2ppm_reshape,'.-b');
plot(t(1:size(train_data_reshape))',train_data_reshape,'*-g');
plot(t(1:size(train_data_long_reshape))',train_data_long_reshape, ...
                    'color',[1 1 0],'marker','o');

plot(t(end-size(forecast_reshape)+1:end)',forecast_reshape,'*-r');
plot(t(end-size(forecast_long_reshape)+1:end)',forecast_long_reshape, ...
                    'color',[1 0.5 0],'marker','o');
title(strcat('Parameters: threshold=',num2str(threshold) , ...
             ', relax\_threshold=', num2str(relax_threshold), ...
             ',\newlinetraining\_percentage=', num2str(training_percentage), ...
             ', max\_iterations=', num2str(max_iterations), ...
             ', beta=', num2str(beta), '\newline \newline', ...
             ' total\_time\_short=', num2str(total_time),...
             ',total\_time\_long=', num2str(total_time_long)));
legend('Original Data','Short Training Data','Long Training Data', ...
    'Shrot Forecast Data','Long Forecast Data','location','northwest');
% Label axes
xlabel('Time');
ylabel('co2ppm');
grid on
hold off

% Visualize the model residuals.
real_values = co2ppm(cutoff_index:end,:);
error = forecast - real_values;
error_reshape = reshape(error',numel(error),1);
train_error_reshape = reshape(train_errors',numel(train_errors),1);

error_long = forecast_long - real_values;
error_long_reshape = reshape(error_long',numel(error_long),1);
train_error_long_reshape = reshape(train_errors_long', ...
                    numel(train_errors_long),1);

% Calculate goodness of fit measures on the testing data.
rmse = mean(sqrt(mean((error.^2))));
rmse_long = mean(sqrt(mean((error_long.^2))));

fprintf('Goodness of fit measures on short training data:\n');
fprintf('------------------------------------------\n');
fprintf('RMSE: %d\n',mean(train_rmse));
fprintf('Goodness of fit measures on short testing data:\n');
fprintf('------------------------------------------\n');
fprintf('RMSE: %d\n',rmse);

fprintf('Goodness of fit measures on long training data:\n');
fprintf('------------------------------------------\n');
fprintf('RMSE: %d\n',mean(train_rmse_long));
fprintf('Goodness of fit measures on long testing data:\n');
fprintf('------------------------------------------\n');
fprintf('RMSE: %d\n',rmse_long);

% Perform the actual plotting operations.
figure_name = 'Fit Residuals Visualization';
figure('Name',figure_name)%,'Visible','off');
hold on
plot(t(13:size(train_error_reshape)+12)',train_error_reshape,'.-g');
plot(t(end-size(error_reshape)+1:end)',error_reshape,'.-r');

plot(t(13:size(train_error_long_reshape)+12)',train_error_long_reshape, ...
                'color',[1 1 0],'marker','o');
plot(t(end-size(error_long_reshape)+1:end)',error_long_reshape, ...
                'color',[1 0.5 0],'marker','o');

title(strcat('Parameters: threshold=',num2str(threshold) , ...
            ', relax\_threshold=', num2str(relax_threshold), ...
             ',\newlinetraining\_percentage=', num2str(training_percentage), ...
             ', max\_iterations=', num2str(max_iterations), ...
             ', beta=', num2str(beta), '\newline \newline', ...
             ' train\_rmse=', num2str(mean(train_rmse)), ...
             ', rmse=', num2str(rmse), '\newline', ...
             ' train\_long\_rmse=', num2str(mean(train_rmse_long)), ...
             ',rmse\_long=', num2str(rmse_long)));
         
legend('Short Training Residuals', 'Short Forecast Residuals', ...
               'Long Training Residuals', 'Long Forecast Residuals',...
               'location','northwest');
% Label axes
xlabel('Time');
ylabel('Residuals');
grid on
hold off