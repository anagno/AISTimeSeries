% This script file demonstrates the fitting of a AIS 
% model on the monthly CO2 dataset.

% Perform clearing operations
%clc
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

training_percentage = 0.9;
cutoff_index = round(size(co2ppm,1) * training_percentage);

training_data = co2ppm(1:cutoff_index-1,:);
testing_data = co2ppm(cutoff_index:end,:);
%parameters = [] ;

%for threshold = 0.001:0.002:0.011
%    for relax_threshold = 0.01:0.5:0.15
%        for training_percentage = 0.65:0.15:0.95
%            for max_iterations = 40:50:140
%                for beta = 0.04:0.15:1

% !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!%
% Set the parameters.
threshold = 0.003;
relax_threshold = 0.1;
max_iterations = 50 ;
beta = 0.05 ;

% Set up the AIS Short model forecasting
[forecast, confidence, antibodies,iter, total_time, antigens, ...
            train_rmse, train_errors]= ...
            AISShortForecasting(training_data, threshold, relax_threshold, ...
            max_iterations, beta,testing_data,[],true);
                            
% !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!%                            

% Visualize the fitted model.
figure_name = 'Fit Visualization';
legend_name = 'Fit';
figure('Name',figure_name)%,'Visible','off');
hold on

train_data_reshape = reshape(training_data',numel(training_data),1);
forecast_reshape = reshape(forecast',numel(forecast),1);
co2ppm_reshape = reshape(co2ppm',numel(co2ppm),1);

plot(t',co2ppm_reshape,'.b');
plot(t(1:size(train_data_reshape))',train_data_reshape,'.g');
plot(t(size(train_data_reshape)+1:end)',forecast_reshape,'.r');
title(strcat('Parameters: threshold=',num2str(threshold) , ...
             ', relax\_threshold=', num2str(relax_threshold), ...
             ',\newlinetraining\_percentage=', num2str(training_percentage), ...
             ', max\_iterations=', num2str(max_iterations), ...
             ', beta=', num2str(beta), '\newline \newline', ...
             ' total\_time=', num2str(total_time)));
legend('Original Data','Training Data','Testing Data','location','northwest');
% Label axes
xlabel('Time');
ylabel('co2ppm');
grid on
hold off
%print(strcat('doc/data/vis_threshold_',num2str(threshold), ...
%             '_relax_threshold_',num2str(relax_threshold), ...
%             '_training_percentage_', num2str(training_percentage), ...
%             '_max_iterations_', num2str(max_iterations), ...
%             '_beta_', num2str(beta),...
%             '.eps'),'-depsc');
%close;

% Visualize the model residuals.
real_values = co2ppm(cutoff_index:end,:);
error = forecast - real_values;
error_reshape = reshape(error',numel(error),1);
train_error_reshape = reshape(train_errors',numel(train_errors),1);


% Calculate goodness of fit measures on the testing data.
rmse = mean(sqrt(mean((error.^2))));
%fprintf('Goodness of fit measures on training data:\n');
%fprintf('------------------------------------------\n');
%fprintf('RMSE: %d\n',mean(train_rmse));

%fprintf('Goodness of fit measures on testing data:\n');
%fprintf('------------------------------------------\n');
%fprintf('RMSE: %d\n',rmse);


% Perform the actual plotting operations.
figure_name = 'Fit Residuals Visualization';
figure('Name',figure_name)%,'Visible','off');
hold on
plot(t(1:size(train_error_reshape))',train_error_reshape,'.g');
plot(t(size(train_data_reshape)+1:end)',error_reshape,'.r');
title(strcat('Parameters: threshold=',num2str(threshold) , ...
            ', relax\_threshold=', num2str(relax_threshold), ...
             ',\newlinetraining\_percentage=', num2str(training_percentage), ...
             ', max\_iterations=', num2str(max_iterations), ...
             ', beta=', num2str(beta), '\newline \newline', ...
             ' train\_rmse=', num2str(mean(train_rmse)), ...
             ', rmse=', num2str(rmse)));
legend('Training Residuals', 'Testing Residuals', 'location','northwest');
% Label axes
xlabel('Time');
ylabel('Residuals');
grid on
hold off
%print(strcat('doc/data/res_threshold_',num2str(threshold), ...
%             '_relax_threshold_',num2str(relax_threshold), ...
%             '_training_percentage_', num2str(training_percentage), ...
%             '_max_iterations_', num2str(max_iterations), ...
%             '_beta_', num2str(beta),...
%             '.eps'),'-depsc');
%close;




%new_row = [threshold relax_threshold training_percentage...
%     max_iterations beta iter total_time mean(rmse) mean(train_rmse) size(antigens,1)]

%parameters = vertcat(parameters,new_row);


%                end
%            end
%        end
%    end
%end              
                
%parameters
