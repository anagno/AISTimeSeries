% This script file demonstrates the fitting of a AIS 
% model on the monthly CO2 dataset.

% Perform clearing operations
clc
clear all
% Load montly CO2 data from corresponding .mat file.
% This operation results in loading variable co2pppm.
load ('MaunaLoaMonthlyCO2.mat');

% Remove zero values.
%co2ppm = co2ppm(co2ppm==0);
% Set the time interval for each observation to be exactly one month.
%dt = 1/12;
% Set the corresponding time limits;
%MinYear = 1958;
%MaxYear = 2014;
% Set the corresponding time interval.
%t_min = MinYear + 3*dt;
%t_max = MaxYear + 2*dt;
%t = [t_min:dt:t_max];
% Scale time.
%t_scaled = (t-mean(t))/std(t);
% Set the proportion of training data.
%training_percentage = 0.90;
% Get the corresponding training and testing data points.
%[x_train,x_test,y_train,y_test] = SplitTrainingTestingData(t_scaled,co2ppm,training_percentage);

% Remove incoplete data
co2ppm = co2ppm([2:end-1],:);
% Set up the AIS Short model forecasting
forecast= AISShortForecasting(co2ppm, 0.0020,2/3)





