% This script file demonstrates the fitting of a AIS 
% model on the monthly CO2 dataset.

% Perform clearing operations
clc
clear all
% Load montly CO2 data from corresponding .mat file.
% This operation results in loading variable co2pppm.
load ('MaunaLoaMonthlyCO2.mat');

% Remove incoplete data
co2ppm = co2ppm([2:end-1],:);
% Set up the AIS Short model forecasting
[forecast, for_ant,iter]= AISShortForecasting(co2ppm, 0.0020,2/3);
real_values = co2ppm(size(co2ppm) - size(forecast):end-1,:)

error = real_values - forecast




