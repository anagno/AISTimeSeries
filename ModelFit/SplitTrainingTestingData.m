function [x_train,x_test,y_train,y_test] = SplitTrainingTestingData(x,y,training_percentage)
% This function splits a given data set into training and testing samples
% given the percentage of training samples that should be used.

% Get the number of observations.
N = max(length(x),length(y));

% Set the corresponding cuttoff index.
cutoff_index = round(N*training_percentage);

% Set the training and testing data for the x-values of the linear model.
if(~isempty(x))
    x_train = x([1:cutoff_index])';
    x_test = x([cutoff_index+1:N])';
else
    x_train = [];
    x_test = [];
end;
% Set the training and testing data for the y-values of the linear model.
if(~isempty(y))
    y_train = y([1:cutoff_index])';
    y_test = y([cutoff_index+1:N])';
else
    y_train = [];
    y_test = [];
end;
end