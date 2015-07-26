function [ forecast, confidence, antibodies, iterations, total_time, ...
    rmse, errors ] = AISLongForecasting( original_data, num_periods, ...
    threshold, relax_threshold, max_iterations, beta, input_antibodies, ...
    diagnostics )
%AISLongForecasting A function that use the AISShortForecasting function to
%make longer forecasts. It uses the produce forecast as input in the AIS
%to continue the forecast precedure.
%   For more information look on the AISShortForecasting function

% INPUT VARIABLES:
% original_data: the array of input data. It should be an array with the
%   lines represinting the a hole period of measumerents.
% num_periods: the number of future periods we want to forecast
% threshold: declares the the cross-reactivity threshold r, for which an
%   antigen is activated
% relax_threshold: if an antigen does not react the threshold is relax
%   according to this particular percentage.
% max_iterations: The maximum number of iterations 
% beta: The shape parameter
% input_antibodies: Antibodies that can be used in the AIS
% forecast_data: The data for which a forecast from the AIS will be done.
%   If it is left the forecast will go on only for one period
% diagnostics: If set to true diagnostics messages will be printed

% OUTPUT VARIABLES:
% forecast: the forecast values 
% confidence: the confidence of the forecast values. If there is not an
%   antibody that reacts to the values then treshold is relaxed and that is
%   shown in the confidence values that is between 0 and 1
% antibodies: the produced antibodies from the AIS
% iterations: the iterations that the AIS needed
% total_time: the total running time of the algorithm
% forecast_antigen: the forecast antigens that were used to produce the
%   forecast
% rmse: the rmse in the original data of the AIS
% errors: the errors in the training data from which the rmse is produced

switch nargin
    case 8
        
    case 7
        diagnostics = false;
    case 6
        input_antibodies = [];
        diagnostics = false;
    case 5
        beta = 0.04;
        input_antibodies = [];
        diagnostics = false;
    case 4
        max_iterations = 50;
        beta = 0.04;
        input_antibodies = [];
        diagnostics = false;
    case 3
        max_iterations = 50;
        beta = 0.04;
        input_antibodies = [];
        diagnostics = false;
    case 2
        relax_threshold = 0.01; 
        max_iterations = 50;
        beta = 0.04;
        input_antibodies = [];
        diagnostics = false;
    otherwise
        error ('Too few or too many arguments were entered');
end

assert(num_periods >= 1,'The number of periods must be at least 1');

period_size = size(original_data,2);

forecast = zeros(num_periods,period_size);
confidence = zeros(num_periods,1);

[forecast(1,:), confidence(1,:), antibodies, iterations_local, total_time_loc,...
    ~, rmse, errors]= AISShortForecasting(original_data, ...
    threshold, relax_threshold, max_iterations, beta, ...
    [],input_antibodies,true, diagnostics);

iterations = iterations_local;
total_time = total_time_loc;
current_iteration = 2;

while(num_periods +1 > current_iteration)
    
    if(diagnostics)
        fprintf('%d remaining. \n',num_periods-current_iteration)
    end
    
    [forecast(current_iteration,:),  confidence(current_iteration,:), ...
        antibodies, iterations_local, total_time_loc, ~,rmse, errors]= ...
        AISShortForecasting(original_data, threshold, relax_threshold, ...
        max_iterations, beta, forecast(current_iteration-1,:), ...
        antibodies, false, diagnostics);
        
    iterations = iterations + iterations_local;
    total_time = total_time + total_time_loc;
    current_iteration = current_iteration + 1;
    
end


end

