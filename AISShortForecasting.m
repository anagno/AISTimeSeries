function [forecast, confidence, antibodies, iterations, total_time,  ...
    forecast_antigen, rmse, errors] = AISShortForecasting( original_data, ...
    threshold, relax_threshold, max_iterations, ...
    beta,forecast_data, input_antibodies, hypermutation, diagnostics )
%AISShortForecasting This is a function for forecasting time series using an
%artificial immune system
%   The function is implemented using the algorithm presented in [1].

%   Bibliography
%   [1] Dudek, Grzegorz. "Artificial immune system for short-term electric 
%   load forecasting." Artificial Intelligence and Soft Computing–ICAISC 
%   2008. Springer Berlin Heidelberg, 2008. 1007-1017.
%
% INPUT VARIABLES:
% original_data: the array of input data. It should be an array with the
%   lines represinting the a hole period of measumerents.
% threshold: declares the the cross-reactivity threshold r, for which an
%   antigen is activated
% relax_threshold: if an antigen does not react the threshold is relax
%   according to this particular percentage.
% max_iterations: The maximum number of iterations 
% beta: The shape parameter
% input_antibodies: Antibodies that can be used in the AIS
% forecast_data: The data for which a forecast from the AIS will be done.
%   If it is left the forecast will go on only for one period
% hypermutation: A parameter that instructs the AIS not to produce new
%   antibodies. Should only be used if the antibodies are provided
% diagnostics: If set to true diagnostics messages will be printed
%
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
    case 9 
        
    case 8
        hypermutation = false;
    case 7
        diagnostics = false;
        hypermutation = false;
    case 6
        input_antibodies = [];
        diagnostics = false;
    case 5
        beta = 0.04;
        input_antibodies = [];
        forecast_data = [];
        diagnostics = false;
        hypermutation = false;
    case 4
        max_iterations = 50;
        beta = 0.04;
        input_antibodies = [];
        forecast_data = [];
        diagnostics = false;
        hypermutation = false;
    case 3
        max_iterations = 50;
        beta = 0.04;
        input_antibodies = [];
        forecast_data = [];
        diagnostics = false;
        hypermutation = false;
    case 2
        relax_threshold = 0.01; 
        max_iterations = 50;
        beta = 0.04;
        input_antibodies = [];
        forecast_data = [];
        diagnostics = false;
        hypermutation = false;
    otherwise
        error ('Too few or too many arguments were entered');
end

% !!! check if data is good with no null variables. !!!
assert(sum(sum(isnan(original_data)))==0,'NaN values in data');

% Start timer for perfomance
time_start = tic;

% Preprocess the data in order to get rid of the time series trend and 
% seasonality, and simplify the model.

period_size = size(original_data,2);
average_period_data = mean(original_data,2);
data = original_data ./ repmat(average_period_data,1,period_size);
data_next = original_data(2:end,:) ./ repmat(average_period_data(1:end-1), ...
                            1,period_size);

% Loading of the training set of antigens
% The whole dataset is divided into two subsets – training one and test 
% one. The first sequences of the time series (typically two thirds of the 
% whole time series) are included in the training set and the latest 
% sequences are included in the test set. Immune memory is trained using 
% the training set, and after learning the model is tested using the 
% test set.

if(size(forecast_data,1) > 0 )
   assert(size(forecast_data,2) == size(original_data,2), ...
      'Periods are of forecast data and training data is not the same');
   train_data = data;
   average_train = average_period_data;
   average_test = mean(forecast_data,2);
   forecast_data = forecast_data ./ repmat(average_test,1,period_size);
   test_data = forecast_data;
else
    train_data = data(1:end,:);
    average_train = average_period_data(1:end,:);
    test_data = data(end,:);
    average_test = average_period_data(end,:);
end

antigens(:,:) = horzcat(train_data(1:size(train_data)-1,:), ...
                                    data_next(:,:));

% Generation of the initial antibody population. An initial antibody 
% population is created by copying all the antigens from the training
% set (antibodies and antigens have the same structure). This way of
% initialization prevents inserting antibodies in empty regions without 
% antigens. Also the user can choose the antibodies which can be used.

if(size(input_antibodies,1) > 0 )
    antibodies = input_antibodies;
else
    antibodies = antigens;
end

% Calculation of the affinity of antibodies for antigens
% see function calculationOfAfinityTable

affinity_table = calculationOfAfinityTable(antigens,antibodies);

% Activated antibody detection and evaluation
% see function affinityCalculation

% If the affinity of the antibody for the antigen is smaller than or
% equal to the cross-reactivity r, it means that the antigen lies in the
% antibody recognition region (the antibody is activated by the antigen)
% If several antigens lie in the antibody recognition region, the error is 
% calculated for each of them. The mean error delta is applied to evaluate 
% the antibody and is minimized in the following iterations of 
% the algorithm.

[delta_table,enabled_antigens] = antibodyEvaluation(antigens,...
                                    antibodies,affinity_table,threshold);
                                
total_enable_antigens = sum(sum(enabled_antigens));

% Do until the stop criterion is reached
iterations = 1;
while(hypermutation == true && iterations <= max_iterations)
    
    % Store the population of antibodies to check if they change
    old_antibodies = antibodies;
    
    % Clonal selection
    % Each antibody cases secreting as many clones as many antigens are
    % in its recognition region. Thus most clones are generated in the
    % dense clusters of antigens.
    
    new_antibodies = zeros(total_enable_antigens, period_size*2);
    new_antibody_count = 1;
           
    for antibody = 1:size(antibodies,1)
        for antigen = 1:size(antigens,1)
            if (enabled_antigens(antigen,antibody))
                
                % Clone hypermutation
                % The main goal of hypermutation is to improve the
                % diversity of the immune system in order to effectively
                % recognize new antigens.
                
                mutated_antibody = mutateAntibody(...
                        antigens(antigen,:), antibodies(antibody,:), ...
                        delta_table(antibody) , beta);
                
                new_antibodies(new_antibody_count,:) = mutated_antibody;
                new_antibody_count = new_antibody_count +1;
            end
        end
    end

    antibodies = vertcat (antibodies, new_antibodies);

    % Antibody affinity calculation
    % see function affinityCalculation and previous comments
    
    affinity_table = calculationOfAfinityTable(antigens,antibodies);
             
    % Activated antibody detection and evaluation
    % see function antibodyDetection and previous comments
    
    [delta_table,enabled_antigens] = antibodyEvaluation(antigens, ...
                            antibodies,affinity_table,threshold);

    % Selection of the best antibodies
    % For each antigen from the training set, the set of antibodies 
    % activated by this antigen is determined. Only one antibody from 
    % this set, with the best evaluation delta, is selected to the next 
    % population. So the clonal expansion, unnecessary in this model, is
    % halted. The maximum number of antibodies in the next population is 
    % equal to the number of antigens, but the real number of antibodies 
    % is usually smaller because the same antibody could be selected by 
    % the several antigens (it depends on the value of the cross-reactivity 
    % threshold r). Outlier, i.e. antigen lying away from other antigens,
    % is represented by the separate antibody.
    
    new_antibodies = zeros(size(antigens));
    
    for antigen = 1:size(antigens,1)
        best_delta = Inf;
        best_antibody = [];
        for antibody = 1:size(antibodies,1)
            if (enabled_antigens(antigen,antibody))
                delta = delta_table(antibody);
                if(delta < best_delta)
                   best_delta = delta;
                   best_antibody = antibodies(antibody,:);
                end
            end
        end
        new_antibodies(antigen,:) = best_antibody;
    end

    antibodies = unique(new_antibodies,'rows');
    
    % Antibody affinity calculation
    % see function affinityCalculation and previous comments
    
    affinity_table = calculationOfAfinityTable(antigens,antibodies);
             
    % Activated antibody detection and evaluation
    % see function antibodyDetection and previous comments
    
    [delta_table,enabled_antigens] = antibodyEvaluation(antigens, ...
                            antibodies,affinity_table,threshold);
                        
    total_enable_antigens = sum(sum(enabled_antigens));
       
    if(total_enable_antigens == 0)
        if(diagnostics)
            fprintf('Total enable antigens are zero. Breaking out. \n')
        end
        break
    end

    if(isequal(old_antibodies,antibodies))
       if(diagnostics)
           fprintf('Antibodies have not evolved. Breaking out. \n')
       end
       break
    end
    
    % counter for the end of while
    iterations = iterations + 1;
    if(diagnostics)
        fprintf('%.1f %% ready (%d out of %d iterations). \n', ...
            ((iterations-1)/max_iterations)*100,iterations-1,max_iterations);
        
        fprintf('There are %d antibodies in this iteration. \n', ...
            size(antibodies,1) );
    end
    
end

% rmse for training data

% conversion of the train data to antigen with zero y-chain
% minus one is due to the last value is use to forecast
train_data_antigen = horzcat(train_data(1:end-1,:), ...
                        zeros(size(train_data(1:end-1,:))));

forecast_train_antigen = zeros(size(train_data_antigen,1),period_size*2);
forecast_train = zeros(size(train_data_antigen,1),period_size);
errors = zeros(size(train_data_antigen,1),period_size);

for antigen = 1:size(train_data_antigen,1)
    omega_set = [];
    new_threshold = threshold;
    while(isempty(omega_set))
        affinity_table = calculationOfAfinityTable(train_data_antigen,antibodies);
        for antibody = 1:size(antibodies)
            if (affinity_table(antigen,antibody) <= new_threshold)
                omega_set = vertcat(omega_set, antibodies(antibody,:));
            end
        end

        if(isempty(omega_set))
            new_threshold = new_threshold + (new_threshold*relax_threshold);
            if(diagnostics)
                fprintf('Relaxing threshold in training data. It should NEVER happen !!! \n')
            end
        end
    end
                           
    forecast_train_antigen(antigen,:)=forecastChain(omega_set, ...
                                train_data_antigen(antigen,:),new_threshold);
                            
    temp_forecast= forecast_train_antigen(antigen,:) .*  ...
        repmat(average_train(antigen),1,period_size*2);

    forecast_train(antigen,:) = temp_forecast(period_size+1:period_size*2);

    errors(antigen,:) = forecast_train (antigen,:) - ...
                (antigens(antigen,period_size+1:period_size*2) .*  ...
                 repmat(average_train(antigen),1,period_size) );
end

rmse = sqrt(mean((errors.^2)));

% Forecast procedure
% After learning the antibodies represent overlapping clusters of
% similar antigens. In the forecast procedure new antigen having only 
% x-chain is presented. The Omega set of antibodies, activated by 
% this antigen, is determined. The y-chains of these antibodies storage 
% average y-chains of antigens from the training set with similar x-chains.
% The y-chain of the input antigen is reconstructed from the y-chains of 
% the antibodies contained in the Omega set.
% Antibodies, closer to the antigen, have the higher influence on the 
% forecast forming. If an antigen is not recognized by antibodies, it 
% means that it represents a new shape of the load curve, not contained 
% in the training set. In this case the cross-reactivity threshold r is 
% consistently being increased until the antigen is recognized by one or
% more antibodies. The level of confidence in the forecast in such a case 
% is low and the forecast should be verified.


% conversion of the test data to antigen with zero y-chain
test_data_antigen = horzcat(test_data,zeros(size(test_data)));

forecast_antigen = zeros(size(test_data_antigen,1),period_size*2);
forecast = zeros(size(test_data_antigen,1),period_size);
confidence = zeros(size(test_data_antigen,1),1);

for antigen = 1:size(test_data_antigen)
    omega_set = [];
    new_threshold = threshold;
    while(isempty(omega_set))
        affinity_table = calculationOfAfinityTable(test_data_antigen,antibodies);
        for antibody = 1:size(antibodies)
            if (affinity_table(antigen,antibody) <= new_threshold)
                omega_set = vertcat(omega_set, antibodies(antibody,:));
            end
        end

        if(isempty(omega_set))
            new_threshold = new_threshold + (new_threshold*relax_threshold);
            if(diagnostics)
                fprintf('Relaxing threshold. \n')
            end
        end
    end
    
    confidence(antigen) = threshold /new_threshold;
    
    forecast_antigen(antigen,:)=forecastChain(omega_set, ...
                            test_data_antigen(antigen,:),new_threshold);

    temp_forecast = forecast_antigen(antigen,:) .* ...
                        repmat(average_test(antigen),1,period_size*2);

    forecast(antigen,:) = temp_forecast(period_size+1:period_size*2);
end

total_time = toc(time_start);

if(diagnostics)
    fprintf('The AIS run for %d minutes and %f seconds\n', ...
        floor(total_time/60),rem(total_time,60));
end

end

function [affinity_table] = calculationOfAfinityTable(antigens,antibodies)
%calculationOfAfinityTable This is a function for calculating the afinity
%table

    affinity_table = zeros(size(antigens,1),size(antibodies,1));

    for antigen = 1:size(antigens,1)
        for antibody = 1:size(antibodies,1)
            affinity_table(antigen,antibody) =  affinityCalculation( ...
                            antigens(antigen,:), antibodies(antibody,:));
        end    
    end

end

function [affinity_meas] = affinityCalculation(antigen, antibody)
%affinityCalculation This is a function for calculating the affinity 
%measure. 
%   The affinity measure is based on the distance between x-chains of
%   antigens and antibodies. The Euclidean distance is used

    % the size of antiges and antibody must match.
    x_size = size(antigen,2)/2;  
    affinity_meas = sqrt(sum(...
        power(antigen(1:x_size) - antibody(1:x_size),2)));
	
end

function [delta_table,enabled_antigens] = antibodyEvaluation( ...
                antigens,antibodies,affinity_table,threshold)
%antibodyEvaluation This is a function for detecting and evaluating
%the afinity table
    delta_table = zeros(size(antibodies,1),1);
    enabled_antigens = false(size(affinity_table));

    for antibody = 1:size(antibodies,1)
        delta = 0;
        for antigen = 1:size(antigens,1)
            if (affinity_table(antigen,antibody) <= threshold)
                delta = delta + forecastErrorCalculation( ...
                    antigens(antigen,:), antibodies(antibody,:));
                enabled_antigens(antigen,antibody) = true;
            end
        end
        % Possible bug
        % There is a chance of 0 enable_antigens. What should i do ???
        delta_table(antibody) = delta/sum(enabled_antigens(:,antibody));        
    end

end

function [forecast_error] = forecastErrorCalculation(antigen, antibody)
%forecastErrorCalculation This is a function for calulating the forecast
%error.
%   The forecast error is based on MAPE, which is traditionally used in
%   STFL models
    
    % the size of antiges and antibody must match.
    y_size = size(antigen,2)/2;   
    forecast_error = (sum(abs((antibody(y_size+1:end) ...
            - antigen(y_size+1:end)) ./ antigen(y_size+1:end)))*100)/y_size;
        
end

function [mutated_antibody] = mutateAntibody(antigen, antibody, delta, beta)
%mutateAntibody This is a function for mutating antibodies. 
%   The hypermutation is realized as follows. Each clone of the antibody 
%   is shifted towards different antigen lying in the recognition 
%   region of this antibody. The bigger the error δ for the given antigen 
%   is, the bigger shift toward this antigen is. This type of 
%   hypermutation produces new antibodies only in the regions covered 
%   by antigens

    % the size of antiges and antibody must match.
    ag_size = size(antigen,2);

    mutated_antibody = antibody(1:end) + ...
        ((2*ones(1,ag_size)./(exp(-normrnd(1,0.1,ag_size,1)'*beta*delta)+1))-1).* ...
        (antigen(1:end) - antibody(1:end));
    
end

function [forecast_antigen] = forecastChain(omega, new_antigen, threshold)
%forecastChain This is a function for forecasting a new antigen 
%   The y-chains of these antibodies storage average y-chains of antigens 
%   from the training set with similar x-chains. The y-chain of the input 
%   antigen is reconstructed from the y-chains of the antibodies contained 
%   in the omega set

    % the size of antiges and antibody must match.
    antibody_size = size(new_antigen,2);

    w = 1- (sqrt(sum(power(omega(:,1:antibody_size/2) - ...
      repmat(new_antigen(1:antibody_size/2), size(omega,1),1),2),2))/threshold);

    forecast= sum( (omega(:,antibody_size/2 + 1:antibody_size) ...
                    .* repmat(w,1,antibody_size/2)),1) ...
            ./ repmat(sum(w,1),1,antibody_size/2);

    forecast_antigen = horzcat( new_antigen(:,1:antibody_size/2),forecast);

end


