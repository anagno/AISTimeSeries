function [ forecast,forecast_antigen,iterations ] = AISShortForecasting( original_data, ... 
    threshold, training_percentage, max_iterations, beta )
%AISShortForecasting This is a function for forecasting time series using an
%artificial immune system
%   The function is implemented using the algorithm presented in [1].

%   Bibliography
%   [1] Dudek, Grzegorz. "Artificial immune system for short-term electric 
%   load forecasting." Artificial Intelligence and Soft Computing–ICAISC 
%   2008. Springer Berlin Heidelberg, 2008. 1007-1017.
%
% INPUT VARIABLES:
% data: the array of input data. It should be an array with the lines 
%   represinting the a hole period of measumerents.
% training_percentage: declares the percentage of the initilial
%   population that will be used for training the AIS
% threshold: declares the the cross-reactivity threshold r, for which an
%   antigen is activated
% max_iterations: The maximum number of iterations 
%
% OUTPUT VARIABLES:
%
%
%

switch nargin
    case 5
        
    case 4
        max_iterations = 50;        
    case 3
        max_iterations = 50;
        beta = 0.04;        
    case 2
        training_percentage = 2/3; 
        max_iterations = 50;
        beta = 0.04;
    otherwise
        error ('Too few or too many arguments were entered');
end

% !!! check if data is good with no null variables. !!!
assert(sum(sum(isnan(original_data)))==0,'NaN values in data');

% Preprocess the data in order to get rid of the time series trend and 
% seasonality, and simplify the model.

period_size = size(original_data,2);
average_period_data = mean(original_data,2);
data = original_data ./ repmat(average_period_data,1,period_size);

% Loading of the training set of antigens
% The whole dataset is divided into two subsets – training one and test 
% one. The first sequences of the time series (typically two thirds of the 
% whole time series) are included in the training set and the latest 
% sequences are included in the test set. Immune memory is trained using 
% the training set, and after learning the model is tested using the 
% test set.

cutoff_index = round(size(data,1) * training_percentage);
train_data = data([1:cutoff_index-1],:);
test_data = data([cutoff_index:end],:);

antigens = zeros(size(train_data,1)-1,period_size*2);

for n = 2:size(train_data,1)
    antigens(n-1,:) = horzcat(train_data(n-1,:),train_data(n,:));
end

if(training_percentage == 1)
   test_data = horzcat(data(end,:),zeros(1,period_size));
else
    test_data = horzcat(test_data,zeros(size(test_data)));
end

% Generation of the initial antibody population. An initial antibody 
% population is created by copping all the antigens from the training 
% set (antibodies and antigens have the same structure). This way of
% initialization prevents inserting antibodies in empty regions without 
% antigens

antibodies = antigens;

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
while(iterations <= max_iterations)
    
    % Store the population of antibodies to check if they change
    old_antibodies = antibodies;
    
    % Clonal selection
    % Each antibody cases secreting as many clones as many antigens are
    % in its recognition region. Thus most clones are generated in the
    % dense clusters of antigens.
    
    %new_antibodies = [];
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
                        
    total_enable_antigens = sum(sum(enabled_antigens));
                        
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
    
    %delta_table = zeros(size(antibodies,1),1);
    %enabled_antigens = false(size(affinity_table));

    for antigen = 1:size(antigens,1)
        best_delta = Inf;
        best_antibody = [];
        for antibody = 1:size(antibodies,1)
            if (enabled_antigens(antigen,antibody))
                %delta = forecastErrorCalculation( ...
                    %antigens(antigen,:), antibodies(antibody,:));
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
                   
    if(total_enable_antigens == 0)
        fprintf('Total enable antigens are zero. Breaking out. \n')
        break
    end

    if(isequal(old_antibodies,antibodies))
       fprintf('Antibodies have not evolved. Breaking out. \n')
       break
    end
    
    % counter for the end of while
    iterations = iterations + 1;
    fprintf('%.1f %% ready (%d out of %d iterations). \n', ...
        ((iterations-1)/max_iterations)*100,iterations-1,max_iterations);

    fprintf('There are %d antibodies in this iteration. \n', ...
        size(antibodies,1) );
    
end

% Forecast procedure
% After learning the antibodies represent overlapping clusters of
% similar antigens. In the forecast procedure new antigen having only 
% x-chain is presented. The Omega set of antibodies, activated by 
% this antigen, is determined. The y-chains of these antibodies storage 
% average y-chains of antigens from the training set with similar x-chains.
% The y-chain of the input antigen is reconstructed from the y-chains of 
% the antibodies contained in the Omega set.

forecast_antigen = [];

affinity_table = calculationOfAfinityTable(test_data,antibodies);

for antigen = 1:size(test_data)
    omega_set = [];
    for antibody = 1:size(antibodies)
        if (affinity_table(antigen,antibody) <= threshold)
            omega_set = vertcat(omega_set, antibodies(antibody,:) );
        end        
    end
    forecast_antigen(antigen,:)=forecastChain(omega_set, ...
                                test_data(antigen,:),threshold);
    if(training_percentage == 1)
        temp_forecast(antigen,:) = forecast_antigen(antigen,:) .*  ...
            repmat(average_period_data(size(average_period_data,1)),1,period_size*2);
    else
        temp_forecast(antigen,:) = forecast_antigen(antigen,:) .* ...
            repmat(average_period_data(cutoff_index + antigen -1),1,period_size*2);
    end
    forecast = temp_forecast (:,period_size+1:period_size*2);
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
                enabled_antigens(antigen,antibody) = 1;
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

function [forecast_antigen] = forecastChain(omega, antibody, threshold)
%mutateAntibody This is a function for forecasting a new antigen 
%   The hypermutation is realized as follows. Each clone of the antibody 
%   is shifted towards different antigen lying in the recognition 
%   region of this antibody. The bigger the error δ for the given antigen 
%   is, the bigger shift toward this antigen is. This type of 
%   hypermutation produces new antibodies only in the regions covered 
%   by antigens

    % the size of antiges and antibody must match.
    antibody_size = size(antibody,2);
    
    forecast_antigen = antibody;
       
    for k = (antibody_size/2):(antibody_size)
        sum_w = 0;
        sum_wy = 0;     
        for antigen = 1:size(omega)
            w = 1 - ((antibody(k) - omega(antigen,k))/threshold);
            sum_w = sum_w + w;
            sum_wy = sum_wy + omega(antigen,k)*w;            
        end
        forecast_antigen(1,k)= (sum_wy)/sum_w;
    end
    
end


