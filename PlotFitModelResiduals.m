function PlotFitModelResiduals(x,y,training_percentage,model_name,fitresult)

% This function visualizes the fitted model residuals on a given set of 
% training and testing data points.

% Get the estimated y-values.
y_est = feval(fitresult,x);
% Split real and estimated y-values into corresponding training and testing
% subsets.
[y_train,y_test,y_est_train,y_est_test] = SplitTrainingTestingData(y,y_est',training_percentage);
% Split x-values into corresponding training and testing subsets.
[x_train,x_test,~,~] = SplitTrainingTestingData(x,[],training_percentage);
% Compute residuals for the training and testings stages.
res_train = y_train - y_est_train;
res_test = y_test - y_est_test;
% Get the minimum and maximum x-values.
x_min = min(x);
x_max = max(x);
% Get the minimum and maximum y-values.
y_min = min(min(res_train),min(res_test));
if(y_min>0)
    y_min = 0;
end;
y_max = max(max(res_train),max(res_test));
% Perform the actual plotting operations.
figure_name = strcat([model_name 'Fit Residuals Visualization']);
figure('Name',figure_name);
hold on
plot(x_train,res_train,'*g');
plot(x_test,res_test,'*r');
plot(x,zeros(1,length(x)),'-b','LineWidth',1.6);
axis([x_min-2 x_max+2 y_min-2 y_max+2]);
legend('Training Residuals','Testing Residuals','Location','NorthWest');
% Label axes
xlabel('Time');
ylabel('Residuals');
grid on
hold off
end

