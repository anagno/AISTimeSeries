function PlotFitModel(x,y,training_percentage,model_name,fitresult)

% This function visualizes the fitted model on a given set of training and
% testing data points.

% Get the corresponding training and testing subsets.
[x_train,x_test,y_train,y_test] = SplitTrainingTestingData(x,y,training_percentage);
% Get the corresponding range for the estimated y-values.
y_est = feval(fitresult,x);
% Get the minimum and maximum x-values.
x_min = min(x);
x_max = max(x);
% Get the minimum and maximum y-values.
y_min = min(min(y),min(y_est));
y_max = max(max(y),max(y_est));
% Perform the actual plotting operations.
figure_name = strcat([model_name 'Fit Visualization']);
legend_name = strcat([model_name,'Fit']);
figure('Name',figure_name);
hold on
plot(x_train,y_train,'*g');
plot(x_test,y_test,'*r');
plot(x,y_est,'-b','LineWidth',1.6);
axis([x_min-1 x_max+1 y_min-1 y_max+1]);
legend('Training Data','Testing Data',legend_name,'Location','NorthWest');
% Label axes
xlabel('Time');
ylabel('co2ppm');
grid on
hold off
end

