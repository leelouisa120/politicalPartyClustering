%% Error of Gaussian Model vs. Number of Data Points
clear all;
close all;

% Seeding Random Number Generator
rng('default');
rng(1);

% Generate Gaussian Data
% p = 10:10:10000;
p = [1000,2000,3000];

for ii = 1:length(p)
    MU1 = [2 -1];
    SIGMA1 = [1 0.6; 0.6 1.5];
    MU2 = [-2 1];
    SIGMA2 = [1 -0.3; -0.3 1.5];
    X = [mvnrnd(MU1,SIGMA1,p(ii));mvnrnd(MU2,SIGMA2,p(ii))];
    
    % Find the Gaussian Fit Model
    obj = gmdistribution.fit(X,2);
    
    % Show Parameter Values
    model_means = obj.mu;
%     model_means = sortrows(model_means,1);
    model_cov = obj.Sigma;
%     model_cov = sort(model_cov,3);

    
    % Mean Error
    Error1 = norm(abs(abs(MU1)-abs(model_means(2,:))));
%     Error = sqrt((MU1(1,1)-model_means(1,1))^2 + ((MU1(1,2)-model_means(1,2))^2));
    Error2 = norm(abs(abs(MU2)-abs(model_means(1,:))));
    MeanError(ii) = Error1 + Error2;
    
    % Covariance Error
    Error3 = norm(SIGMA1 - model_cov(:,:,2));
    Error4 = norm(SIGMA2 - model_cov(:,:,1));
    CovError(ii) = Error3 + Error4;
    
end
plot(p,MeanError);
title('Mean Error');
xlabel('Number of Data Points');
ylabel('Error');

figure;
plot(p,CovError);
title('Covariance Error');
xlabel('Number of Data Points');
ylabel('Error');
