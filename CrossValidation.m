close all;
clear all;

% 2012-2 Issues: Reading Data from Excel Sheet
% num_of_cluster = input('Number of Clusters: ');
% filename = 'ANES 2012 Data.xlsx';
% sheet = 4;
% xlRange = 'E5:G5918';
% data = xlsread(filename,sheet,xlRange);

% 2012-14 Issues: Reading Data from Excel Sheet
num_of_cluster = input('Number of Clusters: ');
filename = 'ANES 2012 Data.xlsx';
sheet = 4;
xlRange = 'E5:S5918';
data = xlsread(filename,sheet,xlRange);

subset = data(:,2:end);    % All columns but first are political issues
partyscale = data(:,1);    % first column is party affiliation

for jj = 1:num_of_cluster
    disp(jj);
    N = size(subset,1);
    K = 5;
    
    for kk = 1:25
        Indices = crossvalind('Kfold', N, K);
        
        for ii = 1:K
            
            train_data = subset(~(Indices == ii),:);
            test_data = subset(Indices == ii,:);
            
            % Find the Gaussian Fit Model
%             'Display','final',
            options = statset('MaxIter',5000,'TolFun',1e-9);
            obj = gmdistribution.fit(train_data,jj,'CovType','diagonal','Options',options,'Regularize',1e-5);
            
            % Show Parameter Values
            model_means = obj.mu;
            model_cov = obj.Sigma;
            
            clear probfunc;
            
            % Calculate Probability Density Function for Each Individual
            for cc = 1:jj
                mu = model_means(cc,:);
                sigma = model_cov(:,:,cc);
                probfunc(:,cc) = mvnpdf(test_data,mu,sigma);
            end
            
            log_of_pdf = log((sum(probfunc,2))./jj);
            LogLikelihood(ii) = sum(log_of_pdf);
        end
        AVG_LogLikelihood_iteration(kk,:) = mean(LogLikelihood);
    end
    AVG_LogLikelihood(jj,:) = mean(AVG_LogLikelihood_iteration);
    numParam = jj + (size(subset,2)*jj);
    aic(jj,:) = aicbic(AVG_LogLikelihood(jj,:),numParam);
end

clusters = (1:num_of_cluster)';
plot(clusters,aic,'o-','LineWidth',3)
title('AIC Cluster Analysis');
xlabel('Number of Clusters');
ylabel('AIC');
