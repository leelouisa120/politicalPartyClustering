clear all;
close all;

dimension = input('Number of Variables: ');

switch(dimension)
    case 2
        % Reading Data from Excel Sheet
        num_of_cluster = input('Number of Clusters: ');
        filename = 'ANES 2012-2008 Data.xlsx';
        sheet = 4;
        xlRange = 'F5:G8240';
        subset = xlsread(filename,sheet,xlRange);
        
        % Finding Correlation Matrix for Data
        [RHO,PVAL] = corr(subset);
        
        % Find the Gaussian Fit Model
        options = statset('Display','final','MaxIter',500,'TolFun',1e-9);
        obj = gmdistribution.fit(subset,num_of_cluster,'CovType','diagonal','Options',options,'Regularize',1e-5);
       
        % Save 2 Parameter Values
        model_means = obj.mu;
        model_cov = obj.Sigma;
        
        % Adding Noise to the repeated rows
        [~,unique_rows,~] = unique(subset,'rows');
        duplicate_rows = setdiff(1:size(subset,1),unique_rows);
        subset(duplicate_rows,:) = jitter(subset(duplicate_rows,:), [], 1);
        
        % Inidicating Party Affiliation through Color
        partysheet = 4;
        Range = 'E5:E8240';
        partyscale = xlsread(filename,partysheet,Range);
        
        for ii = 1:length(partyscale)
            color = mapcolor(partyscale(ii));
            partycolor(ii,:) = color;
        end
        
        % Plotting Data Set
        scatter(subset(:,1),subset(:,2),25,partycolor,'filled');
        hold on;

        % Plot the contour of the Gaussian Model
        h = ezcontour(@(x,y)pdf(obj,[x y]),[0 10],[0 10]);
        h.LineWidth = 2.0;
        
        
    case 8
        % Reading Data from Excel Sheet
        num_of_cluster = input('Number of Clusters: ');
        filename = 'ANES 2012-2008 Data.xlsx';
        sheet = 4;
        xlRange = 'F5:M8240';
        subset = xlsread(filename,sheet,xlRange);
        
        % Finding Correlation Matrix for Data
        [RHO,PVAL] = corr(subset);
        
        % Find the Gaussian Fit Model
        options = statset('Display','final','MaxIter',500,'TolFun',1e-9);
        obj = gmdistribution.fit(subset,num_of_cluster,'CovType','diagonal','Options',options,'Regularize',1e-5);
        
        % Show Parameter Values
        model_means = obj.mu;
        model_cov = obj.Sigma;   
        
        % Plotting cluster means onto PCA axes
        mu = mean(subset); % find mean of data
        
        % Adding Noise to the repeated rows
        [~,unique_rows,~] = unique(subset,'rows');
        duplicate_rows = setdiff(1:size(subset,1),unique_rows);
        subset(duplicate_rows,:) = jitter(subset(duplicate_rows,:), [], 1);
       
        % Indicating Party Affiliation through Color
        partysheet = 4;
        Range = 'E5:E8240';
        partyscale = xlsread(filename,partysheet,Range);
        
        for ii = 1:length(partyscale)
            color = mapcolor(partyscale(ii));
            partycolor(ii,:) = color;
        end

        % Conductin PCA Analysis on the data
        [coefs,score] = pca(subset);

        % Plotting Data with Party Affiliation colors onto PCA axes
        for jj = 1:length(partycolor)
            plot(score(jj,1),score(jj,2),'o','Color',partycolor(jj,:));
            hold on;
        end
        
        title('2012-2008 Principle Component Analysis');
        xlabel('Component 1');
        ylabel('Component 2');
        hold on;

        % find model mean in the centered space
        for kk = 1:size(model_means,1)
        mean_centered(kk,:) = model_means(kk,:) - mu;  
        end
        
        % now mapp centered model_mean to the pca components 
        mean_mapped = mean_centered/coefs';
        plot(mean_mapped(:,1), mean_mapped(:,2), 'kx','LineWidth',4,'MarkerSize',20);
        
end

        
        