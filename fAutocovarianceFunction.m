function [ gamma ] = fAutocovarianceFunction( y, m )
% returns the autocovariance function up to lag m
gamma = fautocovariance(y, 0);
for k = 1:m 
    gamma = [gamma; fautocovariance(y,k)];
end
end

