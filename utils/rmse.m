% function r=rmse(data,estimate)
% % Function to calculate root mean square error from a data vector or matrix 
% % and the corresponding estimates.
% % Usage: r=rmse(data,estimate)
% % Note: data and estimates have to be of same size
% % Example: r=rmse(randn(100,100),randn(100,100));
% 
% % delete records with NaNs in both datasets first
% I = ~isnan(data) & ~isnan(estimate); 
% data = data(I); estimate = estimate(I);
% 
% r=sqrt(sum((data(:)-estimate(:)).^2)/numel(data));

function [ rmse ] = rmse( in, true, use_abs )
%RMSE Summary of this function goes here
%   Detailed explanation goes here

if nargin < 3
    rmse = 100 * norm(in(:) - true(:)) / norm(true(:));
end

if (nargin == 3) && (use_abs == 1)
    rmse = 100 * norm(abs(in(:)) - abs(true(:))) / norm(abs(true(:)));
end

end

