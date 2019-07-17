function [x,w,ind] = systematicResample(x,w)
% Performs a systematic resampling of a given set of samples

N = length(w);
u = ([0:N-1]'+rand(1))/N;
% u = ([0:N-1]+rand(1,N))/N;
wc = cumsum(w);
wc = wc/wc(N);
[~,ind1] = sort([u;wc]);
ind2 = find(ind1<=N);
ind = ind2 - (0:N-1)';
x = x(ind,:);
w = ones(N,1)./N;
end