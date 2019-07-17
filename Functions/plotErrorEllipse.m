function plotErrorEllipse(mean, cov, confidence, color, linewidth)
%INPUT: X is the unnormalized, uncentered data, c is the confidence
%region (i.e., 0.95 implies 95% of the data explained, 5% rejected)

%OUTPUT: e is a vector with elliptical data points for plotting along the
%rows (# rows = size of X)

r = sqrt(chi2inv(confidence,length(mean))); % Get radius
[V, D] = eig(cov); %Eigenvalue decomposition of covariance matrix, V is [m by m] and D is [m by m]
[D, order] = sort(diag(D), 'descend'); %Sorts the eigenvalues into a vector in descending order
D = diag(D); %Converts the vector into a diagonal matrix
V = V(:, order); %Orders the columns of V according to the sort
theta = 0:0.1:2.2*pi; %Angle values
e = [cos(theta); sin(theta)]; %Unit circle points
VV = V*sqrt(D)*r; %Scaled eigenvectors of S
e = bsxfun(@plus, VV*e, mean); %Project circle onto original space

plot(e(1,:), e(2,:), '-', 'color', color, 'LineWidth', linewidth)

end