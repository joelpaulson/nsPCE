function [a,b] = getBoundsClass(myTd, tf, class, P_lb, P_ub)
% Calculates the outer-bounding box for a given class
% class = 1 -- before discontinuity has occurred
% class = 2 -- after discontinuity has occurred

X = myTd.X;
Y = myTd.Y;
n = size(X,2);
if class == 1
    index = find(Y > tf);
elseif class == 2
    index = find(Y <= tf);
end
a = zeros(1,n);
b = zeros(1,n);
for i = 1:n
    a(i) = min(X(index,i))*0.999;
    b(i) = max(X(index,i))*1.001;
    if abs((a(i) - P_lb(i))/P_lb(i)) < 1e-2
        a(i) = P_lb(i);
    end
    if abs((b(i) - P_ub(i))/P_ub(i)) < 1e-2
        b(i) = P_ub(i);
    end
end
end