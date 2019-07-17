function Xin = getSamplesClass(N, myTd, tf, class, sample_type)
% Gets samples belonging to a single class 
% class = 1 -- before discontinuity has occurred
% class = 2 -- after discontinuity has occurred

index_in = [];
factor = 5;
while length(index_in) < N
    X = uq_getSample(myTd.myInput, factor*N, sample_type);
    Y = uq_evalModel(myTd.myPCE, X);
    if class == 1
        index_in = find(Y > tf);
    elseif class == 2
        index_in = find(Y <= tf);
    end
    factor = 2*factor;
end
Xin = X(index_in(1:N),:);
end