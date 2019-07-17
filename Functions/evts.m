function [value,isterminal,direction] = evts(t,y,INFO)
% fprintf('Event time: %d \n', t);
eps = INFO.tolevt;
lexID = INFO.lexID;
nmodel = INFO.nmodel;
bmodel = INFO.b;
lbct = INFO.lbct;
indlb = INFO.indlb;
indub = INFO.indub;
U = INFO.U;
L = INFO.L;
P = INFO.P;
Q = INFO.Q;    
%% Update solutions
[lbx,ubx] = RHS( t,y,INFO );
ct = 0;
total = 0;
for i=1:nmodel
   total = length(bmodel{i}) + total; 
end
value = zeros(total,1);
isterminal = ones(total,1); % stop the integration
direction = -1*ones(total,1); % negative direction

for i=1:nmodel
    b = bmodel{i};
    lb = lbx(i,1:lexID(i));
    ub = ubx(i,1:lexID(i));
    lb(indlb{i}) = [];
    ub(indub{i})=[];
    b(1:length(lb)) = lb;
    b(length(lb)+lbct(i)+1:length(lb)+lbct(i)+length(ub)) = ub;
    x = (L{i}\(P{i}*b));
    x = U{i}\x;
    x = Q{i}*x;
%     x = U{i}\(L{i}\(P{i}*b));
%     x = INFO.B{i}\b;
% Detect when a basic variable crosses zero.
    value(1+ct:length(x)+ct) = x + eps;
    ct = ct + length(x);
end
end


