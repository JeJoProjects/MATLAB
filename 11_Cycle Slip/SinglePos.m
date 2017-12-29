function beta = SinglePos(DPR,relvec,t) 
ind=find(DPR(:,2)==t);
Y=DPR(ind,1);
vec=relvec(ind,1:3);
l=sqrt(diag(vec*vec'));
e1=-diag(1./l)*vec;
A=[e1,ones(length(l),1)];
beta=inv(A'*A)*A'*Y;


end

