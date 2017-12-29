function beta = SingleDiffPhase( DPhiR,relvec )
lambda=0.190294; %%% wavelength of f1
Y=DPhiR(:,1);
vec=relvec(:,1:3);
l=sqrt(diag(vec*vec'));
e1=-diag(1.0./l)*vec;
A=[e1,repmat(lambda*eye(11),20,1)];
beta=inv(A'*A)*A'*Y;
end

