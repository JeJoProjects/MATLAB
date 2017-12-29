function [ beta, GDOP ] = SPP_epoch(PR,sat,epoch,x0 )



ind=find(PR(:,2)==epoch);
sat_epoch=sat(ind,1:3);
PR_epoch=PR(ind,1);

relvec=sat_epoch-repmat(x0,length(ind),1);

rho0=sqrt(diag(relvec*relvec'));

Y=PR_epoch-rho0;

A=[-diag(1.0./rho0)*relvec, ones(length(ind),1)];
    
Niv=inv(A'*A);

GDOP=sqrt(trace(Niv));
beta=Niv*A'*Y;

end

