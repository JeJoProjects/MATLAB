function beta = SD_epoch( Pr1,Pr2,x01,x02,SatPos,epoch )

ind=find(Pr1(:,2)==epoch);
sat_epoch=SatPos(ind,1:3);
PR1_epoch=Pr1(ind,1);
PR2_epoch=Pr2(ind,1);

DeltaRho=PR1_epoch-PR2_epoch;  %%%% single differences of pseudoranges

rv1=sat_epoch-repmat(x01,length(ind),1);
rv2=sat_epoch-repmat(x02,length(ind),1);

rho01=sqrt(diag(rv1*rv1'));
rho02=sqrt(diag(rv2*rv2'));

Y=DeltaRho-(rho01-rho02);

A=[-diag(1.0./rho01)*rv1, ones(length(ind),1)];

beta=A\Y;
end

