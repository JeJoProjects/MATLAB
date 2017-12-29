function y = CRshiftet(t,tau,Deltaf )

global If

y=cos(2*pi*(If+Deltaf)*(t-tau));

end

