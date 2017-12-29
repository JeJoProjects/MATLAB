function [yy,mm,dd,hh,Mn,Sc] = Birthday( YY,MM,DD )
%%%% Age after 10000 years

    jd=G2JD(YY,MM,DD,0.0,0.0,0.0);
    jd=jd+10000.0;
    [yy,mm,dd,hh,Mn,Sc]=JD2G(jd);

end

