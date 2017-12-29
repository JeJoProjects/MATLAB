function JD=G2JD(yy,mm,dd,hh,min,sec)
    if mm >2
        Y=yy;
        M=mm;
    else
        Y=yy-1;
        M=mm+12;
    end
    D=dd;
    H=(hh/24.0+min/1440.0+sec/86400.0);
    A=floor(Y/100.0);
    B=2-A+floor(A/4);
    JD=floor(365.25*(Y+4716.0))+floor(30.6001*(M+1))+D+H+B-1524.5;
end