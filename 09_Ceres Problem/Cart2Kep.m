function el = Cart2Kep( pos,vel,t )
el=zeros(3,length(t));
CR=cross(pos,vel);
el(4,:)=atan2(-CR(2,:),CR(1,:));
el(3,:)=atan2(CR(3,:),sqrt(CR(1,:).^2+CR(2,:).^2));
end

