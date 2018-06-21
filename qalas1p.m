
function[M,Mmeas]=qalas1p(Minit,M0,T1,T2,TR,TE_T2prep,flipAngle,nacq,dt)

[M,Mmeas]=qalas1time(Minit,M0,T1,T2,TR,TE_T2prep,flipAngle,nacq,dt);
for iii=1:100
    [M,Mmeas]=qalas1time(M(end),M0,T1,T2,TR,TE_T2prep,flipAngle,nacq,dt);
    if abs(M(1)-M(end))<=0.0001; break; end;
end

end

function[M,Mmeas]=qalas1time(Minit,M0,T1,T2,TR,TE_T2prep,flipAngle,nacq,dt)

star=(1-exp(-TR./T1))./(1-cosd(flipAngle).*exp(-TR./T1));

% T2 sensitization
M(1)=Minit;        % assume initialization?
M(2)=M(1).*exp(-TE_T2prep./T2);
M(3)=M0.*star-(M0.*star-M(2)).*exp(-dt(3)./(T1.*star));
M(4)=M0-(M0-M(3)).*exp(-dt(4)./T1);

% T1 sensitization
M(5)=-M(4);             % Assume perfect inversion? 100 ms inversion time
M(6)=M0-(M0-M(5)).*exp(-dt(6)./T1);

% Post T1 sens acquisitions
for iii=1:nacq-1
    M(5+2*iii)=M0.*star-(M0.*star-M(4+2*iii)).*exp(-dt(5+2*iii)./(T1.*star));
    M(6+2*iii)=M0-(M0-M(5+2*iii)).*exp(-dt(6+2*iii)./T1);
end

% M=sind(flipAngle).*M;
Mmeas=sind(flipAngle).*[M(2),M(6:2:end-1)];

end