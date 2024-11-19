function [X, Y, Z, Vx, Vy, Vz, Vp, W]=getDryParticleInfoFromData(data, idcase, idRange)
eval(['data=data.data',num2str(idcase),';']);
Z=zeros(size(data,2),length(idRange));
X=zeros(size(Z));
Y=zeros(size(Z));
Vp=zeros(size(Z));
Vz=zeros(size(Z));
Vx=zeros(size(Z));
Vy=zeros(size(Z));
W=zeros(size(Z));

for i=1:size(data,2)%遍历所有步长
    for j=idRange(1):idRange(end)%遍历指定颗粒
        X(i,j-idRange(1)+1)=data{i}.Position(j,1);
        Y(i,j-idRange(1)+1)=data{i}.Position(j,2);
        Z(i,j-idRange(1)+1)=data{i}.Position(j,3);
        V=data{i}.Velocity(j,:);
        Vx(i,j-idRange(1)+1)=V(1);
        Vz(i,j-idRange(1)+1)=V(3);
        Rp=data{i}.Radius(j);
        Vp(i,j-idRange(1)+1)=4/3*pi*Rp^3;
    end
end
end