function [Z, Vx, Vy, Vz, Vp]=getDryParticleInfo(data, idRange)
Velocity=data.Velocity;%all the particle velocities at every t
Position=data.Position;
Radii=data.Radii;

Z=zeros(length(Position),length(idRange));
Vp=zeros(size(Z));
Vz=zeros(size(Z));
Vx=zeros(size(Z));

for i=1:length(Position)%遍历所有步长
    for j=idRange(1):idRange(end)%遍历指定颗粒
        XYZ=str2num(Position{i}{j});
        Z(i,j-idRange(1)+1)=XYZ(3);
        V=str2num(Velocity{i}{j});
        Vx(i,j-idRange(1)+1)=V(1);
        Vy(i,j-idRange(1)+1)=V(2);
        Vz(i,j-idRange(1)+1)=V(3);
        Rp=str2num(Radii{i}{j});
        Vp(i,j-idRange(1)+1)=4/3*pi*Rp^3;
    end
end
end