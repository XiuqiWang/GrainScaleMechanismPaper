function [ezVector, exVector, LdiffVector, VimVector, thetaimVector, thetareVector, EimVector, Smean] = portDataToVectorSal(Par,Vx,Vz,X,Z,Vp,LC,N_inter)
%Arrays containing values from N_inter intervals
ez = cell(1,N_inter);ex=ez;Vim=ez;thetaim=ez;thetare=ez;Eim=ez;Ldiff=ez;
%vectors in the whole 5s
ezVector = [];exVector = [];LdiffVector = [];VimVector=[];thetaimVector=[];thetareVector=[];EimVector=[];
%Loop through each cell, flatten the multi-dimensional array, and concatenate
for i = 1:numel(Par)
    %Restitution coefficient and Particle liquid change
    ID_S = Par{i}{1}; %the id of saltation for particle i
    [flattenedArray_ez,flattenedArray_ex,flattenedArray_L,flattenedArray_Vim,flattenedArray_thetaim,flattenedArray_thetare,flattenedArray_Eim, ID_realS] = CalRestitutionLCdiff(ID_S,Vz(:,i),X(:,i),Z(:,i),Vx(:,i),Vp(1,i),LC(:,i));
    ezVector = [ezVector; flattenedArray_ez];
    exVector = [exVector; flattenedArray_ex];
    EimVector = [EimVector; flattenedArray_Eim];
    LdiffVector = [LdiffVector; flattenedArray_L];
    VimVector = [VimVector; flattenedArray_Vim];
    thetaimVector = [thetaimVector; flattenedArray_thetaim];
    thetareVector = [thetareVector;flattenedArray_thetare];
    
    %determine which interval the values should go, loop within the real
    %saltation ids
if ~isempty(ID_realS)
    for j = 1:length(ID_realS)
    % 找到值的范围
    idx = ceil(ID_realS(j) / (size(X,1)/N_inter)); % 每个cell表示一个(size(X,1)/Ninter个元素的范围
    if idx > N_inter
        idx = N_inter; % 防止溢出
    end
    ez{idx}=[ez{idx}, flattenedArray_ez(j)];
    ex{idx}=[ex{idx}, flattenedArray_ex(j)];
    Ldiff{idx}=[Ldiff{idx}, flattenedArray_L(j)];
    Vim{idx}=[Vim{idx}, flattenedArray_Vim(j)];
    thetaim{idx}=[thetaim{idx}, flattenedArray_thetaim(j)];
    thetare{idx}=[thetare{idx}, flattenedArray_thetare(j)];
    end
end
end
Smean=NaN(6,N_inter);
for i=1:size(Smean,2)
    if ~isempty(ez{i})
    Smean(1,i)=mean(ez{i});
    end
    if ~isempty(ex{i})
    Smean(2,i)=mean(ex{i});
    end
    if ~isempty(Ldiff{i})
    Smean(3,i)=mean(Ldiff{i});
    end
    if ~isempty(Vim{i})
    Smean(4,i)=sum(Vim{i});
    end
    if ~isempty(thetaim{i})
    Smean(5,i)=mean(thetaim{i});
    end
    if ~isempty(thetare{i})
    Smean(6,i)=mean(thetare{i});
    end
end
end

function [ez, ex, dLC, Vim, thetaim, thetare, Eim, ID_realS]=CalRestitutionLCdiff(IDi,Vzi,Xi,Zi,Vxi,Vpi,LCi)
ez=[];ex=[];dLC=[];Vim=[];thetaim=[];thetare=[];Eim=[];ID_realS=[];
if ~isempty(IDi)
    [ex_ini, ez_ini, V1, thetaim_ini, thetare_ini, Eim_ini, dLC_ini, ID_realS_ini]=CalEfromXZ(Xi,Zi,Vxi,Vpi,LCi,IDi);
   %ID_neg = IDi(Vzi(IDi)<0);
%    V1n=sqrt(Vxi(ID_neg).^2+Vzi(ID_neg).^2);
%    V2n=sqrt(Vxi(ID_neg+1).^2+Vzi(ID_neg+1).^2);
%    ex_ini = V2n./V1n;
%    ez_ini = abs(Vzi(ID_neg+1)./Vzi(ID_neg));
%    ID_lto = intersect(IDi(ex_ini<1),IDi(ez_ini<1));%only keep the ids where ez and ex are both below 1
%    V1n_new = sqrt(Vxi(ID_lto).^2+Vzi(ID_lto).^2);
%    V2n_new = sqrt(Vxi(ID_lto+1).^2+Vzi(ID_lto+1).^2);
%    ex_ini = V2n_new./V1n_new;
%    ez_ini = abs(Vzi(ID_lto+1)./Vzi(ID_lto));

   ez=[ez; ez_ini];
   ex=[ex; ex_ini];
   Vim=[Vim; V1];
   Eim=[Eim; Eim_ini];
   thetaim=[thetaim; thetaim_ini];
   thetare=[thetare; thetare_ini];
   dLC=[dLC; dLC_ini];
   ID_realS=[ID_realS; ID_realS_ini];
   
%    ID_pos = IDi(Vzi(IDi)>0);
%    V1p=sqrt(Vxi(ID_pos-1).^2+Vzi(ID_pos-1).^2);
%    V2p=sqrt(Vxi(ID_pos).^2+Vzi(ID_pos).^2);
%    ex_ini2 = V2p./V1p;
%    ez_ini2 = abs(Vzi(ID_pos)./Vzi(ID_pos-1));
%    ID_lto2 = intersect(ID_pos(ex_ini2<1),ID_pos(ez_ini2<1));%only keep the ids where ez and ex are both below 1
%    V1p_new = sqrt(Vxi(ID_lto2-1).^2+Vzi(ID_lto2-1).^2);
%    V2p_new = sqrt(Vxi(ID_lto2).^2+Vzi(ID_lto2).^2);
%    ex_ini2 = V2p_new./V1p_new;
%    ez_ini2 = abs(Vzi(ID_lto2)./Vzi(ID_lto2-1));
%    dLC_ini2 = LCi(ID_lto2)-LCi(ID_lto2-1);
%    theta2 = atand(abs(Vzi(ID_lto2-1))./Vxi(ID_lto2-1));
%    theta2r = atand(abs(Vzi(ID_lto2))./Vxi(ID_lto2));
%    Vim=[Vim; V1p_new];
%    thetaim=[thetaim; theta2];
%    thetare=[thetare; theta2r];
%    ez=[ez; ez_ini2];
%    ex=[ex; ex_ini2];
%    dLC=[dLC; dLC_ini2];
end
end

function [ex_ini, ez_ini, V1, thetaim_ini, thetare_ini, Eim_ini, dLC_ini, ID_realS_ini]=CalEfromXZ(X,Z,Vx,Vp,LC,ID)
dt=0.01;Lx=0.025;
dLC_ini = LC(ID+1)-LC(ID-1);

Vx2 = (X(ID+1)-X(ID))/dt;
Vx2stored = Vx(ID+1);
Index_neg2 = find(Vx2<0 & Vx2stored>0);
Vx2(Index_neg2) = Vx2(Index_neg2)+Lx/dt;
Index_pos2 = find(Vx2>0 & Vx2stored<0);
Vx2(Index_pos2) = Vx2(Index_pos2)-Lx/dt;

Vx1 = (X(ID)-X(ID-1))/dt;
Vx1stored = Vx(ID-1);
Index_neg1 = find(Vx1<0 & Vx1stored>0);
Vx1(Index_neg1) = Vx1(Index_neg1)+Lx/dt;
Index_pos1 = find(Vx1>0 & Vx1stored<0);
Vx1(Index_pos1) = Vx1(Index_pos1)-Lx/dt;

Vz2 = (Z(ID+1)-Z(ID))/dt;
Vz1 = (Z(ID)-Z(ID-1))/dt;

V2n = sqrt(Vx2.^2+Vz2.^2);
V1n = sqrt(Vx1.^2+Vz1.^2);
ex_ini = V2n./V1n;
ez_ini = abs(Vz2./Vz1);

ID_lto = find(ex_ini<1 & ez_ini<1);
ID_realS_ini = ID(ID_lto);
ex_ini = ex_ini(ID_lto);
ez_ini = ez_ini(ID_lto);
dLC_ini = dLC_ini(ID_lto);
V1 = V1n(ID_lto);
Eim_ini = 0.5*Vp*2650*V1.^2;%impact energy
thetaim_ini = atand(abs(Vz1(ID_lto)./Vx1(ID_lto)));
thetare_ini = atand(abs(Vz2(ID_lto)./Vx2(ID_lto)));
end