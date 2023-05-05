clear all
close all
clc


%gridgen
%% Data
param.xl=4000;param.yl=4000;param.r=100*sqrt(2);param.del=param.r/sqrt(2);
g=gridgen(param);
param.alpha=0.5;param.beta=0.1;param.Q=100;param.Tmi=0.1;param.Tma=1;param.Nv=5;param.Tar=10;
num_iter=100;param.mode='targetdel';param.Nx=g.Nx;param.Ny=g.Ny;
ntl=[];
%% Initialisation
for it=1:10
upd.Tau=ones(g.Nx,g.Ny);upd.NTL=0;
upd.UAVlocs=[randi(g.Nx,param.Nv,1) randi(g.Ny,param.Nv,1) ];
%upd.UAVlocs=[10 10;14 10;12 14;12 6;13 9];
upd.UAVlocstheta=randi(8,param.Nv,1);%[1;5;7;3;4];%%[2];%4;6;7;2];
%upd.UAVlocstheta=[1;5;7;3;4];
upd.Tarlocs2=[randi(g.Nx,param.Tar,1) randi(g.Ny,param.Tar,1)];%RandomTargets
upd.Tarlocs=upd.Tarlocs2;
Mat=[];
upd.TF=[];
for i=1:num_iter
    delTau=zeros(param.Nx,param.Ny);
    for j=1:param.Nv
        if i==1
        upd.UAV(j).locsnext=[];
        upd.UAV(j).thetanext=[];
        toplot.UAV(j).locations=[];
        upd.UAV(j).locs=upd.UAVlocs(j,:);
        upd.UAV(j).theta=upd.UAVlocstheta(j);
        end
        %storing the locations
        toplot.UAV(j).locations(i,:)=[upd.UAV(j).locs upd.UAV(j).theta];

         upd.UAV(j).tau=zeros(param.Nx,param.Ny);
        %check for target match
        upd=targetMatch(param,upd,i,j);
        upd=nextstate(param,upd,i,j);
        upd.UAV(j).locs=  upd.UAV(j).locsnext;
        upd.UAV(j).theta=  upd.UAV(j).thetanext;
        delTau=delTau+upd.UAV(j).tau;
       % delTau=upd.UAV(j).tau; 

    end
    upd.Tau=upd.Tau-delTau;
    
    for k1=1:param.Nx
        for k2=1:param.Ny
            if(upd.Tau(k1,k2)<param.Tmi)
                upd.Tau(k1,k2)=param.Tmi;
            elseif(upd.Tau(k1,k2)>param.Tma)
                upd.Tau(k1,k2)=param.Tma;
            else
                upd.Tau(k1,k2)=upd.Tau(k1,k2);
            end
        end
    end
    if size(upd.Tarlocs,1)>=1
        for i2=1:size(upd.Tarlocs,1)
            r2=upd.Tarlocs(i2,1);c2=upd.Tarlocs(i2,2);
               upd.Tau(r2,c2)=1;
        end
    end
    
    if upd.NTL>=1
    for i1 =1:length(upd.NTL)
        r=upd.TF(i1,1);c=upd.TF(i1,2);
        upd.Tau(r,c)=0;
    end
    end
   Mat(:,:,i)= upd.Tau;
    %heatmap(upd.Tau)
end
%end
for i=1:param.Nv
   
    plot(toplot.UAV(i).locations(:,1),toplot.UAV(i).locations(:,2))
    hold on
    grid on
end
 %comet(toplot.UAV(i).locations(:,1),toplot.UAV(i).locations(:,2))
 plot(upd.Tarlocs2(:,1),upd.Tarlocs2(:,2),'or')
hold off


%end
ntl(it)=upd.NTL;
end
ave=sum(ntl)/10
%functions
function y=nextstate(param,upd,it,j)
nx=upd.UAV(j).locs(1);ny=upd.UAV(j).locs(2);theta=upd.UAV(j).theta;
xvec=[1 1 0 -1 -1 -1 0 1]';yvec=[0 1 1 1 0 -1 -1 -1]';
indices=theta-2:theta+2;
indices2=indices;
for i=1:5
    if(indices(i)>8)
        indices(i)=indices(i)-8;
    elseif(indices(i)<1)
        indices(i)=8+indices(i);
    else
        disp('');
    end
end
nn=[nx*ones(5,1)+xvec(indices) ny*ones(5,1)+yvec(indices)];
vec=[1/6 5/24 1/4 5/24 1/6]';
for i=1:length(nn(:,1))
    if(nn(i,1)>param.Nx)
        nn(i,1)=param.Nx;
        vec(i)=0;
    elseif(nn(i,1)<1)
        nn(i,1)=1;
        vec(i)=0;
    else
        nn(i,1)=nn(i,1);
    end
     if(nn(i,2)>param.Ny)
        nn(i,2)=param.Ny;
        vec(i)=0;
    elseif(nn(i,2)<1)
        nn(i,2)=1;
        vec(i)=0;
    else
        nn(i,2)=nn(i,2);
    end
end
nn_cords=nn*param.del;
distance=sqrt((nx*param.del-nn_cords(:,1)).^2 +(ny*param.del-nn_cords(:,2)).^2);


for i=1:length(distance)
    if distance(i)==param.del
        eta(i)=1;
    elseif distance(i)==sqrt(2)*param.del
        eta(i)=1.4;
    else
        eta(i)=0;
    end
end
deltheta=[90 45 0 45 90]';
etaijk=vec.*eta';
tau_p=diag(upd.Tau(nn(:,1),nn(:,2)));
prob_N=(tau_p.^(param.alpha)).*(etaijk.^(param.beta));
prob=prob_N/sum(prob_N);
[~,id]=max(prob);
upd.UAV(j).locsnext=nn(id,:);lol=nn(id,:)
upd.UAV(j).thetanext=indices(id);
%contribution to change in pheromones
%calc distances from curr location to other
xc=nx*param.del;yc=ny*param.del;


for s=1:param.Nx
    for t=1:param.Ny
        D(s,t)=sqrt((xc-s*param.del)^2 +(yc-t*param.del)^2);
        if (D(s,t)>0 && D(s,t)<500)
            upd.UAV(j).tau(s,t)=exp(-0.5*(param.Q/D(s,t))^2)/(2*sqrt(2)*pi);
        elseif(D(s,t)==0)
            upd.UAV(j).tau(s,t)=0.9;
        else
            upd.UAV(j).tau(s,t)=0;
        end
    end
end
y=upd;
end


function y=targetMatch(param,upd,i,j)

if sum(ismember(upd.Tarlocs,upd.UAV(j).locs,'rows'))>0
    row_id=find(ismember(upd.Tarlocs,upd.UAV(j).locs,'rows'));
    upd.NTL=upd.NTL+1;
    upd.TF(upd.NTL,:)=upd.UAV(j).locs;
    if(param.mode=='targetadd')
        upd.Tarlocs(row_id,:)=[randi(param.Nx,1,1) randi(param.Ny,1,1)];
    elseif(param.mode=='targetdel')
        upd.Tarlocs(row_id,:)=[];
    else
        disp('');
    end

end
y=upd;
end



function y=gridgen(p)
y.Nx=p.xl/p.del;y.Ny=p.yl/p.del;
y.xcord=[1:1:y.Nx]*p.del;y.ycord=[1:1:y.Ny]*p.del;
end