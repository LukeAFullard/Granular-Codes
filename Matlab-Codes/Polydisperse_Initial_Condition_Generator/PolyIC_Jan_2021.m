%This code creates a box full of randomised polydisperse particles of a set size distribution and also creates a
%top and bottom plate to shear the box of particles. Outputs a file to use
%with MFiX DEM as an initial condition.

clear all
close all
clc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%      NEW Sept 2019     %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%    Top and Bottom plate stuff   %%%%

dr1 = 0.00133;   % Bottom plate particle diameter.
Hr1 = 0.01;     % Bottom plate height.
rho1 = 1050;     % Bottom plate particle density.
delta1 = 0.5*dr1;  %Bottom plate position randomness.
gap1 = 0.1*dr1;   % Bottom plate X spacing between particles.

dr2 = 0.00133;   % Top plate particle diameter.
Hr2 = 0.05;     % Top plate height.
rho2 = 1050;     % Top plate particle density.
delta2 = 0.5*dr2;  %Top plate position randomness.
gap2 = 0.1*dr1;   % Top plate X spacing between particles.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

xmin=-0.04;   %left limit on x domain
xmax=0.04;   %right limit on x domain
ymin=-0.04;   %lower limit on y doman
ymax=0.04;   %upper limit on y doman
zmin=0;   %front limit on z doman
zmax=0.4;   %back limit on z doman

u=0;  %Initial x velocity of particles.
v=0;   %Initial y velocity of particles.
w=0;    %Initial z velocity of particles.

NumberOf_x_Splits=5;    %How many splits to make in x direction  make sure dimension/split=3xlargest particle
NumberOf_y_Splits=5;    %How many splits to make in y direction
NumberOf_z_Splits=20;    %How many splits to make in z direction
       

Total_Number_Of_Splits=NumberOf_x_Splits*NumberOf_z_Splits*NumberOf_y_Splits;
Number_Of_Random_blocks=10;    %Increase this to increase randomisation, but also increase processing time.

ParticleDiameters=[0.005 0.004 0.003 0.00229 0.00133 0.0005]; %large to small particles left to right   %enter here all particle diameters. Number of entries and order must match number fraction vector;
NumberFraction=[0.02 0.03 0.05 0.1 0.7 0.1];  %enter here the number fraction of each particle size. Number of entries must = number of enteries in particle diameter vector. Sum of entries must = 1
Densities=[1050 1050 1050 1050 1050 1050];    %enter here the density of each particle size in the same order as above
CumulativeNumberFraction=cumsum(NumberFraction);

ApproxTotalSolidsFraction=0.3;  %Set this low to enable the code to converge easily. Then use DEM to settle the particles.

%-----------------------------------------------------------------------
%-----------------------------------------------------------------------
%-----------------------------------------------------------------------
%-----------------------------------------------------------------------
%-----------------------------------------------------------------------
%-----------------------------------------------------------------------

dx=(xmax-xmin)./NumberOf_x_Splits;
dy=(ymax-ymin)./NumberOf_y_Splits;
dz=(zmax-zmin)./NumberOf_z_Splits;


TestTestTest=min(((xmax-xmin)/NumberOf_x_Splits/max(ParticleDiameters)),((zmax-zmin)/NumberOf_z_Splits/max(ParticleDiameters)));
TestTestTest=min(TestTestTest,((ymax-ymin)/NumberOf_y_Splits/max(ParticleDiameters)));

if TestTestTest<3.01
   warning('Likely you have too many splits for the size of the particles. Please consider reducing NumberOf_x/y/z_Splits')
    pause(10)
end


% if sum(NumberFraction)~=1
if abs(sum(NumberFraction)-1)>1e-13
    error('Number fraction must sum to 1')
    return
end
if length(ParticleDiameters)~=length(NumberFraction)
    error('Number of enteries in NumberFraction vector must equal number of entries in ParticleDiameters vector')
    return
end
if length(ParticleDiameters)~=length(Densities)
    error('Number of enteries in Densities vector must equal number of entries in ParticleDiameters vector')
    return
end


x_lim=xmin + (xmax-xmin)./NumberOf_x_Splits;
z_lim=zmin + (zmax-zmin)./NumberOf_z_Splits;
y_lim=ymin + (ymax-ymin)./NumberOf_y_Splits;
SplitVolume=(x_lim-xmin)*(z_lim-zmin)*(y_lim-ymin);

N_split_particles=round(ApproxTotalSolidsFraction.*SplitVolume./(4./3.*pi.*sum(NumberFraction.*((ParticleDiameters./2).^3))));

NumberOfEachParticle=round(N_split_particles.*NumberFraction)
ParticlesVolumes=sum((4/3.*pi.*(ParticleDiameters./2).^3).*NumberOfEachParticle)

tic
for jjj=1:Number_Of_Random_blocks
    
    Npos=1;

    for iii=1:length(NumberFraction)
        for abc=1:NumberOfEachParticle(iii)
            clear diam TestPoint Distances Errortest whilemarker
            
            whilemarker=0;
            while whilemarker==0
                clear diam TestPoint Distances Errortest TestBit wmark2
                wmark2=0;
                diam=ParticleDiameters(iii);
                while wmark2==0
                    clear TestPoint xmoo zmoo ymoo

                    xmoo=((x_lim-diam/2)-(xmin+diam/2)).*rand(1) + (xmin+diam/2);
                    zmoo=((z_lim-diam/2)-(zmin+diam/2)).*rand(1) + (zmin+diam/2);
                    ymoo=((y_lim-diam/2)-(ymin+diam/2)).*rand(1) + (ymin+diam/2);

                        Testpoint=[xmoo,ymoo,zmoo,diam,Densities(iii),u,v,w];
                    if Testpoint(1)>(xmin+diam/2) && Testpoint(1)<(x_lim-diam/2) && Testpoint(3)>(zmin+diam/2) && Testpoint(3)<(z_lim-diam/2) && Testpoint(2)>(ymin+diam/2) && Testpoint(2)<(y_lim-diam/2)
                        wmark2=1;
                    end
                end

                if iii==1 && abc==1
                   ParticlePosition(Npos,:,jjj)=Testpoint; 
                   Npos=Npos+1;
                    whilemarker=1;
                else
                    Distances=sqrt((ParticlePosition(:,1,jjj)-Testpoint(1)).^2 + (ParticlePosition(:,2,jjj)-Testpoint(2)).^2 +(ParticlePosition(:,3,jjj)-Testpoint(3)).^2);
                    TestBit=(Distances-(diam/2+ParticlePosition(:,4,jjj)./2));
                    if isempty(find(TestBit<0))
                    ParticlePosition(Npos,:,jjj)=Testpoint;
                    Npos=Npos+1;
                    [Number_Of_Random_blocks-jjj sum(NumberOfEachParticle)-Npos Npos]
                    whilemarker=1;
                    end
                end
            end
                 
        end
    end



end
toc

display('Please wait...')

BlankSlate=[];
save ('particle_input.dat','BlankSlate','-ASCII')



clear iii jjj abc

for iii=1:NumberOf_x_Splits
    for jjj=1:NumberOf_y_Splits
        for abc=1:NumberOf_z_Splits
            [NumberOf_z_Splits-abc NumberOf_y_Splits-jjj NumberOf_x_Splits-iii]
            clear RRR PF xadd yadd zadd
            RRR=randi(Number_Of_Random_blocks,1);
            PF=ParticlePosition(:,:,RRR);
            xadd=(iii-1)*dx;
            yadd=(jjj-1)*dy;
            zadd=(abc-1)*dz;
            PF(:,1)=PF(:,1)+xadd;
            PF(:,2)=PF(:,2)+yadd;
            PF(:,3)=PF(:,3)+zadd;
            PF(:,4)=PF(:,4)./2;  %radius
            save ('particle_input.dat','PF','-ASCII','-append')
        end
    end
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%      NEW Sept 2019     %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%    Top and Bottom plate stuff   %%%%
%%% Bottom plate
x_bot_Left = xmin + (gap1./2 + dr1./2);x_bot_right = xmax - (gap1./2 + dr1./2);
x_bot = [x_bot_Left:gap1:x_bot_right-gap1 x_bot_right];
z_bot_Left = zmin + (gap1./2 + dr1./2);z_bot_right = zmax - (gap1./2 + dr1./2);
z_bot = [z_bot_Left:gap1:z_bot_right-gap1 z_bot_right];
[X_BOT,Z_BOT]=meshgrid(x_bot,z_bot);
Y_BOT = Hr1 + ((rand(size(X_BOT))-0.5).*2).*delta1;   %%% NOTE: the ((rand(size(X_BOT))-0.5).*2)  part gives us a uniform random variable between -1 and 1. We muliply this by delta1 to give the randomness (defined at beginning of file)


%%% Top plate
x_top_Left = xmin + (gap2./2 + dr2./2);x_top_right = xmax - (gap2./2 + dr2./2);
x_top = [x_top_Left:gap2:x_top_right-gap2 x_top_right];
z_top_Left = zmin + (gap2./2 + dr2./2);z_top_right = zmax - (gap2./2 + dr2./2);
z_top = [z_top_Left:gap2:z_top_right-gap2 z_top_right];
[X_TOP,Z_TOP]=meshgrid(x_top,z_top);
Y_TOP = Hr2 + ((rand(size(X_TOP))-0.5).*2).*delta2;   %%% NOTE: the ((rand(size(X_BOT))-0.5).*2)  part gives us a uniform random variable between -1 and 1. We muliply this by delta1 to give the randomness (defined at beginning of file)


BottomPlate=[X_BOT(:) Y_BOT(:) Z_BOT(:) (ones(size(X_BOT(:)))).*dr1./2   (ones(size(X_BOT(:)))).*rho1./2      (ones(size(X_BOT(:)))).*u     (ones(size(X_BOT(:)))).*v     (ones(size(X_BOT(:)))).*w];
TopPlate=[X_TOP(:) Y_TOP(:) Z_TOP(:) (ones(size(X_TOP(:)))).*dr2./2   (ones(size(X_TOP(:)))).*rho2./2      (ones(size(X_BOT(:)))).*u     (ones(size(X_BOT(:)))).*v     (ones(size(X_BOT(:)))).*w];

display('Please wait, saving particles...')


save ('particle_input.dat','BottomPlate','-ASCII','-append')
save ('particle_input.dat','TopPlate','-ASCII','-append')


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
display('Please wait, plotting example...')
pause(1)

for bbb=1:1
figure
axes('NextPlot', 'Add');
CMap = parula(4);
A = ParticlePosition(:,:,bbb);
A(:,4)=A(:,4)./2;
[X, Y, Z] = sphere;
for k = 1:sum(NumberOfEachParticle)
  XX = X * A(k, 4) + A(k, 1);
  YY = Y * A(k, 4) + A(k, 2);
  ZZ = Z * A(k, 4) + A(k, 3);
  surface(XX, YY, ZZ, 'FaceAlpha', 0.5, 'FaceColor', CMap(4, :), ...
          'EdgeColor', [0.4, 0.4, 0.4]);
end
axis equal
end

display('FINISHED!')





