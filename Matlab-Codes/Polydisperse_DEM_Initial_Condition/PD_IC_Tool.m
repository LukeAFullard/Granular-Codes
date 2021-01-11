%This code creates a box full of randomised polydisperse particles of a set size distribution and also creates a
%top and bottom plate to shear the box of particles. Outputs a file to use
%with MFiX DEM as an initial condition.

clear all
close all
clc


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%      NEW Sept 2019     %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% Use with CAUTION %%%
%%% Use with CAUTION %%%
%%% Use with CAUTION %%%
plot_all=0;  %Set = 1 if you want to plot all particles in domain after running code. 
%Note, this can be very many particles and will crash your computer if using too many.
%IF IN DOUBT, set = 0   %%% Use with CAUTION %%%
%%% Use with CAUTION %%%
%%% Use with CAUTION %%%
%%% Use with CAUTION %%%

%%%%    Top and Bottom plate stuff   %%%%
Large_Particles_Regular_Grid=1;  %Set this = 1 for a regular grid of particles. Otherwise will be random

%%%%   Usually the mesh options below will be needed when you offset the large particles relative to the small ones.
BottomMeshGrid=1;   %Set this = 1 if you want a wall below bottom plate to ensure particles do not get lost.
TopMeshGrid=1; %Set this = 1 if you want a wall above top plate to ensure particles do not get lost.
%--------------------------------------------------------

dr1 = 0.7./1000;   % Bottom plate small particle diameter.
% Hr1 = dr1;     % Bottom plate height.
rho1 = 2501;     % Bottom plate particle density.
delta1 = 0.5*dr1;  %Bottom plate position small particle randomness.
gap1 = 0.01*dr1;   % Bottom plate X spacing between small particles. NOTE: This gap is in addition to the particle diameter, so if  gap1=0 then particles will be spaced exactly touching
        dr1_L = 11./1000;   %Bottom plate large particle diameter
        delta1_L = 0.1*dr1_L; %Bottom plate position large particle randomness.
        N1_L = 60;   %Number of randomly placed large particles (For random placement only)
        Large_Gap_1 = dr1*1.75;   %Min gap between large particles
        OS1_L = 0.0*dr1_L;     %Offset for large particles. Positive offset pushes particles down away from the main flow
            % 24th April contol large particle placement
            OSL1_left = dr1_L./2 + dr1_L./8;  %Controls how close large particles can be placed to left side of domain for Bottom plate 
            OSL1_right = dr1_L./2 + dr1_L./8;  %Controls how close large particles can be placed to right side of domain for Bottom plate 
            OSL1_top = dr1_L./2 + dr1_L./8;  %Controls how close large particles can be placed to top side of domain for Bottom plate 
            OSL1_bottom = dr1_L./2 + dr1_L./8;  %Controls how close large particles can be placed to bottom side of domain for Bottom plate 
            N1_L_x = 7;    %Controls how many large particles in the x-direction, bottom plate
            N1_L_z = 7;    %Controls how many large particles in the z-direction, bottom plate

Hr1 = -dr1_L*0.525 - delta1_L ;     % Bottom plate height.
plot_Bottom=0;   %Set = 1 to plot the bottom plate for checking

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

dr2 = 0.7./1000;   % Top plate small particle diameter
rho2 = 2502;     % Top plate particle density.
delta2 = 0.5*dr2;  %Top plate position small particle randomness.
gap2 = 0.01*dr1;   % Top plate X spacing between small particles.
        dr2_L = 11./1000;   %Top plate large particle diameter
        delta2_L = 0.1*dr2_L; %Top plate position large particle randomness.
        N2_L = 49;   %Number of randomly placed large particles
        Large_Gap_2 = dr2*1.75;   %Min gap between large particles
        OS2_L = 0.0*dr2_L;     %Offset for large particles. Poitive offset pushes particles up away from the main flow
                        % 24th April contol large particle placement
                        OSL2_left = dr2_L./2 + dr2_L./8;  %Controls how close large particles can be placed to left side of domain for top plate 
                        OSL2_right = dr2_L./2 + dr2_L./8;  %Controls how close large particles can be placed to right side of domain for top plate 
                        OSL2_top = dr2_L./2 + dr2_L./8;  %Controls how close large particles can be placed to top side of domain for top plate 
                        OSL2_bottom = dr2_L./2 + dr2_L./8;  %Controls how close large particles can be placed to bottom side of domain for top plate 
                        N2_L_x = 7;    %Controls how many large particles in the x-direction, top plate
                        N2_L_z = 7;    %Controls how many large particles in the z-direction, top plate
Hr2 = 0.12+dr2_L*0.525 + delta2_L;     % Top plate height.
plot_Top=1;   %Set = 1 to plot the top plate for checking
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

xmin=0.0; %+ 2.*dr1;   %left limit on x domain
xmax=0.1;   %right limit on x domain
ymin=dr1*2;   %lower limit on z domain
ymax=0.12;   %uppr limit on z domain
zmin=0;   %front limit on z domain
zmax=0.1;   %back limit on z domain


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%--- 24th April 2020 large particle checks
if Large_Particles_Regular_Grid==1
    N1_L = N1_L_x.*N1_L_z;
    N2_L = N2_L_x.*N2_L_z;
    
     XL_test = xmin  +  OSL1_left;
    XR_test = xmax  -  OSL1_right;
    ZL_test = zmin  +  OSL1_top;
    ZR_test = zmax  -  OSL1_bottom;
%     x_try_test = [XL_test: (XR_test-XL_test)./(N1_L_x-1) :XR_test]
%     z_try_test = [ZL_test: (ZR_test - ZL_test)./(N1_L_z-1) :ZR_test]
    if ((XR_test-XL_test)./(N1_L_x-1)) < dr1_L + Large_Gap_1
        error('You are using too many large particles in the x-directiob of the bottom plate. Please reduce the number')
    end
    if ((ZR_test-ZL_test)./(N1_L_z-1)) < dr1_L + Large_Gap_1
        error('You are using too many large particles in the z-directiob of the in the bottom plate. Please reduce the number')
    end
    
    XL_test = xmin  +  OSL2_left;
    XR_test = xmax  -  OSL2_right;
    ZL_test = zmin  +  OSL2_top;
    ZR_test = zmax  -  OSL2_bottom;
    
    if ((XR_test-XL_test)./(N2_L_x-1)) < dr2_L + Large_Gap_2
        error('You are using too many large particles in the x-directiob of the in the top plate. Please reduce the number')
    end
    if ((ZR_test-ZL_test)./(N2_L_z-1)) < dr2_L + Large_Gap_2
        error('You are using too many large particles in the z-directiob of the in the top plate. Please reduce the number')
    end
    
end
% pause

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




u=0;  %Initial x velocity
v=0;    %Initial y velocity
w=0;       %Initial z velocity

NumberOf_x_Splits=3;    %How many splits to make in x direction  make sure dimension/split=3xlargest particle
NumberOf_y_Splits=3;    %How many splits to make in y direction
NumberOf_z_Splits=3;    %How many splits to make in z direction



Total_Number_Of_Splits=NumberOf_x_Splits*NumberOf_z_Splits*NumberOf_y_Splits;
Number_Of_Random_blocks=5;    %Increase this to increase randomisation, but also increase processing time.

ParticleDiameters=[11.0 10.0 9.0 1.1 1.0 0.9]./1000; %large to small particles left to right   %enter here all particle diameters. Number of entries and order must match number fraction vector;
NumberFraction=[0.0003 0.0006 0.0003 0.2497 0.4994 0.2497];  %enter here the number fraction of each particle size. Number of entries must = number of enteries in particle diameter vector. Sum of entries must = 1
Densities=[2500 2500 2500 2500 2500 2500];    %enter here the density of each particle size in the same order as above
CumulativeNumberFraction=cumsum(NumberFraction);

ApproxTotalSolidsFraction=0.30;  %Set this low to enable the code to converge easily. Then use DEM to settle the particles.


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


if sum(NumberFraction)~=1
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
% tic
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
%                     Testpoint=[(x_lim-xmin).*rand(1) + xmin,(z_lim-zmin).*rand(1) + zmin,(y_lim-ymin).*rand(1) + ymin,diam];
                    xmoo=((x_lim-diam/2)-(xmin+diam/2)).*rand(1) + (xmin+diam/2);
                    zmoo=((z_lim-diam/2)-(zmin+diam/2)).*rand(1) + (zmin+diam/2);
                    ymoo=((y_lim-diam/2)-(ymin+diam/2)).*rand(1) + (ymin+diam/2);
%                     Testpoint=[xmoo,zmoo,ymoo,diam,Densities(iii),u,v,w];
                        Testpoint=[xmoo,ymoo,zmoo,diam,Densities(iii),u,v,w];
                    if Testpoint(1)>(xmin+diam/2) && Testpoint(1)<(x_lim-diam/2) && Testpoint(3)>(zmin+diam/2) && Testpoint(3)<(z_lim-diam/2) && Testpoint(2)>(ymin+diam/2) && Testpoint(2)<(y_lim-diam/2)
                        wmark2=1;
                    end
                end
                %NEED TO CHECK HERE THAT IT IS IN THE WALL BOUNDS TOO
        %         Switch=1;
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
%             pause
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




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%      NEW 20th Nov 2019     %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear numofeach
ABCDE=load('particle_input.dat');
for ooo=1:numel(ParticleDiameters)
   numofeach(ooo)=numel(find(ABCDE(:,4).*2==ParticleDiameters(ooo)));
end


D43 = sum(numofeach.*((ParticleDiameters).^4))./(sum(numofeach.*((ParticleDiameters).^3)))
D32 = sum(numofeach.*((ParticleDiameters).^3))./(sum(numofeach.*((ParticleDiameters).^2)))


 save ('average_particle_diam.dat','D43','D32','ParticleDiameters','numofeach','-ASCII')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%








%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%      NEW April 2020     %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%    Top and Bottom plate stuff   %%%%
%%% Bottom plate
x_bot_Left = xmin + (gap1./8 + dr1./2);x_bot_right = xmax - (gap1./8 + dr1./2);
x_bot = [x_bot_Left:(gap1+dr1):x_bot_right-gap1./2];
z_bot_Left = zmin + (gap1./8 + dr1./2);z_bot_right = zmax - (gap1./8 + dr1./2);
z_bot = [z_bot_Left:(gap1+dr1):z_bot_right-gap1./2];
[X_BOT,Z_BOT]=meshgrid(x_bot,z_bot);
Y_BOT = Hr1 + ((rand(size(X_BOT))-0.5).*2).*delta1;   %%% NOTE: the ((rand(size(X_BOT))-0.5).*2)  part gives us a uniform random variable between -1 and 1. We muliply this by delta1 to give the randomness (defined at beginning of file)
%%% Bottom plate large particles

X_BOT_Large = [];
Y_BOT_Large = [];
Z_BOT_Large = [];

if Large_Particles_Regular_Grid==1
    clear x_try y_try z_try X_TRY Z_TRY XL XR ZL ZR

    
    XL = xmin  +  OSL1_left;
    XR = xmax  -  OSL1_right;
    ZL = zmin  +  OSL1_top;
    ZR = zmax  -  OSL1_bottom;

    x_try = [XL: (XR-XL)./(N1_L_x-1) :XR];
    z_try = [ZL: (ZR - ZL)./(N1_L_z-1) :ZR];

    
    [X_TRY,Z_TRY]=meshgrid(x_try,z_try);
    Y_TRY = Hr1 + ((rand(size(X_TRY))-0.5).*2).*delta1_L - OS1_L;
    X_BOT_Large=X_TRY(:);
    Z_BOT_Large=Z_TRY(:);
    Y_BOT_Large=Y_TRY(:);
else
        for iiii=1:N1_L
            dist=0;
            while dist < dr1_L+Large_Gap_1
                clear x_try y_try z_try dist
                x_try = (xmin+dr1_L./2 + gap1./2) + ((xmax-dr1_L./2 - gap1./2)-(xmin+dr1_L./2  + gap1./2)).*rand(1,1);
                z_try = (zmin+dr1_L./2 + gap1./2) + ((zmax-dr1_L./2 - gap1./2)-(zmin+dr1_L./2  + gap1./2)).*rand(1,1);
                y_try = Hr1 + ((rand(1,1)-0.5).*2).*delta1_L - OS1_L;

                if iiii==1
                    dist=100000;
                else
                    dist = min(sqrt((x_try - X_BOT_Large).^2 + (y_try - Y_BOT_Large).^2 + (z_try - Z_BOT_Large).^2));
                end
            end
        %     iiii
            X_BOT_Large(iiii)=x_try;
            Y_BOT_Large(iiii)=y_try;
            Z_BOT_Large(iiii)=z_try;

        end
end
%%% Bottom plate erase overlapping small particles
X_BOT_FIN=X_BOT(:);
Y_BOT_FIN=Y_BOT(:);
Z_BOT_FIN=Z_BOT(:);
for iiii=1:length(X_BOT_Large)
    clear Dists Dist_Locs ExtraDist epsilony Delete_Locs
    Dists = sqrt((X_BOT_Large(iiii) - X_BOT_FIN).^2 + (Y_BOT_Large(iiii) - Y_BOT_FIN).^2 + (Z_BOT_Large(iiii) - Z_BOT_FIN).^2);
    
    Delete_Locs = find(Dists < ((dr1_L - 1.*dr1)./2 - 0.025*dr1));
    
    Dist_Locs = find(Dists < ((dr1_L + dr1)./2 + 0.*dr1));

    ExtraDist = 1.00025.*abs(Dists(Dist_Locs) - ((dr1_L + dr1)./2 + 0.*dr1));

    epsilony = -sqrt((Dists(Dist_Locs)+ExtraDist).^2 - (X_BOT_Large(iiii) - X_BOT_FIN(Dist_Locs)).^2 - (Z_BOT_Large(iiii) - Z_BOT_FIN(Dist_Locs)).^2) + Y_BOT_Large(iiii);

      
      
    Y_BOT_FIN(Dist_Locs) = epsilony; 
    clear Dists
    Dists = sqrt((X_BOT_Large(iiii) - X_BOT_FIN).^2 + (Y_BOT_Large(iiii) - Y_BOT_FIN).^2 + (Z_BOT_Large(iiii) - Z_BOT_FIN).^2);
    if isempty(find(Dists < ((dr1_L + dr1)./2 + 0.*dr1)))
        
    else
        iiii
        error('Issue with particles overlapping')
    end
      X_BOT_FIN(Delete_Locs) = [];
    Y_BOT_FIN(Delete_Locs) = [];  
    Z_BOT_FIN(Delete_Locs) = [];
end

for iiii=1:length(X_BOT_Large)
    X_BOT_FIN(end+1) = X_BOT_Large(iiii);
    Y_BOT_FIN(end+1) = Y_BOT_Large(iiii);
    Z_BOT_FIN(end+1) = Z_BOT_Large(iiii);
end




%------------------------
%----------------------------------------
%-----------------------------------------------------

%%% Top plate
x_top_Left = xmin + (gap2./8 + dr2./2);x_top_right = xmax - (gap2./8 + dr2./2);
x_top = [x_top_Left:(gap2+dr2):x_top_right-gap2./2];
z_top_Left = zmin + (gap2./8 + dr2./2);z_top_right = zmax - (gap2./8 + dr2./2);
z_top = [z_top_Left:(gap2+dr2):z_top_right-gap2./2 ];
[X_TOP,Z_TOP]=meshgrid(x_top,z_top);
Y_TOP = Hr2 + ((rand(size(X_TOP))-0.5).*2).*delta2;   %%% NOTE: the ((rand(size(X_BOT))-0.5).*2)  part gives us a uniform random variable between -1 and 1. We muliply this by delta1 to give the randomness (defined at beginning of file)
%Top particles large
X_TOP_Large = [];
Y_TOP_Large = [];
Z_TOP_Large = [];

if Large_Particles_Regular_Grid==1
    clear x_try y_try z_try X_TRY Z_TRY XL XR ZL ZR


    XL = xmin  +  OSL2_left;
    XR = xmax  -  OSL2_right;
    ZL = zmin  +  OSL2_top;
    ZR = zmax  -  OSL2_bottom;

    x_try = [XL: (XR-XL)./(N2_L_x-1) :XR];
    z_try = [ZL: (ZR - ZL)./(N2_L_z-1) :ZR];

    

    [X_TRY,Z_TRY]=meshgrid(x_try,z_try);
    Y_TRY = Hr2 + ((rand(size(X_TRY))-0.5).*2).*delta2_L + OS2_L;
    X_TOP_Large=X_TRY(:);
    Z_TOP_Large=Z_TRY(:);
    Y_TOP_Large=Y_TRY(:);
else
        for iiii=1:N2_L
            dist=0;
            while dist < dr2_L+Large_Gap_2
                clear x_try y_try z_try dist
                x_try = (xmin+dr2_L./2) + ((xmax-dr2_L./2)-(xmin+dr2_L./2)).*rand(1,1);
                z_try = (zmin+dr2_L./2) + ((zmax-dr2_L./2)-(zmin+dr2_L./2)).*rand(1,1);
                y_try = Hr2 + ((rand(1,1)-0.5).*2).*delta2_L + OS2_L;

                if iiii==1
                    dist=100000;
                else
                    dist = min(sqrt((x_try - X_TOP_Large).^2 + (y_try - Y_TOP_Large).^2 + (z_try - Z_TOP_Large).^2));
                end
            end
        %     iiii
            X_TOP_Large(iiii)=x_try;
            Y_TOP_Large(iiii)=y_try;
            Z_TOP_Large(iiii)=z_try;

        end
end


%%% top plate erase overlapping small particles
X_TOP_FIN=X_TOP(:);
Y_TOP_FIN=Y_TOP(:);
Z_TOP_FIN=Z_TOP(:);
for iiii=1:length(X_TOP_Large)
    



    clear Dists Dist_Locs ExtraDist epsilony Delete_Locs
    Dists = sqrt((X_TOP_Large(iiii) - X_TOP_FIN).^2 + (Y_TOP_Large(iiii) - Y_TOP_FIN).^2 + (Z_TOP_Large(iiii) - Z_TOP_FIN).^2);
    
    Delete_Locs = find(Dists < ((dr2_L - 1.0*dr2)./2 - 0.025*dr2));
    Dist_Locs = find(Dists < ((dr2_L + dr2)./2));
    ExtraDist = 1.00025.*abs(Dists(Dist_Locs) - ((dr2_L + dr2)./2 + 0.*dr2));
    epsilony = -sqrt((Dists(Dist_Locs)+ExtraDist).^2 - (X_TOP_Large(iiii) - X_TOP_FIN(Dist_Locs)).^2 - (Z_TOP_Large(iiii) - Z_TOP_FIN(Dist_Locs)).^2) + Y_TOP_Large(iiii);

    
     Y_TOP_FIN(Dist_Locs) = Y_TOP_Large(iiii) + abs(Y_TOP_Large(iiii) - epsilony);

    clear Dists
    Dists = sqrt((X_TOP_Large(iiii) - X_TOP_FIN).^2 + (Y_TOP_Large(iiii) - Y_TOP_FIN).^2 + (Z_TOP_Large(iiii) - Z_TOP_FIN).^2);
        if isempty(find(Dists < ((dr2_L + dr2)./2 + 0.*dr2)))
        
        else
            iiii
            error('Issue with particles overlapping')
        end
%         pause
    X_TOP_FIN(Delete_Locs) = [];
    Y_TOP_FIN(Delete_Locs) = [];  
    Z_TOP_FIN(Delete_Locs) = [];

end

for iiii=1:length(X_TOP_Large)
    X_TOP_FIN(end+1) = X_TOP_Large(iiii);
    Y_TOP_FIN(end+1) = Y_TOP_Large(iiii);
    Z_TOP_FIN(end+1) = Z_TOP_Large(iiii);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


BottomPlate=[X_BOT_FIN(:) Y_BOT_FIN(:) Z_BOT_FIN(:) (ones(size(X_BOT_FIN(:)))).*dr1./2   (ones(size(X_BOT_FIN(:)))).*rho1      (ones(size(X_BOT_FIN(:)))).*u     (ones(size(X_BOT_FIN(:)))).*v     (ones(size(X_BOT_FIN(:)))).*w];
for iiii=1:length(X_BOT_Large)
    BottomPlate(end - N1_L + iiii,:)=[X_BOT_FIN(end - N1_L + iiii) Y_BOT_FIN(end - N1_L + iiii) Z_BOT_FIN(end - N1_L + iiii) dr1_L./2   rho1      u     v     w];
end




   TopPlate=[X_TOP_FIN(:) Y_TOP_FIN(:) Z_TOP_FIN(:) (ones(size(X_TOP_FIN(:)))).*dr2./2   (ones(size(X_TOP_FIN(:)))).*rho2      (ones(size(X_TOP_FIN(:)))).*u     (ones(size(X_TOP_FIN(:)))).*v     (ones(size(X_TOP_FIN(:)))).*w];

for iiii=1:length(X_TOP_Large)
        TopPlate(end - N2_L + iiii,:)=[X_TOP_FIN(end - N2_L + iiii) Y_TOP_FIN(end - N2_L + iiii) Z_TOP_FIN(end - N2_L + iiii) dr2_L./2   rho2      u     v     w];
end
disp('write rough plates')
save ('particle_input.dat','TopPlate','-ASCII','-append')
save ('particle_input.dat','BottomPlate','-ASCII','-append')


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


if plot_all==1
        load('particle_input.dat')
        for bbb=1:1
        figure
        axes('NextPlot', 'Add');
        CMap = parula(4);
        A = particle_input;
        % A(:,4)=A(:,4)./2;
        [X, Y, Z] = sphere;
        for k = 1:length(A)
          XX = X * A(k, 4) + A(k, 1);
          YY = Y * A(k, 4) + A(k, 2);
          ZZ = Z * A(k, 4) + A(k, 3);
          surface(XX, YY, ZZ, 'FaceAlpha', 0.5, 'FaceColor', CMap(4, :), ...
                  'EdgeColor', [0.4, 0.4, 0.4]);
        end
        axis equal
        end
end



if plot_Bottom==1
    load('particle_input.dat')
    for bbb=1:1
    figure
    axes('NextPlot', 'Add');
    CMap = parula(4);
    A = particle_input(end-length(BottomPlate)+1:end,:);
    % A(:,4)=A(:,4)./2;
    [X, Y, Z] = sphere;
    for k = 1:length(A)
      XX = X * A(k, 4) + A(k, 1);
      YY = Y * A(k, 4) + A(k, 2);
      ZZ = Z * A(k, 4) + A(k, 3);
      surface(XX, YY, ZZ, 'FaceAlpha', 0.5, 'FaceColor', CMap(4, :), ...
              'EdgeColor', [0.4, 0.4, 0.4]);
    end
    view(0,180)
    axis image
    end
end


if plot_Top==1
    load('particle_input.dat')
    for bbb=1:1
    figure
    axes('NextPlot', 'Add');
    CMap = parula(4);
    A = particle_input(end-length(BottomPlate)-length(TopPlate)+1:end-length(BottomPlate),:);
    % A(:,4)=A(:,4)./2;
    [X, Y, Z] = sphere;
    for k = 1:length(A)
      XX = X * A(k, 4) + A(k, 1);
      YY = Y * A(k, 4) + A(k, 2);
      ZZ = Z * A(k, 4) + A(k, 3);
      surface(XX, YY, ZZ, 'FaceAlpha', 0.5, 'FaceColor', CMap(4, :), ...
              'EdgeColor', [0.4, 0.4, 0.4]);
    end
    view(0,180)
    axis image
    end
end


disp('FINISHED')