%Tissueparameters;
rho=1138; %density kg/m^3
f=1.36*1e6; %frequency
c=1569;    %ultrasound speed
Kn=2*pi*f/c; %wave number
Alpha=12.24;%attenuation coefficient(0.9Np/MHz*m)
k=0.52;     %coefficient of thermal conductivity
Cp=3370;    %specific heat capacity

A=0.03; %tissue depth

%Elementparameters;
SphericalR=0.15;                            %radius of Sphericality
ElementR=0.005;                             %radius of each element
NumEachRing=[20,24,32,36];                  %number of elements of each ring
AngleZ=[15.4,19.75,24.20,28.80]*pi/180;       %angle between each ring and Z axis
AngleXFirst=pi./NumEachRing;                %angle between first element of each ring and X axis
RingNum=4;                                  %number of rings
ElementNum=sum(NumEachRing);                %total number of elements

AngleXelement=[];                           %angle between each element and X axis
AngleZelement=[];                           %angle between each element and Z axis
for i=1:RingNum
    AngleXelement=[AngleXelement AngleXFirst(i):2*pi/NumEachRing(i):2*pi-AngleXFirst(i)];
    AngleZelement=[AngleZelement AngleZ(i)*ones(1,NumEachRing(i))];
end

%Focus
FocusX=0;
FocusY=0;
FocusZ=0;

%simulation area and step
Xmax = 0.015;
Xstep = 0.001;

Ymax = 0.015;
Ystep = 0.001;

Zmax = 0.015; 
Zstep = 0.001;

h=Xstep;% length between adjacent points

%Discretized point source

ElementRNum=10;  %consider each element as 10*10 point sources
ElementAngleNum=10; 

DiscretizedR=ElementR/ElementRNum-ElementR/ElementRNum/2:ElementR/ElementRNum:ElementR-ElementR/ElementRNum/2;      %Radius and angle of each point source
DiscretizedAngle=2*pi/ElementAngleNum:2*pi/ElementAngleNum:2*pi;

CenterX=DiscretizedR'*cos(DiscretizedAngle);         %x,y,z coordinates of the point sources
CenterY=DiscretizedR'*sin(DiscretizedAngle);
CenterZ=zeros(ElementRNum,ElementAngleNum);

OuterR=ElementR/ElementRNum:ElementR/ElementRNum:ElementR;        %outer radius of each point source
InnerR=0:ElementR/ElementRNum:ElementR-ElementR/ElementRNum;      %inner radius of each point source
DS=(OuterR.^2-InnerR.^2)*pi/ElementAngleNum;                      %area of each point source
DS=repmat(DS',1,ElementAngleNum);

Xrow=reshape(CenterX,1,ElementRNum*ElementAngleNum);
Yrow=reshape(CenterY,1,ElementRNum*ElementAngleNum);
Zrow=reshape(CenterZ,1,ElementRNum*ElementAngleNum);
XYZrows=[Xrow;Yrow;Zrow];              %Discritized coordinates of one element


%coordinate transformation
for element=1:ElementNum
    SphericalX=SphericalR*sin(AngleZelement(element))*cos(AngleXelement(element));
    SphericalY=SphericalR*sin(AngleZelement(element))*sin(AngleXelement(element));
    SphericalZ=SphericalR*cos(AngleZelement(element));
    SphericalP=[SphericalX;SphericalY;SphericalZ];                  %coordinate of each element on hemispherical
    
    alpha=acos(SphericalZ/sqrt(SphericalX*SphericalX+SphericalY*SphericalY+SphericalZ*SphericalZ));               %angle between Z-axis and radius
    ProjectionR=sqrt(SphericalX*SphericalX+SphericalY*SphericalY);                                               %projection radius of each element
    beta=acos(SphericalX/ProjectionR);                                                                           %angle between radius and X-axis on the projection plane
    if SphericalY<0
        beta=2*pi-beta;
    end
    P1=[cos(alpha), 0, sin(alpha);              %align Z-axis of both coordinates
        0,        1,    0;
     -sin(alpha), 0, cos(alpha)];
 
    P2=[cos(beta),-sin(beta),0;                 %align X-axis of both coordinates
      sin(beta),cos(beta), 0;
        0,      0,       1];
    P0=P2*P1*XYZrows;                         %coordinates of discretized point source of each element
    SphericalP=repmat(SphericalP,1,size(P0,2));
    Pnew=SphericalP+P0;
    
    X(element,:) = Pnew(1,:);
    Y(element,:) = Pnew(2,:);
    Z(element,:) = Pnew(3,:);  
end

%Calculate Forward transfer matrix
U=ones(ElementNum,1);                     %unit normal particle velocity
Stepoints = 31;
[FieldX,FieldY,FieldZ] = meshgrid(-Xmax:Xstep:Xmax,-Ymax:Ystep:Ymax,-Zmax:Zstep:Zmax);
[d1,d2,d3] = size(FieldX);

PointNum = d1*d2*d3;            %all points for simulation
NumOStep = PointNum/Stepoints;

FieldX = reshape(FieldX, PointNum,1);
FieldY = reshape(FieldY, PointNum,1);
FieldZ = reshape(FieldZ, PointNum,1);

DS = repmat(reshape(DS,1,size(DS,1)*size(DS,2)),ElementNum*Stepoints,1);
P = [];

for i = 1:NumOStep
 
    FieldXi = FieldX(Stepoints*(i-1)+1:Stepoints*i);
    FieldYi = FieldY(Stepoints*(i-1)+1:Stepoints*i);
    FieldZi = FieldZ(Stepoints*(i-1)+1:Stepoints*i);
    
    TransducerX = repmat(X,Stepoints,1);
    TransducerY = repmat(Y,Stepoints,1);
    TransducerZ = repmat(Z,Stepoints,1);
    
    SoundFieldX = repmat(FieldXi',size(X,1)*size(X,2),1);
    SoundFieldX = reshape(SoundFieldX,size(X,2),Stepoints*size(X,1))';
    
    SoundFieldY = repmat(FieldYi',size(Y,1)*size(Y,2),1);
    SoundFieldY = reshape(SoundFieldY,size(Y,2),Stepoints*size(Y,1))';

    SoundFieldZ = repmat(FieldZi',size(Z,1)*size(Z,2),1);
    SoundFieldZ = reshape(SoundFieldZ,size(Z,2),Stepoints*size(Z,1))';
   
    
% Forward transfer matrix
    r = sqrt( (TransducerX - SoundFieldX).^2 + (TransducerY - SoundFieldY).^2 + (TransducerZ - SoundFieldZ).^2 );
    H_temp = exp( -(j * Kn) .* r) .* DS ./ r;
    H=( j * rho* f ) .*reshape(sum(H_temp,2),size(X,1),Stepoints).';
  
% Calculate the pressure distribution in the ultrasound field
    P_temp = H * U;
    P = [P;P_temp]; 
end

fprintf('\tOK\n');
P3D = P(1:PointNum);
P3D = reshape(P3D,d1,d2,d3);
P3Dabs = abs(P3D);


I3D =exp(-2*Alpha*A).*P3Dabs.^2 / ( 2 * rho * c );

I=sum(sum(I3D(:,:,16))); %sum of intensity on focal plane
PowerC=Xstep*Xstep*I;
I3Dnew=I3D*200/PowerC;

[xi,zi]=meshgrid(-Xmax:Xstep:Xmax,-Ymax:Ystep:Ymax);
surf(xi,zi,reshape(I3Dnew(16,:,:),[31 31])); %plot ultrasound intensity on XZ plane

Q3D = I3Dnew * 2 * Alpha;

fprintf('\tQover\n');

% ambient temperature:
T0 = 25;	% (degrees C)

% continuous sonication
t_i = 0;	% initial sonication duration (s)

% periodic sonications
n_c = 40;		% number of pulse cycles
D = 60;		% duty cycle (%)
t_p = 0.5;	% pulse cycle period (s)

% cooling period
t_c = 0;	% cool-off duration (s)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%% Computed parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% determine proper timestep:
if(n_c==0)
  dt = min([0.1,t_i/5]);
elseif(t_i==0)
  dt = min(0.1,0.01*D*t_p/5);
else
  dt = min([0.1,0.01*D*t_p/5,t_i/5]);
end

% determine number of integration steps at each stage:
N_i = round(t_i/dt);
if(n_c==0)
  N_p=0;
else
  N_p = round(t_p/dt);
end
N_c = round(t_c/dt);
N = N_i+n_c*N_p+N_c;		% total number of integration steps
%T = dt*(N-1);			% total simulation duration

% build input vector, which catalogs the HIFU beam on/off operation 
if(N_p~=0)
  pulse = zeros(1,N_p);
  for n=1:round(0.01*D*N_p);
    pulse(n) = 1;
  end
  pulses = pulse;
  for m=2:n_c
    pulses = [pulses,pulse];
  end
end
if(N_c~=0)
  cooloff = zeros(1,N_c);
end 
if(N_i~=0)
  initial = ones(1,N_i);
end
for n=1:N_i
  input(n) = initial(n);
end
for n=N_i+1:N_i+n_c*N_p
  input(n) = pulses(n-N_i);
end
for n=N_i+n_c*N_p+1:N
    input(n) = cooloff(n-N_i-n_c*N_p);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%% Computed thermal dodse %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


N100=0;
T=T0*ones(d1,d2,d3);
D=zeros(PointNum,1);
Dmat=zeros(d1,d2,d3);
%Calculate temperature rise with FDTD method
for n=1:N
    T(2:d1-1,2:d2-1,2:d3-1)=T(2:d1-1,2:d2-1,2:d3-1)+dt*k*(T(3:d1,2:d2-1,2:d3-1)+T(2:d1-1,3:d2,2:d3-1)+T(2:d1-1,2:d2-1,3:d3)...
                            +T(1:d1-2,2:d2-1,2:d3-1)+T(2:d1-1,1:d2-2,2:d3-1)+T(2:d1-1,2:d2-1,1:d3-2)-6*T(2:d1-1,2:d2-1,2:d3-1))/Cp/rho/h/h...
                           +dt*input(n)*Q3D(2:d1-1,2:d2-1,2:d3-1)/rho/Cp;
    MaxT(n)=max(max(max(T)));
    if(MaxT(n)>=100)
        if(N100==0)
        N100=n*dt;
        end
    end
    T_temp=reshape(T,PointNum,1);
    D=D+equivalent_time(T_temp,PointNum,T0);
    
    

end
figure
n=1:N;
plot(n*dt,MaxT); %plot temperature vs time

figure   %plot thermal dose on XZ plane
[x,z]=meshgrid(-Xmax:Xstep:Xmax);
Dmat=dt*reshape(D,d1,d2,d3)/60;
Dmat1=reshape(Dmat(:,16,:),[31,31]);
contour(100*z,100*x,Dmat1,[240 240],'r'); %separate larger than 240 and lower with red line
xlabel('z (cm)')
ylabel('x (cm)')
title('Thermal Dose (CEM43C)')
grid on

fprintf('\tover\n');
  

