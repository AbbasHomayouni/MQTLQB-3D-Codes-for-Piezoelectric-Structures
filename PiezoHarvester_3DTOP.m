% 3DTOPIEZO_ENERGY-HARVESTING // Abbas Homayouni-Amlashi et al. 2024
clc ;clear ;close all
%% GENERAL DEFINITIONS
La = 3e-2 ; Wa = 1.0e-2 ; Ha = 0.1e-2 ;  % Pieozoelectric geometrical dimensions (length, width, height) (m)
Lp = 3e-2 ; Wp = 1.0e-2 ; Hp = 0.1e-2 ;  % Passive material geometrical dimensions (length, width, height) (m)
nelx = 3*33 ; nely = 33 ; nelz = 4 ; % Number of elements in each direction
penalKuu = 3; penalKup = 6; penalKpp = 4 ; penalPol = 1 ; % Penalization factors
EL_T = 1; % Element type 1: trilinear, 2- quadratic
volfrac = 0.5; % Volume fraction
Max_loop = 120; % Maximum number of Iteration
pol_dir = 'z'; % Piezoelectric polarization direction
move = 0.2; % Optimization variable update move
ft = 2; % 1= Density filter, 2&3= projection with eta and beta as parameters
ftBC = 'N';
rmin = 1.6; % Filter radius
eta = 0.5; % Threshold
beta = 2; % Sharpness factor
penalCnt = {60,5,10,0.2}; % Continuation scheme on penalKuu {istart, maxPar, isteps, deltaPar}
betaCnt = {60,60,10,2}; % Continuation scheme on beta {istart, maxPar, isteps, deltaPar}
Dir=0; % Direction of force x=2; y=1; z=0;
omega = 150; % Excitation frequency (Hz)
wj = 0.1; % Objective function weigthing factor
Mass = 30; % Mass of attachement (gram)
%% MATERIAL PROPERTIES (PZT 4)
ro_p = 7500; % Density of piezoelectric material (kg/m^3)
C_p = [1.3900    0.7784    0.7428    0.0000    0.0000    0.0000
0.7784    1.3900    0.7428    0.0000    0.0000    0.0000
0.7428    0.7428    1.1541    0.0000    0.0000    0.0000
0.0000    0.0000    0.0000    0.2564    0.0000    0.0000
0.0000    0.0000    0.0000    0.0000    0.2564    0.0000
0.0000    0.0000    0.0000    0.0000    0.0000    0.3058]*1.0e+11; % Piezoelectric stiffness tensor
e   = [0.0000    0.0000    0.0000   0.0000  12.7179    0.0000
0.0000    0.0000    0.0000  12.7179   0.0000    0.0000
-5.2028   -5.2028   15.0804   0.0000   0.0000    0.0000]; % Piezoelectric coupling matrix
Ep  = [0.6746    0.0000    0.0000
0.0000    0.6746    0.0000
0.0000    0.0000    0.5867]*1.0e-08; % Piezoelectric permittivity matrix
ro_s = 2710; % Density of Passive material (kg/m^3)
EE = 70e9; % Young modulus of elasticity
NU = 0.30 ; % Poisson ratio
C_s = (EE/((1+NU)*(1-2*NU)))*[1-NU NU NU 0 0 0 ; NU 1-NU NU 0 0 0 ; NU NU 1- NU 0 0 0 ;
0 0 0 (1-2*NU)/2 0 0 ; 0 0 0 0 (1- 2*NU)/2 0 ; 0 0 0 0 0 (1- 2*NU)/2]; % Stiffness matrix for passive material
%% PREPARE FINITE ELEMENT ANALYSIS
[C_p_1,e_1,Ep_1] = Matrix_Rotation(C_p,e,Ep,pol_dir);
[kuu,kup,kpp,m_p,ndofPZT,EL_NN,TOPNODS,BOTNODS,FRNODS,BAKNODS,LEFNODS,RTNODS] = FEM(La,Wa,Ha,nelz,nelx,nely,C_p_1,e_1,Ep_1,ro_p,EL_T); % Pizeoelectric elemental matrices
[ks,~,~,m_s,~,~,~,~,~,~,~,~] = FEM(Lp,Wp,Hp,max(1,nelz-2),nelx,nely,C_s,zeros(3,6),zeros(3,3),ro_s,EL_T); % Passive material elemental matrices
k0 = max(abs(kuu(:)));beta0 = max(kpp(:));alpha0 = max(kup(:)); M0 = max(m_p(:));  % Normalization Factors
kuu = kuu/k0;ks = ks/k0;kup = kup/alpha0;kpp = kpp/beta0;gamma = (k0*beta0)/(alpha0^2); % Application of normalization I
m_p = m_p/M0; m_s = m_s/M0; omega = M0*(omega*2*pi)^2/k0; % Application of normalization II
kuu_LT = kuu(tril(true(size(kuu)))); % Vector of lower triangular matrix
kpp_LT = kpp(tril(true(size(kpp)))); % Vector of lower triangular element of piezoelectric dielectric stifness matrix
mp_LT = m_p(tril(true(size(m_p)))); % Vector of lower triangular element of piezoelectric mass matrix
ms_LT = m_s(tril(true(size(m_s)))); % Vector of lower triangular element of passive material mass matrix
ks_LT = ks(tril(true(size(ks)))); % Vector of lower triangular matrix
ndof = 3*ndofPZT; % Mechanical degrees of freedom
nele = nelx*nelz*nely; % Number of elements
ElNum = reshape(1:nele,nelz,nelx,nely); % Element indexing
% Building connectivity matrix
NNlinear=(nelz+1)*(nelx+1)*(nely+1);
edg1=reshape(1:NNlinear,nelz+1,nelx+1,nely+1);
edg2=reshape(NNlinear+1:NNlinear+(nely+1)*(nelx+1)*nelz,nelz,nelx+1,nely+1);
edg3=reshape(NNlinear+(nely+1)*(nelx+1)*nelz+1:NNlinear+(nely+1)*(nelx+1)*nelz+(nelz+1)*(nely+1)*nelx,nelz+1,nelx,nely+1);
edg4=reshape(NNlinear+(nely+1)*(nelx+1)*nelz+(nelz+1)*(nely+1)*nelx+1:...
NNlinear+(nely+1)*(nelx+1)*nelz+(nelz+1)*(nely+1)*nelx+nely*(nelz+1)*(nelx+1),nelz+1,nelx+1,nely);
n=0;
for i=1:nely
for j=1:nelx
for k = 1:nelz
n=n+1;
EDG1=edg1([k,k+1],[j,j+1],[i,i+1]);
EDG2=edg2(k,[j,j+1],[i,i+1]);
EDG3=edg3([k,k+1],j,[i,i+1]);
EDG4=edg4([k,k+1],[j,j+1],i);
ED(n,:)=[EDG1(:);EDG2(:);EDG3(:);EDG4(:)]';
end
end
end
EDM = ED(:,[2,4,3,1,6,8,7,5,14,10,13,9,20,19,17,18,16,12,15,11]);
edofMatPZT = EDM(:,1:EL_NN); % Electrical connectivitry matrix
edofMat(:,3:3:3*EL_NN)=3*edofMatPZT;edofMat(:,2:3:3*EL_NN)=3*edofMatPZT-1;edofMat(:,1:3:3*EL_NN)=3*edofMatPZT-2; % Mechanical connectivitry matrix
[sI,sII] = deal([]);
for j = 1:3*EL_NN
sI = cat(2,sI,j:3*EL_NN);
sII = cat(2,sII,repmat(j,1,3*EL_NN-j+1));
end
[iK,jK] = deal(edofMat(:,sI)',edofMat(:,sII)');
Iar = sort([iK(:),jK(:)],2,'descend'); clear iK jK % Assembly indexing (stiffness matrix)
[sI,sII] = deal([]);
for j = 1 : EL_NN
sI = cat(2,sI,1:3*EL_NN);
sII = cat(2,sII,repmat(j,1,3*EL_NN));
end
[iKup,jKup] = deal(edofMat(:,sI)',edofMatPZT(:,sII)');
Iar_up = [iKup(:),jKup(:)]; clear iKup jKup; % Assembly indexing for piezoelectric coupling matrix
[sI,sII] = deal([]);
for j = 1:EL_NN
sI = cat(2,sI,j:EL_NN);
sII = cat(2,sII,repmat(j,1,EL_NN-j+1));
end
[iKp,jKp] = deal(edofMatPZT(:,sI)',edofMatPZT(:,sII)');
Iar_p = sort([iKp(:),jKp(:)],2,'descend'); clear iKp jKp % Assembly indexing for piezoelectric dielectric stiffness matrix
%% ACTIVE & PASSIVE DOMAINS
Passive_el=ElNum([2:nelz-1],:,:); Passive_el=Passive_el(:); % Definition of passive elements
Active_el=setdiff(1:nele,Passive_el); Active_el=Active_el(:);% Definition of active elements
%% DEFINITION OF BOUNDARY CONDITION
DE = ElNum(:,1,:); DE=DE(:); % Desired element for left clamped side
DNN = LEFNODS; % Desired node numbers (elemental left nodes)
fixeddof = edofMat(DE,[3*DNN-1,3*DNN-2]); fixeddof = fixeddof(:); % Fix mechanical DOFs
freedofs = setdiff(1:ndof,fixeddof); lf = length(freedofs); % Free mechanical DOFs
%% DEFINITION OF ELECTRODES
PE1 = edofMatPZT(ElNum(1,:,:),TOPNODS); PE1 = PE1(:);
PE2 = edofMatPZT(ElNum(1,:,:),BOTNODS); PE2 = PE2(:);
PE3 = edofMatPZT(ElNum(nelz,:,:),TOPNODS); PE3 = PE3(:);
PE4 = edofMatPZT(ElNum(nelz,:,:),BOTNODS); PE4 = PE4(:);
en = [PE1;PE2;PE3;PE4]; en = unique(en(:)); % Equipotential nodes
pn = edofMatPZT(Passive_el,:); pn = unique(setdiff(pn(:),en(:))); % Nodes of passive elements
fn = setdiff(1:ndofPZT,[en;pn]); fn = unique(fn(:)); % FreeNodes
Nelec=2; % Number of potential electrodes
B = sparse(ndofPZT,Nelec); % Creation of null Bolean matrix
B(PE1',1)=1; B(PE4',2)=1;  % Creation of Bolean matrix
Up=zeros(ndofPZT,1);ADJ1 = zeros(ndof,1); ADJ2 = zeros(ndof,1); % Creation of null displacement vector
%% FORCE DEFINITION
nf = 1; % Number of forces
F = sparse(ndof,nf); % Definition of null vector for the force
DE = ElNum(:,1,:);DE=DE(:); % Desired element for appliation of force
DNN = LEFNODS; % Desired elemental node number
Fe=edofMat(DE,3*DNN-Dir); Fe = Fe(:); % Desired mechanical degrees of freedom
F(Fe,1) = +1; % Amplitude of the force
Ftot = [F(freedofs,:);zeros(length(fn),nf);zeros(Nelec,nf)];
%% DEFINITION OF ATTACHMENT MASS
sMass=zeros(nele,1);
sMass(ElNum([2:nelz-1],nelx,ceil(0.4*nely):ceil(0.6*nely)))=1; % Distribution of mass
le = Lp/nelx; we = Wp/nely; he = Ha; % Dimension of each element
ro_M = Mass*1e-3/(le*we*he)/length(find(sMass)); % Density of heavy elements
sMMass = (ro_M/ro_p)*mp_LT(:).*sMass';
sMMass = reshape(  sMMass, length( mp_LT(:) ) * nele, 1 );
M_Att = sparse(Iar( :, 1 ),Iar( :, 2 ),sMMass(:)); % Mass matrix containing only the attachement mass
%% SOLID & VOID DOMAINS
VOID = []; % Definition of void elements
SOLID = []; % Definition of solid elements
NVS = setdiff(1:nele,union(VOID(:),SOLID(:))); NVS = NVS(:); % Definition of non-void and non-solid elements
%% FILTER INITIALIZATION [Ferrari & Sigmund 2020]
if ftBC == 'N', bcF = 'symmetric'; else, bcF = 0; end
prj = @(v,eta,beta) (tanh(beta*eta)+tanh(beta*(v(:)-eta)))./(tanh(beta*eta)+tanh(beta*(1-eta)));  % projection
deta = @(v,eta,beta) - beta * csch( beta ) .* sech( beta * ( v( : ) - eta ) ).^2 .*sinh( v( : ) * beta ) .* sinh( ( 1 - v( : ) ) * beta );  % projection eta-derivative
dprj = @(v,eta,beta) beta*(1-tanh(beta*(v-eta)).^2)./(tanh(beta*eta)+tanh(beta*(1-eta)));% proj. x-derivative
cnt = @(v,vCnt,l) v+(l>=vCnt{1}).*(v<vCnt{2}).*(mod(l,vCnt{3})==0).*vCnt{4};
[dy,dz,dx] = meshgrid(-ceil(rmin)+1:ceil(rmin)-1,-ceil(rmin)+1:ceil(rmin)-1,-ceil(rmin)+1:ceil(rmin)-1 );
h = max( 0, rmin - sqrt( dx.^2 + dy.^2 + dz.^2 )); % Conv. kernel
Hs = imfilter( ones( nelz, nelx, nely ), h, bcF ); dHs = Hs; % Matrix of weights (filter)
%% INITIALIZE ITERATION
x = repmat(volfrac,nelz,nelx,nely); x(VOID) = 0; x(SOLID) = 1; xPhys = x;  % Initial guess for the densities
pol = repmat(0.5,[nelz,nelx,nely]); % Initial values for polarization
loop = 0;
Density_change = 1;
E0 = 1; Emin = 1e-9;
e0 = 1; eMin = 1e-9;
eps0 = 1; epsMin = 1e-9;
as = []; % Initialize asymptotes
dv0 = ones(nelz,nelx,nely); % Volume sensitivity
penalratio_up = penalKup/penalKuu; penalratio_pp = penalKpp/penalKuu;  % Penalty ratios for continuation scheme
xold1 = [x(:);pol(:)]; % Vector of variables for previous iteration
xold2 = [x(:);pol(:)]; % Vector of variables for 2nd previous iteration
%% OPTIMIZATION ITERATIONS
while  loop < Max_loop; tic
loop = loop+1;
%% COMPUTE PHYSICAL DENSITY FIELD (AND ETA IF PROJECT.) [Ferrari & Sigmund 2020]
xTilde  = imfilter( reshape( x, nelz, nelx, nely ), h, bcF ) ./ Hs; xPhys(NVS) = xTilde(NVS);      % filtered field
if ft > 1    % Compute optimal eta* with Newton
f = ( mean( prj( xPhys, eta, beta ) ) - volfrac )  * (ft == 3);  % Function (volume)
while abs( f ) > 1e-6   % Newton process for finding opt. eta
eta = eta - f / mean( deta( xPhys, eta, beta ) );
f = mean( prj( xPhys, eta, beta ) ) - volfrac;
end
dHs = Hs ./ reshape( dprj( xPhys, eta, beta ), nelz, nelx, nely );   % Sensitivity modification
xPhys = prj( xPhys, eta, beta );     % Projected (physical) field
end
%% FE-ANALYSIS
xPhys = reshape(xPhys,nelz,nelx,nely);
sM = ones(length(mp_LT( : )),1).*(Emin+xPhys(:)'.^penalKuu*(E0-Emin));
sM(:,Active_el)= mp_LT( : ).* sM(:,Active_el);
sM(:,Passive_el)= ms_LT( : ).* sM(:,Passive_el);
sM = reshape(  sM, length( mp_LT(:) ) * nele, 1 );
sKuu = ones(length(kuu_LT( : )),1).*(Emin+xPhys(:)'.^penalKuu*(E0-Emin));
sKuu(:,Active_el)= kuu_LT( : ).* sKuu(:,Active_el);
sKuu(:,Passive_el)= ks_LT( : ).* sKuu(:,Passive_el);
sKuu = reshape(sKuu, length(kuu_LT(:)) * nele, 1 );
sKup = ones(length(kup(:)),1)*(eMin+xPhys(:)'.^penalKup*(e0-eMin).*((2*pol(:)-1)'.^penalPol));
sKup(:,Active_el)=kup(:).* sKup(:,Active_el);
sKup(:,Passive_el)=zeros(size(kup(:))).*sKup(:,Passive_el);
sKpp = ones(length(kpp_LT(:)),1)*(epsMin+xPhys(:)'.^penalKpp*(eps0-epsMin));
sKpp(:,Active_el)= kpp_LT(:).* sKpp(:,Active_el);
sKpp(:,Passive_el)= zeros(size(kpp_LT(:))).*sKpp(:,Passive_el);
M= sparse(Iar( : , 1 ),Iar( : , 2 ),sM); % Global Mass matrix
Kuu = sparse(Iar( : , 1 ),Iar( : , 2 ),sKuu); KuuM = Kuu-(M)*omega; % Global stiffness matrix
Kuu = Kuu-(M+M_Att)*omega;Kuu = Kuu+Kuu'-diag(diag(Kuu)); % Global dynamic stiffness matrix
Kup = sparse( Iar_up( :, 1 ), Iar_up ( :, 2 ),sKup(:)); % Global piezoelectric coupling matrix
Kpp = sparse( Iar_p( :, 1 ), Iar_p ( :, 2 ),sKpp(:)); Kpp =Kpp+Kpp'-diag(diag(Kpp)); % Global piezoelectric permittivity matrix
Ktot = [Kuu(freedofs,freedofs),Kup(freedofs,fn),Kup(freedofs,en)*B(en,:);
Kup(freedofs,fn)',-gamma*Kpp(fn,fn),-gamma*Kpp(fn,en)*B(en,:);
B(en,:)'*Kup(freedofs,en)',-gamma*B(en,:)'*Kpp(fn,en)',-gamma*B(en,:)'*Kpp(en,en)*B(en,:)];
U =  ( Ktot \ Ftot); % Response of the system
Up(fn,:) = U(lf+1:lf+length(fn),:); Up(en,:) = B(en,:)*U(lf+length(fn)+1:end,:);
%% OBJECTIVE FUNCTION AND SENSITIVITY ANALYSIS
lambda1 = zeros (ndof,nf);lambda2 = zeros (ndof,nf);mu1 = zeros (ndofPZT,nf);mu2 = zeros (ndofPZT,nf);
ADJ1 =  Ktot\[-Kuu(freedofs,freedofs)*U(1:lf,:);zeros(length(fn),nf);zeros(Nelec,nf)]; % First adjoint vector
lambda1(freedofs,:) = ADJ1(1:lf,:); mu1(fn,:) = ADJ1(lf+1:lf+length(fn),:);mu1(en,:) = B(en,:)*ADJ1(lf+length(fn)+1:end,:);
ADJ2 =  Ktot\[zeros(lf,nf);-gamma*Kpp(fn,fn)*Up(fn,:)-gamma*Kpp(fn,en)*Up(en,:);-gamma*B(en,:)'*Kpp(en,en)*Up(en,:)-gamma*B(en,:)'*Kpp(en,fn)*Up(fn,:)]; % Second adjoint vector
lambda2(freedofs,:) = ADJ2(1:lf,:); mu2(fn,:) = ADJ2(lf+1:lf+length(fn),:);mu2(en,:) = B(en,:)*ADJ2(lf+length(fn)+1:end,:);
Uu_i = zeros(ndof,1); Wm = 0 ; We = 0 ; dc = zeros(nelz,nelx,nely) ; dp = zeros(nelz,nelx,nely);
for i = 1:nf
Uu_i(freedofs,1) = U(1:lf,i);Up_i = Up(:,i);
lambda1_i = lambda1(:,i); lambda2_i = lambda2(:,i);
mu1_i = mu1(:,i);mu2_i = mu2(:,i);
Wm = Wm+Uu_i'*KuuM*Uu_i; % Mechanical energy
We = We+Up_i'*Kpp*Up_i*gamma; % Electrical energy
dcKuuE(Active_el,:) = wj*((((1/2)*Uu_i(edofMat(Active_el,:)) + lambda1_i(edofMat(Active_el,:)))*kuu).*Uu_i(edofMat(Active_el,:)))-(1-wj)*((lambda2_i(edofMat(Active_el,:))*kuu).*Uu_i(edofMat(Active_el,:)));
dcKuuE(Passive_el,:) = wj*((((1/2)*Uu_i(edofMat(Passive_el,:)) + lambda1_i(edofMat(Passive_el,:)))*ks).*Uu_i(edofMat(Passive_el,:)))-(1-wj)*((lambda2_i(edofMat(Passive_el,:))*ks).*Uu_i(edofMat(Passive_el,:)));
dcKupE(Active_el,:) = wj*((lambda1_i(edofMat(Active_el,:))*kup).*Up_i(edofMatPZT(Active_el,:)) + ((Uu_i(edofMat(Active_el,:)))*kup).*mu1_i(edofMatPZT(Active_el,:)))-(1-wj)*((lambda2_i(edofMat(Active_el,:))*kup).*Up_i(edofMatPZT(Active_el,:)) + ((Uu_i(edofMat(Active_el,:)))*kup).*mu2_i(edofMatPZT(Active_el,:)));
dcKupE(Passive_el,:) = 0;
dcKppE(Active_el,:) = wj*((-mu1_i(edofMatPZT(Active_el,:))*kpp).*Up_i(edofMatPZT(Active_el,:)))-(1-wj)*((1/2)*(Up_i(edofMatPZT(Active_el,:))*kpp).*Up_i(edofMatPZT(Active_el,:)) - (mu2_i(edofMatPZT(Active_el,:))*kpp).*Up_i(edofMatPZT(Active_el,:)));
dcKppE(Passive_el,:) = 0;
dcME(Active_el,:) = wj*((((1/2)*Uu_i(edofMat(Active_el,:)) + lambda1_i(edofMat(Active_el,:)))*(-m_p*omega)).*Uu_i(edofMat(Active_el,:)))-(1-wj)*((lambda2_i(edofMat(Active_el,:))*(-m_p*omega)).*Uu_i(edofMat(Active_el,:)));
dcME(Passive_el,:) = wj*((((1/2)*Uu_i(edofMat(Passive_el,:)) + lambda1_i(edofMat(Passive_el,:)))*(-m_s*omega)).*Uu_i(edofMat(Passive_el,:)))-(1-wj)*((lambda2_i(edofMat(Passive_el,:))*(-m_s*omega)).*Uu_i(edofMat(Passive_el,:)));
dcKuu = reshape(full(sum(dcKuuE,2)),[nelz,nelx,nely]);
dcKup = reshape(full(sum(dcKupE,2)),[nelz,nelx,nely]);
dcKpp = gamma*reshape(full(sum(dcKppE,2)),[nelz,nelx,nely]);
dcM = reshape(full(sum(dcME,2)),[nelz,nelx,nely]);
dc = dc + penalKuu*(E0-Emin)*xPhys.^(penalKuu-1).*dcKuu+penalKup*(e0-eMin)*xPhys.^(penalKup-1).*dcKup.*((2*pol-1).^(penalPol))+penalKpp*(eps0-epsMin)*xPhys.^(penalKpp-1).*dcKpp+dcM; % Density variable sensitivity
dp = dp + (e0-eMin)*2*penalPol*((2*pol-1).^(penalPol-1)).*xPhys.^penalKup.*dcKup;   % Polarization variable sensitivity
end
c = wj*Wm-(1-wj)*We; % Objective function
dc = imfilter( reshape( dc, nelz, nelx, nely ) ./ dHs, h, bcF ); % Filter objective sensitivity
dv = imfilter( reshape( dv0, nelz, nelx, nely ) ./ dHs, h, bcF ); % Filter compliance sensitivity
%% UPDATING OPTIMIZATION VARIABLES (Ferrari & Sigmund 2020)
[Xupdate,as ,lmid ]= ocUpdate (loop , [x(:);pol(:)], [dc(:);dp(:)] ,[sum(xPhys(:))/(volfrac*nele) - 1],[dv(:)' / (volfrac*nele),0*pol(:)']' ,[move ,0.7 ,1.2] ,xold1 ,xold2 ,as , beta );
xnew = Xupdate(1:nele,1);xnew(VOID)=0;xnew(SOLID)=0;% Vector of updated density variables
Density_change = max(abs(xnew(:)-x(:)));
xold2 = xold1(:);xold1 = [x(:);pol(:)];
pol = reshape(Xupdate(nele+1:2*nele,1),nelz,nelx,nely); % Vector of updated polarization variables
x (NVS)= xnew(NVS);
%% CONTINIUATION SCHEME ON PENALIZATION FACTORS & BETA
[penalKuu ,beta] = deal(cnt(penalKuu ,penalCnt,loop),cnt(beta,betaCnt,loop));
penalKup=penalKuu*penalratio_up; penalKpp=penalKuu*penalratio_pp;
%% PRESENTATION OF RESULTS
fprintf(' It:%2.0i Time:%3.2fs Obj:%3.4e Wm.:%3.4e We.:%3.4e Vol:%3.3f ch:%3.3f\n ',loop,toc,c,Wm,We,mean(xPhys(:)),Density_change);
Display(xPhys,pol,nelz,nelx,nely,Active_el,Passive_el,nele,ElNum)
end
%% PLOT DEFORMATION (ELEMENTAL)
Uu = zeros(ndof,1);
Uu(freedofs,1) = U(1:lf,1);
Uu(edofMat(:,[3:3:24])) = Uu(edofMat(:,3:3:24))- Uu(edofMat(1,3));
Deformation(Uu,xPhys,nelz,nelx,nely,edofMat,ElNum) % Plot the deformation

%% PLOT DEFORMATION (ELEMENTAL) // ABBAS HOMAYOUNI-AMLASHI 2024
function Deformation(Uu,xPhys,nelz,nelx,nely,edofMat,ElNum)
figure (2); ax = gca;
AMP = 5/ max(abs(Uu(:))); % Amplification & normalization of deformation
view( [ 50, 20 ] ); % View angle
face = [ 1 2 3 4 ; 2 6 7 3 ; 4 3 7 8 ; 1 5 8 4 ; 1 2 6 5 ; 5 6 7 8 ];
xPhys=reshape(xPhys,nelz,nelx,nely);
for elz = 1:nelz
for elx = 1:nelx
for ely = 1:nely
if xPhys(elz,elx,ely)>0.5
Ue = AMP*Uu(edofMat(ElNum(elz,elx,ely),:));
ly = -(ely-nely)-1; lx = elx-1; lz = -(elz-nelz)-1;
xx_box = [lx lx+1 lx+1 lx lx lx+1 lx+1 lx]';
yy_box = [ly ly ly ly ly+1 ly+1 ly+1 ly+1]';
zz_box = [lz lz lz+1 lz+1 lz lz lz+1 lz+1]';
patch('Faces',face,'Vertices',[xx_box,yy_box,zz_box],'FaceColor','none')
xx = [Ue(1,1)+lx Ue(4,1)+lx+1 Ue(7,1)+lx+1 Ue(10,1)+lx Ue(13,1)+lx Ue(16,1)+lx+1 Ue(19,1)+lx+1 Ue(22,1)+lx]';
yy = [Ue(2,1)+ly+1 Ue(5,1)+ly+1 Ue(8,1)+ly+1 Ue(11,1)+ly+1 Ue(14,1)+ly Ue(17,1)+ly Ue(20,1)+ly Ue(23,1)+ly]';
zz = [Ue(3,1)+lz Ue(6,1)+lz Ue(9,1)+lz+1 Ue(12,1)+lz+1 Ue(15,1)+lz Ue(18,1)+lz Ue(21,1)+lz+1 Ue(24,1)+lz+1]';
Dis_R = max(abs(Uu(edofMat(ElNum(elz,elx,ely),:)))/max(abs(Uu(:))));
patch('Faces',face,'Vertices',[xx,yy,zz],'FaceColor',[Dis_R,0.4,1-Dis_R])
end
end
end
end
box on; ax.BoxStyle = 'full'; axis equal; set(gca,'XTick',[], 'YTick', [],'ZTick', []);axis off; drawnow; end

%% PRESENTATION OF RESULTS (2D & 3D) // ABBAS HOMAYOUNI-AMLASHI 2024
function Display(xPhys,pol,nelz,nelx,nely,Active_el,Passive_el,nele,ElNum)
figure (1)
if nelz<=6;AX= subplot (4, nelz, [1:2*nelz]); cla();
else; AX=subplot (1, 2, 1); cla(); end
Xactive=zeros(nele,1);Xactive(Active_el,1)=xPhys(Active_el);
Xpassive=zeros(nele,1);Xpassive(Passive_el,1)=xPhys(Passive_el);
isovals_active = shiftdim( flipud(reshape( Xactive, nelz, nelx, nely )), 1);
isovals_active = smooth3( isovals_active, 'box', 1 );
patch(isosurface(isovals_active, .5),'FaceColor','m','EdgeColor','none');
patch(isocaps(isovals_active, .5),'FaceColor','m','EdgeColor','none');
isovals_Passive = shiftdim( flipud(reshape( Xpassive, nelz, nelx, nely )), 1);
isovals_Passive = smooth3( isovals_Passive, 'box', 1 );
patch(isosurface(isovals_Passive, .5),'FaceColor',[0 1 1],'EdgeColor','none');
patch(isocaps(isovals_Passive, .5),'FaceColor',[0 1 1],'EdgeColor','none');
isovals = shiftdim( reshape( ones(nele,1), nelz, nelx, nely ), 1);
isovals = smooth3( isovals, 'box', 1 );
patch(isosurface(isovals, .5),'FaceColor','none','EdgeColor','none');
patch(isocaps(isovals, .5),'FaceColor','none','EdgeColor','none');
set(AX,'XTick',[], 'YTick', [],'ZTick', []);
title('Density');view( [ 120, 30 ] ); axis equal tight; camlight; box on;AX.BoxStyle = 'full';drawnow
if nelz>6; n=0;
for i=1:nelz; for j=1:nelx; for k=1:nely
if  xPhys (i,j,k)>0.9 && ismember(ElNum(i,j,k),Active_el)==1;
n=n+1; Coordinate(n,:)=0.01*[i,j,k];
COl(n,:) = [0.5-0.5*pol(i,j,k),0,0.5+0.5*pol(i,j,k)];
end
end; end; end
if exist('Coordinate')
AX2=subplot (1, 2, 2);
pcshow([Coordinate(:,2) Coordinate(:,3) Coordinate(:,1)],COl,'MarkerSize',120,'BackgroundColor',[1,1,1])
set(AX2,'XTick',[], 'YTick', [],'ZTick', []);title('Polarization');view( [ 30, 30 ] );axis equal tight;box on;AX2.BoxStyle = 'full';drawnow
end;end
if nelz<=6
for NL=2*nelz+1:3*nelz
XX(:,:)=xPhys(NL-2*nelz,:,:);
ax(NL)=subplot ( 4 , nelz , NL );
imagesc(1-XX(:,:)') ;colormap(ax(NL),gray) ;
set ( ax(NL) , 'XTick' , [ ] , 'YTick' , [ ] , 'XTicklabel' , [ ] ,...
'YTicklabel' , [ ] , 'xcolor' , 'w' , 'ycolor' , 'w')
xlabel ( sprintf ( 'Layer Number = %.0f' , NL-2*nelz ) , 'Color' , 'k')
title('Densities');
axis equal ; axis tight ; caxis([0 1]); drawnow ; hold on
end
for NL=3*nelz+1:4*nelz
XX(:,:)=(xPhys(NL-3*nelz,:,:).*(pol(NL-3*nelz,:,:).*2-1));
ax(NL)=subplot ( 4 , nelz , NL );
imagesc(XX(:,:)') ;colormap(ax(NL),jet) ;
set ( ax(NL) , 'XTick' , [ ] , 'YTick' , [ ] , 'XTicklabel' , [ ] ,...
'YTicklabel' , [ ] , 'xcolor' , 'w' , 'ycolor' , 'w')
xlabel ( sprintf ( 'Layer Number = %.0f' , NL-3*nelz ) , 'Color' , 'k')
title('Polarization');
axis equal ; axis tight ;caxis([-1 1]); drawnow ; hold on
end
else;end
end

%% TRANSFORMATION OF TENSOR CONSTANTS FOR ANISOTROPIC MATERIALS // ABBAS HOMAYOUNI-AMLASHI 2024
function [Cnew,enew,Epnew]=Matrix_Rotation(C,e,Ep,pol_dir)
if pol_dir == 'y' % For Polarization in the y direction
alpha= 0;beta=pi/2;gamma=0;
elseif pol_dir == 'z' % For Polarization in the z direction
alpha= 0;beta=0;gamma=0;
elseif pol_dir == 'x' % For Polarization in the x direction
alpha= 0;beta=pi/2;gamma=pi/2;
end
xi1=cos(gamma)*cos(alpha)-cos(beta)*sin(alpha)*sin(gamma);
xi2=-sin(gamma)*cos(alpha)-cos(beta)*sin(alpha)*cos(gamma);
xi3=sin(beta)*sin(alpha);
theta1=cos(gamma)*sin(alpha)+cos(beta)*cos(alpha)*sin(gamma);
theta2=-sin(gamma)*sin(alpha)+cos(beta)*cos(alpha)*cos(gamma);
theta3=-sin(beta)*cos(alpha);
psi1=sin(gamma)*sin(beta);psi2=cos(gamma)*sin(beta);psi3=cos(beta);
L=[xi1,theta1,psi1;xi2,theta2,psi2;xi3, theta3, psi3];
Z=[ xi1^2,theta1^2,psi1^2,2*theta1*psi1,2*psi1*xi1,2*xi1*theta1;
xi2^2,theta2^2,psi2^2,2*theta2*psi2,2*psi2*xi2,2*xi2*theta2;
xi3^2,theta3^2,psi3^2,2*theta3*psi3,2*psi3*xi3,2*xi3*theta3;
xi2*xi3,theta2*theta3,psi2*psi3,theta2*psi3+theta3*psi2,psi2*xi3+psi3*xi2,xi2*theta3+xi3*theta2;
xi3*xi1,theta3*theta1,psi3*psi1,theta1*psi3+theta3*psi1,psi1*xi3+psi3*xi1,xi1*theta3+xi3*theta1;
xi1*xi2,theta1*theta2,psi1*psi2,theta1*psi2+theta2*psi1,psi1*xi2+psi2*xi1,xi1*theta2+xi2*theta1];
Cnew=Z*C*Z'; % New stiffness matrix
enew=L*e*Z'; % New coupling matrix
Epnew=L*Ep*L'; % New permittivity matrix
end

%% OCUpdate Algorithm (F. Ferrari et al. 2021)
function [x,as , lmid ]= ocUpdate (loop ,xT ,dg0 ,g1 ,dg1 ,ocPar ,xOld ,xOld1 ,as , beta )
[xU,xL] = deal ( min(xT+ ocPar (1) ,1) , max (xT-ocPar (1) ,0));
if loop <2.5 || beta > 4
as = xT +[ -0.5 ,0.5].*( xU -xL) ./( beta +1) ;
else
tmp = (xT - xOld ) .*( xOld - xOld1 );
gm = ones ( length (xT) ,1);
[gm(tmp >0) , gm(tmp <0) ] = deal ( ocPar (3) ,ocPar (2) );
as = xT + gm .* [-( xOld -as (: ,1)) ,(as (: ,2) -xOld )];
end
xL = max (0.9* as (: ,1) +0.1* xT ,xL); % adaptive lower bound
xU = min (0.9* as (: ,2) +0.1* xT ,xU); % adaptive upper bound
% ----- split (+) and (-) parts of the objective and constraint derivatives
p0_0 = (dg0 >0) .* dg0 ; q0_0 = (dg0 <0) .* dg0 ;
p1_0 = (dg1 >0) .* dg1 ; q1_0 = (dg1 <0) .* dg1 ;
[p0 ,q0] = deal ( p0_0 .*( as (: ,2) -xT).^2 , - q0_0 .*( xT -as (: ,1)) .^2) ;
[p1 ,q1] = deal ( p1_0 .*( as (: ,2) -xT).^2 , - q1_0 .*( xT -as (: ,1)) .^2) ;
% ---------------------- define the primal projection map and dual function
primalProj = @(lm) min (xU ,max (xL ,( sqrt (p0+lm*p1).* as (: ,1)+ sqrt (q0+lm*q1).* as (: ,2)) ...
./( sqrt (p0+lm*p1)+ sqrt (q0+lm*q1))));
psiDual = @(lm) g1-((as(:,2)-xT)'*p1_0-(xT-as(: ,1))'*q1_0)+sum(p1./(max(as(:,2)-primalProj(lm),1e-12))+q1./( max(primalProj(lm)-as (:,1),1e-12)));
% ----------------------- compute the Lagrange multiplier through bisection
lmUp = 1e6; x = xT; lmid = -1;
if psiDual ( 0 ) * psiDual ( lmUp ) < 0 % check if LM is within the interval
lmid = fzero ( psiDual , [ 0, lmUp ] );
x = primalProj ( lmid ); % update desing variables
elseif psiDual (0) < 0 % constraint cannot be active
lmid =0; x= primalProj ( lmid );
elseif psiDual ( lmUp ) > 0 % constraint cannot be fulfilled
lmid = lmUp ; x= primalProj ( lmid );
end
end

%% FINITE ELEMENT MATRICES // ABBAS HOMAYOUNI-AMLASHI 2024
function [kuu,kup,kpp,m,ndofPZT,EL_NN,TOPNODS,BOTNODS,FRNODS,BAKNODS,LEFNODS,RTNODS] = FEM(L,W,H,nelz,nelx,nely,C,e,Ep,ro,EL_T)
le = L/nelx; we = W/nely; he=H; % Element geometry
if EL_T == 1 % Trilinear elements
g = 1/sqrt(3);
GP = [-g -g -g;g -g -g;g g -g;-g g -g;-g -g g;g -g g;g g g;-g g g]; % Gauss quadrature points
J = [ le/2, 0, 0; 0, we/2, 0; 0, 0, he/2]; % Jacobian matrix
detJ = (he*le*we)/8; % Determinant of Jacobian matrix
kuu = 0 ; kup = 0 ; kpp=0 ; m = 0; % Initial values for piezoelectric matrices
for ii=1:8
s=GP(ii,1);t=GP(ii,2);u=GP(ii,3); % s,t,u (natural coordinates)
N1 = (1/8)*(1-s)*(1+t)*(1-u);N2 = (1/8)*(1+s)*(1+t)*(1-u);N3 = (1/8)*(1+s)*(1+t)*(1+u);N4 = (1/8)*(1-s)*(1+t)*(1+u);
N5 = (1/8)*(1-s)*(1-t)*(1-u);N6 = (1/8)*(1+s)*(1-t)*(1-u);N7 = (1/8)*(1+s)*(1-t)*(1+u);N8 = (1/8)*(1-s)*(1-t)*(1+u);
N=[N1,0,0,N2,0,0,N3,0,0,N4,0,0,N5,0,0,N6,0,0,N7,0,0,N8,0,0;
0,N1,0,0,N2,0,0,N3,0,0,N4,0,0,N5,0,0,N6,0,0,N7,0,0,N8,0;
0,0,N1,0,0,N2,0,0,N3,0,0,N4,0,0,N5,0,0,N6,0,0,N7,0,0,N8]; % Matrix of interpolation functions
DN =[ ((t + 1)*(u - 1))/8, -((t + 1)*(u - 1))/8, ((t + 1)*(u + 1))/8, -((t + 1)*(u + 1))/8, -((t - 1)*(u - 1))/8, ((t - 1)*(u - 1))/8, -((t - 1)*(u + 1))/8, ((t - 1)*(u + 1))/8;
(s/8 - 1/8)*(u - 1), -(s/8 + 1/8)*(u - 1), (s/8 + 1/8)*(u + 1), -(s/8 - 1/8)*(u + 1), -(s/8 - 1/8)*(u - 1), (s/8 + 1/8)*(u - 1), -(s/8 + 1/8)*(u + 1), (s/8 - 1/8)*(u + 1);
(s/8 - 1/8)*(t + 1), -(s/8 + 1/8)*(t + 1), (s/8 + 1/8)*(t + 1), -(s/8 - 1/8)*(t + 1), -(s/8 - 1/8)*(t - 1), (s/8 + 1/8)*(t - 1), -(s/8 + 1/8)*(t - 1), (s/8 - 1/8)*(t - 1)]; % dN/d(s,t,u)
Bphi=J\DN; % Piezo Gradient interpolation matrix (Potential to electrical field matrix)
Bu(1,1:3:24)=Bphi(1,:);Bu(2,2:3:24)=Bphi(2,:);Bu(3,3:3:24)=Bphi(3,:);
Bu(6,1:3:24)=Bphi(2,:);Bu(6,2:3:24)=Bphi(1,:);
Bu(4,2:3:24)=Bphi(3,:);Bu(4,3:3:24)=Bphi(2,:);
Bu(5,1:3:24)=Bphi(3,:);Bu(5,3:3:24)=Bphi(1,:); % Strain-displacement matrix
kuu = kuu + transpose(Bu)*C*Bu*detJ; % Stiffness matrix
kup = kup + Bu'*e'*Bphi*detJ;  % Piezoelectric coupling matrix
kpp = kpp + Bphi'*Ep*Bphi*detJ; % Dielectric stiffness matrix
m = m+detJ*ro*(N'*N); % mass matrix
end
EL_NN = 8; % Elemental node numbers
ndofPZT = (nelx+1)*(nelz+1)*(nely+1); % Total electrical degrees of freedom
TOPNODS=[4,3,7,8]; BOTNODS=[1,2,5,6];
FRNODS=[5,6,7,8];  BAKNODS=[1,2,3,4];
LEFNODS=[1,4,5,8]; RTNODS=[2,3,6,7];
elseif EL_T == 2 % Quadratic elements
GPW = [-sqrt(3/5),5/9;0,8/9;sqrt(3/5),5/9]; % Gauss quadrature points and weights
x1=0; y1=we; z1=0; x2=le; y2=we; z2=0; x3=le; y3=we; z3=he; x4=0; y4=we; z4=he;
x5=0; y5=0; z5=0; x6=le; y6=0; z6=0; x7=le; y7=0; z7=he; x8=0; y8=0; z8=he;
x9=(x1+x2)/2;y9=(y1+y2)/2;z9=(z1+z2)/2;
x10=(x2+x3)/2;y10=(y2+y3)/2;z10=(z2+z3)/2;
x11=(x3+x4)/2;y11=(y3+y4)/2;z11=(z3+z4)/2;
x12=(x1+x4)/2;y12=(y1+y4)/2;z12=(z1+z4)/2;
x13=(x2+x6)/2;y13=(y2+y6)/2;z13=(z2+z6)/2;
x14=(x3+x7)/2;y14=(y3+y7)/2;z14=(z3+z7)/2;
x15=(x4+x8)/2;y15=(y4+y8)/2;z15=(z4+z8)/2;
x16=(x1+x5)/2;y16=(y1+y5)/2;z16=(z1+z5)/2;
x17=(x5+x6)/2;y17=(y5+y6)/2;z17=(z5+z6)/2;
x18=(x6+x7)/2;y18=(y6+y7)/2;z18=(z6+z7)/2;
x19=(x8+x7)/2;y19=(y8+y7)/2;z19=(z8+z7)/2;
x20=(x8+x5)/2;y20=(y8+y5)/2;z20=(z8+z5)/2;
xs=le/2;xt=0;xu=0;ys=0;yt=we/2;yu=0;zs=0;zt=0;zu=he/2;
J = [xs ys zs; xt yt zt; xu yu zu]; % Jacobian matrix
detJ = xs*(yt*zu - zt*yu) - ys*(xt*zu - zt*xu) + zs*(xt*yu - yt*xu); % Determinant of Jacobian matrix
kuu = 0 ; kup = 0 ; kpp=0 ; m = 0; % Initial values for piezoelectric matrices
for i=1:3
s = GPW(i,1);
for j=1:3
t = GPW(j,1);
for k=1:3
u = GPW(k,1);
NV(1)=(1-s)*(1+t)*(1-u)*(-s+t-u-2)/8; NV(2)=(1+s)*(1+t)*(1-u)*(s+t-u-2)/8;
NV(3)=(1+s)*(1+t)*(1+u)*(s+t+u-2)/8;NV(4)=(1-s)*(1+t)*(1+u)*(-s+t+u-2)/8;
NV(5)=(1-s)*(1-t)*(1-u)*(-s-t-u-2)/8;NV(6)=(1+s)*(1-t)*(1-u)*(s-t-u-2)/8;
NV(7)=(1+s)*(1-t)*(1+u)*(s-t+u-2)/8;NV(8)=(1-s)*(1-t)*(1+u)*(-s-t+u-2)/8;
NV(9)=(1+t)*(1-u)*(1-s^2)/4;  NV(10)=(1+s)*(1+t)*(1-u^2)/4;
NV(11)=(1+t)*(1+u)*(1-s^2)/4; NV(12)=(1-s)*(1+t)*(1-u^2)/4;
NV(13)=(1+s)*(1-u)*(1-t^2)/4; NV(14)=(1+s)*(1+u)*(1-t^2)/4;
NV(15)=(1-s)*(1+u)*(1-t^2)/4;NV(16)=(1-s)*(1-u)*(1-t^2)/4;
NV(17)=(1-t)*(1-u)*(1-s^2)/4;NV(18)=(1+s)*(1-t)*(1-u^2)/4;
NV(19)=(1-t)*(1+u)*(1-s^2)/4;NV(20)=(1-s)*(1-t)*(1-u^2)/4;
N=zeros(3,60);N(1,1:3:60)=NV(:)';N(2,2:3:60)=NV(:)';N(3,3:3:60)=NV(:)';
DN=[-((t + 1)*(u - 1)*(2*s - t + u + 1))/8, -((t + 1)*(u - 1)*(2*s + t - u - 1))/8, ((t + 1)*(u + 1)*(2*s + t + u - 1))/8, ((t + 1)*(u + 1)*(2*s - t - u + 1))/8, ((t - 1)*(u - 1)*(2*s + t + u + 1))/8, -((t - 1)*(u - 1)*(t - 2*s + u + 1))/8, -((t - 1)*(u + 1)*(2*s - t + u - 1))/8, -((t - 1)*(u + 1)*(2*s + t - u + 1))/8, (s*(t + 1)*(u - 1))/2, -((u^2 - 1)*(t + 1))/4, -(s*(t + 1)*(u + 1))/2, ((u^2 - 1)*(t + 1))/4, ((t^2 - 1)*(u - 1))/4, -((t^2 - 1)*(u + 1))/4, ((t^2 - 1)*(u + 1))/4, -((t^2 - 1)*(u - 1))/4, -(s*(t - 1)*(u - 1))/2, ((u^2 - 1)*(t - 1))/4, (s*(t - 1)*(u + 1))/2, -((u^2 - 1)*(t - 1))/4;
-((s - 1)*(u - 1)*(s - 2*t + u + 1))/8, -((s + 1)*(u - 1)*(s + 2*t - u - 1))/8, ((s + 1)*(u + 1)*(s + 2*t + u - 1))/8, ((s - 1)*(u + 1)*(s - 2*t - u + 1))/8, ((s - 1)*(u - 1)*(s + 2*t + u + 1))/8, -((s + 1)*(u - 1)*(2*t - s + u + 1))/8, -((s + 1)*(u + 1)*(s - 2*t + u - 1))/8, -((s - 1)*(u + 1)*(s + 2*t - u + 1))/8, ((s^2 - 1)*(u - 1))/4, -((u^2 - 1)*(s + 1))/4, -((s^2 - 1)*(u + 1))/4, ((u^2 - 1)*(s - 1))/4, (t*(s + 1)*(u - 1))/2, -(t*(s + 1)*(u + 1))/2, (t*(s - 1)*(u + 1))/2, -(t*(s - 1)*(u - 1))/2, -((s^2 - 1)*(u - 1))/4, ((u^2 - 1)*(s + 1))/4, ((s^2 - 1)*(u + 1))/4, -((u^2 - 1)*(s - 1))/4;
-((s - 1)*(t + 1)*(s - t + 2*u + 1))/8, -((s + 1)*(t + 1)*(s + t - 2*u - 1))/8, ((s + 1)*(t + 1)*(s + t + 2*u - 1))/8, ((s - 1)*(t + 1)*(s - t - 2*u + 1))/8, ((s - 1)*(t - 1)*(s + t + 2*u + 1))/8, -((s + 1)*(t - 1)*(t - s + 2*u + 1))/8, -((s + 1)*(t - 1)*(s - t + 2*u - 1))/8, -((s - 1)*(t - 1)*(s + t - 2*u + 1))/8, ((s^2 - 1)*(t + 1))/4, -(u*(s + 1)*(t + 1))/2, -((s^2 - 1)*(t + 1))/4, (u*(s - 1)*(t + 1))/2, ((t^2 - 1)*(s + 1))/4, -((t^2 - 1)*(s + 1))/4, ((t^2 - 1)*(s - 1))/4, -((t^2 - 1)*(s - 1))/4, -((s^2 - 1)*(t - 1))/4, (u*(s + 1)*(t - 1))/2, ((s^2 - 1)*(t - 1))/4, -(u*(s - 1)*(t - 1))/2];
Bphi= double(J\DN);
B=zeros(6,60);
B(1,1:3:60)= Bphi(1,:);B(2,2:3:60)= Bphi(2,:);B(3,3:3:60)= Bphi(3,:);
B(6,1:3:60)= Bphi(2,:);B(6,2:3:60)= Bphi(1,:);
B(4,2:3:60)= Bphi(3,:);B(4,3:3:60)= Bphi(2,:);
B(5,1:3:60)= Bphi(3,:);B(5,3:3:60)= Bphi(1,:);
kuu = kuu + GPW(i,2)*GPW(j,2)*GPW(k,2)*transpose(B)*C*B*detJ;
kup = kup + GPW(i,2)*GPW(j,2)*GPW(k,2)*B'*e'*Bphi*detJ;
kpp = kpp+ GPW(i,2)*GPW(j,2)*GPW(k,2)*Bphi'*Ep*Bphi*detJ;
m = m + GPW(i,2)*GPW(j,2)*GPW(k,2)* detJ*ro*(N'*N);
end
end
end
EL_NN = 20; % Elemental node numbers
ndofPZT =((2*nelz+1)*(2*nelx+1)-nelx*nelz)*(nely+1)+(nelx+1)*(nelz+1)*nely; % Electrical degrees of freedom
TOPNODS=[4,3,7,8,11,14,15,19]; BOTNODS=[1,2,5,6,9,13,16,17];
FRNODS=[5,6,7,8,17,18,19,20];  BAKNODS=[1,2,3,4,9,10,11,12];
LEFNODS=[1,4,5,8,12,15,16,20]; RTNODS=[2,3,6,7,10,13,14,18];
end
end

