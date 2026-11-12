%% Setup surface, nodes, and radial kernel
% Red blood cell (RBC) surface
r0=3.39; c0=0.81/r0; c2=7.83/r0; c4=-4.39/r0; a=3.91/r0;
rbc = @(la,th) [a*cos(la).*cos(th) a*sin(la).*cos(th) ...
   0.5*sin(th).*(c0+c2*(cos(th)).^2+c4*(cos(th)).^4)];
syms la th; % Compute the surface normal vectors symbolically
nr = inline(cross(diff(rbc(la,th),la),diff(rbc(la,th),th)));
% Use minimum energy nodes projected to the sphere
x=load('me63.4096'); x=x(:,1:3); N=size(x,1);
[la,th]=cart2sph(x(:,1),x(:,2),x(:,3)); x=rbc(la,th); 
nr=nr(la,th); nr=nr./repmat(sqrt(sum(nr.^2,2)),[1 3]);
% Radial kernel
ep=4;  % Shape parameter
phi=@(r2) 1./sqrt(1+ep^2*r2);         % IMQ
dphi=@(r2) -ep^2./sqrt(1+ep^2*r2).^3; % Derivative of IMQ over r

%%  Compute the Surface Laplacian - not specific to RBC surface
xij=repmat(x(:,1),[1 N]); xij=xij-xij.'; nxi=repmat(nr(:,1),[1 N]);
yij=repmat(x(:,2),[1 N]); yij=yij-yij.'; nyi=repmat(nr(:,2),[1 N]);
zij=repmat(x(:,3),[1 N]); zij=zij-zij.'; nzi=repmat(nr(:,3),[1 N]);
r2=xij.^2 + yij.^2 + zij.^2; A=dphi(r2);
DPx=((1-nxi.^2).*xij - nxi.*nyi.*yij - nxi.*nzi.*zij).*A;
DPy=(-nxi.*nyi.*xij + (1-nyi.^2).*yij - nyi.*nzi.*zij).*A;
DPz=(-nxi.*nzi.*xij - nyi.*nzi.*yij + (1-nzi.^2).*zij).*A;
A=chol(phi(r2)); DPx=(DPx/A)/A.'; DPy=(DPy/A)/A.'; DPz=(DPz/A)/A.';
Lap=DPx*DPx+DPy*DPy+DPz*DPz;  % Surface Laplacian

%% Turing spot pattern example using SBDF2
del=0.0045; d=0.516; tau1=0.02; tau2=0.2; alp=0.899; bet=-0.91; gam=-alp;
ufun=@(u,v) alp*u.*(1-tau1*v.^2)+v.*(1-tau2*u);
vfun=@(u,v) bet*v.*(1+alp*tau1/bet*u.*v)+u.*(gam+tau2*v);
% Initial condition
stream=RandStream('mrg32k3a','seed',7122005); u=rand(stream,N,2)-0.5; 
id=find(abs(x(:,3))>=0.1); u(id,:)=0; v=u(:,2); u=u(:,1);
% Implicit systems for SBDF2
dt=0.05; tfinal=400;
Du=1.5*eye(N)-del*d*dt*(Lap); [Lu,Uu,pu]=lu(Du,'vector');
Dv=1.5*eye(N)-del*dt*(Lap);   [Lv,Uv,pv]=lu(Dv,'vector');
optsu.UT=true; optsl.LT=true;
% One step of backward Euler to bootstrap SBDF2
rhsu=u+dt*ufun(u,v); rhsv=v+dt*vfun(u,v);
u0=u; u=(eye(N)-del*d*dt*Lap)\rhsu;
v0=v; v=(eye(N)-del*dt*Lap)\rhsv;
for j=1:tfinal/dt
   rhsu=2*u-0.5*u0+dt*(2*ufun(u,v)-ufun(u0,v0)); % SBDF2
   rhsv=2*v-0.5*v0+dt*(2*vfun(u,v)-vfun(u0,v0));
   u0=u; u=linsolve(Uu,linsolve(Lu,rhsu(pu),optsl),optsu);
   v0=v; v=linsolve(Uv,linsolve(Lv,rhsv(pv),optsl),optsu);
end

%% Interpolate the solution to a grid on the surface using RBFs and plot
sz=[101 201]; M = prod(sz);  % surface grid parameters
[ll,tt]=meshgrid(linspace(-pi,pi,sz(2)),linspace(-pi/2,pi/2,sz(1)));
xx=rbc(ll(:),tt(:));
% Interpolate to the grid
re2=(repmat(xx(:,1),[1 N])-repmat(x(:,1).',[M 1])).^2;
re2=re2+(repmat(xx(:,2),[1 N])-repmat(x(:,2).',[M 1])).^2;
re2=re2+(repmat(xx(:,3),[1 N])-repmat(x(:,3).',[M 1])).^2;
uu=reshape(phi(re2)*(A\(A.'\u)),sz); 
% Plot the results
yy=reshape(xx(:,2),sz); zz=reshape(xx(:,3),sz); xx=reshape(xx(:,1),sz);
surf(xx,yy,zz,uu); shading interp; daspect([1 1 1]); axis tight;