function q = chern(lam,xi0,kFa,eps0,Nc,ksteps)
% Returns the Chern number q

del = 1;

xvec = 1:Nc;
yvec = 1:Nc;
[xij,yij] = meshgrid(xvec,yvec);
rij = bsxfun(@hypot,xij,yij);
invr = rij.^-1;

Nm = 1 + lam/sqrt(1+lam^2);
Np = 1 - lam/sqrt(1+lam^2);
kFm = kFa*(sqrt(1+lam^2) + lam);
kFp = kFa*(sqrt(1+lam^2) - lam);

I1m = @(rij) Nm/del*real( besselj(0,(kFm+1i/xi0)*rij) + 1i*StruveH0((kFm+1i/xi0)*rij) );
I1p = @(rij) Np/del*real( besselj(0,(kFp+1i/xi0)*rij) + 1i*StruveH0((kFp+1i/xi0)*rij) );
I0m = @(rij) Nm*real( 1i*besselj(1,(kFm+1i/xi0)*rij) + StruveH_1((kFm+1i/xi0)*rij) );
I0p = @(rij) Np*real( 1i*besselj(1,(kFp+1i/xi0)*rij) + StruveH_1((kFp+1i/xi0)*rij) );

A0 = -eps0; %the element at the origin
A = del^2/2*(I1m(rij) + I1p(rij));
Ax = del^2/2*(I1m(xvec) + I1p(xvec));
Ay = del^2/2*(I1m(yvec.') + I1p(yvec.'));
B = del/2*(I0m(rij) - I0p(rij)).*xij.*invr;
Bx = del/2*(I0m(xvec) - I0p(xvec)); % y=0 and x assumed to be positive
C = del/2*(I0m(rij) - I0p(rij)).*yij.*invr;
Cy = del/2*(I0m(yvec.') - I0p(yvec.')); % x=0 and y assumed to be positive



% Fourier summation. The matrix product produces a meshgrid of Fourier sums
% for the vectors kx and ky.

kx = linspace(0,pi,ksteps);
ky = linspace(0,pi,ksteps);
dk = (pi-0)/(ksteps-1);

coskx = cos(xvec.'*kx);
cosky = cos(ky.'*yvec);
sinkx = sin(xvec.'*kx);
sinky = sin(ky.'*yvec);

dz = A0 + 2*ones(ksteps,1)*Ax*coskx + 2*cosky*Ay*ones(1,ksteps) ...
    + 4*cosky*A*coskx;
dy = 2*ones(ksteps,1)*Bx*sinkx + 4*cosky*B*sinkx;
dx = 2*sinky*Cy*ones(1,ksteps) + 4*sinky*C*coskx;

dxdkx = -4*sinky*(C.*xij)*sinkx;
dxdky = 2*cosky*(Cy.*yvec.')*ones(1,ksteps) + 4*cosky*(C.*yij)*coskx;
dydkx = 2*ones(ksteps,1)*(Bx.*xvec)*coskx + 4*cosky*(B.*xij)*coskx;
dydky = -4*sinky*(B.*yij)*sinkx;
dzdkx = -2*ones(ksteps,1)*(Ax.*xvec)*sinkx - 4*cosky*(A.*xij)*sinkx;
dzdky = -2*sinky*(Ay.*yvec.')*ones(1,ksteps) - 4*sinky*(A.*yij)*coskx;

qk = ( dx.*(dydkx.*dzdky - dzdkx.*dydky) + dy.*(dzdkx.*dxdky - dxdkx.*dzdky) ...
    + dz.*(dxdkx.*dydky - dydkx.*dxdky) ).*(dz.^2 + dy.^2 + dx.^2).^(-3/2);

q = sum(sum( qk./pi ))*dk^2;

