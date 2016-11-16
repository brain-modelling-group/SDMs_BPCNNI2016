function SDMs_examples_bistability()
% Illustration of bistability in a simple oscillator model with noise.
%
% Model of an oscillator in polar coords, with additive noise:
%
%   dr = (-r^5 + lambda*r^3 + beta*r)*dt + sigma*dW
%   dphi = theta*dt
%
% Generates Figure 4 of Roberts, Friston, Breakspear (2016) 

% Incorporates code written by Matthew Aburn and Frank Freyer

rng(1) % random seed used in the paper

figure('position',[560   528-350   560   420+350],'paperpositionmode','auto')
axs=zeros(1,6);


% supercritical bifurcation diagram
lambda=-8;
bb=linspace(-15,15,5000);
r_stab_f=@(beta) sqrt(2*lambda+2*sqrt(lambda.^2+4*beta))/2;
r_stab=r_stab_f(bb);
r_stab(imag(r_stab)~=0)=NaN;
axs(1)=axes('position',[0.1089    0.8000    0.7750    0.1577]);
plot([min(bb),0],[0 0],'k',...
     [0,max(bb)],[0 0],'k:',...
     bb,r_stab,'g')
ylim([-0.3 2])
xlabel('\beta')
ylabel('r')


% subcritical bifurcation diagram, with bistability
lambda=8;
bb=linspace(-20,10,10000);
r_unstab_f=@(beta) sqrt(2*lambda-2*sqrt(lambda.^2+4*beta))/2;
r_stab_f=@(beta) sqrt(2*lambda+2*sqrt(lambda.^2+4*beta))/2;
r_unstab=r_unstab_f(bb);
r_stab=r_stab_f(bb);
r_unstab(imag(r_unstab)~=0)=NaN;
r_stab(imag(r_stab)~=0)=NaN;
axs(2)=axes('position',[0.1089    0.5700    0.7750    0.1577]);
plot([min(bb),0],[0 0],'k',...
     [0,max(bb)],[0 0],'k:',...
     bb,r_unstab,'r',...
     bb,r_stab,'g',...
     -14*[1 1],[0 r_stab_f(-14)],'b',...
     -10*[1 1],[0 r_stab_f(-10)],'b')
ylim([-0.3 3.5])
xlabel('\beta')
ylabel('r')
text(-15,0.6,'D')
text(-9.5,0.6,'F')


% example time series

% low-amplitude case
% parameters
theta=20;
beta=-14;
lambda=8;
sigma=1.2;

% initial condition
y0=[0;0];

% integrate SDE using the Heun algorithm
trange=0:0.01:50;
sol = Heun(@f, @g, trange, y0, false);

% extract the values of interest from the final results
t=sol.x;
r=sol.y(2,:);

% radii of separatrix and stable orbit, for plotting
sep=sqrt(2*lambda-2*sqrt(lambda.^2+4*beta))/2;
stab=sqrt(2*lambda+2*sqrt(lambda.^2+4*beta))/2;

% time series in cartesian coords
axs(3)=axes('position',[0.1089    0.3000    0.3793    0.1819]);
plot(t,r.*cos(theta*t),'b')
xlabel('t (s)'), ylabel('x')
ylim([-3 3]), xlim([0 max(t)])

% trajectory in phase space
axs(4)=axes('position',[0.5650    0.3000    0.3875    0.1819]);
plot(r.*cos(theta*t),r.*sin(theta*t),'b',...
     sep.*cos(theta*t),sep.*sin(theta*t),'r',...
     stab.*cos(theta*t),stab.*sin(theta*t),'g')
axis square
xlim([-3 3])
ylim([-3 3])
xlabel('x'), ylabel('y')


% high-amplitude case
% parameters
theta=20;
lambda=8;
beta=-10;
sigma=1.2;

% initial condition
y0=[0;0];

% integrate SDE using the Heun algorithm
trange=0:0.01:50; sol = Heun(@f, @g, trange, y0, false);

% extract the values of interest from the final results
t=sol.x;
r=sol.y(2,:);

% radii of separatrix and stable orbit, for plotting
sep=sqrt(2*lambda-2*sqrt(lambda.^2+4*beta))/2;
stab=sqrt(2*lambda+2*sqrt(lambda.^2+4*beta))/2;

% time series in cartesian coords
axs(5)=axes('position',[0.1089    0.0500    0.3793    0.1819]);
plot(t,r.*cos(theta*t),'b')
ylim([-3 3]), xlim([0 max(t)])
xlabel('t (s)'), ylabel('x')

% trajectory in phase space
axs(6)=axes('position',[0.5650    0.0500    0.3875    0.1819]); %#ok<NASGU>
plot(r.*cos(theta*t),r.*sin(theta*t),'b',sep.*cos(theta*t),sep.*sin(theta*t),'r',stab.*cos(theta*t),stab.*sin(theta*t),'g')
axis square
xlim([-3 3])
ylim([-3 3])
xlabel('x'), ylabel('y')


% plot labeling
letters=upper({'a','b','c','d','e','f'});
ax2=axes('position',[0 0 1 1],'visible','off');
textopts={'parent',ax2,'units','normalized','fontweight','bold'};
text(0.02,0.97,letters(1),textopts{:});
text(0.02,0.74,letters(2),textopts{:});
text(0.02,0.50,letters(3),textopts{:});
text(0.55,0.50,letters(4),textopts{:});
text(0.02,0.25,letters(5),textopts{:});
text(0.55,0.25,letters(6),textopts{:});


% the bistable Hopf SDE

    % deterministic term 
    function ret = f(t, y)
        ret = [theta; -y(2)^5+lambda*y(2)^3+beta*y(2)];
    end
    
    % noise term
    function ret = g(t, y)
        ret  = [0; sigma]; % noise only in r
    end

end