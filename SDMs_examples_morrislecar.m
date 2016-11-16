function SDMs_examples_morrislecar()
% Example time series for the Morris-Lecar model, with:
% A,B: no noise
% C,D: additive noise
% E,F: multiplicative noise
%
% Generates Figure 1 of Roberts, Friston, Breakspear (2016)

% Incorporates code written by Stewart Heitmann and Matthew Aburn


rng(5) % random seed used in the paper

% constants (See Breakspear & Jirsa, pg 29).
cm = 1;        % membrane capacitance per unit area (nF/cm^2)

EL = -80;      % equilibrium potential of passive leaky membrane (mV)
gL = 8;        % conductance of passive leaky membrane (mS/cm^2)

ENa = 60;      % equilibrium potential of Na ion channel (mV)
gNa = 20;      % maximal conductance of Na ion channel (mS/cm^2)
Vm  = -20;     % threshold potential of the Na activation gate (mV)
km  = 15;      % variance of the Na activation gate threshold potential
max_m = 1;

EK_hi = -90;      % equilibrium potential of K+ ion channel (mV)
gK_hi = 10;       % maximal conductance of K+ ion channel (mS/cm^2)
Vn_hi = -25;      % threshold potential of the K activation gate (mV)

EK = EK_hi;
gK = gK_hi;
Vn = Vn_hi;
kn = 5;        % variance of the K activation gate threshold potential
max_n = 1;
tau_n = 1;


figure('position',[560   528-200   625   420+200],'paperpositionmode','auto');
axs=zeros(1,6);


% no noise
inj=20;
[t,V,n]=simulate(inj);
axs(1)=opsubplot(3,2,1); plot(t,V,'b'), ylabel('V (mV)'), xlabel('t (ms)')
axs(2)=opsubplot(3,2,2); plotbyorientation(n,V), ylabel('V (mV)'), xlabel('n')

% additive noise
inj=20;
sigma=10; % noise amplitude
[t,V,n]=simulate_sde(inj,sigma);
axs(3)=opsubplot(3,2,3); plot(t,V,'b'), ylabel('V (mV)'), xlabel('t (ms)')
axs(4)=opsubplot(3,2,4); plotbyorientation(n,V), ylabel('V (mV)'), xlabel('n')

% multiplicative noise
inj=20;
sigma=0.5; % noise amplitude
[t,V,n]=simulate_sde(inj,sigma,2);
axs(5)=opsubplot(3,2,5); plot(t,V,'b'), ylabel('V (mV)'), xlabel('t (ms)')
axs(6)=opsubplot(3,2,6); plotbyorientation(n,V), ylabel('V (mV)'), xlabel('n')


set(axs(1:2:5),'xlim',[0 70],'ylim',[-120 30])
set(axs(2:2:6),'xlim',[0 0.65],'ylim',[-120 30],'box','on')

% plot labeling
letters=upper({'a','b','c','d','e','f'});
for j=1:6
    set(axs(j),'position',get(axs(j),'position')+[0.03 0 0 0])
    text(-0.27,1.05,letters{j},'parent',axs(j),'units','normalized','fontweight','bold')
end


% noise-free (ODE) case

% This function uses an ODE solver to evaluate the temporal evolution of
% the model. The "inject" parameter controls the amplitude of
% the injection current. The timing of the injection pulses are hard-coded
% into this function by piecewise application of the ODE solver.
    function [t,V,n] = simulate(inject)
        
        % initial values.
        V0 = -66.5;                    % resting potential
        n0 = sigmoid(max_n,V0,Vn,kn);   % steady-state prob of open K activation gate when V=V0
        
        % Here we apply the ODE solver piecewise for each injection pulse.
        % ie = Injection current per unit surface area of membrane (see odefun below)
        ie = 0;      [ t1, X1] = ode15s(@odefun,[  0   5],[V0 n0]);          % no injection current
        ie = inject; [ t2, X2] = ode15s(@odefun,[  5  70], X1(end,:));       % injection current
        
        % collate the results
        t = [t1 ; t2];
        X = [X1 ; X2];
        
        % extract the values of interest from the final results
        V = X(:,1);     % membrane potential
        n = X(:,2);     % activation variable of K conductance channel
        
        % reconstruct m from V (remember that m is m_infinity in this model)
        m = sigmoid(max_m,V,Vm,km);
        
        
        % This function is repeatedly called by the ODE solver to evaluated the
        % evolution of the planar model over time.
        function dXdt = odefun(t,X)
            
            % extract our variables from parameter X
            V = X(1);      % membrane potential
            n = X(2);      % activation variable of K conductance channel
            
            % The planar model assumes the Na activation channel
            % instantaneously assumes its steady state value (m infinity)
            m_inf = sigmoid(max_m,V,Vm,km);
            
            % determine the transmembrane current
            im = gL*(V-EL) + gK*n*(V-EK) + gNa*m_inf*(V-ENa);
            
            % determine the change in membrane potential
            dVdt = (ie - im)/cm;
            
            % determine the change (dn/dt) in the (slow) K activation gate
            n_inf = sigmoid(max_n,V,Vn,kn);
            dndt = (n_inf - n)/tau_n;
            
            % return the results
            dXdt = [dVdt; dndt];
        end
    end


% noise cases
    function [t,V,n]=simulate_sde(inject,sigma,whichnoise,tend)
        % simulate Morris-Lecar SDEs
        %
        % whichnoise=1 ==> additive noise
        % whichnoise=2 ==> multiplicative noise
        
        if nargin<3
            whichnoise=1;
        end
        if nargin<4
            tend=70;
        end
        
        
        % initial conditions
        V0 = -66.5;                    % resting potential
        n0 = sigmoid(max_n,V0,Vn,kn);   % steady-state prob of open K activation gate when V=V0
        y0 = [V0 n0]';
        
        noisefun=@g;
        if whichnoise==2
            noisefun=@g_mult;
        end
        
        % Integrate SDE using Heun algorithm
        ton=5;
        dt=0.01;
        ie = 0; trange=0:dt:ton; sol = Heun(@f, noisefun, trange, y0, false); X1=sol.y;
        ie = inject; trange=ton:dt:tend; sol = Heun(@f, noisefun, trange, X1(:,end), false); X2=sol.y;
        X = [X1' ; X2(:,2:end)'];
        % extract the values of interest from the final results
        V = X(:,1);     % membrane potential
        n = X(:,2);     % activation variable of K conductance channel
        t=0:dt:tend;
        
        % deterministic part
        function ret = f(t, y)
            ret=sdefun_f(t,y);
        end
        
        % additive noise
        function ret = g(t, y)
            ret  = [sigma; 0]; % noise only in V
        end
        
        % multiplicative noise
        function ret = g_mult(t, y)
            ret  = [sigma*y(1); 0]; % noise only in V
        end
        
        % This function is repeatedly called by the DE solver to evaluate the
        % evolution of the planar model over time.
        function dXdt = sdefun_f(t,X)
            
            % extract our variables from parameter X
            V = X(1);      % membrane potential
            n = X(2);      % activation variable of K conductance channel
            
            % The planar model assumes the Na activation channel
            % instantaneously assumes its steady state value (m infinity)
            m_inf = sigmoid(max_m,V,Vm,km);
            
            % determine the transmembrane current
            im = gL*(V-EL) + gK*n*(V-EK) + gNa*m_inf*(V-ENa);
            
            % determine the change in membrane potential
            dVdt = (ie - im)/cm;
            
            % determine the change (dn/dt) in the (slow) K activation gate
            n_inf = sigmoid(max_n,V,Vn,kn);
            dndt = (n_inf - n)/tau_n;
            
            % return the results
            dXdt = [dVdt; dndt];
        end
     
    end


% sigmoid function as used to generate steady-state values of m,n
    function y = sigmoid(ymax,V,Vion,kion)
        y = ymax ./ ( 1 + exp((Vion-V)/kion) );
    end

end




% -------------------------------------------------------------------------
function plotbyorientation(x,y)
c = atan2(diff(y)/abs(max(diff(y))),diff(x)/abs(max(diff(x))));
x=x(1:end-1); y=y(1:end-1); z=zeros(size(x));
surface([x(:),x(:)],[y(:),y(:)],[z(:),z(:)],[c(:),c(:)], ...
    'EdgeColor','flat','FaceColor','none');
colormap(hsv(numel(x)))
end


function out=opsubplot(nrows,ncols,ind)
%  h=opsubplot(nrows,ncols,[ind])
% Create subplots with maximal outerpositions in a grid nrows-by-ncols,
% starting from the top left.
%
% If specifying index ind (like in subplot), give just that axis. If not,
% create all nrows*ncols axes and return their handles (row-wise) in a
% vector.

if nargin<3
    ind=1:nrows*ncols;
end

naxes=length(ind);
h=zeros(naxes,1);

for j=1:naxes
    h(j)=createonesubplot(nrows,ncols,ind(j));
end

if nargout>0
    out=h;
end

end

function h=createonesubplot(nrows,ncols,ind)
% create one subplot with maximal outerposition within its grid
[col,row]=ind2sub([ncols nrows],ind); % ind2sub works column-wise
h=axes('outerposition',[(col-1)/ncols (nrows-row)/nrows 1/ncols 1/nrows]);
end