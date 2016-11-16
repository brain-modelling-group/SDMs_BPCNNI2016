function SDMs_examples_NMM()
% Example time series for the neural mass model, with:
% A,B: no noise
% C,D: additive noise
% E,F: multiplicative noise
%
% Generates Figure 2 of Roberts, Friston, Breakspear (2016)

% Incorporates code written by Stewart Heitmann and Matthew Aburn


rng(4) % random seed used in the paper

% NMM params
V1 = -0.01; V2 = 0.15; V3 = 0; V4 = 0.3; V5 = 0; V7 = 0; V9 = 0.3; V8 = 0.15;
gCa = 1; gK = 2.0; gL = 0.5; gNa = 6.7;
VK = -0.7; VL = -0.5; I = 0.3; b = 0.1; phi = 0.7; VNa = 0.53; VCa = 1;
ani = 0.4; aei = 2; aie = 2; aee = 0.36; ane = 1; rnmda = 0.25;
V6 = 0.65;


figure('position',[560   528-250   625   420+250],'paperpositionmode','auto');
axs=zeros(1,6);


% no noise
nse = 0; [t,v,w,z]=simulate_NMM(0);
axs(1)=opsubplot(3,2,1); plot(t,v,'b'), ylabel('V'), xlabel('t (ms)')
axs(2)=opsubplot(3,2,2); plot3byorientation(w,v,z), xlabel('W'), ylabel('V'), zlabel('Z')

% additive noise
I=0.32; nse=0.01; [t,v,w,z]=simulate_NMM(1);
axs(3)=opsubplot(3,2,3); plot(t,v,'b'), ylabel('V'), xlabel('t (ms)')
axs(4)=opsubplot(3,2,4); plot3byorientation(w,v,z), xlabel('W'), ylabel('V'), zlabel('Z')

% multiplicative noise
I=0.32; nse=0.053; [t,v,w,z]=simulate_NMM(2);
axs(5)=opsubplot(3,2,5); plot(t,v,'b'), ylabel('V'), xlabel('t (ms)')
axs(6)=opsubplot(3,2,6); plot3byorientation(w,v,z), xlabel('W'), ylabel('V'), zlabel('Z')


set(axs(2:2:6),'ylim',[-0.5 0.5])
set(axs(2),'zlim',[0.07 0.22])
set(axs([4 6]),'zlim',[0.1 0.25])
for j=2:2:6, view(axs(j),100,-8), end

% plot labeling
letters=upper({'a','b','c','d','e','f'});
textopts={'units','normalized','fontweight','bold'};
for j=1:2:5
    set(axs(j),'position',get(axs(j),'position')+[0.02 0 0 0])
    text(-0.25,1.05,letters{j},'parent',axs(j),textopts{:})
end
for j=2:2:6
    text(-0.25,1.05,letters{j},'parent',axs(j),textopts{:})
end


    function [t,v,w,z]=simulate_NMM(whichnoise)
        % simulate NMM SDEs
        %
        % whichnoise=0 ==> no noise
        % whichnoise=1 ==> additive noise
        % whichnoise=2 ==> multiplicative noise
        
        if nargin<1
            whichnoise=0;
        end
        
        % initial conditions
        y0 = [-0.2; 0.3; 0.12];
        
        noisefun=@(t,y) zeros(3,1);
        if whichnoise==1
            noisefun=@g;
        elseif whichnoise==2
            noisefun=@g_mult;
        end
        
        % Integrate SDE using Heun algorithm
        trange=0:0.1:2000; sol = Heun(@f, noisefun, trange, y0, false);
        
        % extract the values of interest from the final results
        t=sol.x;
        v=sol.y(1,:);     % membrane potential
        w=sol.y(2,:);     %
        z=sol.y(3,:);

        % deterministic part
        function ret = f(t, y)
            ret=nrlmass_dde(t,y);
        end
        
        % additive noise
        function ret = g(t, y)
            ret  = [ane.*nse; 0; 0]; % noise only in V
        end
        
        % multiplicative noise
        function ret = g_mult(t, y)
            ret  = [ane.*nse*y(1); 0; 0]; % noise only in V
        end
        
        function out=nrlmass_dde(t,y)
            % init 'out' and vectors for the variable types
            out = zeros(3,1);
            v = y(1);
            w = y(2);
            z = y(3);
            
            % calculate the inter-node inputs using transfer function
            gainv = gain(v,V5,V6,0.5);

            out(1) = -(gCa+rnmda.*aee.*gainv).*gain(v,V1,V2,0.5).*(v-VCa) ...
                -gK.*w.*(v-VK) ...
                -gL.*(v-VL) ...
                -(gNa.*gain(v,V9,V8,0.5)+aee.*gainv).*(v-VNa) ...
                +ane.*I ...
                -gain(z,V7,V6,aie.*0.5).*z;
            
            out(2) = phi.*(gain(v,V3,V4,0.5)-w);
            
            out(3) = b.*(ani.*I+gain(v,V5,V6,aei.*0.5).*v);
        end
        
        function f=gain(VAR,C1,C2,C3)
            %nonlinear gain function for coupled ODE neural lump
            f=C3*(1+tanh((VAR-C1)./C2));
        end
    end
end



% -------------------------------------------------------------------------
function plot3byorientation(x,y,z)
% by 2d orientation in the x-y plane
c = atan2(diff(y)/abs(max(diff(y))),diff(x)/abs(max(diff(x))));
C=[c(:) c(:)];

x=x(1:end-1); y=y(1:end-1); z=z(1:end-1);
surface([x(:),x(:)],[y(:),y(:)],[z(:),z(:)], C, ...
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

