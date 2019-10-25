function [C] = Draw_39_Map(DEF_S,DEF_R)
% DRAW_39_MAP    Draw the DEF map for the 39 bus system
%
% Inputs:
% 1) DEF_S       Sending end DEF data
% 2) DEF_R       Receiving end DEF data
%
% Outputs:
% 1) C           Colorbar handle
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
initpsat;                           % Initialize PSAT global variables
datafile = 'NE39';                  % Test case data file
runpsat(datafile,'data');           % Initialize datafile
Settings.freq = 60;                 % Change System Freq from default to 60
runpsat('pf');                      % Run power flow

% Process DEF data
DEF_S_slope = abs(DEF_S);    % Slope of DEF integral
DEF_S_sign  = sign(DEF_S);   % Direction of DEF integral
DEF_R_slope = abs(DEF_R);
DEF_R_sign  = sign(DEF_R);

% Get Useful Power Flow Results
Pgs = Bus.Pg;
Qgs = Bus.Qg;
Pls = Bus.Pl;
Qls = Bus.Ql;

% Gens and Load Buses
PQs  = 1:29;
Gens = 30:39;

% Re/Active Power Injection
Sl = sqrt(Pls.^2 + Qls.^2);
Sg = sqrt(Pgs.^2 + Qgs.^2);
Si = Sg - Sl;

% Determine facecolor
max_Si = max(Si);
min_Si = min(Si);
max_Si = max(abs(max_Si),abs(min_Si));

% Color Map (256 points)
p     = 128;

% Original Color Map
% % v     = fliplr(linspace(0.25,1,p)) .';
% % v1    = [ones(size(v)); v];
% % v2    = [flipud(v); v];
% % v3    = [flipud(v); ones(size(v))];
% % mymap = horzcat(v1,v2,v3);

% Load nice color map from (https://www.kennethmoreland.com/color-advice/)
load('Color_Warm')
mymap = Color_Warm2;

% Assign facecolor index based on where it falls in this range
color_vec = zeros(length(Si),3);
for ii = 1:length(Si);
    Si_val = Si(ii);
    frac   = Si_val/max_Si;
    val    = round(frac*p);
    val    = p - val; 
    if val < 0                 % This is a load
        color_vec(ii,:) = mymap(val,:);
    else
        color_vec(ii,:) = mymap(2*p-val,:);
    end
end

% Take the Slope and Normalize such that the largest width is 7
DEF_S_slope = 7*DEF_S_slope/max(DEF_S_slope);
DEF_R_slope = 7*DEF_R_slope/max(DEF_R_slope);

% Set minimum width
DEF_S_slope(DEF_S_slope<0.25) = 0.25;
DEF_R_slope(DEF_R_slope<0.25) = 0.25;

%% Draw System Map
crd_Data = csvread('39_Bus_Coords.csv',1);
bus = crd_Data(:,1);
X   = crd_Data(:,2);
Y   = crd_Data(:,3);

% Call Line connections
fr_bus = Line.con(:,1);
to_bus = Line.con(:,2);
n      = length(bus);

%% Prepare the Figure
clf; hold on;
set(gcf,'Color','w'); % set bg to white
axis off;

%% Draw the Branches
for ii = 1:length(fr_bus)
    
    % Average Flow strength and direction
    flow1 = DEF_S_sign(ii);
    flow2 = DEF_R_sign(ii);
    f     = fr_bus(ii);
    t     = to_bus(ii);
    x     = [X(f) X(t)];
    y     = [Y(f) Y(t)];
    
    % We want branch size to be based on DEF slope (average both sides)
    width = (DEF_S_slope(ii) + DEF_R_slope(ii));
    
    % In drawing the line, we need to know which direction the DEF is
    % flowing. If flow == 1, it is flowing "f" => "t", otherwise "t" => "f"
    drawline(x,y,width);
end

% Draw Buses (Equisized Rectangles)
for ii = PQs
    x = X(ii);
    y = Y(ii);
    h = 20;
    w = 30;
    
    % Choose facecolor based on apparent power injection
    facecolor = color_vec(ii,:);
    
    % Call Function
    drawrectangle(x,y,h,w,facecolor)
end

% Draw Now!
drawnow

% Draw the Generators
for ii = Gens
    gx = X(ii);
    gy = Y(ii);
    diameter = 60;
    % Choose facecolor based on apparent power injection
    facecolor = color_vec(ii,:);
    
    % Cal Function
    drawcircle(gx,gy,diameter,facecolor);
end

% Draw Arrows
DEF_S_slope = 10*DEF_S_slope;
DEF_R_slope = 10*DEF_R_slope;
scale_fac   = max([DEF_S_slope; DEF_R_slope]);
DEF_S_slope = 35*DEF_S_slope/scale_fac;
DEF_R_slope = 35*DEF_R_slope/scale_fac;
for ii = 1:length(fr_bus)
    
    % Average Flow strength and direction
    flow1 = DEF_S_sign(ii);
    flow2 = DEF_R_sign(ii);
    f     = fr_bus(ii);
    t     = to_bus(ii);
    x     = [X(f) X(t)];
    y     = [Y(f) Y(t)];
    width1 = DEF_S_slope(ii);
    width2 = DEF_R_slope(ii);
    drawarrows(x,y,width1,width2,flow1,flow2);
end

% Draw Now!
drawnow

% Add a colorbar with tick labels
colormap(mymap)
tick = [1, (500/4), 250, (3*500/4) 500];
C  = colorbar('Location', 'east','Ticks',tick,'YTickLabel',{'-8.3','-4.15','0','4.15','8.3',},'AxisLocation','out');
ap = C.Position;
C.Label.String      = 'Apparent Power Injection';
C.Label.Interpreter = 'latex';
C.Label.FontSize = 13;

% Move the colorbar
set(C,'Position',[ap(1)-0.125 ap(2) ap(3) ap(4)])

%% Show Bus Number
for ii=1:n
    fs = 9;
    text(X(ii),Y(ii),num2str(bus(ii,1)),'HorizontalAlignment','center','FontSize',fs,'FontWeight','bold');
end

% Move x
mv = 100;

% Receive
xr = [760  799.1]+5;
yr = [920   920];
dg = [0 0.6 0];
arrow_width2 = 25;
davinci('arrow','X',xr+mv,'Y',yr,'Head.Width',arrow_width2,'Head.Length',arrow_width2,'Head.Sweep',5,'Color',dg,'Shaft.Width',5);

% Sending
xs = [760  799.1]+5;
ys = [960   960];
dg = [0, 0.4470, 0.7410];
arrow_width1 = 25;
davinci('arrow','X',xs+mv,'Y',ys,'Head.Width',arrow_width1,'Head.Length',arrow_width1,'Head.Sweep',5,'Color',dg,'Shaft.Width',5);

% Add legend text
text(540+mv,960,'${\rm Sending \; End \;} P^{\star}$','Interpreter','latex','fontsize',13)
text(517+mv,920,'${\rm Receiving \; End \;} P^{\star}$','Interpreter','latex','fontsize',13)

w = 150;
h = 40;
rectangle('Position',[612 903 2*w 2*h],'LineWidth' ,0.5,'EdgeColor','black')

% Draw Now!
drawnow


%% Keep the x and y proportional
axis equal;
caxis([0 500]);

return

% Draw a line
function drawline(x,y,width)
dx = x(2) - x(1);
dy = y(2) - y(1);
w  = .0001;
ortho = [-dy dx]/sqrt(dx^2 + dy^2);
p = [[x(1) y(1)] + ortho*w;
     [x(1) y(1)] - ortho*w;
     [x(2) y(2)] - ortho*w;
     [x(2) y(2)] + ortho*w;];

% First draw a line
line(p(:,1),p(:,2),'LineWidth',width,'color','k');

function drawarrows(x,y,width1,width2,flow1,flow2)
xs = x;
ys = y;
xr = x;
yr = y;

% Test Direction (S)
if flow1 ~= 1 % DEF is "f" to "t"
    temp  = xs(1);
    xs(1) = xs(2);
    xs(2) = temp;
    
    temp  = ys(1);
    ys(1) = ys(2);
    ys(2) = temp;
    
    % Differences
    dxs = xs(2) - xs(1);
    dys = ys(2) - ys(1);
    
    xs = [xs(1) (3.8*dxs/5)+xs(1)];
    ys = [ys(1) (3.8*dys/5)+ys(1)];
    
else
    % Differences
    dxs = xs(2) - xs(1);
    dys = ys(2) - ys(1);
    
    xs = [xs(1) (1.7*dxs/5)+xs(1)];
    ys = [ys(1) (1.7*dys/5)+ys(1)];
end

% Test Direction (R)
if flow2 ~= 1 % DEF is "f" to "t"
    temp  = xr(1);
    xr(1) = xr(2);
    xr(2) = temp;
    
    temp  = yr(1);
    yr(1) = yr(2);
    yr(2) = temp;
    
    % Differences
    dxr = xr(2) - xr(1);
    dyr = yr(2) - yr(1);
    
    xr = [xr(1) (1.7*dxr/5)+xr(1)];
    yr = [yr(1) (1.7*dyr/5)+yr(1)];
else
    % Differences
    dxr = xr(2) - xr(1);
    dyr = yr(2) - yr(1);

    xr = [xr(1) (3.8*dxr/5)+xr(1)];
    yr = [yr(1) (3.8*dyr/5)+yr(1)];
end


% Scale arrow width at a slower rate
arrow_width1 = width1;
arrow_width2 = width2;

if arrow_width1 < 10
    arrow_width1 = 10;
end

if arrow_width2 < 10
    arrow_width2 = 10;
end

% Draw first arrow (S)
dg = [0, 0.4470, 0.7410];
davinci('arrow','X',xs,'Y',ys,'Head.Width',arrow_width1,'Head.Length',arrow_width1,'Head.Sweep',5,'Color',dg,'Shaft.Width',0.0001);

% Draw second arrow (R)
dg = [0 0.6 0];
davinci('arrow','X',xr,'Y',yr,'Head.Width',arrow_width2,'Head.Length',arrow_width2,'Head.Sweep',5,'Color',dg,'Shaft.Width',0.0001);

% Draw a Circle
function drawcircle(x,y,diameter,facecolor)
angles = (0:pi/(12*4):2*pi)';
r = diameter/2;
dy = sin(angles)*r;
dx = cos(angles)*r;
patch( x+dx, y+dy, 'white','LineWidth' ,0.5,'EdgeColor','black','FaceColor',facecolor);

% Draw a Rectangle
function drawrectangle(x,y,h,w,facecolor)
rectangle('Position',[x-w y-h 2*w 2*h],'LineWidth' ,0.5,'EdgeColor','black','FaceColor',facecolor,'Curvature',[0.5,0.5])

