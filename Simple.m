%% Info
% Toerental: 1500 RPM
% SOA van 4.2º voor TDC
% Resolutie van 0.2º CA
% Data voor 69 cycles (maximale van de Smetec, de OGO gensets kunnen in principe “onbeperkt” aan)
% 
%% init
clearvars; clc;close all;
%% Units
mm      = 1e-3;dm=0.1;
bara    = 1e5;
MJ      = 1e6;
kWhr    = 1000*3600;
volperc = 0.01; % Emissions are in volume percentages
ppm     = 1e-6; % Some are in ppm (also a volume- not a mass-fraction)
g       = 1e-3;
s       = 1;
%% Load NASA maybe you need it at some point?
% Global (for the Nasa database in case you wish to use it).
global Runiv
Runiv = 8.314;
[SpS,El]        = myload(fullfile('NasaThermalDatabase.mat'),{'Diesel','O2','N2','CO2','H2O'});
%% Engine geom data (check if these are correct)
Cyl.Bore                = 104*mm;
Cyl.Stroke              = 85*mm;
Cyl.CompressionRatio    = 21.5;
Cyl.ConRod              = 136.5*mm;
Cyl.TDCangle            = 0;
% -- Valve closing events can sometimes be seen in fast oscillations in the pressure signal (due
% to the impact when the Valve hits its seat).
CaIVO = -355;
CaIVC = -135;
CaEVO = 149;
CaEVC = -344;
CaSOI = -3.2;
% Write a function [V] = CylinderVolume(Ca,Cyl) that will give you Volume
% for the given Cyl geometry. If you can do that you can create pV-diagrams
%% Fuel and data file
Fuel     = FuelMix();                                   % edit the table in Functions/FuelMix.m
FullName = fullfile('Data','ExampleDataSet.txt');       % measurement file for this fuel
dataIn          = table2array(readtable(FullName));
[Nrows,Ncols]   = size(dataIn);                    % Determine size of array
NdatapointsperCycle = 720/0.2;                     % Nrows is a multitude of NdatapointsperCycle
Ncycles         = Nrows/NdatapointsperCycle;       % This must be an integer. If not checkwhat is going on
Ca              = reshape(dataIn(:,1),[],Ncycles); % Both p and Ca are now matrices of size (NCa,Ncycles)
MeanCa          = mean(Ca,2);
p               = reshape(dataIn(:,2),[],Ncycles)*bara; % type 'help reshape' in the command window if you want to know what it does (reshape is a Matlab buit-in command
Meanp           = mean(p,2);
%% Plotting 
f1=figure(1);
set(f1,'Position',[ 200 800 1200 400]);             % Just a size I like. Your choice
pp = plot(MeanCa,Meanp/bara,'LineWidth',1);                 % Plots the whole matrix
xlabel('Ca');ylabel('p [bar]');                     % Always add axis labels
xlim([-360 360]);ylim([0 50]);                      % Matter of taste
iselect = 10;                                    % Plot cycle 10 again in the same plot to emphasize it. Just to show how to access individual cycles.
line(Ca(:,iselect),p(:,iselect)/bara,'LineWidth',2,'Color','r');
YLIM = ylim;
% Add some extras to the plot
line([CaIVC CaIVC],YLIM,'LineWidth',1,'Color','b'); % Plot a vertical line at IVC. Just for reference not a particular reason.
line([CaEVO CaEVO],YLIM,'LineWidth',1,'Color','r'); % Plot a vertical line at EVO. Just for reference not a particular reason.
set(gca,'XTick',[-360:60:360],'XGrid','on','YGrid','on');        % I like specific axis labels. Matter of taste
title('All cycles in one plot.')

%% pV-diagram
V = CylinderVolume(Ca(:,iselect),Cyl);
f2 = figure(2);
set(f2,'Position',[ 200 400 600 800]);              % Just a size I like. Your choice
subplot(2,1,1)
plot(V/dm^3,p(:,iselect)/bara);
xlabel('V [dm^3]');ylabel('p [bar]');               % Always add axis labels
xlim([0 0.8]);ylim([0.5 50]);                      % Matter of taste
set(gca,'XTick',[0:0.1:0.8],'XGrid','on','YGrid','on');        % I like specific axis labels. Matter of taste
title({'pV-diagram, cycle 10'})
subplot(2,1,2)
loglog(V/dm^3,p(:,iselect)/bara);
xlabel('V [dm^3]');ylabel('p [bar]');               % Always add axis labels
xlim([0.02 0.8]);ylim([0.5 50]);
set(gca,'XTick',[0.02 0.05 0.1 0.2 0.5 0.8],...
    'YTick',[0.5 1 2 5 10 20 50],'XGrid','on','YGrid','on');        % I like specific axis labels. Matter of taste
title('pV-diagram (log-log), cycle 10')
%% CA-V diagram
f3 = figure(3);
plot(Ca(:,iselect), V/dm^3)
xlabel('Crank Angle [deg]');ylabel('Volume [bar]');
title('Crank angle vs. Volume');




%% ===================== Engine and fuel results =====================
% --- Inputs (change per experiment)
RPM        = 1500;      % engine speed [rpm]
mdot_fuel  = NaN;       % measured fuel mass flow [kg/s]; NaN = estimate from heat release
gammaC     = 1.35;      % constant ratio of specific heats (first estimate)

% --- Volume and mean cycle
Vall  = CylinderVolume(Ca, Cyl);             % volume for all cycles [m^3]
Vd    = pi/4*Cyl.Bore^2*Cyl.Stroke;          % displacement volume [m^3]
Cam   = Ca(:,1);                             % crank angle of one cycle [deg]
pm    = mean(p, 2);                          % mean pressure over all cycles [Pa]
Vm    = CylinderVolume(Cam, Cyl);            % volume of one cycle [m^3]

% --- Work, IMEP and power
Wcyc = zeros(1, Ncycles);
for i = 1:Ncycles
    Wcyc(i) = trapz(Vall(:,i), p(:,i));      % closed integral of p dV [J]
end
W     = mean(Wcyc);   dW    = std(Wcyc);     % work per cycle and its scatter [J]
IMEP  = W/Vd;         dIMEP = dW/Vd;         % [Pa]
P     = W*RPM/120;                           % indicated power, 4-stroke [W]

% --- Apparent heat release (constant gamma), filtered mean cycle
pf    = smoothdata(pm, 'sgolay', 21);        % Savitzky-Golay, 21 points = 4.2 deg
dCa   = Cam(2) - Cam(1);
dpdCa = gradient(pf, dCa);                   % [Pa/deg]
dVdCa = gradient(Vm, dCa);                   % [m^3/deg]
aROHR = gammaC/(gammaC-1)*pf.*dVdCa + 1/(gammaC-1)*Vm.*dpdCa;   % [J/deg]

iw    = Cam >= CaIVC & Cam <= CaEVO;         % closed part of the cycle (valves shut)
Caw   = Cam(iw);
aHR   = cumtrapz(Caw, aROHR(iw));            % cumulative heat release [J]
[~, imax]      = max(aHR);
[aHRmin, imin] = min(aHR(1:imax));
Qhr   = aHR(imax) - aHRmin;                  % heat released by combustion [J]
frac  = (aHR - aHRmin)/Qhr;                  % burned fraction [-]
CA10  = Caw(find(frac(imin:imax) >= 0.10, 1) + imin - 1);
CA50  = Caw(find(frac(imin:imax) >= 0.50, 1) + imin - 1);
CA90  = Caw(find(frac(imin:imax) >= 0.90, 1) + imin - 1);

% --- Fuel
if isnan(mdot_fuel)
    mf = Qhr/Fuel.LHV;                       % estimate: all released heat comes from the fuel
    mdot_fuel = mf*RPM/120;
    fprintf('No fuel flow given: fuel mass estimated from heat release (lower bound).\n');
else
    mf = mdot_fuel*120/RPM;                  % fuel mass per cycle [kg]
end
Qfuel     = mf*Fuel.LHV;                     % chemical energy of the fuel per cycle [J]
eta       = W/Qfuel;                         % indicated efficiency [-]
Qloss     = Qfuel - W;                       % energy not converted into work [J]
mair_st   = mf*Fuel.AFRst;                   % stoichiometric air per cycle [kg]
mCO2      = mf*Fuel.wC*44.01/12.011;         % CO2 from complete combustion per cycle [kg]
mdot_CO2  = mCO2*RPM/120;                    % [kg/s]
isfc      = (mdot_fuel/P)/(g/kWhr);          % indicated specific fuel consumption [g/kWh]
isCO2     = (mdot_CO2/P)/(g/kWhr);           % indicated specific CO2 [g/kWh]

% --- Print
fprintf('\n=========== Results: %s ===========\n', Fuel.Name);
fprintf('Work per cycle        W      = %8.1f +/- %.1f J\n', W, dW);
fprintf('IMEP                         = %8.2f +/- %.2f bar\n', IMEP/bara, dIMEP/bara);
fprintf('Indicated power       P      = %8.2f kW\n', P/1000);
fprintf('Heat released         Qhr    = %8.1f J\n', Qhr);
fprintf('Fuel energy           Qfuel  = %8.1f J\n', Qfuel);
fprintf('Energy not to work    Qloss  = %8.1f J\n', Qloss);
fprintf('Fuel mass per cycle   mf     = %8.2f mg\n', mf*1e6);
fprintf('Stoich. air per cycle        = %8.1f mg\n', mair_st*1e6);
fprintf('Indicated efficiency  eta    = %8.1f %%\n', eta*100);
fprintf('isfc                         = %8.0f g/kWh\n', isfc);
fprintf('isCO2                        = %8.0f g/kWh\n', isCO2);
fprintf('CA10 / CA50 / CA90           = %.1f / %.1f / %.1f deg\n', CA10, CA50, CA90);


degtime = 1/(6*RPM)*1000;                     % duration of 1 deg crank angle [ms]
fprintf('Fuel use                     = %8.2f kg/h = %.2f L/h\n', mdot_fuel*3600, mdot_fuel*3600/Fuel.rho*1000);
fprintf('CO2 emitted                  = %8.2f kg/h\n', mdot_CO2*3600);
fprintf('Fuel per cycle               = %8.1f mg\n', mf*1e6);
fprintf('Of 100 J of fuel energy      = %8.0f J becomes work, %.0f J is lost\n', 100*eta, 100*(1-eta));
fprintf('Ignition delay (SOI to CA10) = %8.1f deg = %.2f ms\n', CA10 - CaSOI, (CA10 - CaSOI)*degtime);
fprintf('Burn duration (CA10 to CA90) = %8.1f deg = %.2f ms\n', CA90 - CA10, (CA90 - CA10)*degtime);
fprintf('CA50 after TDC               = %8.1f deg (handbook: should be closer to TDC)\n', CA50);

%% Functions
function [V] = CylinderVolume(Ca,Cyl)

% This function provides the cylinder volume as function of 
% Ca : Crankangle [degrees]
% Cyl :  a struct containing
%   Cyl.Stroke              : Stroke [m]
%   Cyl.Bore                : Bore [m]
%   Cyl.ConRod              : Connecting Rod length
%   Cyl.CompressionRatio    : Compession Ratio
%   Cyl.TDCangle            : Angle associated with the Top Dead Center

B   = Cyl.Bore;
S   = Cyl.Stroke;
cr  = Cyl.CompressionRatio;
r   = S/2;
l   = Cyl.ConRod;
%-------------------------------------------------------------------------------------------------------
Ap    = pi/4*B^2;
Vd    = Ap*S;
Vc    = Vd/(cr-1);
theta = Ca - Cyl.TDCangle;
V     = Vc + Ap*( r + l - ( r*cosd(theta) + sqrt(l^2 - r^2*sind(theta).^2) ) );
end
function Fuel = FuelMix()
% FuelMix  Properties of the fuel used in the experiment.
%   Edit ONLY the table below: set the volume fraction of each component.
%   - Fractions must add up to 1 (e.g. 0.5 B7 + 0.5 HVO for a 50/50 blend).
%   - To add a component, add a row. Components at 0 may keep NaN values.
%   Blend properties are mass-weighted averages (ideal mixing).
%   Sources: Project Handbook 4CBLB10 Table 2.1, group SSA 1 reports,
%   IEA-AMF (HVO composition). See the comment at the end of each row.

%   Name       VolFrac  LHV      rho      AFRst   wC      wH      wO
%              [-]      [MJ/kg]  [kg/m3]  [-]     [-]     [-]     [-]
% --- Fuels available in the lab (List of Available Fuels 4CBLB10)
T = {
    'B7',        1.0,   43.00,   836.1,   14.50,  0.858,  0.136,  0.006   % handbook Table 2.1
    'GTL',       0.0,   44.00,   777.1,   14.98,  0.847,  0.153,  0.000   % handbook Table 2.1
    'HVO',       0.0,   44.00,   782.5,   14.93,  0.850,  0.150,  0.000   % LHV, rho: HVO SSA (KIT reFuels, rho 765-800); C/H/O: IEA-AMF; AFRst calculated
    % --- Not available in the lab (for comparison only)
    'GTLB30',    0.0,   42.00,   808.3,   14.24,  0.824,  0.143,  0.032   % handbook Table 2.1
    'FAME',      0.0,   38.27,   NaN,     12.59,  0.770,  0.122,  0.108   % LHV: FAME SSA (Tesfa, RME B100); C/H/O: methyl oleate C19H36O2; rho not found
    'OME3-5',    0.0,   19.80,   1055.5,  6.38,   0.457,  0.090,  0.453   % LHV, rho, O: OMEx SSA (Holzer); C/H estimated from OME4; AFRst calculated
    'Propanol',  0.0,   30.68,   800.0,   10.35,  0.600,  0.134,  0.266   % LHV, rho: higher alcohols SSA; C/H/O, AFRst from C3H8O
    'Butanol',   0.0,   33.11,   804.0,   11.19,  0.648,  0.136,  0.216   % LHV, rho: higher alcohols SSA; C/H/O, AFRst from C4H10O
    'Pentanol',  0.0,   34.76,   812.0,   11.77,  0.681,  0.137,  0.181   % LHV, rho: higher alcohols SSA; C/H/O, AFRst from C5H12O
    'Hexanol',   0.0,   35.97,   816.0,   12.18,  0.705,  0.138,  0.157   % LHV, rho: higher alcohols SSA; C/H/O, AFRst from C6H14O
    'Octanol',   0.0,   37.59,   826.0,   12.74,  0.738,  0.139,  0.123   % LHV, rho: higher alcohols SSA; C/H/O, AFRst from C8H18O
    'Decanol',   0.0,   38.65,   826.0,   13.10,  0.759,  0.140,  0.101   % LHV, rho: higher alcohols SSA; C/H/O, AFRst from C10H22O
    };

names = T(:,1);
x     = cell2mat(T(:,2));            
P     = cell2mat(T(:,3:8));          % [LHV rho AFRst wC wH wO]

% Checks
if abs(sum(x) - 1) > 1e-6
    error('FuelMix:sum', 'Volume fractions add up to %.3f, not 1.', sum(x));
end
used = x > 0;
missing = used & any(isnan(P), 2);
if any(missing)
    error('FuelMix:missing', 'Properties missing for: %s', strjoin(names(missing)', ', '));
end

% Mixing (ideal): volume fractions -> mass fractions
m = x(used) .* P(used,2);            % kg of each component per m^3 of blend
y = m / sum(m);                      % mass fractions

parts = {};
for k = find(used)'
    parts{end+1} = sprintf('%s %g%%', names{k}, 100*x(k)); %#ok<AGROW>
end
Fuel.Name   = strjoin(parts, ' + ');
Fuel.rho    = sum(m);                % [kg/m^3]
Fuel.LHV    = (y' * P(used,1)) * 1e6;% [J/kg]
Fuel.AFRst  = y' * P(used,3);        % [-]
Fuel.wC     = y' * P(used,4);        % carbon mass fraction [-]
Fuel.wH     = y' * P(used,5);        % hydrogen mass fraction [-]
Fuel.wO     = y' * P(used,6);        % oxygen mass fraction [-]

% Consistency check: AFRst from the composition
AFRcalc = 31.998/0.2314 * (Fuel.wC/12.011 + Fuel.wH/(4*1.008) - Fuel.wO/31.998);
if abs(AFRcalc - Fuel.AFRst) > 0.3
    warning('FuelMix:AFR', 'AFRst from table (%.2f) and from composition (%.2f) differ.', Fuel.AFRst, AFRcalc);
end
end
function [SpS,El] = myload(filename,species)
% This function allows us to easily extract data from Nasa files
% filename : The name of the data file being loaded
% species  : Values being extracted from the file
load(filename);
isp = myfind({Sp.Name},species);
SpS = Sp(isp);
for i=1:length(El)
    El(i).Mass = El(i).Mass*1e-3;
end
end
%% Nasa Functions
function [Cp,CpU] = CpNasa(T,Sp)
%CpNasa(T,Sp):: computes heat capacity of species at Temp T
% 
%   Input: T, temperature
%          Sp, thermal database entry. So-called NASA polynomials
global Runiv
Cp=zeros(size(T));CpU=Cp;
if (isempty(Runiv))
    fprintf('[CpNasa] Assign global Runiv\n');
    return
end
if (isnan(Sp.Ts))
    Tl=T;
    a=Sp.Pol(1,:);
    CpU=a(1)+a(2).*Tl+a(3).*Tl.^2+a(4).*Tl.^3+a(5).*Tl.^4;        % Formula 5.4 of lecture notes
else
    ilow = (T <= Sp.Ts);
    Tl=T(ilow);
    a=Sp.Pol(1,:);
    CpU(ilow)=a(1)+a(2).*Tl+a(3).*Tl.^2+a(4).*Tl.^3+a(5).*Tl.^4;        % Formula 5.4 of lecture notes

    ihigh = (T > Sp.Ts);
    Tl=T(ihigh);
    a=Sp.Pol(2,:);
    CpU(ihigh)=a(1)+a(2).*Tl+a(3).*Tl.^2+a(4).*Tl.^3+a(5).*Tl.^4;        % Formula 5.4 of lecture notes
end
Cp=CpU.*Runiv/Sp.Mass;
end
function [Cv,CvU] = CvNasa(T,Sp)

%CvNasa(T,Sp):: computes heat capacity of species at Temp T
% 
%   Input: T, temperature
%          Sp, thermal database entry. So-called NASA polynomials
global Runiv
if (isempty(Runiv))
    fprintf('[CvNasa] Assign global Runiv\n');
    return
end
[~,CpU]=CpNasa(T,Sp);
CvU = CpU-1;
Cv=CvU.*Runiv/Sp.Mass;
end
function [H,HU] = HNasa(T,Sp)
%HNasa(T,Sp):: computes internal specific enthalpy of species at Temp T
% 
%   Input: T, temperature
%          Sp, thermal database entry. So-called NASA polynomials
global Runiv
H=zeros(size(T));HU=H;
if (isempty(Runiv))
    fprintf('[HNasa] Assign global Runiv\n');
    return
end
if (isnan(Sp.Ts))
    Tl=T;
    a=Sp.Pol(1,:);
    HU=a(6) + a(1).*Tl+(a(2)/2).*Tl.^2+(a(3)/3).*Tl.^3+(a(4)/4).*Tl.^4+(a(5)/5).*Tl.^5;        % Formula 5.4 of lecture notes    
else
    ilow = (T <= Sp.Ts);
    Tl=T(ilow);
    a=Sp.Pol(1,:);
    HU(ilow)=a(6) + a(1).*Tl+(a(2)/2).*Tl.^2+(a(3)/3).*Tl.^3+(a(4)/4).*Tl.^4+(a(5)/5).*Tl.^5;        % Formula 5.4 of lecture notes

    ihigh = (T > Sp.Ts);
    Tl=T(ihigh);
    a=Sp.Pol(2,:);
    HU(ihigh)=a(6) + a(1).*Tl+(a(2)/2).*Tl.^2+(a(3)/3).*Tl.^3+(a(4)/4).*Tl.^4+(a(5)/5).*Tl.^5;        % Formula 5.4 of lecture notes
end
H=HU.*Runiv/Sp.Mass;
end
function [ip]=myfind(S,str)
n=length(S);
m=length(str);
for j=1:m
    i=1;
    while ~strcmp(S(i),str(j)) & i < n
        i=i+1;
    end;
    if (strcmp(S(i),str(j)))
        %   disp(['Found ' str '=' S(i)]);
        ip(j)=i;
    else
        disp([str(j) ' NOT found']);
        ip(j)=-1;
    end;
end
end
function [S,SU] = SNasa(T,Sp)
%SNasa(T,Sp):: computes internal specific entropy of species at Temp T
% 
%   Input: T, temperature
%          Sp, thermal database entry. So-called NASA polynomials
global Runiv
S=zeros(size(T));SU=S;
if (isempty(Runiv))
    fprintf('[SNasa] Assign global Runiv\n');
    return
end
if (isnan(Sp.Ts))
    Tl=T;
    a=Sp.Pol(1,:);
    SU=  a(7) + a(1).*log(Tl) +(a(2)/1).*Tl.^1+(a(3)/2).*Tl.^2+(a(4)/3).*Tl.^3+(a(5)/4).*Tl.^4;        % Formula 5.4 of lecture notes
else
    ilow = (T <= Sp.Ts);
    Tl=T(ilow);
    a=Sp.Pol(1,:);
%     SU(ilow)=a(7) + a(1).*log(Tl) + (a(2)/1) +(a(3)/2).*Tl.^1+(a(4)/3).*Tl.^2+(a(5)/4).*Tl.^3;        % Formula 5.4 of lecture notes
    SU(ilow)= a(7) + a(1).*log(Tl) +(a(2)/1).*Tl.^1+(a(3)/2).*Tl.^2+(a(4)/3).*Tl.^3+(a(5)/4).*Tl.^4;        % Formula 5.4 of lecture notes
    
    ihigh = (T > Sp.Ts);
    Tl=T(ihigh);
    a=Sp.Pol(2,:);
%     SU(ihigh)=a(7) + a(1).*log(Tl) + (a(2)/1) +(a(3)/2).*Tl.^1+(a(4)/3).*Tl.^2+(a(5)/4).*Tl.^3;        % Formula 5.4 of lecture notes
    SU(ihigh)= a(7) + a(1).*log(Tl) +(a(2)/1).*Tl.^1+(a(3)/2).*Tl.^2+(a(4)/3).*Tl.^3+(a(5)/4).*Tl.^4;        % Formula 5.4 of lecture notes
end
S=SU.*Runiv/Sp.Mass;
end
function [Cp,Cv,H,E] = ThermoMix(Y,T,Sp)
global Runiv
%
% Y = squeeze(Yf);
Ma=0;
for i=1:length(Sp)
    Hi(i) = HNasa(T,Sp(i));
    Cpi(i)= CpNasa(T,Sp(i));
    Ma    = Ma+Y(i)/Sp(i).Mass;
end
Ma = 1/Ma;
H  = Y*Hi';  
E  = H-Runiv/Ma*T;
Cp = Y*Cpi';
Cv = Cp- Runiv/Ma;
end
function [U,UU] = UNasa(T,Sp)
%CvNasa(T,Sp):: computes heat capacity of species at Temp T
% 
%   Input: T, temperature
%          Sp, thermal database entry. So-called NASA polynomials
global Runiv
if (isempty(Runiv))
    fprintf('[UNasa] Assign global Runiv\n');
    return
end
[~,HU]=HNasa(T,Sp);
UU = HU-T;
U=UU.*Runiv/Sp.Mass;
end
