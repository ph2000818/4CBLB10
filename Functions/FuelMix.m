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