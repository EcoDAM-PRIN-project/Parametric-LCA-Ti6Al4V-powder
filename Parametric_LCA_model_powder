function optimizeLCA()
    % Optimizes an LCA process for Ti6Al4V powder production
    % Copyright (C) 2025 - Christian Spreafico, Baris Ördek

    % This program is free software: you can redistribute it and/or modify
    % it under the terms of the GNU General Public License as published by
    % the Free Software Foundation, either version 3 of the License, or
    % (at your option) any later version.

    % This program is distributed in the hope that it will be useful,
    % but WITHOUT ANY WARRANTY; without even the implied warranty of
    % MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    % GNU General Public License for more details.

    % You should have received a copy of the GNU General Public License
    % along with this program.  If not, see <https://www.gnu.org/licenses/>.
       
    % ================== CONTROL PARAMETERS AND CONSTANTS ==================
    params = struct(...
        'm_finalPowder', input('Enter the mass of the produced powder [kg]: '), ...                 
        'impact_category', input('Enter the impact category to minimize (TA, GW, FET, MET, TET, FF, ME, HTPc, HTPnc, IR, LO, SO, OD, PMF, HOF, EOF, WC): ', 's'), ...
        'country', input('Enter the country (EU, CN): ', 's'), ...                    
        'd_target', input('Enter the target diameter of the powder [µm]: '), ...                    
        'alpha_SM', 0.53, ...                   % Mass ratio melted ilmenite / titanium slag
        'alpha_CR', 2.38, ...                   % Mass ratio TiO2 / TiCl4
        'alpha_RD', 0.25, ...                   % Mass ratio TiCl4 / Ti sponge
        'R', 8.314, ...                         % Universal gas constant [J/mol K]
        'M_argon', 0.039948, ...                % Argon molar mass [kg/mol]
        'T_ATargon', 293, ...                   % Argon temperature [K]
        'Vm_ATargon', 0.03, ...                 % Argon volume [m3]
        'mm_SMpetrolPitch', 0.040, ...          % Specific mass of petrol pitch [kg/kg]
        'mm_SMpetrolCoke', 0.072, ...           % Specific mass of petrol coke [kg/kg]
        'mm_SMrawCoal', 0.0013, ...             % Specific mass of raw coal [kg/kg]
        'mm_SMcrudeOil', 0.081, ...             % Specific mass of crude oil [kg/kg]
        'mm_SMgraphite', 0.0096, ...            % Specific mass of grafite [kg/kg]
        'mm_SMsodiumOleate', 0.00037, ...       % Specific mass of sodium oleate [kg/kg]
        'Enm_CR_ng', 0.0066, ...                % Specific energy of chlorination and reduction [kWh/kg]
        'mm_CRfreshWater', 0.098, ...           % Specific mass of fresh water [kg/kg]
        'mm_CRpetrolCoke', 0.145, ...           % Specific mass of petrol coke [kg/kg]
        'mm_CRsodiumHydrox', 1.004, ...         % Specific mass of sodium hydroxide [kg/kg]
        'mm_CRchlorine', 0.36, ...              % Specific mass of chlorine gas [kg/kg]
        'mm_CRrawCoal', 0.002, ...              % Specific mass of raw oil [kg/kg]
        'mm_CRcrudeOil', 0.126, ...             % Specific mass of crude oil [kg/kg]
        'Enm_RD_el', 1.22, ...                  % Specific energy of reduction and distillation [kWh/kg]
        'mm_RDmagnesium', 0.01, ...             % Specific mass of magnesium [kg/kg]
        'Enm_CS_el', 0.02295, ...               % Specific energy of compaction and sintering [kWh/kg]
        'Enm_PS_el', 0.16, ...                  % Specific energy of powder sieving [kWh/kg]
        'Delta_ATargon', 0.7, ...               % Mass of recycled argon
        'kd', 5e-3, ...                         % Empirical constant
        'nhu_m', 5.0e-5, ...                    % Melted titanium viscosity [Pa·s]
        'nhu_g', 2.125e-5, ...                  % Argon viscosity [Pa·s]
        'We', 50, ...                           % Weber number
        'beta_TiO2Min', 0.75, ...               % Constraint: min TiO2 content in Ti slag
        'beta_TiO2Max', 0.9, ...                % Constraint: max TiO2 content in Ti slag
        'p_ATargonMin', 5.5, ...                % Constraint: min argon pressure [MPa]
        'p_ATargonMax', 7, ...                  % Constraint: max argon pressure [MPa]
        'phi_ATelectrodeMin', 0.05, ...         % Constraint: min electrode diameter [m]
        'phi_ATelectrodeMax', 0.10, ...         % Constraint: max electrode diameter [m]
        % Impact coefficients from LCA database
    );

    % ================== INDEPENDENT VARIABLES ==================
    x0 = [0.06, 6, 0.86];      % Initial values:  phi, p_ATargon, beta_TiO2
    lb = [params.phi_ATelectrodeMin, params.p_ATargonMin, params.beta_TiO2Min];    % Lower bounds
    ub = [params.phi_ATelectrodeMax, params.p_ATargonMax, params.beta_TiO2Max];   % Upper bounds

    % ================== INPUT VALIDATION ==================
    try
        validateattributes(x0, {'numeric'}, {'real', 'vector', 'numel', 3, 'finite'}, 'optimizeLCA', 'x0');
        validateattributes(lb, {'numeric'}, {'real', 'vector', 'numel', 3, 'finite'}, 'optimizeLCA', 'lb');
        validateattributes(ub, {'numeric'}, {'real', 'vector', 'numel', 3, 'finite'}, 'optimizeLCA', 'ub');

        for i = 1:length(x0)
            validateattributes(lb(i), {'numeric'}, {'real', 'scalar', 'finite', '<=', x0(i)}, 'optimizeLCA', sprintf('lb(%d)', i));
            validateattributes(ub(i), {'numeric'}, {'real', 'scalar', 'finite', '>=', x0(i)}, 'optimizeLCA', sprintf('ub(%d)', i));
        end

    catch ME
        error('optimizeLCA:InvalidInput', 'Input error: %s', ME.message);
    end

    % ================== OPTIMIZATION ==================
    options = optimoptions('fmincon', 'Display', 'iter', 'MaxFunctionEvaluations', 10000, 'MaxIterations', 5000, 'Algorithm', 'sqp', 'ConstraintTolerance', 1e-6);
    [x_opt, fval, exitflag, output] = fmincon(@(x) objectiveFunction(x, params), x0, [], [], [], [], lb, ub, @(x) constraintFunction(x, params), options);
    [~, I_impact_breakdown] = objectiveFunction(x_opt, params);

    % ================== OUTPUT RESULTS ==================
    var_names = {'phi_ATelectrode (m)', 'p_ATargon (MPa)', 'beta_TiO2'};
    disp('-----------------------------------------------------');
    disp('                 OPTIMIZATION RESULTS                  ');
    disp('-----------------------------------------------------');
    for i = 1:length(x_opt)
        fprintf('%s: %.5f\n', var_names{i}, x_opt(i));
    end

    % ================== CALCULATE DEPENDENT VARIABLES ==================
    dependent_vars = calculateDependentVariables(x_opt, params);

    % ================== OUTPUT DEPENDENT VARIABLES ==================
    fprintf('mm_ATargon (kg): %.5f\n', dependent_vars.mm_ATargon);
    fprintf('m_ilmenite (kg): %.5f\n', dependent_vars.m_ilmenite);
    fprintf('m_TiO2 (kg): %.5f\n', dependent_vars.m_TiO2);
    fprintf('m_TiCl4 (kg): %.5f\n', dependent_vars.m_TiCl4);
    fprintf('m_TiSponge (kg): %.5f\n', dependent_vars.m_TiSponge);
    fprintf('m_Ti64ingot (kg): %.5f\n', dependent_vars.m_Ti64ingot);
    fprintf('m_atomizedPowder (kg): %.5f\n', dependent_vars.m_atomizedPowder);
    fprintf('m_wastePowder (kg): %.5f\n', dependent_vars.m_wastePowder);
    fprintf('En_SM_el (kWh): %.5f\n', dependent_vars.m_ilmenite*dependent_vars.Enm_SM_el);
    fprintf('En_CR_el (kWh): %.5f\n', dependent_vars.m_ilmenite*dependent_vars.Enm_CR_el);
    fprintf('En_RE_el (kWh): %.5f\n', dependent_vars.m_Ti64ingot*dependent_vars.Enm_RE_el);
    fprintf('En_ATmelt_el (kWh): %.5f\n', dependent_vars.m_Ti64ingot*dependent_vars.Enm_ATmelt_el);
    fprintf('m_argon (kg): %.5f\n', dependent_vars.m_atomizedPowder*dependent_vars.mm_ATargon*(1-params.Delta_ATargon));
    disp('---------------------------------');
    disp('Environmental Impact Breakdown:');
    fprintf('Mineral Extraction: %.5f\n', I_impact_breakdown.MI);
    fprintf('Smelting: %.5f\n', I_impact_breakdown.SM);
    fprintf('Chlorination: %.5f\n', I_impact_breakdown.CR);
    fprintf('Reduction: %.5f\n', I_impact_breakdown.RD);
    fprintf('Compaction and sintering: %.5f\n', I_impact_breakdown.CS);
    fprintf('Remelting: %.5f\n', I_impact_breakdown.RE);
    fprintf('Atomization: %.5f\n', I_impact_breakdown.AT);
    fprintf('Powder Sieving: %.5f\n', I_impact_breakdown.PS);
    fprintf('Objective function value: %.5f\n', fval);
end

function [I_tot, I_impact_breakdown] = objectiveFunction(x, params)
    %   Minimize this function to optimize the process.

    % Extract independent variables
    phi_ATelectrode = x(1);
    p_ATargon = x(2);
    beta_TiO2 = x(3);

    % Calculate dependent variables
    dependent_vars = calculateDependentVariables(x, params);

    switch params.impact_category
        case 'TA'
            Ic_Ilmenite = params.Ic_Ilmenite_GLO_TA;
            Ic_PetroleumCoke = params.Ic_PetroleumCoke_GLO_TA;
            Ic_SodiumOleate = params.Ic_SodiumOleate_ROW_TA; 
            Ic_SodiumHydroxide = params.Ic_SodiumHydroxide_GLO_TA;
            Ic_Magnesium = params.Ic_Magnesium_GLO_TA;
            Ic_Aluminium = params.Ic_Aluminium_GLO_TA;
            Ic_Vanadium = params.Ic_Vanadium_CN_TA;
            Ic_NaturalGas = params.Ic_NaturalGas_GLO_TA;
            if strcmp(params.country, 'EU')
                Ic_PetroleumPitch = params.Ic_PetroleumPitch_Europe_TA;
                Ic_RawCoal = params.Ic_HardCoal_Europe_TA;
                Ic_CrudeOil = params.Ic_CrudeOil_Europe_TA;
                Ic_Graphite = params.Ic_Graphite_RoW_TA;
                Ic_FreshWater = params.Ic_TapWater_Europe_TA;
                Ic_Chlorine = params.Ic_Chlorine_Europe_TA;
                Ic_Argon = params.Ic_Argon_Europe_TA;
                Ic_Electricity = params.Ic_Electricity_Europe_TA;
            end
            if strcmp(params.country, 'CN')
                Ic_PetroleumPitch = params.Ic_PetroleumPitch_RoW_TA;
                Ic_RawCoal = params.Ic_HardCoal_CN_TA;
                Ic_CrudeOil = params.Ic_CrudeOil_RoW_TA;
                Ic_Graphite = params.Ic_Graphite_CN_TA;
                Ic_FreshWater = params.Ic_TapWater_RoW_TA;
                Ic_Chlorine = params.Ic_Chlorine_RoW_TA;
                Ic_Argon = params.Ic_Argon_RoW_TA;
                Ic_Electricity = params.Ic_Electricity_CN_TA;
            end
        case 'GW'
            Ic_Ilmenite = params.Ic_Ilmenite_GLO_GW;
            Ic_PetroleumCoke = params.Ic_PetroleumCoke_GLO_GW;
            Ic_SodiumOleate = params.Ic_SodiumOleate_ROW_GW; 
            Ic_SodiumHydroxide = params.Ic_SodiumHydroxide_GLO_GW;
            Ic_Magnesium = params.Ic_Magnesium_GLO_GW;
            Ic_Aluminium = params.Ic_Aluminium_GLO_GW;
            Ic_Vanadium = params.Ic_Vanadium_CN_GW;
            Ic_NaturalGas = params.Ic_NaturalGas_GLO_GW;
            if strcmp(params.country, 'EU')
                Ic_PetroleumPitch = params.Ic_PetroleumPitch_Europe_GW;
                Ic_RawCoal = params.Ic_HardCoal_Europe_GW;
                Ic_CrudeOil = params.Ic_CrudeOil_Europe_GW;
                Ic_Graphite = params.Ic_Graphite_RoW_GW;
                Ic_FreshWater = params.Ic_TapWater_Europe_GW;
                Ic_Chlorine = params.Ic_Chlorine_Europe_GW;
                Ic_Argon = params.Ic_Argon_Europe_GW;
                Ic_Electricity = params.Ic_Electricity_Europe_GW;
            end
            if strcmp(params.country, 'CN')
                Ic_PetroleumPitch = params.Ic_PetroleumPitch_RoW_GW;
                Ic_RawCoal = params.Ic_HardCoal_CN_GW;
                Ic_CrudeOil = params.Ic_CrudeOil_RoW_GW;
                Ic_Graphite = params.Ic_Graphite_CN_GW;
                Ic_FreshWater = params.Ic_TapWater_RoW_GW;
                Ic_Chlorine = params.Ic_Chlorine_RoW_GW;
                Ic_Argon = params.Ic_Argon_RoW_GW;
                Ic_Electricity = params.Ic_Electricity_CN_GW;
            end
        case 'FET'
            Ic_Ilmenite = params.Ic_Ilmenite_GLO_FET;
            Ic_PetroleumCoke = params.Ic_PetroleumCoke_GLO_FET;
            Ic_SodiumOleate = params.Ic_SodiumOleate_ROW_FET; 
            Ic_SodiumHydroxide = params.Ic_SodiumHydroxide_GLO_FET;
            Ic_Magnesium = params.Ic_Magnesium_GLO_FET;
            Ic_Aluminium = params.Ic_Aluminium_GLO_FET;
            Ic_Vanadium = params.Ic_Vanadium_CN_FET;
            Ic_NaturalGas = params.Ic_NaturalGas_GLO_FET;
            if strcmp(params.country, 'EU')
                Ic_PetroleumPitch = params.Ic_PetroleumPitch_Europe_FET;
                Ic_RawCoal = params.Ic_HardCoal_Europe_FET;
                Ic_CrudeOil = params.Ic_CrudeOil_Europe_FET;
                Ic_Graphite = params.Ic_Graphite_RoW_FET;
                Ic_FreshWater = params.Ic_TapWater_Europe_FET;
                Ic_Chlorine = params.Ic_Chlorine_Europe_FET;
                Ic_Argon = params.Ic_Argon_Europe_FET;
                Ic_Electricity = params.Ic_Electricity_Europe_FET;
            end
            if strcmp(params.country, 'CN')
                Ic_PetroleumPitch = params.Ic_PetroleumPitch_RoW_FET;
                Ic_RawCoal = params.Ic_HardCoal_CN_FET;
                Ic_CrudeOil = params.Ic_CrudeOil_RoW_FET;
                Ic_Graphite = params.Ic_Graphite_CN_FET;
                Ic_FreshWater = params.Ic_TapWater_RoW_FET;
                Ic_Chlorine = params.Ic_Chlorine_RoW_FET;
                Ic_Argon = params.Ic_Argon_RoW_FET;
                Ic_Electricity = params.Ic_Electricity_CN_FET;
            end
        case 'MET'
            Ic_Ilmenite = params.Ic_Ilmenite_GLO_MET;
            Ic_PetroleumCoke = params.Ic_PetroleumCoke_GLO_MET;
            Ic_SodiumOleate = params.Ic_SodiumOleate_ROW_MET; 
            Ic_SodiumHydroxide = params.Ic_SodiumHydroxide_GLO_MET;
            Ic_Magnesium = params.Ic_Magnesium_GLO_MET;
            Ic_Aluminium = params.Ic_Aluminium_GLO_MET;
            Ic_Vanadium = params.Ic_Vanadium_CN_MET;
            Ic_NaturalGas = params.Ic_NaturalGas_GLO_MET;
            if strcmp(params.country, 'EU')
                Ic_PetroleumPitch = params.Ic_PetroleumPitch_Europe_MET;
                Ic_RawCoal = params.Ic_HardCoal_Europe_MET;
                Ic_CrudeOil = params.Ic_CrudeOil_Europe_MET;
                Ic_Graphite = params.Ic_Graphite_RoW_MET;
                Ic_FreshWater = params.Ic_TapWater_Europe_MET;
                Ic_Chlorine = params.Ic_Chlorine_Europe_MET;
                Ic_Argon = params.Ic_Argon_Europe_MET;
                Ic_Electricity = params.Ic_Electricity_Europe_MET;
            end
            if strcmp(params.country, 'CN')
                Ic_PetroleumPitch = params.Ic_PetroleumPitch_RoW_MET;
                Ic_RawCoal = params.Ic_HardCoal_CN_MET;
                Ic_CrudeOil = params.Ic_CrudeOil_RoW_MET;
                Ic_Graphite = params.Ic_Graphite_CN_MET;
                Ic_FreshWater = params.Ic_TapWater_RoW_MET;
                Ic_Chlorine = params.Ic_Chlorine_RoW_MET;
                Ic_Argon = params.Ic_Argon_RoW_MET;
                Ic_Electricity = params.Ic_Electricity_CN_MET;
            end
        case 'TET'
            Ic_Ilmenite = params.Ic_Ilmenite_GLO_TET;
            Ic_PetroleumCoke = params.Ic_PetroleumCoke_GLO_TET;
            Ic_SodiumOleate = params.Ic_SodiumOleate_ROW_TET; 
            Ic_SodiumHydroxide = params.Ic_SodiumHydroxide_GLO_TET;
            Ic_Magnesium = params.Ic_Magnesium_GLO_TET;
            Ic_Aluminium = params.Ic_Aluminium_GLO_TET;
            Ic_Vanadium = params.Ic_Vanadium_CN_TET;
            Ic_NaturalGas = params.Ic_NaturalGas_GLO_TET;
            if strcmp(params.country, 'EU')
                Ic_PetroleumPitch = params.Ic_PetroleumPitch_Europe_TET;
                Ic_RawCoal = params.Ic_HardCoal_Europe_TET;
                Ic_CrudeOil = params.Ic_CrudeOil_Europe_TET;
                Ic_Graphite = params.Ic_Graphite_RoW_TET;
                Ic_FreshWater = params.Ic_TapWater_Europe_TET;
                Ic_Chlorine = params.Ic_Chlorine_Europe_TET;
                Ic_Argon = params.Ic_Argon_Europe_TET;
                Ic_Electricity = params.Ic_Electricity_Europe_TET;
            end
            if strcmp(params.country, 'CN')
                Ic_PetroleumPitch = params.Ic_PetroleumPitch_RoW_TET;
                Ic_RawCoal = params.Ic_HardCoal_CN_TET;
                Ic_CrudeOil = params.Ic_CrudeOil_RoW_TET;
                Ic_Graphite = params.Ic_Graphite_CN_TET;
                Ic_FreshWater = params.Ic_TapWater_RoW_TET;
                Ic_Chlorine = params.Ic_Chlorine_RoW_TET;
                Ic_Argon = params.Ic_Argon_RoW_TET;
                Ic_Electricity = params.Ic_Electricity_CN_TET;
            end
        case 'FF'
            Ic_Ilmenite = params.Ic_Ilmenite_GLO_FF;
            Ic_PetroleumCoke = params.Ic_PetroleumCoke_GLO_FF;
            Ic_SodiumOleate = params.Ic_SodiumOleate_ROW_FF; 
            Ic_SodiumHydroxide = params.Ic_SodiumHydroxide_GLO_FF;
            Ic_Magnesium = params.Ic_Magnesium_GLO_FF;
            Ic_Aluminium = params.Ic_Aluminium_GLO_FF;
            Ic_Vanadium = params.Ic_Vanadium_CN_FF;
            Ic_NaturalGas = params.Ic_NaturalGas_GLO_FF;
            if strcmp(params.country, 'EU')
                Ic_PetroleumPitch = params.Ic_PetroleumPitch_Europe_FF;
                Ic_RawCoal = params.Ic_HardCoal_Europe_FF;
                Ic_CrudeOil = params.Ic_CrudeOil_Europe_FF;
                Ic_Graphite = params.Ic_Graphite_RoW_FF;
                Ic_FreshWater = params.Ic_TapWater_Europe_FF;
                Ic_Chlorine = params.Ic_Chlorine_Europe_FF;
                Ic_Argon = params.Ic_Argon_Europe_FF;
                Ic_Electricity = params.Ic_Electricity_Europe_FF;
            end
            if strcmp(params.country, 'CN')
                Ic_PetroleumPitch = params.Ic_PetroleumPitch_RoW_FF;
                Ic_RawCoal = params.Ic_HardCoal_CN_FF;
                Ic_CrudeOil = params.Ic_CrudeOil_RoW_FF;
                Ic_Graphite = params.Ic_Graphite_CN_FF;
                Ic_FreshWater = params.Ic_TapWater_RoW_FF;
                Ic_Chlorine = params.Ic_Chlorine_RoW_FF;
                Ic_Argon = params.Ic_Argon_RoW_FF;
                Ic_Electricity = params.Ic_Electricity_CN_FF;
            end
        case 'FE'
            Ic_Ilmenite = params.Ic_Ilmenite_GLO_FE;
            Ic_PetroleumCoke = params.Ic_PetroleumCoke_GLO_FE;
            Ic_SodiumOleate = params.Ic_SodiumOleate_ROW_FE; 
            Ic_SodiumHydroxide = params.Ic_SodiumHydroxide_GLO_FE;
            Ic_Magnesium = params.Ic_Magnesium_GLO_FE;
            Ic_Aluminium = params.Ic_Aluminium_GLO_FE;
            Ic_Vanadium = params.Ic_Vanadium_CN_FE;
            Ic_NaturalGas = params.Ic_NaturalGas_GLO_FE;
            if strcmp(params.country, 'EU')
                Ic_PetroleumPitch = params.Ic_PetroleumPitch_Europe_FE;
                Ic_RawCoal = params.Ic_HardCoal_Europe_FE;
                Ic_CrudeOil = params.Ic_CrudeOil_Europe_FE;
                Ic_Graphite = params.Ic_Graphite_RoW_FE;
                Ic_FreshWater = params.Ic_TapWater_Europe_FE;
                Ic_Chlorine = params.Ic_Chlorine_Europe_FE;
                Ic_Argon = params.Ic_Argon_Europe_FE;
                Ic_Electricity = params.Ic_Electricity_Europe_FE;
            end
            if strcmp(params.country, 'CN')
                Ic_PetroleumPitch = params.Ic_PetroleumPitch_RoW_FE;
                Ic_RawCoal = params.Ic_HardCoal_CN_FE;
                Ic_CrudeOil = params.Ic_CrudeOil_RoW_FE;
                Ic_Graphite = params.Ic_Graphite_CN_FE;
                Ic_FreshWater = params.Ic_TapWater_RoW_FE;
                Ic_Chlorine = params.Ic_Chlorine_RoW_FE;
                Ic_Argon = params.Ic_Argon_RoW_FE;
                Ic_Electricity = params.Ic_Electricity_CN_FE;
            end
        case 'ME'
            Ic_Ilmenite = params.Ic_Ilmenite_GLO_ME;
            Ic_PetroleumCoke = params.Ic_PetroleumCoke_GLO_ME;
            Ic_SodiumOleate = params.Ic_SodiumOleate_ROW_ME; 
            Ic_SodiumHydroxide = params.Ic_SodiumHydroxide_GLO_ME;
            Ic_Magnesium = params.Ic_Magnesium_GLO_ME;
            Ic_Aluminium = params.Ic_Aluminium_GLO_ME;
            Ic_Vanadium = params.Ic_Vanadium_CN_ME;
            Ic_NaturalGas = params.Ic_NaturalGas_GLO_ME;
            if strcmp(params.country, 'EU')
                Ic_PetroleumPitch = params.Ic_PetroleumPitch_Europe_ME;
                Ic_RawCoal = params.Ic_HardCoal_Europe_ME;
                Ic_CrudeOil = params.Ic_CrudeOil_Europe_ME;
                Ic_Graphite = params.Ic_Graphite_RoW_ME;
                Ic_FreshWater = params.Ic_TapWater_Europe_ME;
                Ic_Chlorine = params.Ic_Chlorine_Europe_ME;
                Ic_Argon = params.Ic_Argon_Europe_ME;
                Ic_Electricity = params.Ic_Electricity_Europe_ME;
            end
            if strcmp(params.country, 'CN')
                Ic_PetroleumPitch = params.Ic_PetroleumPitch_RoW_ME;
                Ic_RawCoal = params.Ic_HardCoal_CN_ME;
                Ic_CrudeOil = params.Ic_CrudeOil_RoW_ME;
                Ic_Graphite = params.Ic_Graphite_CN_ME;
                Ic_FreshWater = params.Ic_TapWater_RoW_ME;
                Ic_Chlorine = params.Ic_Chlorine_RoW_ME;
                Ic_Argon = params.Ic_Argon_RoW_ME;
                Ic_Electricity = params.Ic_Electricity_CN_ME;
            end
        case 'HTPc'
            Ic_Ilmenite = params.Ic_Ilmenite_GLO_HTPc;
            Ic_PetroleumCoke = params.Ic_PetroleumCoke_GLO_HTPc;
            Ic_SodiumOleate = params.Ic_SodiumOleate_ROW_HTPc; 
            Ic_SodiumHydroxide = params.Ic_SodiumHydroxide_GLO_HTPc;
            Ic_Magnesium = params.Ic_Magnesium_GLO_HTPc;
            Ic_Aluminium = params.Ic_Aluminium_GLO_HTPc;
            Ic_Vanadium = params.Ic_Vanadium_CN_HTPc;
            Ic_NaturalGas = params.Ic_NaturalGas_GLO_HTPc;
            if strcmp(params.country, 'EU')
                Ic_PetroleumPitch = params.Ic_PetroleumPitch_Europe_HTPc;
                Ic_RawCoal = params.Ic_HardCoal_Europe_HTPc;
                Ic_CrudeOil = params.Ic_CrudeOil_Europe_HTPc;
                Ic_Graphite = params.Ic_Graphite_RoW_HTPc;
                Ic_FreshWater = params.Ic_TapWater_Europe_HTPc;
                Ic_Chlorine = params.Ic_Chlorine_Europe_HTPc;
                Ic_Argon = params.Ic_Argon_Europe_HTPc;
                Ic_Electricity = params.Ic_Electricity_Europe_HTPc;
            end
            if strcmp(params.country, 'CN')
                Ic_PetroleumPitch = params.Ic_PetroleumPitch_RoW_HTPc;
                Ic_RawCoal = params.Ic_HardCoal_CN_HTPc;
                Ic_CrudeOil = params.Ic_CrudeOil_RoW_HTPc;
                Ic_Graphite = params.Ic_Graphite_CN_HTPc;
                Ic_FreshWater = params.Ic_TapWater_RoW_HTPc;
                Ic_Chlorine = params.Ic_Chlorine_RoW_HTPc;
                Ic_Argon = params.Ic_Argon_RoW_HTPc;
                Ic_Electricity = params.Ic_Electricity_CN_HTPc;
            end
        case 'HTPnc'
            Ic_Ilmenite = params.Ic_Ilmenite_GLO_HTPnc;
            Ic_PetroleumCoke = params.Ic_PetroleumCoke_GLO_HTPnc;
            Ic_SodiumOleate = params.Ic_SodiumOleate_ROW_HTPnc; 
            Ic_SodiumHydroxide = params.Ic_SodiumHydroxide_GLO_HTPnc;
            Ic_Magnesium = params.Ic_Magnesium_GLO_HTPnc;
            Ic_Aluminium = params.Ic_Aluminium_GLO_HTPnc;
            Ic_Vanadium = params.Ic_Vanadium_CN_HTPnc;
            Ic_NaturalGas = params.Ic_NaturalGas_GLO_HTPnc;
            if strcmp(params.country, 'EU')
                Ic_PetroleumPitch = params.Ic_PetroleumPitch_Europe_HTPnc;
                Ic_RawCoal = params.Ic_HardCoal_Europe_HTPnc;
                Ic_CrudeOil = params.Ic_CrudeOil_Europe_HTPnc;
                Ic_Graphite = params.Ic_Graphite_RoW_HTPnc;
                Ic_FreshWater = params.Ic_TapWater_Europe_HTPnc;
                Ic_Chlorine = params.Ic_Chlorine_Europe_HTPnc;
                Ic_Argon = params.Ic_Argon_Europe_HTPnc;
                Ic_Electricity = params.Ic_Electricity_Europe_HTPnc;
            end
            if strcmp(params.country, 'CN')
                Ic_PetroleumPitch = params.Ic_PetroleumPitch_RoW_HTPnc;
                Ic_RawCoal = params.Ic_HardCoal_CN_HTPnc;
                Ic_CrudeOil = params.Ic_CrudeOil_RoW_HTPnc;
                Ic_Graphite = params.Ic_Graphite_CN_HTPnc;
                Ic_FreshWater = params.Ic_TapWater_RoW_HTPnc;
                Ic_Chlorine = params.Ic_Chlorine_RoW_HTPnc;
                Ic_Argon = params.Ic_Argon_RoW_HTPnc;
                Ic_Electricity = params.Ic_Electricity_CN_HTPnc;
            end
        case 'IR'
            Ic_Ilmenite = params.Ic_Ilmenite_GLO_IR;
            Ic_PetroleumCoke = params.Ic_PetroleumCoke_GLO_IR;
            Ic_SodiumOleate = params.Ic_SodiumOleate_ROW_IR; 
            Ic_SodiumHydroxide = params.Ic_SodiumHydroxide_GLO_IR;
            Ic_Magnesium = params.Ic_Magnesium_GLO_IR;
            Ic_Aluminium = params.Ic_Aluminium_GLO_IR;
            Ic_Vanadium = params.Ic_Vanadium_CN_IR;
            Ic_NaturalGas = params.Ic_NaturalGas_GLO_IR;
            if strcmp(params.country, 'EU')
                Ic_PetroleumPitch = params.Ic_PetroleumPitch_Europe_IR;
                Ic_RawCoal = params.Ic_HardCoal_Europe_IR;
                Ic_CrudeOil = params.Ic_CrudeOil_Europe_IR;
                Ic_Graphite = params.Ic_Graphite_RoW_IR;
                Ic_FreshWater = params.Ic_TapWater_Europe_IR;
                Ic_Chlorine = params.Ic_Chlorine_Europe_IR;
                Ic_Argon = params.Ic_Argon_Europe_IR;
                Ic_Electricity = params.Ic_Electricity_Europe_IR;
            end
            if strcmp(params.country, 'CN')
                Ic_PetroleumPitch = params.Ic_PetroleumPitch_RoW_IR;
                Ic_RawCoal = params.Ic_HardCoal_CN_IR;
                Ic_CrudeOil = params.Ic_CrudeOil_RoW_IR;
                Ic_Graphite = params.Ic_Graphite_CN_IR;
                Ic_FreshWater = params.Ic_TapWater_RoW_IR;
                Ic_Chlorine = params.Ic_Chlorine_RoW_IR;
                Ic_Argon = params.Ic_Argon_RoW_IR;
                Ic_Electricity = params.Ic_Electricity_CN_IR;
            end
        case 'LO'
            Ic_Ilmenite = params.Ic_Ilmenite_GLO_LO;
            Ic_PetroleumCoke = params.Ic_PetroleumCoke_GLO_LO;
            Ic_SodiumOleate = params.Ic_SodiumOleate_ROW_LO; 
            Ic_SodiumHydroxide = params.Ic_SodiumHydroxide_GLO_LO;
            Ic_Magnesium = params.Ic_Magnesium_GLO_LO;
            Ic_Aluminium = params.Ic_Aluminium_GLO_LO;
            Ic_Vanadium = params.Ic_Vanadium_CN_LO;
            Ic_NaturalGas = params.Ic_NaturalGas_GLO_LO;
            if strcmp(params.country, 'EU')
                Ic_PetroleumPitch = params.Ic_PetroleumPitch_Europe_LO;
                Ic_RawCoal = params.Ic_HardCoal_Europe_LO;
                Ic_CrudeOil = params.Ic_CrudeOil_Europe_LO;
                Ic_Graphite = params.Ic_Graphite_RoW_LO;
                Ic_FreshWater = params.Ic_TapWater_Europe_LO;
                Ic_Chlorine = params.Ic_Chlorine_Europe_LO;
                Ic_Argon = params.Ic_Argon_Europe_LO;
                Ic_Electricity = params.Ic_Electricity_Europe_LO;
            end
            if strcmp(params.country, 'CN')
                Ic_PetroleumPitch = params.Ic_PetroleumPitch_RoW_LO;
                Ic_RawCoal = params.Ic_HardCoal_CN_LO;
                Ic_CrudeOil = params.Ic_CrudeOil_RoW_LO;
                Ic_Graphite = params.Ic_Graphite_CN_LO;
                Ic_FreshWater = params.Ic_TapWater_RoW_LO;
                Ic_Chlorine = params.Ic_Chlorine_RoW_LO;
                Ic_Argon = params.Ic_Argon_RoW_LO;
                Ic_Electricity = params.Ic_Electricity_CN_LO;
            end
        case 'SO'
            Ic_Ilmenite = params.Ic_Ilmenite_GLO_SO;
            Ic_PetroleumCoke = params.Ic_PetroleumCoke_GLO_SO;
            Ic_SodiumOleate = params.Ic_SodiumOleate_ROW_SO; 
            Ic_SodiumHydroxide = params.Ic_SodiumHydroxide_GLO_SO;
            Ic_Magnesium = params.Ic_Magnesium_GLO_SO;
            Ic_Aluminium = params.Ic_Aluminium_GLO_SO;
            Ic_Vanadium = params.Ic_Vanadium_CN_SO;
            Ic_NaturalGas = params.Ic_NaturalGas_GLO_SO;
            if strcmp(params.country, 'EU')
                Ic_PetroleumPitch = params.Ic_PetroleumPitch_Europe_SO;
                Ic_RawCoal = params.Ic_HardCoal_Europe_SO;
                Ic_CrudeOil = params.Ic_CrudeOil_Europe_SO;
                Ic_Graphite = params.Ic_Graphite_RoW_SO;
                Ic_FreshWater = params.Ic_TapWater_Europe_SO;
                Ic_Chlorine = params.Ic_Chlorine_Europe_SO;
                Ic_Argon = params.Ic_Argon_Europe_SO;
                Ic_Electricity = params.Ic_Electricity_Europe_SO;
            end
            if strcmp(params.country, 'CN')
                Ic_PetroleumPitch = params.Ic_PetroleumPitch_RoW_SO;
                Ic_RawCoal = params.Ic_HardCoal_CN_SO;
                Ic_CrudeOil = params.Ic_CrudeOil_RoW_SO;
                Ic_Graphite = params.Ic_Graphite_CN_SO;
                Ic_FreshWater = params.Ic_TapWater_RoW_SO;
                Ic_Chlorine = params.Ic_Chlorine_RoW_SO;
                Ic_Argon = params.Ic_Argon_RoW_SO;
                Ic_Electricity = params.Ic_Electricity_CN_SO;
            end
        case 'OD'
            Ic_Ilmenite = params.Ic_Ilmenite_GLO_OD;
            Ic_PetroleumCoke = params.Ic_PetroleumCoke_GLO_OD;
            Ic_SodiumOleate = params.Ic_SodiumOleate_ROW_OD; 
            Ic_SodiumHydroxide = params.Ic_SodiumHydroxide_GLO_OD;
            Ic_Magnesium = params.Ic_Magnesium_GLO_OD;
            Ic_Aluminium = params.Ic_Aluminium_GLO_OD;
            Ic_Vanadium = params.Ic_Vanadium_CN_OD;
            Ic_NaturalGas = params.Ic_NaturalGas_GLO_OD;
            if strcmp(params.country, 'EU')
                Ic_PetroleumPitch = params.Ic_PetroleumPitch_Europe_OD;
                Ic_RawCoal = params.Ic_HardCoal_Europe_OD;
                Ic_CrudeOil = params.Ic_CrudeOil_Europe_OD;
                Ic_Graphite = params.Ic_Graphite_RoW_OD;
                Ic_FreshWater = params.Ic_TapWater_Europe_OD;
                Ic_Chlorine = params.Ic_Chlorine_Europe_OD;
                Ic_Argon = params.Ic_Argon_Europe_OD;
                Ic_Electricity = params.Ic_Electricity_Europe_OD;
            end
            if strcmp(params.country, 'CN')
                Ic_PetroleumPitch = params.Ic_PetroleumPitch_RoW_OD;
                Ic_RawCoal = params.Ic_HardCoal_CN_OD;
                Ic_CrudeOil = params.Ic_CrudeOil_RoW_OD;
                Ic_Graphite = params.Ic_Graphite_CN_OD;
                Ic_FreshWater = params.Ic_TapWater_RoW_OD;
                Ic_Chlorine = params.Ic_Chlorine_RoW_OD;
                Ic_Argon = params.Ic_Argon_RoW_OD;
                Ic_Electricity = params.Ic_Electricity_CN_OD;
            end
        case 'PMF'
            Ic_Ilmenite = params.Ic_Ilmenite_GLO_PMF;
            Ic_PetroleumCoke = params.Ic_PetroleumCoke_GLO_PMF;
            Ic_SodiumOleate = params.Ic_SodiumOleate_ROW_PMF; 
            Ic_SodiumHydroxide = params.Ic_SodiumHydroxide_GLO_PMF;
            Ic_Magnesium = params.Ic_Magnesium_GLO_PMF;
            Ic_Aluminium = params.Ic_Aluminium_GLO_PMF;
            Ic_Vanadium = params.Ic_Vanadium_CN_PMF;
            Ic_NaturalGas = params.Ic_NaturalGas_GLO_PMF;
            if strcmp(params.country, 'EU')
                Ic_PetroleumPitch = params.Ic_PetroleumPitch_Europe_PMF;
                Ic_RawCoal = params.Ic_HardCoal_Europe_PMF;
                Ic_CrudeOil = params.Ic_CrudeOil_Europe_PMF;
                Ic_Graphite = params.Ic_Graphite_RoW_PMF;
                Ic_FreshWater = params.Ic_TapWater_Europe_PMF;
                Ic_Chlorine = params.Ic_Chlorine_Europe_PMF;
                Ic_Argon = params.Ic_Argon_Europe_PMF;
                Ic_Electricity = params.Ic_Electricity_Europe_PMF;
            end
            if strcmp(params.country, 'CN')
                Ic_PetroleumPitch = params.Ic_PetroleumPitch_RoW_PMF;
                Ic_RawCoal = params.Ic_HardCoal_CN_PMF;
                Ic_CrudeOil = params.Ic_CrudeOil_RoW_PMF;
                Ic_Graphite = params.Ic_Graphite_CN_PMF;
                Ic_FreshWater = params.Ic_TapWater_RoW_PMF;
                Ic_Chlorine = params.Ic_Chlorine_RoW_PMF;
                Ic_Argon = params.Ic_Argon_RoW_PMF;
                Ic_Electricity = params.Ic_Electricity_CN_PMF;
            end
        case 'HOF'
            Ic_Ilmenite = params.Ic_Ilmenite_GLO_HOF;
            Ic_PetroleumCoke = params.Ic_PetroleumCoke_GLO_HOF;
            Ic_SodiumOleate = params.Ic_SodiumOleate_ROW_HOF; 
            Ic_SodiumHydroxide = params.Ic_SodiumHydroxide_GLO_HOF;
            Ic_Magnesium = params.Ic_Magnesium_GLO_HOF;
            Ic_Aluminium = params.Ic_Aluminium_GLO_HOF;
            Ic_Vanadium = params.Ic_Vanadium_CN_HOF;
            Ic_NaturalGas = params.Ic_NaturalGas_GLO_HOF;
            if strcmp(params.country, 'EU')
                Ic_PetroleumPitch = params.Ic_PetroleumPitch_Europe_HOF;
                Ic_RawCoal = params.Ic_HardCoal_Europe_HOF;
                Ic_CrudeOil = params.Ic_CrudeOil_Europe_HOF;
                Ic_Graphite = params.Ic_Graphite_RoW_HOF;
                Ic_FreshWater = params.Ic_TapWater_Europe_HOF;
                Ic_Chlorine = params.Ic_Chlorine_Europe_HOF;
                Ic_Argon = params.Ic_Argon_Europe_HOF;
                Ic_Electricity = params.Ic_Electricity_Europe_HOF;
            end
            if strcmp(params.country, 'CN')
                Ic_PetroleumPitch = params.Ic_PetroleumPitch_RoW_HOF;
                Ic_RawCoal = params.Ic_HardCoal_CN_HOF;
                Ic_CrudeOil = params.Ic_CrudeOil_RoW_HOF;
                Ic_Graphite = params.Ic_Graphite_CN_HOF;
                Ic_FreshWater = params.Ic_TapWater_RoW_HOF;
                Ic_Chlorine = params.Ic_Chlorine_RoW_HOF;
                Ic_Argon = params.Ic_Argon_RoW_HOF;
                Ic_Electricity = params.Ic_Electricity_CN_HOF;
            end
        case 'EOF'
            Ic_Ilmenite = params.Ic_Ilmenite_GLO_EOF;
            Ic_PetroleumCoke = params.Ic_PetroleumCoke_GLO_EOF;
            Ic_SodiumOleate = params.Ic_SodiumOleate_ROW_EOF; 
            Ic_SodiumHydroxide = params.Ic_SodiumHydroxide_GLO_EOF;
            Ic_Magnesium = params.Ic_Magnesium_GLO_EOF;
            Ic_Aluminium = params.Ic_Aluminium_GLO_EOF;
            Ic_Vanadium = params.Ic_Vanadium_CN_EOF;
            Ic_NaturalGas = params.Ic_NaturalGas_GLO_EOF;
            if strcmp(params.country, 'EU')
                Ic_PetroleumPitch = params.Ic_PetroleumPitch_Europe_EOF;
                Ic_RawCoal = params.Ic_HardCoal_Europe_EOF;
                Ic_CrudeOil = params.Ic_CrudeOil_Europe_EOF;
                Ic_Graphite = params.Ic_Graphite_RoW_EOF;
                Ic_FreshWater = params.Ic_TapWater_Europe_EOF;
                Ic_Chlorine = params.Ic_Chlorine_Europe_EOF;
                Ic_Argon = params.Ic_Argon_Europe_EOF;
                Ic_Electricity = params.Ic_Electricity_Europe_EOF;
            end
            if strcmp(params.country, 'CN')
                Ic_PetroleumPitch = params.Ic_PetroleumPitch_RoW_EOF;
                Ic_RawCoal = params.Ic_HardCoal_CN_EOF;
                Ic_CrudeOil = params.Ic_CrudeOil_RoW_EOF;
                Ic_Graphite = params.Ic_Graphite_CN_EOF;
                Ic_FreshWater = params.Ic_TapWater_RoW_EOF;
                Ic_Chlorine = params.Ic_Chlorine_RoW_EOF;
                Ic_Argon = params.Ic_Argon_RoW_EOF;
                Ic_Electricity = params.Ic_Electricity_CN_EOF;
            end
        case 'WC'
            Ic_Ilmenite = params.Ic_Ilmenite_GLO_WC;
            Ic_PetroleumCoke = params.Ic_PetroleumCoke_GLO_WC;
            Ic_SodiumOleate = params.Ic_SodiumOleate_ROW_WC; 
            Ic_SodiumHydroxide = params.Ic_SodiumHydroxide_GLO_WC;
            Ic_Magnesium = params.Ic_Magnesium_GLO_WC;
            Ic_Aluminium = params.Ic_Aluminium_GLO_WC;
            Ic_Vanadium = params.Ic_Vanadium_CN_WC;
            Ic_NaturalGas = params.Ic_NaturalGas_GLO_WC;
            if strcmp(params.country, 'EU')
                Ic_PetroleumPitch = params.Ic_PetroleumPitch_Europe_WC;
                Ic_RawCoal = params.Ic_HardCoal_Europe_WC;
                Ic_CrudeOil = params.Ic_CrudeOil_Europe_WC;
                Ic_Graphite = params.Ic_Graphite_RoW_WC;
                Ic_FreshWater = params.Ic_TapWater_Europe_WC;
                Ic_Chlorine = params.Ic_Chlorine_Europe_WC;
                Ic_Argon = params.Ic_Argon_Europe_WC;
                Ic_Electricity = params.Ic_Electricity_Europe_WC;
            end
            if strcmp(params.country, 'CN')
                Ic_PetroleumPitch = params.Ic_PetroleumPitch_RoW_WC;
                Ic_RawCoal = params.Ic_HardCoal_CN_WC;
                Ic_CrudeOil = params.Ic_CrudeOil_RoW_WC;
                Ic_Graphite = params.Ic_Graphite_CN_WC;
                Ic_FreshWater = params.Ic_TapWater_RoW_WC;
                Ic_Chlorine = params.Ic_Chlorine_RoW_WC;
                Ic_Argon = params.Ic_Argon_RoW_WC;
                Ic_Electricity = params.Ic_Electricity_CN_WC;
            end
    end
    
    % Calculate environmental impacts
    I_MI = dependent_vars.m_ilmenite * Ic_Ilmenite;
    I_SM = dependent_vars.m_ilmenite * (dependent_vars.Enm_SM_el * Ic_Electricity + params.mm_SMpetrolPitch * Ic_PetroleumPitch + params.mm_SMpetrolCoke * Ic_PetroleumCoke + params.mm_SMrawCoal * Ic_RawCoal + params.mm_SMcrudeOil * Ic_CrudeOil + params.mm_SMgraphite * Ic_Graphite + params.mm_SMsodiumOleate * Ic_SodiumOleate);
    I_CR = dependent_vars.m_TiSlag * (dependent_vars.Enm_CR_el * Ic_Electricity + params.Enm_CR_ng * Ic_NaturalGas + params.mm_CRfreshWater * Ic_FreshWater + params.mm_CRpetrolCoke * Ic_PetroleumCoke + params.mm_CRsodiumHydrox * Ic_SodiumHydroxide + params.mm_CRchlorine * Ic_Chlorine + params.mm_CRrawCoal * Ic_RawCoal + params.mm_CRcrudeOil * Ic_CrudeOil);
    I_RD = dependent_vars.m_TiCl4 * (params.Enm_RD_el * Ic_Electricity + params.mm_RDmagnesium * Ic_Magnesium);
    I_CS = dependent_vars.m_Ti64ingot * params.Enm_CS_el * Ic_Electricity + dependent_vars.m_aluminum * Ic_Aluminium + dependent_vars.m_vanadium * Ic_Vanadium;
    I_RE = dependent_vars.m_Ti64ingot * dependent_vars.Enm_RE_el * Ic_Electricity;
    I_AT = dependent_vars.m_Ti64ingot * ((dependent_vars.Enm_ATmelt_el + dependent_vars.Enm_ATcompr_el) * Ic_Electricity + (1-params.Delta_ATargon) * dependent_vars.mm_ATargon * Ic_Argon);
    I_PS = dependent_vars.m_atomizedPowder * params.Enm_PS_el * Ic_Electricity;

    % Store results in a struct
    I_impact_breakdown = struct(...
        'MI', I_MI, ...
        'SM', I_SM, ...
        'CR', I_CR, ...
        'RD', I_RD, ...
        'CS', I_CS, ...
        'RE', I_RE, ...
        'AT', I_AT, ...
        'PS', I_PS ...
    );

    I_tot = I_MI + I_SM + I_CR + I_RD + I_CS + I_RE + I_AT + I_PS;
end

function [c, ceq] = constraintFunction(x, params)
    % Defines the constraints for optimization

    % Extract variables
    phi_ATelectrode = x(1);
    p_ATargon = x(2);
    beta_TiO2 = x(3);

    mm_ATargon = 448.82 * exp(-30.61*phi_ATelectrode);
    factor = (params.nhu_m / (params.nhu_g * params.We)) * (1 + (1/mm_ATargon));
    d50 = params.kd * phi_ATelectrode * sqrt(factor) * 1e6;
    d50 = d50-40*(p_ATargon-5.5);
    eta = 0.8 * exp(-0.02 * abs(d50-params.d_target));  
    m_atomizedPowder = params.m_finalPowder / eta;
    m_wastePowder = m_atomizedPowder - params.m_finalPowder;
    m_Ti64ingot = m_atomizedPowder;
    m_TiSponge = (m_Ti64ingot - m_wastePowder) / (1 + 0.0638 + 0.0426);
    m_aluminum = 0.0638 * m_TiSponge;
    m_vanadium = 0.0426 * m_TiSponge;
    m_TiCl4 = m_TiSponge / params.alpha_RD;
    m_TiO2 = m_TiCl4 / params.alpha_CR;
    m_TiSlag = m_TiO2 / beta_TiO2;
    m_ilmenite = m_TiSlag / params.alpha_SM;

    Enm_SM_el = 21.14 * beta_TiO2^2 - 23.16 * beta_TiO2 + 7.41;
    Enm_CR_el = -7.63 * beta_TiO2^2 + 4.65 * beta_TiO2 + 3.43;
    Enm_RE_el = 81.67 * phi_ATelectrode^2 - 26.95 * phi_ATelectrode + 3.09;
    Enm_ATmelt_el = 4364.7 * phi_ATelectrode^2 - 401.13 * phi_ATelectrode + 12.407;
    Enm_ATcompr_el = (mm_ATargon/params.M_argon)*(3.44*10^-6)*params.T_ATargon*(((p_ATargon/0.1)^0.4)-1);

    c = [
        phi_ATelectrode - params.phi_ATelectrodeMax;
        params.phi_ATelectrodeMin - phi_ATelectrode;
        p_ATargon - params.p_ATargonMax;
        params.p_ATargonMin - p_ATargon;
        -m_wastePowder;
        -m_TiSponge;                         % m_TiSponge >= 0
        -m_TiCl4;                            % m_TiCl4 >= 0
        -m_TiSlag;                           % m_TiSlag >= 0
        -m_ilmenite;                         %m_ilmenite >=0
        params.beta_TiO2Min - beta_TiO2;     % beta_TiO2 >= beta_TiO2Min;
    ];

    ceq = [
           m_TiCl4 / params.alpha_CR - m_TiO2;
           m_ilmenite - m_TiSlag / params.alpha_SM;
           m_TiCl4 - m_TiSponge / params.alpha_RD;
    ];
end

function dependent_vars = calculateDependentVariables(x_opt, params)
    % Calculates dependent variables from the optimization results.

    % Extract the optimized independent variables
    phi_ATelectrode = x_opt(1);
    p_ATargon = x_opt(2);
    beta_TiO2 = x_opt(3);

    mm_ATargon = 448.82 * exp(-30.61*phi_ATelectrode);
    factor = (params.nhu_m / (params.nhu_g * params.We)) * (1 + (1/mm_ATargon));
    d50 = params.kd * phi_ATelectrode * sqrt(factor) * 1e6;
    d50 = d50-40*(p_ATargon-5.5);
    eta = 0.8 * exp(-0.02 * abs(d50-params.d_target));  
    m_atomizedPowder = params.m_finalPowder / eta;
    m_wastePowder = m_atomizedPowder - params.m_finalPowder;
    m_Ti64ingot = m_atomizedPowder;
    m_TiSponge = (m_Ti64ingot - m_wastePowder) / (1 + 0.0638 + 0.0426);
    m_aluminum = 0.0638 * m_TiSponge;
    m_vanadium = 0.0426 * m_TiSponge;
    m_TiCl4 = m_TiSponge / params.alpha_RD;
    m_TiO2 = m_TiCl4 / params.alpha_CR;
    m_TiSlag = m_TiO2 / beta_TiO2;
    m_ilmenite = m_TiSlag / params.alpha_SM;

    Enm_SM_el = 21.14 * beta_TiO2^2 - 23.16 * beta_TiO2 + 7.41;
    Enm_CR_el = -7.63 * beta_TiO2^2 + 4.65 * beta_TiO2 + 3.43;
    Enm_RE_el = 81.67 * phi_ATelectrode^2 - 26.95 * phi_ATelectrode + 3.09;
    Enm_ATmelt_el = 4364.7 * phi_ATelectrode^2 - 401.13 * phi_ATelectrode + 12.407;
    Enm_ATcompr_el = (mm_ATargon/params.M_argon)*(3.44*10^-6)*params.T_ATargon*(((p_ATargon/0.1)^0.4)-1);

    dependent_vars = struct(...
        'm_TiO2', m_TiO2, ...
        'm_TiSlag', m_TiSlag, ...
        'm_TiCl4', m_TiCl4, ...
        'm_TiSponge', m_TiSponge, ...
        'm_Ti64ingot', m_Ti64ingot, ...
        'm_atomizedPowder', m_atomizedPowder, ...
        'mm_ATargon', mm_ATargon, ...
        'd50', d50, ...
        'eta', eta, ...
        'm_wastePowder', m_wastePowder, ...
        'm_ilmenite', m_ilmenite, ...
        'm_aluminum', m_aluminum, ...
        'm_vanadium', m_vanadium, ...
        'Enm_SM_el', Enm_SM_el, ... 
        'Enm_CR_el', Enm_CR_el, ...
        'Enm_RE_el', Enm_RE_el, ...
        'Enm_ATmelt_el', Enm_ATmelt_el, ...
        'Enm_ATcompr_el', Enm_ATcompr_el ... 
    );
end
