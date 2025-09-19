classdef LorenzSystemSim < handle
    % Lorenz system simulation.
    %       Wiki: https://en.wikipedia.org/wiki/Lorenz_system
    %       Original Article: https://journals.ametsoc.org/view/journals/atsc/20/2/1520-0469_1963_020_0130_dnf_2_0_co_2.xml
    %       Tau sims: https://arxiv.org/pdf/2207.00521.pdf
    %       Rho exp: https://arxiv.org/pdf/2207.00521.pdf
    %       Rho sin: https://pubs.aip.org/aip/cha/article/31/3/033149/342213/Using-machine-learning-to-predict-statistical
    % Terminology: 
    %       Time span or time range: duration of time over which the simulation runs.
    %       Time step or delta time: incremental time interval the simulation progresses.
    %       Initial conditions (IC): initial values of the dependent variables and derivatives
    % USAGE:
    %       Initialize the constructor with:
    %           Needed: sigma, beta, rho, rho_0, rho_1, tau, tspan, xinit
    %           Optional: non_st, noise
    %       Call the Run method.
    %       Call the save to CSV method.
    %       See run_lorenz_sim script for an usage example!
    %
    % In: system parameters:
    %           sigma: Prandtl number
    %           beta: Beta parameter > 0
    %           rho: Rayleigh number
    %           rho_0, rho_1, tau: For non-stationarity by using time-dependent rho
    %               Additional info: https://arxiv.org/pdf/2207.00521.pdf
    %           tspan: Time span.
    %           xinit: Initial conditions for x, y and z.
    % In: optional [default 1]:
    %           non_st: To use time-dependent rho param
    %           noise: To add white noise (gaussian) to the functions,    
    %
    % Out: run_simple(NO ARGS):
    %          t, obj(tspan,3): t and x, y, z values over the tspan.
    %      run_extended_random_ic (nr sim,  range_x, range_y, range_z):
    %          Initial conditions for x,y,z random in ranges.
    %          2D array: nr_sim, nr_subsim, t, x, y and z.
    %      run_extended_random_ic_function(nr_sim, mode)
    %          x,y,z randomly initialized from function values.
    %          functions: sin, cos, exp, time-dependent.     
    %          2D array: nr_sim, nr_subsim, t, x, y and z.
    %      run_extended_rho_time_dependent (nr_sim, x_range, y_range, z_range, rho_0_range, rho_0_step)
    %           rho_0 in range [rho_0_range(1) rho_0_range(2)].
    %           Initial conditions for x,y,z random in range [0 1].
    %           rho_0_range, rho_0_step optional, default will run only once.
    %           2D array: nr_sim, nr_subsim, t, x, y, z, rho_0 value.
    %      run_extended_rho_dynamic(rho_min, rho_max)
    %          rho in a range, used for plots. 
    %          Each sim init random IC for x,y and z.
    %
    % CSV: writeToCSVFile (data, file)
    %       data: obtained by running a sim mode (see above).
    %       file: desired CSV file name.
    %   
    %  ! All Simulations output a 2D array containing:
    %       Simulation nr, Subsimulation nr, t, x, y and z.
    %       For time-dependent rho there will be an extra column!!!
    %       Each simulation contains tspan time-steps (observations).

    properties(Access = private)
        % Add white noise to x,y,z.
        noise = 1;
        noise_mean = 0;
        noise_std = 1;
        noise_min = -1;
        noise_max = 1;
        
        % Nonlinear time-dependent rho param for non-stationarity (T/F)!
        rho_non_st = 0;
        
        % Step change rho
        step_rho = 0;

        % Params defaults.
        sigma, beta, rho;
        rho_values = 0;
        rho_change = [];

        % Params for time-dep rho
        rho_0, rho_1, tau;
        gama, mu;

        % Time span and time vars
        tspan
        delta_rho = 0.01

        % Functions values
        x, y, z, t;
        t_limits = 0

        % Initial parameter values
        xinit

        % Run Mode
        mode

        % Counter for continuous generation
        inner_sim_counter = 0

        % Others
        step, temp_rho, rho_column = [];
        counter = 1;
        segment_lengths = 0;
    end
    
    methods (Access = public)

        function obj = LorenzSystemSim(varargin)
            if (nargin < 7)
                throw(MException('Id:id','Arguments Error \n sigma, beta, rho, rho_0, rho_1, tau, tspan, xinit, non_st, noise, gama, mu'));
            elseif (nargin == 7)
                % sigma, beta, rho, tspan, xinit, noise, step_rho);
                obj.sigma = varargin{1};
                obj.beta = varargin{2};
                obj.rho = varargin{3};
                obj.tspan = varargin{4};
                obj.xinit = varargin{5};
                obj.noise = varargin{6};
                obj.step_rho = varargin{7};
            elseif (nargin >= 8 && nargin <= 13)
                obj.sigma = varargin{1};
                obj.beta = varargin{2};
                obj.rho = varargin{3};
                obj.rho_0 = varargin{4};
                obj.rho_1 = varargin{5};
                obj.tau = varargin{6};
                obj.tspan = varargin{7};
                obj.xinit = varargin{8};
                if (nargin == 13)
                    obj.rho_non_st = varargin{9};
                    obj.noise = varargin{10};
                    obj.step_rho = varargin{11};
                    obj.gama = varargin{12};
                    obj.mu = varargin{13};
                end
            else
                throw(MException('Id:id','Arguments Error, will be extended in future'));
            end
        end
               
        function [t, outputArg] = run_simple(obj)

            [t,xf] = ode45(@(t,xf) obj.lorenz(t, xf), obj.tspan, obj.xinit);

            obj.x = xf(:,1);
            obj.y = xf(:,2);
            obj.z = xf(:,3);
            obj.t = t;

            outputArg(:,1) = obj.x;
            outputArg(:,2) = obj.y;
            outputArg(:,3) = obj.z;
            
            obj.mode = 1;
        end

        function out = run_extended_random_ic(obj, nr_sim, range_x, range_y, range_z)

            for i = 1:nr_sim

                % Generate random initial conditions
%                 x = (range_x(2) - range_x(1)).*rand + range_x(1);
%                 y = (range_y(2) - range_y(1)).*rand + range_y(1);
%                 z = (range_z(2) - range_z(1)).*rand + range_z(1);

                % Generate random initial conditions
                x = randn;
                y = randn;
                z = randn;

                obj.xinit = [x,y,z];

                [t,xf] = ode45(@(t,xf) obj.lorenz(t, xf), obj.tspan, obj.xinit);
                
                simulation_nr = zeros(size(xf(:,1),1),1) + i;

                sim(:,1) = simulation_nr;
                sim(:,2) = simulation_nr;
                sim(:,3) = t;
                sim(:,4) = xf(:,1);
                sim(:,5) = xf(:,2);
                sim(:,6) = xf(:,3);
                
                if (i==1)
                    out = sim;
                else
                    out = [out; sim];
                end

            end
            obj.mode = 2;
         end
        
        function out = run_extended_random_ic_rho_intervals(obj, nr_sim, range_x, range_y, range_z, rho_range, decreasing)
            % Rho 25, 20
            j = 1;
            for r = rho_range(1) : rho_range(2) : rho_range(3)
                obj.rho = r;

                % Generate training and testing data in that range. 
                for i = 1:nr_sim
                
                    % Generate random initial conditions
                    x = (range_x(2) - range_x(1)).*rand + range_x(1);
                    y = (range_y(2) - range_y(1)).*rand + range_y(1);
                    z = (range_z(2) - range_z(1)).*rand + range_z(1);

%                     x = randn;
%                     y = randn;
%                     z = randn;
            
                     obj.xinit = [x,y,z];
            
                     [t,xf] = ode45(@(t,xf) obj.lorenz(t, xf), obj.tspan, obj.xinit);
            
                    simulation_nr = zeros(size(xf(:,1),1),1) + i;
            
                    sim(:,1) = simulation_nr;
                    sim(:,2) = simulation_nr;
                    sim(:,3) = t;
                    sim(:,4) = xf(:,1);
                    sim(:,5) = xf(:,2);
                    sim(:,6) = xf(:,3);
                    sim(:,7) = obj.rho;
            
                    if (i == 1 && j == 1)
                        out = sim;
                    else
                        out = [out; sim];
                    end
            
                end
                j = j + 1;
            end
            obj.mode = 7;
            j = 1;
         end
        
        function out = run_extended_random_ic_rho_intervals_cont(obj, nr_sim, range_x, range_y,...
                                                                    range_z, rho_range, decreasing, ...
                                                                    random_shuffle)
            j = 1; 
            
            if decreasing == true
                rho_values_inc = rho_range(1):rho_range(2):rho_range(3);
                rho_values_dec = rho_values_inc(end:-1:1);
                obj.rho_values = [rho_values_inc(1:end-1) rho_values_dec];
            else
                obj.rho_values = rho_range(1):rho_range(2):rho_range(3);
            end

            num_segments = length(obj.rho_values);
            obj.segment_lengths = int32(size(obj.tspan,2) / num_segments);
            
            % obj.t_limits = linspace(obj.tspan(1), obj.tspan(end), num_segments);
            obj.t_limits = round(1 : obj.segment_lengths : size(obj.tspan,2));
            obj.t_limits = obj.tspan( obj.t_limits);

            % Generate training and testing data in that range. 
            for i = 1:nr_sim
                if random_shuffle
                    obj.rho_values = obj.rho_values(randperm(length(obj.rho_values)));
                end

                obj.counter = 1;
                obj.inner_sim_counter = 2;
                obj.rho_column = [];
                obj.rho = obj.rho_values(1);
                obj.rho_change = [];

                if i == 1
                    % Generate random initial conditions
                    x = (range_x(2) - range_x(1)).*rand + range_x(1);
                    y = (range_y(2) - range_y(1)).*rand + range_y(1);
                    z = (range_z(2) - range_z(1)).*rand + range_z(1);
                else
                    x = sim(end,4);
                    y = sim(end,5);
                    z = sim(end,6);

                end
        
                obj.xinit = [x,y,z];
                
                [t,xf] = ode45(@(t,xf) obj.lorenz_cont(t, xf), obj.tspan, obj.xinit);
        
                simulation_nr = zeros(size(xf(:,1),1),1) + i;
                                    
                obj.rho_change = [obj.rho_change, t(end)];
                edges = [1, obj.rho_change];  
                segment_idx = discretize(t, edges);
                rho_column = obj.rho_values(segment_idx);
                
                sim(:,1) = simulation_nr;
                sim(:,2) = simulation_nr;
                sim(:,3) = t;
                sim(:,4) = xf(:,1);
                sim(:,5) = xf(:,2);
                sim(:,6) = xf(:,3);
                sim(:,7) = rho_column;
        
                if (i == 1 && j == 1)
                    out = sim;
                else
                    out = [out; sim];
                end
        
            end
            j = j + 1;
            obj.mode = 7;
          end

        function out = run_extended_random_ic_function(obj, nr_sim, mode)
                outputArg = {};
                for i = 1:nr_sim
                    points = 200;
                    % Generate random initial conditions
                    if (mode == "sin")
                        rands = linspace(1, 2*pi, points);
                        rands = sin(rands);
                    elseif(mode == "cos")
                        rands = linspace(1, 2*pi, points);
                        rands = cos(rands);
                    elseif(mode == "add_others")
                        % TODO: Add othhers if needed.
                        rands = linspace(1, 2*pi, points);
                        rands = exp(rands);
                    end

                    simulation = {};

                    for j = 1:length(rands)
                        
                        x = rands(j);
                        y = rands(j);
                        z = rands(j);
                        
                        obj.xinit = [x,y,z];
        
                        [t,xf] = ode45(@(t,xf) obj.lorenz(t, xf), obj.tspan, obj.xinit);

                        simulation_nr = zeros(size(xf(:,1),1),1) + i;
                        subsim_nr = zeros(size(xf(:,1),1),1) + j;

                        sim(:,1) = simulation_nr;
                        sim(:,2) = subsim_nr;
                        sim(:,3) = t;
                        sim(:,4) = xf(:,1);
                        sim(:,5) = xf(:,2);
                        sim(:,6) = xf(:,3);
                        
                        if (i == 1 && j == 1)
                            out = sim;
                        else
                            out = [out; sim];
                        end
                        %simulation{j,1} = sim;

                    end
                    %outputArg{i,1} = sim;
                end

                obj.mode = 3;
            end

        function out = run_extended_rho_time_dependent(varargin)
            if (nargin <= 4)
                throw(MException('Id:id','Arguments Error '));
            elseif (nargin > 4)
                obj = varargin{1};
                nr_sim = varargin{2};
                range_x = varargin{3};
                range_y = varargin{4};
                range_z = varargin{5};
                if (nargin == 7)
                    rho_0_range = varargin{6};
                    rho_0_step = varargin{7};
                else
                    rho_0_range = [1 1];
                    rho_0_step = 1;
                end
            end
            
            if (obj.step_rho ~= 0)
                 throw(MException('Id:id','Disable rho_step arg for time dependent rho values ! '));
            end
            if (obj.rho_non_st == 0)
                 throw(MException('Id:id','Set rho_non_st to valid values [1, 2, ...]. '));
            end
            j = 1;
            for rho_0 = rho_0_range(1):rho_0_step:rho_0_range(2) 
                obj.rho_0 = rho_0;
                
                for i = 1:nr_sim
                    
                    % Generate random initial conditions
                    x = (range_x(2) - range_x(1)).*randn + range_x(1);
                    y = (range_y(2) - range_y(1)).*randn + range_y(1);
                    z = (range_z(2) - range_z(1)).*randn + range_z(1);
                    
                    obj.xinit = [x,y,z];
    
                    [t,xf] = ode45(@(t,xf) obj.lorenz(t, xf), obj.tspan, obj.xinit);
                    
                    simulation_nr = zeros(size(xf(:,1),1),1) + i;
                    subsim_nr = zeros(size(xf(:,1),1),1) + j;

                    sim(:,1) = subsim_nr;
                    sim(:,2) = simulation_nr;
                    sim(:,3) = t;
                    sim(:,4) = xf(:,1);
                    sim(:,5) = xf(:,2);
                    sim(:,6) = xf(:,3);
                    sim(:,7) = rho_0;

                    if (j == 1 && i == 1)
                        out = sim;
                    else
                        out = [out; sim];
                    end
                end
                j = j + 1;
            end
            obj.mode = 4;
        end
    
        function out = run_extended_step_change_rho(obj, nr_sim, range_x, range_y, range_z)

            if (obj.step_rho == 0)
                 throw(MException('Id:id','step_rho arg must be 1 and rho_range must valid [min, step, max]! '));
            end
            if (obj.rho_non_st ~= 0)
                 throw(MException('Id:id','rho values can not be step changes and time dependent. \n Set rho_non_st to 0. '));
            end

            for i = 1:nr_sim
                % Generate random initial conditions
                x = (range_x(2) - range_x(1)).*randn + range_x(1);
                y = (range_y(2) - range_y(1)).*randn + range_y(1);
                z = (range_z(2) - range_z(1)).*randn + range_z(1);

                obj.xinit = [x,y,z];
                
                [t,xf] = ode45(@(t,xf) obj.lorenz(t, xf), obj.tspan, obj.xinit);
                
                simulation_nr = zeros(size(xf(:,1),1),1) + i;

                sim(:,1) = simulation_nr;
                sim(:,2) = simulation_nr;
                sim(:,3) = t;
                sim(:,4) = xf(:,1);
                sim(:,5) = xf(:,2);
                sim(:,6) = xf(:,3);
                
                if (i==1)
                    out = sim;
                else
                    out = [out; sim];
                end

            end
            obj.mode = 6;
        end
        
        function outputArg = run_extended_rho_dynamic(obj, rho_min, rho_max)
            
            % Used for ploting rho vs z maxima.

            % Using machine learning to predict statistical properties of
            % non-stationary dynamical processes: System climate,regime
            % transitions, and the effect of stochasticity
            % outputArg = {};
            
            z_m = 0;
            idx = 1;
            for rho = rho_min : obj.delta_rho : rho_max
                
               obj.rho = rho;
               for j = 1:1000

                     x = randn;
                     y = randn;
                     z = randn;
                    
                     obj.xinit = [x,y,z];

                    [t,xf] = ode45(@(t,xf) obj.lorenz(t, xf), obj.tspan, obj.xinit);
                    ex = 0;    
                    for i=1:size(xf,1)
                        x_t = xf(i,1);
                        y_t = xf(i,2);
                        z_t = xf(i,3);

                        ex(i) = (x_t)^2 *(obj.rho - z_t) ...
                             + obj.sigma * (y_t)^2 ... 
                             - x_t*y_t*(obj.sigma + obj.beta + 1) ...
                             - obj.beta^2*z_t;
                        if (ex(i) < 0)
                            % if (round(x_t*y_t) ==  round(obj.beta * z_t))
                            if (x_t*y_t < (obj.beta * z_t) + 0.2 && x_t*y_t > (obj.beta * z_t) - 0.2 )
                                z_m(1,end) = z_t;
                                z_m(2,end+1) = rho;
                            end
                        end
                    end  
               end
            end
            outputArg = z_m;  
            obj.mode = 5;
        end
        
        function [outputArg2, outputArg] = getOutputValues(obj)

            outputArg = zeros(length(obj.x), 3);

            outputArg(:,1) = obj.x;
            outputArg(:,2) = obj.y;
            outputArg(:,3) = obj.z;
            
            outputArg2 = obj.t;
        end

        function out = getSegmenLength(obj)
            out = obj.segment_lengths;
        end

        function out = getRhoCols(obj)
            out = obj.rho_column;
        end

        function writeToCSVFile(obj, data, file)
            switch obj.mode
                case 1
                    % Simple Sim 
                    % CSV: Simulation, T, X, Y, Z - 1 simulation.
                    sim = zeros(size(data,1),1);
                    writematrix([sim, sim,  data], file);
                case 2
                    % Random Initial Conditions
                    % CSV: Simulation, T, X, Y, Z - N Simulations.
                    writematrix(data, file);
                case 3
                    % Function for Initial Conditions
                    % CSV: Simulation, T, X, Y, Z - N Simulations.
                    writematrix(data, file);
                case 4
                    % Time-dependent Rho values
                    writematrix(data, file);
                case 5
                    % Rho in fixed range for z_m computation (Plots Only)
                    writematrix(data, file);
                case 6
                    % Rho step change 
                    writematrix(data, file);
                case 7
                    % Rho step change 
                    writematrix(data, file);
            end
    
        end

    end

    methods (Access = private)

        function dxdt = lorenz(obj, t, xf)
    
            % Functions
            x = xf(1);
            y = xf(2);
            z = xf(3);
            
            dxdt = zeros(size(x));
                   
            % Random white 33
            if (obj.noise == 1)
                eps = obj.noise_mean + obj.noise_std * randn(1, 3);
            elseif (obj.noise == 2)
                eps = (obj.noise_max - obj.noise_min) * rand(1, 3) + obj.noise_min;
            else
                eps = [0; 0; 0];
            end
            
            % Nonlinear time-dependant poarameter rho
            if (obj.rho_non_st == 1 && obj.step_rho == 0)   
                rho = obj.rho_0  + obj.rho_1 * exp(t/obj.tau);
            elseif (obj.rho_non_st == 2 && obj.step_rho == 0)
                rho = obj.rho_0 + obj.rho_1 * sin(obj.gama * t) + (obj.mu * t);
            elseif (obj.rho_non_st == 0 && obj.step_rho == 0)
                rho = obj.rho;
            end
                      
            % Step change rho

            if (obj.step_rho == 1 && obj.rho_non_st == 0)

               ranges_a = [1, 5.999; 30, 35.999; 54, 59.999; 78, 83.999];
               ranges_b = [6, 11.999; 36, 41.999; 60, 65.999; 84, 89.999];  
               ranges_c = [12, 17.999; 42, 47.999; 66, 71.999; 90, 95.999];
               ranges_d = [24, 29.999; 48, 53.999; 72, 77.999; 96, 101];

               is_in_a = any(t >= ranges_a(:, 1) & t <= ranges_a(:, 2));
               is_in_b = any(t >= ranges_b(:, 1) & t <= ranges_b(:, 2));
               is_in_c = any(t >= ranges_c(:, 1) & t <= ranges_c(:, 2));
               is_in_d = any(t >= ranges_d(:, 1) & t <= ranges_d(:, 2));

               if is_in_a
                   rho = obj.rho;  
               elseif is_in_b
                   rho = obj.rho + 0.15 * obj.rho;   
               elseif is_in_c
                   rho = obj.rho;  
               elseif is_in_d
                   rho = obj.rho - 0.15 * obj.rho;
               end
               temp_rho = rho;
            end

            % ODE's
            dxdt(1) = obj.sigma*(y - x) + eps(1);
            dxdt(2) = x*(rho - z) - y + eps(2);
            dxdt(3) = x*y - obj.beta * z + eps(3);
            
            dxdt = dxdt';
        end
    
        function dxdt = lorenz_cont(obj, t, xf)
    
            % Functions
            x = xf(1);
            y = xf(2);
            z = xf(3);
            
            dxdt = zeros(size(x));
                   
            % Random white 33
            if (obj.noise == 1)
                eps = obj.noise_mean + obj.noise_std * randn(1, 3);
            elseif (obj.noise == 2)
                eps = (obj.noise_max - obj.noise_min) * rand(1, 3) + obj.noise_min;
            else
                eps = [0; 0; 0];
            end

            % based on t switch rho value in the same simulation
            % Determine rho based on time interval

            % idx = discretize(t, obj.t_limits);
            % current_rho = obj.rho_values(idx);
            % obj.rho = current_rho;

            if (obj.inner_sim_counter <= size(obj.rho_values,2) && t >= obj.t_limits(obj.inner_sim_counter) && ...
                    t < obj.tspan(end))
                 

                rho = obj.rho_values(obj.inner_sim_counter);
                obj.rho = rho;
                obj.inner_sim_counter = obj.inner_sim_counter + 1;
                obj.rho_change(end+1) = t;

            end

            % obj.rho_column(end+1) = obj.rho;
            
            % ODE's
            dxdt(1) = obj.sigma*(y - x) + eps(1);
            dxdt(2) = x*(obj.rho - z) - y + eps(2);
            dxdt(3) = x*y - obj.beta * z + eps(3);
            
            dxdt = dxdt';
        end

    end
end

