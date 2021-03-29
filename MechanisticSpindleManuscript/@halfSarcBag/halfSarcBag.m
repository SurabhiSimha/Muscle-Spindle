classdef halfSarcBag < handle
    
    properties
        % These properties can be accessed from the driver script
        
        % INITIAL STATE PARAMETERS %
        cmd_length = 1300;
        hs_length = 1300;   % the initial length of the half-sarcomere in nm
        hs_force;           % the stress (in N m^(-2)) in the half-sarcomere
        f_overlap;
        f_bound;
        f_activated = 0.0;
        cb_force;
        passive_force;

        % DISTRIBUTION BIN PARAMETERS %
        bin_min = -20;      % min x value for myosin distributions in nm
        bin_max = 20;       % max x value for myosin distributions in nm
        bin_width = 0.5;    % width of bins for myosin distributions in nm
        x_bins;             % array of x_bin values
        no_of_x_bins;       % no of x_bins
        
        % INDIVIDUAL CROSSBRIDGE PARAMETERS %
        k_cb = 0.001;       % Cross-bridge stiffness in N m^-1
        power_stroke = 2.5;   % Cross-bridge power-stroke in nm
        
        % PARAMETERS RELATED TO FORWARD AND REVERSE RATES %
        f_parameters = 5e2;
        g_parameters = [0.3 0.01 4000 50];
                            

        thick_filament_length = 815;
                            % length of thick filaments in nm
        thin_filament_length = 1120;
                            % length of thin filaments in nm
        bare_zone_length = 80;
                            % length of thick filament bare zone in nm
        k_falloff = 0.002;  % defines how f_overlap falls hs_length shortens
                            % below optimal
                            
        f;                  % forward rates
        g;                  % reverse rates
        
        bin_pops;           % number of heads bound in each bin
        no_detached;        % number of heads not attached to a binding site
        
        max_rate = 5000;    % clip f or g values above this value
        
        compliance_factor = 0.50;
                            % proportion of delta_hsl that the
                            % cb distribution is moved
                            
        cb_number_density = 6.9e16;
                            % number of cbs in a half-sarcomere with a
                            % cross-sectional area of 1 m^2
                            
        hsl_slack = 1050;   % slack length of half-sarcomere in nm
        k_passive = 100;   % passive stiffness of half-sarcomere in
                            % N m^-2 nm^-1

    end
    
    methods
        
        % BUILD halfSarcBag OBJECT %
        function obj = halfSarcBag(varargin)
            
            % Set up x_bins
            obj.x_bins = obj.bin_min:obj.bin_width:obj.bin_max;
            obj.no_of_x_bins = numel(obj.x_bins);
            
            % Set up rates
            %%% D --> A rate (symmetric) %%%
            obj.f = zeros(size(obj.x_bins)); %Preallocate
            obj.f = obj.f_parameters(1) * obj.bin_width * ...
                exp(-obj.k_cb*10*((obj.x_bins).^2)/(1e18*1.381e-23*288));
            
            %%% A --> D rate (asymmetric) %%%
            obj.g = zeros(size(obj.x_bins)); %Preallocate
            obj.g(obj.x_bins<-6) = obj.g_parameters(1) + ...
                 abs(0.2*((obj.x_bins(obj.x_bins<-6)+6).^3));
            obj.g(obj.x_bins>=-3) = obj.g_parameters(1) + ...
                 0.3*((obj.x_bins(obj.x_bins>=-3)+3).^3);
            obj.g = obj.g + 1;
            
            % Limit max values
            obj.f(obj.f>obj.max_rate) = obj.max_rate;
            obj.g(obj.g>obj.max_rate) = obj.max_rate;
            
            % Initialize bins
            obj.bin_pops = zeros(obj.no_of_x_bins,1);
            
        end
        
        % Other methods
        function update_filamentOverlap(obj)
            
            x_no_overlap = obj.hs_length - obj.thick_filament_length; %Length of thin filament w/0 overlapping thick filament
            x_overlap = obj.thin_filament_length - x_no_overlap; %
            max_x_overlap = obj.thick_filament_length -  ...
                obj.bare_zone_length; %Region of thick filament containing myosin heads
            
            if (x_overlap<0) %This is impossible in a sarcomere
                obj.f_overlap=0;
            end
            
            if ((x_overlap>0)&&(x_overlap<=max_x_overlap)) %Operating range of half sarcomere
                obj.f_overlap = x_overlap/max_x_overlap;
            end
            
            if (x_overlap>max_x_overlap) %This is impossible in a sarcomere
                obj.f_overlap=1;
            end
            
            if (obj.hs_length<obj.thin_filament_length) %This doesn't happen unless hsl < thin filament length
                obj.f_overlap = 1.0 + obj.k_falloff * ...
                    (obj.hs_length - obj.thin_filament_length);
                if (obj.f_overlap < 0)
                    obj.f_overlap = 0;
                end
            end            
        end
        
        
        function update_fracBound(obj)
            obj.f_bound = sum(obj.bin_pops);
        end
        
        function update_thinFilament(obj)
                        
             obj.f_activated = obj.f_activated;
            
        end
        
        
        
        
        function evolve_cbDist(obj,time_step)
            
            % Construct a vector y where
            % y(1) is the number of cbs in the detached state
            % y(2 ... no_of_x_bins+1) is the number of cbs in bins 1 to no_of_x_bins
            
%             obj.no_detached = max([0 obj.f_overlap-obj.f_bound]);
            obj.no_detached = max([0 obj.f_activated-obj.f_bound]);
%             obj.no_detached = max([0
%             min(obj.f_activated,obj.f_overlap)-obj.f_bound]); %USE THIS
%             TO ESTIMATE change in force with filament overlap
            y = [obj.no_detached ; obj.bin_pops];
            
            % Evolve the system
            [~,y_new] = ode23(@derivs,time_step*[0 0.5 1],y,[]);
            
            % Update the bound cross-bridges
            obj.bin_pops = y_new(end,2:end)';
            obj.no_detached = obj.no_detached;
            % Nested function
            function dy = derivs(~,y)
                n=numel(y);
                obj.bin_pops = y(2:end)';
                dy(1) = sum(obj.g .* obj.bin_pops) - obj.no_detached * sum(obj.f);
                dy(2:n) = (obj.f * obj.no_detached) - obj.g .* obj.bin_pops;
                dy=dy';
            end
        end
        
        
        function shift_cbDist(obj,delta_x)
            % Adjust for filament compliance
            delta_x = delta_x * obj.compliance_factor;
            % Shift populations by interpolation
            interp_positions = obj.x_bins - delta_x;
            obj.bin_pops = interp1(obj.x_bins,obj.bin_pops,interp_positions, ...
                'linear',0)';
        end
        

        function calcForces(obj)
            obj.cb_force = obj.cb_number_density * obj.k_cb * 1e-9 * ...
                sum((obj.x_bins + obj.power_stroke) .* obj.bin_pops');
            obj.passive_force = obj.k_passive * (obj.hs_length - obj.hsl_slack);
            obj.hs_force = obj.cb_force + obj.passive_force;
        end

        
        
        function forwardStep(obj,time_step,delta_hsl,delta_cdl,delta_f_activated,evolve,shift)
            %This function uses the methods for updating the half-sarcomere
            %every time step
                        
            obj.cmd_length = obj.cmd_length + delta_cdl;
            obj.hs_length = obj.hs_length + delta_hsl;
            
%             obj.Ca = obj.Ca + delta_Ca;
            obj.f_activated = obj.f_activated + delta_f_activated;


            
            % Change cb distribution based on cycling kinetics
            if evolve == 1
                                % Calculate the change in f_activated
                obj.update_filamentOverlap();
                obj.update_fracBound();
                obj.update_thinFilament();
                obj.evolve_cbDist(time_step);

            end
            
            % Shift cb distribution based on imposed movements
            % Also, perform calculations that are dependent on hs_length
            if shift == 1
                obj.shift_cbDist(delta_hsl);
            end
                            
            % Calculate forces
            obj.calcForces;
            
        end
        
    end
end      
            
            
            
            
        
        