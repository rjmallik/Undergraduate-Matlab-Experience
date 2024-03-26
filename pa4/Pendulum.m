classdef Pendulum

    properties
        % Distance from the pivot to the center of mass.
        length_to_COM = 1.0 % meters
        mass = 1.0 % kilograms
        
        % Acceleration due to gravity.
        g = -9.8 % m/s^2

        % Coefficient of velocity dampening.
        dampening = 0.10;

        moment_of_inertia = 0.5;

        % Coefficients to determine how noisey the measurements are.
        theta_noise_scale
        omega_noise_scale
    end

    methods
        function obj = Pendulum(varargin)
            % Constructor

            default_moment_of_inertia = 0.5;
            default_theta_noise_scale = 0.01;
            default_omega_noise_scale = 0.01;

            p = inputParser;
            validScalarNonnegativeNum = @(x) isnumeric(x) && isscalar(x) && (x >= 0);
            addParameter(p,'theta_noise_scale',default_theta_noise_scale,validScalarNonnegativeNum);
            addParameter(p,'omega_noise_scale',default_omega_noise_scale,validScalarNonnegativeNum);
            % addParameter(p,'dampening',default_dampening,validScalarNonnegativeNum);
            % addParameter(p,'gravity',-9.8);
            addParameter(p, 'moment_of_inertia', default_moment_of_inertia, validScalarNonnegativeNum)
            addParameter(p,'location', 'earth');

            parse(p,varargin{:});
            
            % Set the object properties.
            % obj.dampening = p.Results.dampening;
            obj.theta_noise_scale = p.Results.theta_noise_scale;
            obj.omega_noise_scale = p.Results.omega_noise_scale;
            % obj.g = p.Results.gravity;
            
            % if ischar(p.Results.gravity) || isstring(p.Results.gravity)
            switch lower(p.Results.location)
                case 'sun'
                    obj.g = -274; % m/s^2
                    obj.dampening = 6; 
                case 'mercury'
                    obj.g = -3.7; % m/s^2
                    obj.dampening = 1; 
                case 'venus'
                    obj.g = -8.87; % m/s^2
                    obj.dampening = 1.2;
                case 'moon'
                    obj.g = -1.62; % m/s^2
                    obj.dampening = 0.0001; 
                case 'earth'
                    obj.g = -9.8; % m/s^2
                    obj.dampening = 0.1;  
                case 'mars'
                    obj.g = -1.62; % m/s^2
                    obj.dampening = 0.01; 
                case 'jupiter'
                    obj.g = -24.79; % m/s^2
                    obj.dampening = 2.5; 
                case 'saturn'
                    obj.g = -10.44; % m/s^2
                    obj.dampening = 0.05; 
                case 'slaturn' % Home of Slammy the Slug of Slanta Cruz.
                    obj.g = -pi; % m/s^2
                    obj.dampening = 0.0; 
                case 'uranus'
                    obj.g = -8.69; % m/s^2
                    obj.dampening = 0.6; 
                case 'neptune'
                    obj.g = -11.15; % m/s^2
                    obj.dampening = 0.3; 
                case 'pluto'
                    error('Sorry, Pluto is not a planet :''(')
                otherwise
                    error(['Planet name not recognized: ', gravity])
            end
            % else
            %     error(['Unrecognized type for ''location'': ', p.Results.location])
            % end
        end

        function [thetas, omegas] = simulateAndMeasure(this, t_grid, x0)
            % Compute the evolution of theta and omega.
            [~,x] = ode45(@(t,x) this.flow_map(x), t_grid, x0);

            % Generate random measurement noise.
            theta_noise = this.theta_noise_scale*randn([size(x,1),1]);
            omega_noise = this.omega_noise_scale*randn([size(x,1),1]);

            % Add noise to output.
            thetas = x(:,1) + theta_noise;
            omegas = x(:,2) + omega_noise;
        end

        function plot(this, theta)
            hold on
            % Plot some empty points to make the 
            axis equal
            xlim(2.1*this.length_to_COM*[-1, 1])
            ylim(2.1*this.length_to_COM*[-1, 1])
            % plot(-2.1*this.length_to_COM, -2.1*this.length_to_COM, 'w')
            % plot( 2.1*this.length_to_COM,  2.1*this.length_to_COM, 'w')

            % Create vector pointing to the center of mass.
            center_of_mass = this.length_to_COM*[sin(theta);-cos(theta)];

            % Plot the pendulum arm as a line.
            plot([0, 2*center_of_mass(1)], [0, 2*center_of_mass(2)], 'b')
            hold on
            % Plot the center of mass.
            plot(center_of_mass(1), center_of_mass(2), 'bo', 'MarkerSize', 5, 'MarkerFaceColor', 'b')

            lgd = legend('Pendulum Arm', 'Center of mass');
            lgd.Location = 'eastoutside';
            xlabel('x')
            ylabel('y')
            title('Pendulum')
        end
        
        function x_dot = flow_map(this, x)
            % moment_of_inertia = 0.54321 * this.mass * (2*this.length_to_COM);
        
            I = this.moment_of_inertia;
        
            theta = x(1);
            omega = x(2);
            theta_dot = omega;
            torque = this.length_to_COM*sin(theta)*this.mass*this.g;
            omega_dot = torque/I  - this.dampening*omega;
            x_dot = [theta_dot; omega_dot];
        end

    end

end