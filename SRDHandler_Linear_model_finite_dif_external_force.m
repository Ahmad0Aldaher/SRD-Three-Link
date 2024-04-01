classdef SRDHandler_Linear_model_finite_dif_external_force < SRDHandler
    properties
        dof_state_space_robot;
        dof_configuration_space_robot;
        dof_control;
        dof_Constraint;
        dof_External;
        
        dof_z;
        dof_Constraint_independent;
        
        last_update_q;
        last_update_v;
        last_update_u;
        
        Handler_dynamics_generalized_coordinates_model;
        Handler_dynamics_Linearized_Model;
        Handler_Constraints_Model;
        Handler_ExtenalForce_Model;
        
        Handler_State;
        Handler_Controller;
        
        finite_dif_step_z;
        finite_dif_step_zeta;
        finite_dif_step_u;
        finite_dif_step_w;

        LinearizationType;
        TemporalType;
        
        % x = [q, v]
        %dz/dt = An*z + Ar*zeta + B*u + c
        An;
        Ar;
        B;
        D
        c;
        
        N;
        R;
        
        tol = 10^(-9);
    end
    methods
        
        function obj = SRDHandler_Linear_model_finite_dif_external_force(varargin)
            Parser = inputParser;
            Parser.FunctionName = 'SRDHandler_Linear_model_finite_dif_external_force';
            Parser.addOptional('Handler_dynamics_generalized_coordinates_model', []);
            Parser.addOptional('Handler_Constraints_Model', []);
            Parser.addOptional('Handler_ExtenalForce_Model', []);
            Parser.addOptional('Handler_State',        []);
            Parser.addOptional('Handler_Controller',   []);
            Parser.addOptional('finite_dif_step_z',    0.0001);
            Parser.addOptional('finite_dif_step_zeta', 0.0001);
            Parser.addOptional('finite_dif_step_u',    0.0001);
            Parser.addOptional('finite_dif_step_w',    0.0001);
            Parser.parse(varargin{:});
            
            obj.Handler_dynamics_generalized_coordinates_model = Parser.Results.Handler_dynamics_generalized_coordinates_model;
            obj.Handler_Constraints_Model                      = Parser.Results.Handler_Constraints_Model;
            obj.Handler_ExtenalForce_Model                     = Parser.Results.Handler_ExtenalForce_Model;
            obj.Handler_State                                  = Parser.Results.Handler_State;
            obj.Handler_Controller                             = Parser.Results.Handler_Controller;
            
            if ~strcmp(obj.Handler_dynamics_generalized_coordinates_model.type, "function caller")
                warning("Handler_dynamics_generalized_coordinates_model appears to be of a wrong type")
            end
            
            obj.LinearizationType             = "Finite-difference constrained orthogonal";
            obj.TemporalType                  = "ContiniousTime";
            obj.dof_control                   = obj.Handler_dynamics_generalized_coordinates_model.dof_control;
            obj.dof_configuration_space_robot = obj.Handler_dynamics_generalized_coordinates_model.dof_configuration_space_robot;
            obj.dof_state_space_robot         = 2 * obj.dof_configuration_space_robot;
            obj.dof_Constraint                = obj.Handler_Constraints_Model.dof_Constraint;
            obj.dof_External                  = obj.Handler_ExtenalForce_Model.dof_Constraint;  
            
            
            obj.finite_dif_step_z                = Parser.Results.finite_dif_step_z;
            obj.finite_dif_step_zeta             = Parser.Results.finite_dif_step_zeta;
            obj.finite_dif_step_u                = Parser.Results.finite_dif_step_u;
            obj.finite_dif_step_w                = Parser.Results.finite_dif_step_w;

            
            %implementing serialization for arbitrary cell arrays of handlers seems to
            %be more pain than it is worth
            obj.SerializationPrepNeeded = true;
            obj.PreSerializationPrepFunction = @PreSerializationPrepFunction;
            function PreSerializationPrepFunction(~)
                error('do not attempt to save this function; create a new one on the fly instead')
            end
        end
        
        function Update(obj, ~)
            dof = obj.dof_configuration_space_robot;
            dof_ctrl = obj.dof_control;
            k = obj.dof_Constraint;
            nw = obj.dof_External;
            n = obj.dof_state_space_robot;
            
            squized = obj.Handler_State.get_position_velocity_acceleration();
            q = squized(:, 1);
            v = squized(:, 2);
            u = obj.Handler_Controller.u;
            
            x = [q; v];
            
            J  = obj.Handler_Constraints_Model.get_Jacobian(q);
            dJ = obj.Handler_Constraints_Model.get_Jacobian_derivative(q, v);
            G = [J,  zeros(k, dof); 
                 dJ, J];
            
            G_svd = svd_suit(G);
            obj.N = G_svd.null;
            obj.R = G_svd.row_space;
            obj.dof_Constraint_independent = G_svd.rank;
            obj.dof_z = n - obj.dof_Constraint_independent; 
            
            if isscalar(obj.finite_dif_step_z) 
                obj.finite_dif_step_z = eye(obj.dof_z) * obj.finite_dif_step_z;
            end
            if isscalar(obj.finite_dif_step_zeta) 
                obj.finite_dif_step_zeta = eye(obj.dof_Constraint_independent) * obj.finite_dif_step_zeta;
            end
            if isscalar(obj.finite_dif_step_u) 
                obj.finite_dif_step_u = eye(obj.dof_control) * obj.finite_dif_step_u;
            end
            if isscalar(obj.finite_dif_step_w) 
                obj.finite_dif_step_w = eye(obj.dof_External) * obj.finite_dif_step_w;
            end
            
            w0 = zeros(nw, 1);
            obj.c  = obj.N' * obj.get_acceleration(x, u, w0);

            
            dz_array_N = zeros(obj.dof_z, obj.dof_z);
            for i = 1:obj.dof_z
                xi = x + obj.N * obj.finite_dif_step_z(:, i);
                dzi = obj.N' * obj.get_acceleration(xi, u, w0);
                dz_array_N(:, i) = dzi - obj.c;
            end
            obj.An = dz_array_N * pinv(obj.finite_dif_step_z, obj.tol);
            
            
            dz_array_R = zeros(obj.dof_z, obj.dof_Constraint_independent);
            for i = 1:obj.dof_Constraint_independent
                xi = x + obj.R * obj.finite_dif_step_zeta(:, i);
                dzi = obj.N' * obj.get_acceleration(xi, u, w0);
                dz_array_R(:, i) = dzi - obj.c;
            end
            obj.Ar = dz_array_R * pinv(obj.finite_dif_step_zeta, obj.tol);
            
            
            b_array = zeros(obj.dof_z, dof_ctrl);  
            for j = 1:dof_ctrl
                uj = u + obj.finite_dif_step_u(:, j);
                dzj = obj.N' * obj.get_acceleration(x, uj, w0);
                b_array(:, j) = dzj - obj.c;
            end
            obj.B = b_array * pinv(obj.finite_dif_step_u, obj.tol);

            
            w_array = zeros(obj.dof_z, nw);  
            for ii = 1:nw
                wii = w0 + obj.finite_dif_step_w(:, ii);
                dzii = obj.N' * obj.get_acceleration(x, u, wii);
                w_array(:, ii) = dzii - obj.c;
            end
            obj.D = w_array * pinv(obj.finite_dif_step_w, obj.tol);

            
            obj.last_update_q = q;
            obj.last_update_v = v;
            obj.last_update_u = u;
        end
        
        
        function [a] = get_acceleration(obj, x, u, w)
            k = obj.Handler_Constraints_Model.dof_Constraint;
            n = obj.Handler_dynamics_generalized_coordinates_model.dof_configuration_space_robot;
            %m = obj.Handler_dynamics_generalized_coordinates_model.dof_control;
            
            q = x(1:n);
            v = x((n+1):end);
            
            H = obj.Handler_dynamics_generalized_coordinates_model.get_joint_space_inertia_matrix(q);
            bias = obj.Handler_dynamics_generalized_coordinates_model.get_bias_vector(q, v);
            T = obj.Handler_dynamics_generalized_coordinates_model.get_control_map(q);
            J = obj.Handler_Constraints_Model.get_Jacobian(q);
            dJ = obj.Handler_Constraints_Model.get_Jacobian_derivative(q, v);

            J_ext = obj.Handler_ExtenalForce_Model.get_Jacobian(q);
            %J_ext * w - generalized external force.
            
            M = [H, -J'; J, zeros(k, k)];
            vec = pinv(M) * [(T*u - bias + J_ext'*w); -dJ*v];
            a = [v; vec(1:n)];
            
        end
        
        
        
        function A = get_A(obj, ~, ~, ~, ~)
            A = obj.A;
        end
        function B = get_B(obj, ~, ~, ~)
            B = obj.B;
        end
    end
end