close all; clear; clear classes; clc;


Handler_IK_Solution = SRD_get('Handler_IK_Solution');

qva = Handler_IK_Solution.get_position_velocity_acceleration(0);

Handler_State = SRDHandler_State(...
    'InitialPosition', qva(:, 1) + 0*randn(size(qva, 1), 1), ...
    'InitialVelocity', qva(:, 2) + 0*randn(size(qva, 1), 1));
Handler_State_StateSpace = SRDHandler_StateConverter_GenCoord2StateSpace(...
    'Handler_State', Handler_State);
% %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Handler_dynamics_generalized_coordinates_model = SRD_get('Handler_dynamics_generalized_coordinates_model');
% Handler_dynamics_Linearized_Model = SRD_get('Handler_dynamics_Linearized_Model');
Handler_Constraints_Model  = SRD_get('Handler_Constraints_Model');
Handler_ExtenalForce_Model = SRD_get('Handler_ExtenalForce_Model');
 

Handler_dynamics_GC_model_evaluator = SRDHandler_dynamics_GC_model_evaluator(...
    'Handler_dynamics_generalized_coordinates_model', Handler_dynamics_generalized_coordinates_model, ...
    'Handler_State', Handler_State, ...
    'UsePinv', true);


%% Linear model finite difference A and B
Handler_Linear_model_orthogonal = SRDHandler_Linear_model_finite_dif_external_force(...
     'Handler_dynamics_generalized_coordinates_model', Handler_dynamics_generalized_coordinates_model, ...
    'Handler_Constraints_Model', Handler_Constraints_Model, ...
    'Handler_ExtenalForce_Model', Handler_ExtenalForce_Model, ...
    'Handler_State', Handler_State, ...
    'Handler_Controller', [], ...
    'finite_dif_step_z', 0.0001, 'finite_dif_step_zeta', 0.0001, 'finite_dif_step_u',0.00001);


%%
Handler_Constraints_Model.Handler_dynamics_generalized_coordinates_model = ...
    Handler_dynamics_generalized_coordinates_model;

C = [1,0,0,0,0,0;
     0,1,0,0,0,0;
     0,0,1,0,0,0;
     0,0,0,1,1,0];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %

dt = 0.001;
tf = Handler_IK_Solution.TimeExpiration;
% tf = 1.0;

Handler_Time = SRDHandler_Time('TimeLog', 0:dt:tf);

% %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Handler_ObserverState = SRDHandler_StateSpace(...
    'InitialState', [qva(:, 1); ...
                     qva(:, 2)]);
Handler_ObserverState_genCoord = SRDHandler_StateConverter_StateSpace2GenCoord(...
    'Handler_StateSpace', Handler_ObserverState); 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %

Handler_MeasuredOutput = SRDHandler_LinearMeasuredOutput(...
    'Handler_State_StateSpace', Handler_State_StateSpace, ...
    'C', C);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %

% %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Handler_Desired_State = SRD_get_handler__desired_state(...
    'Handler_ControlInput', Handler_IK_Solution, ...
    'Handler_Time',         Handler_Time);

Handler_InverseDynamics = SRD_get_handler__InverseDynamicsConstrained_QR(...
    'Handler_ControlInput', Handler_Desired_State, ...
    'Handler_Constraints_Model', Handler_Constraints_Model, ...
    'Handler_dynamics_generalized_coordinates_model', Handler_dynamics_GC_model_evaluator, ...
    'Handler_Time', Handler_Time);

Handler_nominal_trajectory_state = SRDHandler_nominal_trajectory_state(...
    'Handler_Time', Handler_Time, ...
    'Handler_ControlInput_qva', Handler_IK_Solution, ...
    'Handler_InverseDynamics',  Handler_InverseDynamics, ...
    'function_Original_Model',  Handler_Constraints_Model.get_FirstOrderSystem_qv_handle, ...
    'ToEvaluateOriginalFunction', true);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %

dummy_controller.u = zeros(Handler_Linear_model_orthogonal.dof_control, 1);
Handler_Linear_model_implicit.Handler_Controller = dummy_controller;
Handler_Linear_model_orthogonal.Handler_Controller = dummy_controller;
% % Handler_dynamics_Linear_model_evaluator.Handler_Controller = dummy_controller;



%%
Handler_Updater = SRDHandler_Updater({...
    Handler_Desired_State, ...
    Handler_State_StateSpace, ...
    Handler_ObserverState_genCoord, ...
    Handler_MeasuredOutput, ...
    Handler_dynamics_GC_model_evaluator,...
    Handler_InverseDynamics, ...
    Handler_nominal_trajectory_state, ...
    Handler_Linear_model_orthogonal ...
    });

Handler_Updater.Update();

%%




Ar = Handler_Linear_model_orthogonal.An;
Arho = Handler_Linear_model_orthogonal.Ar;

disp("Ar")
disp(Ar)
disp("eig(Ar) =")
disp(eig(Ar))

disp("Arho")
disp(Arho)

B = Handler_Linear_model_orthogonal.B;

disp("B")
disp(B)

D = Handler_Linear_model_orthogonal.D;

disp("D")
disp(D)

save('Ar','Ar')
save('Arho','Arho')
save('B','B')
save('D','D')

%%
% x = [q, v]
k = Handler_Constraints_Model.dof_Constraint;
n = Handler_dynamics_generalized_coordinates_model.dof_configuration_space_robot;
m = Handler_dynamics_generalized_coordinates_model.dof_control;

q = Handler_State.q;    v = Handler_State.v;    u = zeros(m, 1);

J = Handler_Constraints_Model.get_Jacobian(q);
dJ = Handler_Constraints_Model.get_Jacobian_derivative(q, v);


G = [ J, zeros(size(J));
     dJ, J];
N = null(G);    


save('J','J')
save('dJ','dJ')

%%
rq0=Handler_Constraints_Model.get_Constraint(qva(:, 1));


