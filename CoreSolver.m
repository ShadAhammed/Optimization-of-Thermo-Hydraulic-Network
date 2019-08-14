%%Initialization of the System

SystemData = LoadSystemData();

%Initialize Common Variables

MeasureMode = 1;        %Measurement device mode is "real"
SystemMode = 1;         %System mode is not "non-leakage"
LossMode = 1;           %Loss Mode is lossless

PrintLevel = 1;
CharacteristicsStatus = -1;
ParameterStatus = -1;
MaterialPropertyStatus = 0;

NodeScaleVal = 1;
PipeScaleVal = 1;
MeshScaleVal = 1;
MeasScaleVal = 1;


K = SystemData.SystemCharacteristics.K; %Total number of nodes
L = SystemData.SystemCharacteristics.L; %Total number of pipes

DevNum = size(SystemData.MeasureCharacteristics.M_in,1)+size(SystemData.MeasureCharacteristics.M_out_m,1)+size(SystemData.MeasureCharacteristics.M_out_e,1)+size(SystemData.MeasureCharacteristics.M_P,1);

%Calculating number of variables 'n' and number of constraints 'm'

if LossMode == 1
    n = 2*K+L;
    m = size(SystemData.MeasureCharacteristics.M_in,1)+size(SystemData.MeasureCharacteristics.M_out_m,1) + size(SystemData.MeasureCharacteristics.M_out_e,1) + size(SystemData.MeasureCharacteristics.M_P,1) + 1;
else n = 3*K+L;
    m = size(SystemData.MeasureCharacteristics.M_in,1)+size(SystemData.MeasureCharacteristics.M_out_m,1) + size(SystemData.MeasureCharacteristics.M_out_e,1) + size(SystemData.MeasureCharacteristics.M_P,1) + 1;
end

Int_eye_L = eye(L,L);

% Calculating valve device matrix

SystemData.MeasureCharacteristics.V_m = diag(SystemData.MeasureCharacteristics.V);

%Updating the coupling matrices based on valve states

SystemData.SystemCharacteristics.tildeM = SystemData.SystemCharacteristics.M*(Int_eye_L - SystemData.MeasureCharacteristics.V_m);
SystemData.SystemCharacteristics.tildeM_T = (Int_eye_L - SystemData.MeasureCharacteristics.V_m)*SystemData.SystemCharacteristics.M_T;
SystemData.SystemCharacteristics.tildeabs_M_T = (Int_eye_L - SystemData.MeasureCharacteristics.V_m)*SystemData.SystemCharacteristics.abs_M_T;

%Initializing system state variable based on pressure and mass flows

SystemData.SystemCharacteristics.x = InitStateVector(n,K,L,SystemData.SystemCharacteristics.p,SystemData.SystemCharacteristics.m_i,SystemData.SystemCharacteristics.m_a,SystemData.SystemCharacteristics.T,LossMode);

%Initializing measured state variable based on pressure,mass flows & pumps

SystemData.MeasureCharacteristics.x = InitMeasVector(n,K,SystemData.MeasureCharacteristics.p,L,SystemData.MeasureCharacteristics.m_i,SystemData.SystemCharacteristics.K_FlowReturn,SystemData.MeasureCharacteristics.in,SystemData.MeasureCharacteristics.out_m,SystemData.MeasureCharacteristics.out_e,SystemData.MeasureCharacteristics.P_dev,SystemData.MeasureCharacteristics.P,SystemData.MeasureCharacteristics.T,LossMode);

% Initializing optimization variable X by passing system state variables

X = InitOptVar(SystemData.SystemCharacteristics.x,n);

%Initializing system constraints based on number of coupled devices

SystemData.MeasureCharacteristics.S = InitCon(SystemData.MeasureCharacteristics.M_in,SystemData.MeasureCharacteristics.M_out_m,SystemData.MeasureCharacteristics.M_out_e,SystemData.MeasureCharacteristics.M_P,K,L,n);

%Calculating box constraints for Optimisation varibles

[XL,XU] = InitBoxCon(K,L,n,SystemData.MeasureCharacteristics.p,SystemData.MeasureCharacteristics.p_dev,SystemData.SystemCharacteristics.K_Flow,SystemData.SystemCharacteristics.x,SystemData.MeasureCharacteristics.dp,SystemData.MeasureCharacteristics.m_i_dev,SystemData.MeasureCharacteristics.V,SystemData.MeasureCharacteristics.m_i,SystemData.MeasureCharacteristics.dm_i,SystemData.SystemCharacteristics.K_FlowReturn,SystemData.MeasureCharacteristics.in_dev,SystemData.MeasureCharacteristics.in,SystemData.MeasureCharacteristics.out_m_dev,SystemData.MeasureCharacteristics.out_e_dev,SystemData.MeasureCharacteristics.P_dev,SystemData.MeasureCharacteristics.out_m,SystemData.MeasureCharacteristics.out_e,SystemData.MeasureCharacteristics.P_m,SystemData.MeasureCharacteristics.P,SystemData.MeasureCharacteristics.din,SystemData.MeasureCharacteristics.dout_m,SystemData.MeasureCharacteristics.dout_e,SystemData.MeasureCharacteristics.dP,SystemData.MeasureCharacteristics.T,SystemData.MeasureCharacteristics.T_dev,SystemData.MeasureCharacteristics.dT,MeasureMode,SystemMode,LossMode);


%%Update the system properties

% Updating density of the nodes

SystemData = UpdateDensity(LossMode,K,SystemData,SystemData.SystemCharacteristics.x);

%Updating Dynamic Viscosity of the nodes

SystemData = UpdateDynViscosity( LossMode,K,SystemData,SystemData.SystemCharacteristics.x );

%Updating Heat Capacity of the nodes
if LossMode == 2
    SystemData = UpdateHeatCapacity( LossMode, SystemData,K );
end
%Updating Geo pressure of the network

SystemData  = UpdateGeoPressure( SystemData,K,L );

%Updating Geometric Resistance of the pipes

SystemData = UpdateGeometricResistance( SystemData,L );

%Updating Pipe Friction

SystemData = UpdatePipeFriction( SystemData,K,L);

%Updating Pipe Resistance

SystemData  = UpdatePipeResistance( SystemData,K,L );

%%Initialize System and Measure Characteristics

%Initializing Measure Chararcteristics

SystemData = InitMeasChar( LossMode,SystemData,K,L,n );

%Initializing System Characteristics

SystemData  = InitSysChar( SystemData,K,L,n);

%Evaluting the initial Objective Function

Obj = EvalObj(X,SystemData.SystemCharacteristics.S_T,SystemData.SystemCharacteristics.x_g,SystemData.SystemCharacteristics.S,SystemData.MeasureCharacteristics.x,SystemData.MeasureCharacteristics.X2);
Initial_Obj = Obj;

%%Evaltualing the linear eqaulity constraints

%Creating the final System Constraint Matrix
SystemData.MeasureCharacteristics.S = [SystemData.SystemCharacteristics.S;SystemData.MeasureCharacteristics.S];

G = EvalCon(X,SystemData.MeasureCharacteristics.S,m);

%%Optimization with fmincon

options = optimoptions('fmincon','Display','iter','Algorithm','interior-point');
[X,fval,ef,output,lambda] = fmincon(@(X) EvalObj(X,SystemData.SystemCharacteristics.S_T,SystemData.SystemCharacteristics.x_g,SystemData.SystemCharacteristics.S,SystemData.MeasureCharacteristics.x,SystemData.MeasureCharacteristics.X2),X,[],[],SystemData.MeasureCharacteristics.S,zeros(size(SystemData.MeasureCharacteristics.S,1),1),XL,XU,[],options);

% Updating state vector

SystemData.SystemCharacteristics.x = X;

%Updating state variable

SystemData.SystemCharacteristics.p = X(1:K);
SystemData.SystemCharacteristics.m_i = X(K+1:K+L);
SystemData.SystemCharacteristics.m_a = X(K+L+1:n);

%Updating Geodetic pressure drop & pressure drop

SystemData  = UpdateGeoPressure( SystemData,K,L );
SystemData.SystemCharacteristics.Dp_g = SystemData.SystemCharacteristics.x_g(K+1:K+L);
SystemData.SystemCharacteristics.Dp = SystemData.SystemCharacteristics.R*SystemData.SystemCharacteristics.m_i;

%%Calculation of pipe friction coefficient

%Calculating Reynold number to decide the flow pattern

Reynold_Number = CalcRe(SystemData.SystemCharacteristics.m_i,SystemData.SystemCharacteristics.d,SystemData.SystemCharacteristics.mu_L);
for i = 1:1:L
    if Reynold_Number(i) > 4000
        fprintf('The flow is turbulent in Pipe %d\n',i);
    end
end

%Calculating Lambda based using Darcy friction factor formulae

SystemData.SystemCharacteristics.lambda = CalcLambda(SystemData.SystemCharacteristics.m_i,SystemData.SystemCharacteristics.d,SystemData.SystemCharacteristics.k,SystemData.SystemCharacteristics.mu_L);
for i = 1:1:L
    disp(['Pipe friction coefficient in Pipe ' num2str(i) '::' num2str(SystemData.SystemCharacteristics.lambda(i))]);
end

%%Calculating sum of Outer Mass Flow

%Calulating outer mass flow from optimised state variables

m_aSum = Calcm_aSum(K,L,SystemData.SystemCharacteristics.x);

%Calculating sum of outer mass flows for flow network

inSum= CalcinSum(SystemData.SystemCharacteristics.K_Flow,SystemData.MeasureCharacteristics.in_dev,X,K,L);

%Updating Objective value

SystemData.SystemCharacteristics.e = fval;
Int_e = (transpose(SystemData.SystemCharacteristics.m_i)*SystemData.SystemCharacteristics.tildeM_T + transpose(SystemData.SystemCharacteristics.m_a))*(SystemData.SystemCharacteristics.tildeM*SystemData.SystemCharacteristics.m_i + SystemData.SystemCharacteristics.m_a);
SystemData.SystemCharacteristics.e_K = Int_e(1,1);
Int_e = (transpose(SystemData.SystemCharacteristics.p)*SystemData.SystemCharacteristics.tildeM + transpose(SystemData.SystemCharacteristics.Dp_g) + transpose(SystemData.SystemCharacteristics.m_i)*SystemData.SystemCharacteristics.tildeR)*(SystemData.SystemCharacteristics.tildeM_T*SystemData.SystemCharacteristics.p + SystemData.SystemCharacteristics.Dp_g + SystemData.SystemCharacteristics.tildeR*SystemData.SystemCharacteristics.m_i);
SystemData.SystemCharacteristics.e_L = Int_e(1,1);
Int_e = (transpose(SystemData.SystemCharacteristics.x) - transpose(SystemData.MeasureCharacteristics.x))*SystemData.MeasureCharacteristics.X2*(SystemData.SystemCharacteristics.x - SystemData.MeasureCharacteristics.x);
SystemData.SystemCharacteristics.e_m = Int_e(1,1);

% Generate new Temporary Data

TempData = SystoTemp(SystemData);

%%Displaying optimised results

disp('Internal State Vector after optimization:  ');
disp(SystemData.SystemCharacteristics.x);
disp(['Sum of Outer mass flow: ' num2str(m_aSum)]);
disp(['Initial Obj: ' num2str(Initial_Obj)]);
disp (['Final Objective:  ' num2str(fval)]);
disp(['Total iterations:  ' num2str(output.iterations)]);
disp(['Optimization completed with exit flag:  ' num2str(ef)]);
disp('Checking Leakage Information.....');
LeakageMass;
for i = 1:1:K
    if EAProbability(i)>0
        fprintf('The Leakage is in Exclusion Area %d\n',i);
    end
end
for i = 1:1:K
    if m_a_ZeroNodes(i)>0
        fprintf('Leakage is in node %d\n', i);
    end
end
fprintf('Leakage mass: %d\n', m_a_sum);