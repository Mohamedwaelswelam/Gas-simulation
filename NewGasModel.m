function NewGasModel()


close
clc


%% DATA
Model_Compressor=1;


[Branch_Connectivity,Branch_Length,Branch_Diameter,Branch_Efficiency,...
    Node_Demand,Compressor_Location,Compressor_Efficiency,...
    Compressor_Costs,Compressor_Ratio,Standard_Temp,Standard_Pressure,...
    Average_Temperature,Gas_Gravity,Gas_Comp,Compressor_Temp,...
    Gas_Heat,Compressor_iComp]=GetGasNetworkData();
Slack=1;
SlackP=1000;%psia


Pressure=[%psia
    1000
    500
    600
    400
    500
    450
    300
    670
    900
    800
    700
    500
    600
    550
    ];


% Branch_Connectivity=Branch_Connectivity(1:5,:);Branch_Connectivity(5,:)=[3 5]
% Compressor_Location=[4;5];
% Compressor_Ratio=Compressor_Ratio(1:2,:);
% Node_Demand=zeros(5,1);Node_Demand(4:5)=[215;148];






%% Calculations
Ng=max(max(Branch_Connectivity));
Nbg=length(Branch_Connectivity(:,1));
Ncg=length(Compressor_Location(:,1));




Node_aux_Full=zeros(Ng,5);
%Nodes where demand from compressor would be added
Node_aux_Full(Branch_Connectivity(Compressor_Location,1),1)=1:Ncg;
%Location of additional demand to be added to those nodes
Node_aux_Full(Branch_Connectivity(Compressor_Location,1),2)=...
    Branch_Connectivity(Compressor_Location,2);
%Nodes with pressure affected by compressors
Node_aux_Full(Branch_Connectivity(Compressor_Location,2),3:4)=...
    Compressor_Ratio;
%Node to be considered;
Node_aux_Full(Branch_Connectivity(Compressor_Location(...
    Compressor_Ratio(:,2)~=0),2),5)=Branch_Connectivity(...
    Compressor_Location(Compressor_Ratio(:,2)~=0),1);


Pressure(Slack)=SlackP;
ConvFlg=Inf;
xcou=0;
while ConvFlg > 0.0001
    xcou=xcou+1;
    %Get flows at each node and thorugh each pipe
    [New_Branch_Flow,New_Node_Flow,New_CompressorIn,Pipe_Location]=GetGasFlows(Pressure,...
        Ng,Nbg,Ncg,Compressor_Location,Branch_Length,Branch_Diameter,...
        Branch_Connectivity,Branch_Efficiency,Standard_Temp,...
        Standard_Pressure,Gas_Gravity,Average_Temperature,...
        Gas_Comp,Compressor_Temp,Compressor_Efficiency,Gas_Heat,...
        Compressor_iComp,Compressor_Costs,Model_Compressor,...
        Compressor_Ratio,Node_aux_Full);


    
    %Get difference between expected and simulated nodal flows
    aux=Branch_Connectivity(Compressor_Location,:);
    aux_Demand=Node_Demand;
    aux_Demand(aux(:,1)) = aux_Demand(aux(:,1)) + aux_Demand(aux(:,2));
    aux_Demand(aux(:,2)) = 0;


    aux=Branch_Connectivity(Compressor_Location,:);
    New_Node_Flow(aux(:,1)) = New_Node_Flow(aux(:,1)) + New_Node_Flow(aux(:,2));
    New_Node_Flow(aux(:,2)) = 0;
    New_DF = New_Node_Flow+aux_Demand;
    New_DF(aux(:,1))=New_DF(aux(:,1)) + New_CompressorIn;
    
    
    Real_Nodes = 1:Ng;
    Real_Nodes = Real_Nodes(Node_aux_Full(:,4)==0);
    
    srn=length(Real_Nodes);
    aux=find(Real_Nodes==Slack);
    if aux==1
        Real_Nodes_NoS=Real_Nodes(2:srn);
    elseif Slack==srn
        Real_Nodes_NoS=Real_Nodes(1:srn-1);
    else
        Real_Nodes_NoS=Real_Nodes([1:aux-1 aux+1:srn]);
    end


    %Build Jacobian
    Jg=GetGasJacobian(Ng,Slack,Pipe_Location,Branch_Connectivity,Pressure,...
        Branch_Efficiency,Branch_Length,Branch_Diameter,Standard_Temp,...
        Standard_Pressure,Gas_Gravity,Average_Temperature,Gas_Comp,Ncg,...
        Compressor_Location,New_Branch_Flow,New_CompressorIn,...
        Compressor_iComp,Gas_Heat,Compressor_Temp,Compressor_Efficiency,...
        Compressor_Costs,Compressor_Ratio,Model_Compressor,Real_Nodes,...
        Node_aux_Full);


    
    %Calculate pressure updates
    Pressure_Correction=Jg\New_DF(Real_Nodes_NoS);
    
    %Update convergence criteria
    ConvFlg=max(abs(Pressure_Correction));
    Pressure(Real_Nodes_NoS)= Pressure(Real_Nodes_NoS) + Pressure_Correction;


    if xcou>100
        error('Focing termination');
    end
end


%     Display pressures
fprintf('Recommended pressures\n');
for x1=1:Ng
    fprintf('P:%2.0f:  %f\n',x1,Pressure(x1))
end
fprintf('Maximum error: %.4f\n',ConvFlg);
fprintf('Iterations: %.d\n',xcou);
% tstaux=input('');




%% Get gas flows
function [New_Branch_Flow,New_Node_Flow,New_CompressorIn,Pipe_Location]=GetGasFlows(...
    Pressure,Ng,Nbg,Ncg,Compressor_Location,Branch_Length,...
    Branch_Diameter,Branch_Connectivity,Branch_Efficiency,Standard_Temp,...
    Standard_Pressure,Gas_Gravity,Average_Temperature,Gas_Comp,...
    Compressor_Temp,Compressor_Efficiency,Gas_Heat,Compressor_iComp,...
    Compressor_Costs,Model_Compressor,Compressor_Ratio,Node_aux_Full)


ax1=ones(Nbg,1);ax1(Compressor_Location)=0;
ax2=1:Nbg;
Pipe_Location=ax2(ax1==1);


%Calculate flow throughout pipes
New_Branch_Flow=zeros(Nbg,1);
New_Node_Flow=zeros(Ng,1);
for x1=Pipe_Location
    if x1==1
        clc
        fprintf('Gas flows:')
    end
    L=Branch_Length(x1);%Miles
    D=Branch_Diameter(x1);%inches
    
    %Friction factor
    F=0.128/D^(1/3);
    
    
    %Get pressures
    Bc = Branch_Connectivity(x1,:);
    p=[0 0];
    for x2=1:2
        switch Node_aux_Full(Bc(x2),4)
            case 0
                p(x2) = Pressure((Bc(x2)));
            case 1
                p(x2) = Pressure(Node_aux_Full(Bc(x2),5)) * Node_aux_Full(Bc(x2),3);
            case 2
                p(x2) = Node_aux_Full(Bc(x2),3);
        end
    end
    
%     p1=Pressure(Branch_Connectivity(x1,1));%psia
%     p2=Pressure(Branch_Connectivity(x1,2));%psia
    
    Pipe_Eff=Branch_Efficiency(x1);
%     Pipe_Eff=1;
    
    %Get flow throughout the pipe
    ax= SignP(p(1),p(2)) * Pipe_Eff * 6.4774 * Standard_Temp * D^2.5 / ...
        Standard_Pressure / sqrt(F*Gas_Gravity*L*Average_Temperature*...
        Gas_Comp) * sqrt( SignP(p(1),p(2))*(p(1)^2-p(2)^2) );
%     ax= SignP(p(1),p(2)) * Pipe_Eff * 6.4774 * Average_Temperature * D^2.5 / ...
%         Standard_Pressure / sqrt(F*Gas_Gravity*L*Standard_Temp*...
%         Gas_Comp) * sqrt( SignP(p(1),p(2))*(p(1)^2-p(2)^2) );
    New_Branch_Flow(x1)=ax*24/1000000;%ft^3/h --> MMCFD
    
    fprintf('%2.0f-%2.0f)%10.4f - %10.4f: %10.4f\n',Bc, p(1),p(2),New_Branch_Flow(x1))
%     if x1 ==1
%         Pipe_Eff * 6.4774 * Standard_Temp * D^2.5 / ...
%         Standard_Pressure / sqrt(F*Gas_Gravity*L*Average_Temperature*...
%         Gas_Comp) * sqrt( SignP(p(1),p(2))*(p(1)^2-p(2)^2) )*24/1000000
%         
%         0.9 * 6.4774 * 492.3000 * 19.6^2.5 / ...
%         14.5038 / sqrt(0.0475*0.6*80.50*520*...
%         0.87) * sqrt( (1000^2-676.65^2) )*24/1000000
% 
%         p(2)
%     end


    %Estimate nodal flows
    New_Node_Flow(Branch_Connectivity(x1,1)) = ...
        New_Node_Flow(Branch_Connectivity(x1,1)) + New_Branch_Flow(x1);
    New_Node_Flow(Branch_Connectivity(x1,2)) = ...
        New_Node_Flow(Branch_Connectivity(x1,2)) - New_Branch_Flow(x1);
end
fprintf('\n');
%Calculate flows throughout compressors
New_CompressorIn=zeros(Ncg,1);
if Model_Compressor==1
    xc=0;
    for x1 = Compressor_Location'
        xc=xc+1;
        p1=Pressure(Branch_Connectivity(x1,1));%psia
        if Compressor_Ratio(xc,2)==1
            p2=p1*Compressor_Ratio(xc,1);
        else
            p2=Compressor_Ratio(xc,1);
        end


        B = 0.08530992*Compressor_Temp/Compressor_Efficiency(xc)*(Gas_Heat/(Gas_Heat-1));
        
        f=-New_Node_Flow(Branch_Connectivity(x1,1))/24*1000000;%


        HP = B * f * ( (p2/p1)^(Compressor_iComp*((Gas_Heat-1)/Gas_Heat)) - 1);
        
        Ca=Compressor_Costs(xc,1);
        Cb=Compressor_Costs(xc,2);
        Cc=Compressor_Costs(xc,3);


        New_CompressorIn(xc)=(Ca+Cb*HP+Cc*HP^2)*24/1000000;%ft^3/h --> MMCFD
        
        New_Branch_Flow(x1) = f/1000000*24 - New_CompressorIn(xc);
        
        %Estimate nodal flows
        New_Node_Flow(Branch_Connectivity(x1,1)) = ...
            New_Node_Flow(Branch_Connectivity(x1,1)) + New_Branch_Flow(x1);
        New_Node_Flow(Branch_Connectivity(x1,2)) = ...
            New_Node_Flow(Branch_Connectivity(x1,2)) - New_Branch_Flow(x1);
    end
end






%% Get gas Jacobian
function Jg=GetGasJacobian(Ng,Slack,Pipe_Location,Branch_Connectivity,...
    Pressure,Branch_Efficiency,Branch_Length,Branch_Diameter,...
    Standard_Temp,Standard_Pressure,Gas_Gravity,Average_Temperature,...
    Gas_Comp,Ncg,Compressor_Location,New_Branch_Flow,New_CompressorIn,...
    Compressor_iComp,Gas_Heat,Compressor_Temp,Compressor_Efficiency,...
    Compressor_Costs,Compressor_Ratio,Model_Compressor,Real_Nodes,...
    Node_aux_Full)


Ng=length(Real_Nodes);
x1=1;
while Slack~=Real_Nodes(x1)
    x1=x1+1;
end
if x1==1
    Xng=Real_Nodes(2:Ng);
elseif x2==Ng
    Xng=Real_Nodes(1:Ng-1);
else
    Xng=Real_Nodes([1:Slack-1 Slack+1:Ng]);
end


Jg=zeros(Ng-1,Ng-1);
Dp=0.000001;




%Auxiliar for assigning values to the Jacobian
xJaux=zeros(max(Xng),1);
xJaux(Xng)=1:Ng-1;


%Get Jacobian elements corresponding to flow through the pipes
for xp = Pipe_Location%For each pipe
    
    %Get pipe efficiency
    Pipe_Eff=Branch_Efficiency(xp);
    
    %Get pipe length
    L=Branch_Length(xp);%Miles
    
    %Get pipe internal diameter
    D=Branch_Diameter(xp);%inches
    
    %Get friction factor
    F=0.128/D^(1/3);
    
    %Identify relevant nodes and get pressures
    Bc = Branch_Connectivity(xp,:);
    
    p=[0 0];
    Maux = [1 1];
    Eaux=[1 1];
    for x1=1:2
        switch Node_aux_Full(Bc(x1),4)
            %The pressure is based on this node
            case 0
                %Get pressure
                p(x1) = Pressure((Bc(x1)));
                
                
                %Pressure is a function of a different node
            case 1
                %Add auxiliar
                Maux(x1) = Node_aux_Full(Bc(x1),3);
                
                %Change node location
                Bc(x1)=Node_aux_Full(Bc(x1),5);
                
                %Get pressure
                p(x1) = Pressure(Bc(x1)) * Maux(x1);
                
                
                %Is pressure fixed at this node?
            case 2
                %Get pressure
                p(x1) = Node_aux_Full(Bc(x1),3);
                
                %Change node location
                Bc(x1)=Node_aux_Full(Bc(x1),5);
                
                %Do not differentiate this side
                Eaux(x1)=0;
        end
    end
    
    
%     xp
%     Bc
%     error('Just stop');


    %Differentiate input and output sides
    aux = [1 -1];
    ax1 = Pipe_Eff * 6.4774 * Standard_Temp * D^2.5 / Standard_Pressure ...
        / sqrt( F * Gas_Gravity * L * Average_Temperature * Gas_Comp) ...
        *24/1000000;
    ax2 = sqrt( SignP(p(1),p(2)) * (p(1)^2 - p(2)^2) );
    for x1=1:2
%         fprintf('First ');
%         Node_aux_Full(Bc(1),4)
        %Differentiate if the pressure is not fixed and it is not the slack
        if Node_aux_Full(Bc(x1),4)~=2 && xJaux(Bc(x1))~=0 && Eaux(x1)==1
%             fprintf(' Second');
            Jval = aux(x1) * ax1/ax2 * Maux(x1) * p(x1) ;
            
            %Store results
            ax3=[-1 1];
            for x2=1:2
                if xJaux(Bc(x2))~=0
% fprintf('Df%d_%dp%d: %f --> J(%d,%d)\n',Bc,Bc(x1),Jval* ax3(x2),xJaux(Bc(x2)),xJaux(Bc(x1)))
% error('Just stop');
                    Jg(xJaux(Bc(x2)),xJaux(Bc(x1))) = ...
                        Jg(xJaux(Bc(x2)),xJaux(Bc(x1))) + Jval * ax3(x2);
                end
            end
        end
%         fprintf('\n');
    end
%     Jg
%     error('Just stop');
end






%Get Jacobian elements corresponding to compressors
if Model_Compressor==1
    for xc=1:Ncg
        Bc = [
            Branch_Connectivity(Compressor_Location(xc),1)
            Branch_Connectivity(Compressor_Location(xc),2)
            ];
        switch Compressor_Ratio(xc,2)
            %Is the pressure a function of p1
            case 1
                aux = Compressor_Ratio(xc,1);
                Jval = ( 0.08530992 * Compressor_Temp * (Gas_Heat/(Gas_Heat-1))...
                    / Compressor_Efficiency(xc) * ( aux^ ( Compressor_iComp * ...
                    (Gas_Heat-1) / Gas_Heat )-1 ) ) * Compressor_Costs(xc,2) * ...
                    -sum(Jg(1:xJaux(Bc(1))-1,xJaux(Bc(1))));
                
            case 2
                %Use numeric integration for the time being
                %Look for pipes connected to this node
                f=[0 0];Dp=[0 0.00001];
                Ps=Pressure;
                for xf=1:2
                    Ps(Bc(1))=Ps(Bc(1))+Dp(xf);
                    for xp = Pipe_Location%For each pipe
                        %Identify relevant nodes
                        Bcp = Branch_Connectivity(xp,:);
                        if Bcp(1)==Bc(1) || Bcp(2)==Bc(1)
                            p=[0 0];
                            for x1=1:2
                                switch Node_aux_Full(Bcp(x1),4)
                                    case 0
                                        p(x1) = Ps((Bcp(x1)));
                                    case 1
                                        p(x1) = Ps(Node_aux_Full(Bcp(x1),5)) * Node_aux_Full(Bcp(x1),3);
                                    case 2
                                        p(x1) = Node_aux_Full(Bcp(x1),3);
                                end
                            end
                            Pipe_Eff=Branch_Efficiency(xp);
                            L=Branch_Length(xp);%Miles
                            D=Branch_Diameter(xp);%inches
                            F=0.128/D^(1/3);
                            
                            
                            aux = SignP(p(1),p(2)) * Pipe_Eff * 6.4774 * ...
                                Standard_Temp * D^2.5 / Standard_Pressure / ...
                                sqrt( F * Gas_Gravity * L * ...
                                Average_Temperature * Gas_Comp) *24/1000000*...
                                sqrt( SignP(p(1),p(2)) * (p(1)^2 - p(2)^2) );
                            
                            if Bcp(1)==Bc(2)
                                f(xf)=f(xf)+aux;
                            else
                                f(xf)=f(xf)-aux;
                            end
                        end
                    end
                end                                


                p=[Pressure(Bc(1)) Compressor_Ratio(xc,1)];
                f(1)= 0.08530992 * Compressor_Temp * (Gas_Heat/(Gas_Heat-1))...
                    / Compressor_Efficiency(xc) * ( (p(2)/p(1))^ ( Compressor_iComp * ...
                    (Gas_Heat-1) / Gas_Heat )-1 ) * Compressor_Costs(xc,2) * f(1);
                
                p=[Pressure(Bc(1))+Dp(2) Compressor_Ratio(xc,1)];
                f(2)= 0.08530992 * Compressor_Temp * (Gas_Heat/(Gas_Heat-1))...
                    / Compressor_Efficiency(xc) * ( (p(2)/p(1))^ ( Compressor_iComp * ...
                    (Gas_Heat-1) / Gas_Heat )-1 ) * Compressor_Costs(xc,2) * f(2);
                Jval = -(f(2)-f(1))/Dp(2);   
        end
        


% fprintf('Df%d_%dp%d: %f --> J(%d,%d)\n',Bc,Bc(1),Jval* ax3(x2),xJaux(Bc(1)),xJaux(Bc(1)))        
%         fprintf('%d-%d/%d: %f \n',Bc,Bc(1),Jval)
%         error('Just stop');
        
        %Add results to the Jacobian
        Jg(xJaux(Bc(1)),xJaux(Bc(1))) = ...
            Jg(xJaux(Bc(1)),xJaux(Bc(1))) - Jval;
    end
end
%  error('Just stop');




%% Colebrook equation
function f=Colebrook(err,Re)
f=0.01;
dt=0.00000001;


faux=[Inf f];
x1=0;
while abs(faux(1)-faux(2))>dt
    x1=x1+1;
    %Assessing the function
    ax1=-2*log(err/3.7+2.51/Re/sqrt(f))-1/sqrt(f);
        
    %Numeric differentiation
    ax2=-2*log(err/3.7+2.51/Re/sqrt(f+dt))-1/sqrt(f+dt);
    ax2=(ax2-ax1)/dt;
    
    %Avoid using negative values
    f=max(f-ax1/ax2,dt);
    
    faux=[faux(2) f];
    
    if x1==100
        fprintf('A solution for the Colebrook equation could not \n');
        fprintf('be found after %d itertions\n',x1);
        fprintf('err:%.4f; Re:%.4f\n',err,Re)
        error('Stopping program');
    end
end


%% Sign function
function no=SignP(a,b)


if a>=b
    no=1;
else
    no=-1;
end


%% Network data
function [Branch_Connectivity,Branch_Length,Branch_Diameter,Branch_Efficiency,...
    Node_Demand,Compressor_Location,Compressor_Efficiency,...
    Compressor_Costs,Compressor_Ratio,Standard_Temp,Standard_Pressure,...
    Average_Temperature,Gas_Gravity,Gas_Comp,Compressor_Temp,...
    Gas_Heat,Compressor_iComp]=GetGasNetworkData()


Branch_Connectivity=[
    1 2   %1
    1 3   %2
    2 3   %3
    2 4   %4
    4 5   %5
    3 6   %6
    6 7   %7
    5 8   %8
    8 9   %9
    7 10  %10
    10 11 %11
    9 12  %12
    11 13 %13
    12 13 %14
    12 14 %15
    13 14 %16
    ];


Branch_Length=[%Miles
    80.5 %1
    88.3 %2
    55.9 %3
    61.1 %4
    0    %5
    67.9 %6
    0    %7
    93.5 %8
    0    %9
    79.7 %10
    0    %11
    73.5 %12
    87.9 %13
    86.6 %14
    79.7 %15
    83.5 %16
    ];

Branch_Diameter=[%inch
    19.6 %1
    19.6 %2
    19.6 %3
    19.6 %4
    0    %5
    19.6 %6
    0    %7
    19.6 %8
    0    %9
    16.7 %10
    0    %11
    16.7 %12
    16.7 %13
    16.7 %14
    16.7 %15
    16.7 %16
    ];

Branch_Efficiency=[
    0.9  %1
    0.9  %2
    0.9  %3
    0.9  %4
    0.9  %5
    0.9  %6
    0.9  %7
    0.85 %8
    0.85 %9
    0.9  %10
    0.9  %11
    0.85 %12
    0.85 %13
    0.9  %14
    0.9  %15
    0.85 %16
    ];
Node_Demand=[
    0   %1 
    30  %2
    90  %3
    0   %4
    0   %5
    0   %6
    0   %7
    0   %8
    0   %9
    0   %10
    0   %11
    110 %12
    40  %13
    90  %14
    ];


Compressor_Location=[%Node Efficiency Ca,Cb,Cc
    5
    7
    9
    11
    ];
Compressor_Efficiency=[
    0.83
    0.84
    0.83
    0.84
    ];
Compressor_Costs=[
    0 0.2e-3 0
    0 0.2e-3 0
    0 0.2e-3 0
    0 0.2e-3 0
    ];
Compressor_Ratio=[
    1.6 1
    1.8 1
    1000 2
    1100 2
    ];


%Standard temperature
Standard_Temp=492.3;%R
%Standard pressure
Standard_Pressure=14.5038;%psia
%Average gas temperature
Average_Temperature=520;%R
%Gas specific gravity
Gas_Gravity=0.6;%no unit
%Gas compressibility factor
Gas_Comp=0.87;%9987;%No units
%Compressor suction temperature
Compressor_Temp=520;
%Gas specific heat ratio
Gas_Heat=1.3049;
%Gas compressibility factor at compressor inlet
Compressor_iComp=0.9987;


