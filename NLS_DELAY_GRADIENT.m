%GAIN OPTIMIZATION FOR AN INDUSTRIAL PID CONTROL
%PERFORMANCE INDEX J = x(N)Hx(N) + sum(xQx + uRu)
%NONLINEAR SYSTEM WITH INPUT DELAY
%29/05/23

clear all
clc
format long

%ITERATIONS FOR THE GRADIENT DESCENT METHOD
pasosMetodo = 5;
%SYSTEM STAGES
pasosSimulacion = 50000;
%INITIAL CONDITIONS FOR THE SYSTEM
ciEstado = 0.25;
ciControl = 0;
%GAIN INITIALIZATION Z-N:
kp(1) = 24.43;
Td(1) = 25;
Ti(1) = 100;
%SAMPLE TIME AND SYSTEM PARAMETERS OBTAINED BY MEANS OF THE LRS ALGORITHM
Ts = 0.2;
aGorro = 1.000146377925391;
bGorro = 3.665347874142847e-05;
cGorro = -9.359303656000177e-04;
dGorro = 7.952239373253265e-04;
%ALPHA INITIALIZATION
alfa1(1) = kp(1) + (kp(1)*Td(1))/Ts + (kp(1)*Ts)/(Ti(1));
alfa2(1) = kp(1) + (2*kp(1)*Td(1))/Ts;
alfa3(1) = (kp(1)*Td(1))/Ts;
%H,Q,R FOR THE PERFORMANCE INDEX
H=10;
Q=10;
R=220;
tiempo=0;
%EPSILON VALUE FOR THE GRADIENT STEP
epsilonGorro = 0.0000000146;
M = 50/Ts;

%PRINCIPAL CYCLE FOR THE GRADIENT DESCENT METHOD
contadorMetodo = 0;
for indexMetodo = 1:1:pasosMetodo     
  contadorMetodo = contadorMetodo + 1;
  %VECTOR OF PARAMETERS
  vectorAlfa = [alfa1(indexMetodo); alfa2(indexMetodo); alfa3(indexMetodo)];
  %ALPHA VALUES FOR ANALYSIS
  valoresAlfa(:,indexMetodo) = vectorAlfa;
  
        %APROXIMATE RESPONSE 
        %CONTROL VECTOR
        %IF t<0, u(t) = 0;
        for index = 1:1:pasosSimulacion
            
            indiceA = index - 1 - M;
            indiceB = index - 2 - M;
            indiceC = index - 3 - M;
   
            if index <= M
               
               x(index) = ciEstado;
               
               if index-1<=0
                   chiK(index,:) = [-ciEstado ciEstado -ciEstado];
                   u(index) = ciControl + chiK(index,:)*vectorAlfa;
               end
       
               if (index-1>0) && (index-1<=(M-1))
                   if index<3
                       chiK(index,:) = [-x(index) x(index-1) -ciEstado];
                       u(index) = u(index-1) + chiK(index,:)*vectorAlfa;
                   else 
                       chiK(index,:) = [-x(index) x(index-1) -x(index-2)];
                       u(index) = u(index-1) + chiK(index,:)*vectorAlfa;
                   end
               end
            end
            
            if indiceA == 0 
              
               x(index) = aGorro*x(index-1) + cGorro*x(index-1)*x(index-1) + dGorro*x(index-1)*x(index-1)*x(index-1) ...
                   + bGorro*ciControl + bGorro*[-ciEstado ciEstado -ciEstado]*vectorAlfa;
               
               chiK(index,:) = [-x(index) x(index-1) -x(index-2)];
               u(index) = u(index-1) + chiK(index,:)*vectorAlfa;  
            end
  
            if indiceB == 0
               
               x(index) = aGorro*x(index-1) + cGorro*x(index-1)*x(index-1) + dGorro*x(index-1)*x(index-1)*x(index-1)...
                   + bGorro*ciControl + bGorro*chiK(index-1-M,:)*vectorAlfa;
               
               chiK(index,:) = [-x(index) x(index-1) -x(index-2)];
               u(index) = u(index-1) + chiK(index,:)*vectorAlfa;
            end
  
            if indiceC == 0
               
               x(index) = aGorro*x(index-1) + cGorro*x(index-1)*x(index-1) + dGorro*x(index-1)*x(index-1)*x(index-1)...
                   + bGorro*u(index-2-M) + bGorro*chiK(index-1-M,:)*vectorAlfa;
               
               chiK(index,:) = [-x(index) x(index-1) -x(index-2)];
               u(index) = u(index-1) + chiK(index,:)*vectorAlfa;
            end
  
            if indiceC >= 1
               
               x(index) = aGorro*x(index-1) + cGorro*x(index-1)*x(index-1) + dGorro*x(index-1)*x(index-1)*x(index-1)...
                   + bGorro*u(index-2-M) + bGorro*chiK(index-1-M,:)*vectorAlfa;
               
               chiK(index,:) = [-x(index) x(index-1) -x(index-2)];
               u(index) = u(index-1) + chiK(index,:)*vectorAlfa;
            end
            
            
            if index-1 < 1
                tiempo(index) = 0;
            end
            if index-1 >= 1
                tiempo(index) = tiempo(index-1) + Ts;
            end
            
        end
        
        %COMPUTATION OF THE COST
        sumaCosto(1) = ciEstado*Q*ciEstado + u(1)*R*u(1);
        for indexSuma = 2:1:pasosSimulacion-1
            sumaCosto(indexSuma) = x(indexSuma)*Q*x(indexSuma) + u(indexSuma)*R*u(indexSuma) + sumaCosto(indexSuma-1);
        end
        jota(indexMetodo) = x(pasosSimulacion)*H*x(pasosSimulacion) + sumaCosto(indexSuma);

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %COMPUTATION OF THE VALUES OF THE GRADIENT
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %WE START WITH FalfaK=bGorro*chi(k-1-M)
        for indiceFalfaK=1:pasosSimulacion
            if indiceFalfaK <= (M+1)
                FalfaK(indiceFalfaK,:) = bGorro*[-ciEstado, ciEstado, -ciEstado]';
            else
                FalfaK(indiceFalfaK,:) = bGorro*chiK(indiceFalfaK-1-M,:)';
            end    
        end
        
        %WE CONTINUE WITH Lalfa(k):
        
        for indiceLalfaK=1:pasosSimulacion
            if indiceLalfaK <= (M+2)
                if indiceLalfaK == 1
                    LalfaK(indiceLalfaK,:) = 2*aGorro*bGorro*Q*x(1)*[-ciEstado, ciEstado, -ciEstado]'...
                    + 2*bGorro*cGorro*Q*x(1)*x(1)*[-ciEstado, ciEstado, -ciEstado]'...
                    + 2*bGorro*dGorro*Q*x(1)*x(1)*x(1)*[-ciEstado, ciEstado, -ciEstado]'...
                    + 2*bGorro*bGorro*Q*ciControl*[-ciEstado, ciEstado, -ciEstado]'...
                    + 2*bGorro*bGorro*Q*([-ciEstado, ciEstado, -ciEstado]')*[-ciEstado, ciEstado, -ciEstado]*vectorAlfa ...
                    + 2*ciControl*R*chiK(indiceLalfaK,:)'...
                    + 2*R*(chiK(indiceLalfaK,:)')*chiK(indiceLalfaK,:)*vectorAlfa;
                else
                    LalfaK(indiceLalfaK,:) = 2*aGorro*bGorro*Q*x(indiceLalfaK-1)*[-ciEstado, ciEstado, -ciEstado]'...
                    + 2*bGorro*cGorro*Q*x(indiceLalfaK-1)*x(indiceLalfaK-1)*[-ciEstado, ciEstado, -ciEstado]'...
                    + 2*bGorro*dGorro*Q*x(indiceLalfaK-1)*x(indiceLalfaK-1)*x(indiceLalfaK-1)*[-ciEstado, ciEstado, -ciEstado]'...
                    + 2*bGorro*bGorro*Q*u(indiceLalfaK-1)*[-ciEstado, ciEstado, -ciEstado]'...
                    + 2*bGorro*bGorro*Q*([-ciEstado, ciEstado, -ciEstado]')*[-ciEstado, ciEstado, -ciEstado]*vectorAlfa ...
                    + 2*R*u(indiceLalfaK-1)*chiK(indiceLalfaK,:)'...
                    + 2*R*(chiK(indiceLalfaK,:)')*chiK(indiceLalfaK,:)*vectorAlfa;
                end
            else
                LalfaK(indiceLalfaK,:) = 2*aGorro*bGorro*Q*x(indiceLalfaK-1)*chiK(indiceLalfaK-1-M,:)'...
                    + 2*bGorro*cGorro*Q*x(indiceLalfaK-1)*x(indiceLalfaK-1)*chiK(indiceLalfaK-1-M,:)'...
                    + 2*bGorro*dGorro*Q*x(indiceLalfaK-1)*x(indiceLalfaK-1)*x(indiceLalfaK-1)*chiK(indiceLalfaK-1-M,:)'...
                    + 2*bGorro*bGorro*Q*u(indiceLalfaK-1)*chiK(indiceLalfaK-1-M,:)'...
                    + 2*bGorro*bGorro*Q*(chiK(indiceLalfaK-1-M,:)')*chiK(indiceLalfaK-1-M,:)*vectorAlfa ...
                    + 2*R*u(indiceLalfaK-1)*chiK(indiceLalfaK,:)'...
                    + 2*R*(chiK(indiceLalfaK,:)')*chiK(indiceLalfaK,:)*vectorAlfa;
            end    
        end
        
        %WE CREATE THE VECTOR Vx(k+1) WITH k=N-1 TO k=0 
        %FOR Vx(N) = 2x(N)H;
        VxKmasUno(pasosSimulacion) = 2*x(pasosSimulacion)*H;
        for indexGradiente = pasosSimulacion-1:-1:1
          VxKmasUno(indexGradiente) = 2*x(indexGradiente)*Q + VxKmasUno(indexGradiente+1)*(aGorro...
              + 2*cGorro*x(indexGradiente) + 3*dGorro*x(indexGradiente)*x(indexGradiente) );
        end
        
        %IN A SIMILAR WAY WE COMPUTE Valpha VECTOR 
        %WE START WITH Valpha(N):
        ValfaK(pasosSimulacion,:) = 2*aGorro*bGorro*H*x(pasosSimulacion-1)*(chiK(pasosSimulacion-1-M)')...
            + 2*bGorro*cGorro*H*x(pasosSimulacion-1)*x(pasosSimulacion-1)*chiK(pasosSimulacion-1-M,:)' ...
            + 2*bGorro*dGorro*H*x(pasosSimulacion-1)*x(pasosSimulacion-1)*x(pasosSimulacion-1)*chiK(pasosSimulacion-1-M,:)' ...
            + 2*bGorro*bGorro*u(pasosSimulacion-2-M)*H*chiK(pasosSimulacion-1-M,:)'...
            + 2*bGorro*bGorro*chiK(pasosSimulacion-1-M,:)'*chiK(pasosSimulacion-1-M,:)*vectorAlfa*H;
        
        for indexGradienteAlfa = pasosSimulacion-1:-1:1%empieza en N-1
          ValfaK(indexGradienteAlfa,:) = LalfaK(indexGradienteAlfa,:) + ValfaK(indexGradienteAlfa+1,:) ...
              + VxKmasUno(indexGradienteAlfa+1)*FalfaK(indexGradienteAlfa,:);
        end
        
        vectorGradienteAlfa = ValfaK(1);
        %WE STORE THE VALUES FOR ANALYSIS
        valoresGradienteAlfa(:,indexMetodo)=vectorGradienteAlfa';        
        
        %NEW ALPHA COMPUTATION
        nuevasAlfa = vectorAlfa - epsilonGorro*vectorGradienteAlfa;
        valoresNuevasAlfa(:,indexMetodo)=nuevasAlfa;
        
        alfa1(indexMetodo+1)= nuevasAlfa(1);
        alfa2(indexMetodo+1)= nuevasAlfa(2);
        alfa3(indexMetodo+1)= nuevasAlfa(3);
         
        %NEW VALUES FOR THE PID GAINS
        kpTd = alfa3(indexMetodo+1)*Ts;
        kp(indexMetodo+1) = alfa2(indexMetodo+1) - ((2*kpTd)/Ts);
        Td(indexMetodo+1) = (alfa3(indexMetodo+1)*Ts)/kp(indexMetodo+1);
        Ti(indexMetodo+1) = (kp(indexMetodo+1)*Ts)/( alfa1(indexMetodo+1) - kp(indexMetodo+1) - (kp(indexMetodo+1)*Td(indexMetodo+1))/(Ts) );
        
        if indexMetodo ==1
           figure(indexMetodo)
           plot(tiempo,x)
   
           grid on
           hold on
        end
         
        if mod(indexMetodo,1)==0
           %figure(indexMetodo)
           plot(tiempo,x)
           
           grid on
           hold on
        end
end
figure(indexMetodo+1)
grid on
hold on
plot(1:1:length(jota), jota)
