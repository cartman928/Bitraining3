%% Initialize Parameters

clc
clear

alpha = 0;    %coefficient for block fading model
beta = 0.8^2;  %attenuation loss from non-direct antennas
w1 = 1;     %weight of MSE
w2 = 1;
w3 = 1;
n0 = 10^(-2);    %noise variance

iternums = 1:20; % number of iterations
N_Realizations = 10;

C1 = zeros(N_Realizations, length(iternums));
C2 = zeros(N_Realizations, length(iternums));
C3 = zeros(N_Realizations, length(iternums));

%% Start Loop
for Realization = 1 : N_Realizations
    Realization
        
    %% Random Channels
    H11 = (randn(2,2)+1i*randn(2,2))/sqrt(2);
    H22 = (randn(2,2)+1i*randn(2,2))/sqrt(2);
    H33 = (randn(2,2)+1i*randn(2,2))/sqrt(2);
    
    H12 = (randn(2,2)+1i*randn(2,2))/sqrt(2/beta); 
    H13 = (randn(2,2)+1i*randn(2,2))/sqrt(2/beta); 
    H21 = (randn(2,2)+1i*randn(2,2))/sqrt(2/beta);
    H23 = (randn(2,2)+1i*randn(2,2))/sqrt(2/beta); 
    H31 = (randn(2,2)+1i*randn(2,2))/sqrt(2/beta); 
    H32 = (randn(2,2)+1i*randn(2,2))/sqrt(2/beta);
    
    
    %% one iteration per block
    g1 = rand(2, 1) + 1i*rand(2, 1);    
    g2 = rand(2, 1) + 1i*rand(2, 1);
    g3 = rand(2, 1) + 1i*rand(2, 1);
    g1/norm(g1);
    g2/norm(g2);
    g3/norm(g3);
 
    v11 = zeros(2, 1); 
    v12 = zeros(2, 1);
    v13 = zeros(2, 1); 
    v21 = zeros(2, 1); 
    v22 = zeros(2, 1);
    v23 = zeros(2, 1); 
    v31 = zeros(2, 1); 
    v32 = zeros(2, 1); 
    v33 = zeros(2, 1); 
    
    for numiters = 1:length(iternums)
        
        %% bi-directional training
            %%Backward Training: sudo-LS Algorithm
            %abs(error)<2*10^(-4)
            for k1 = 1 : 500
            Loops = [Realization numiters k1] 
            v11_o = v11;
            v12_o = v12;
            v13_o = v13;
            v21_o = v21;
            v22_o = v22;
            v23_o = v23;
            v31_o = v31;
            v32_o = v32;
            v33_o = v33;
            
            %{
            z1 = rand;
            if z1 <= (1/3)
                [v11, v12, v13, lambda1] = S_LS_User1_Brutal(H11, H12, H13, H21, H22, H23, H31, H32, H33, g1, g2, g3, v21, v22, v23, v31, v32, v33, n0, w1, w2, w3);
                z11 = rand;
                if z11 >= 0.5
                [v21, v22, v23, lambda2] = S_LS_User2_Brutal(H11, H12, H13, H21, H22, H23, H31, H32, H33, g1, g2, g3, v11, v12, v13, v31, v32, v33, n0, w1, w2, w3);
                [v31, v32, v33, lambda3] = S_LS_User3_Brutal(H11, H12, H13, H21, H22, H23, H31, H32, H33, g1, g2, g3, v11, v12, v13, v21, v22, v23, n0, w1, w2, w3);
                else
                [v31, v32, v33, lambda3] = S_LS_User3_Brutal(H11, H12, H13, H21, H22, H23, H31, H32, H33, g1, g2, g3, v11, v12, v13, v21, v22, v23, n0, w1, w2, w3);
                [v21, v22, v23, lambda2] = S_LS_User2_Brutal(H11, H12, H13, H21, H22, H23, H31, H32, H33, g1, g2, g3, v11, v12, v13, v31, v32, v33, n0, w1, w2, w3);
                end
            elseif rand >= (2/3)
                [v21, v22, v23, lambda2] = S_LS_User2_Brutal(H11, H12, H13, H21, H22, H23, H31, H32, H33, g1, g2, g3, v11, v12, v13, v31, v32, v33, n0, w1, w2, w3);
                z12 = rand;
                if z12 >= 0.5
                [v11, v12, v13, lambda1] = S_LS_User1_Brutal(H11, H12, H13, H21, H22, H23, H31, H32, H33, g1, g2, g3, v21, v22, v23, v31, v32, v33, n0, w1, w2, w3);
                [v31, v32, v33, lambda3] = S_LS_User3_Brutal(H11, H12, H13, H21, H22, H23, H31, H32, H33, g1, g2, g3, v11, v12, v13, v21, v22, v23, n0, w1, w2, w3);
                else
                [v31, v32, v33, lambda3] = S_LS_User3_Brutal(H11, H12, H13, H21, H22, H23, H31, H32, H33, g1, g2, g3, v11, v12, v13, v21, v22, v23, n0, w1, w2, w3);
                [v11, v12, v13, lambda1] = S_LS_User1_Brutal(H11, H12, H13, H21, H22, H23, H31, H32, H33, g1, g2, g3, v21, v22, v23, v31, v32, v33, n0, w1, w2, w3);
                end
            else
                [v31, v32, v33, lambda3] = S_LS_User3_Brutal(H11, H12, H13, H21, H22, H23, H31, H32, H33, g1, g2, g3, v11, v12, v13, v21, v22, v23, n0, w1, w2, w3);
                z13 = rand;
                if z13 >= 0.5
                [v11, v12, v13, lambda1] = S_LS_User1_Brutal(H11, H12, H13, H21, H22, H23, H31, H32, H33, g1, g2, g3, v21, v22, v23, v31, v32, v33, n0, w1, w2, w3);
                [v21, v22, v23, lambda2] = S_LS_User2_Brutal(H11, H12, H13, H21, H22, H23, H31, H32, H33, g1, g2, g3, v11, v12, v13, v31, v32, v33, n0, w1, w2, w3);
                else
                [v21, v22, v23, lambda2] = S_LS_User2_Brutal(H11, H12, H13, H21, H22, H23, H31, H32, H33, g1, g2, g3, v11, v12, v13, v31, v32, v33, n0, w1, w2, w3);
                [v11, v12, v13, lambda1] = S_LS_User1_Brutal(H11, H12, H13, H21, H22, H23, H31, H32, H33, g1, g2, g3, v21, v22, v23, v31, v32, v33, n0, w1, w2, w3);
                end
            end
            %}

            [v11, v12, v13, lambda1] = S_LS_User1_Brutal(H11, H12, H13, H21, H22, H23, H31, H32, H33, g1, g2, g3, v21, v22, v23, v31, v32, v33, n0, w1, w2, w3);
            [v21, v22, v23, lambda2] = S_LS_User2_Brutal(H11, H12, H13, H21, H22, H23, H31, H32, H33, g1, g2, g3, v11, v12, v13, v31, v32, v33, n0, w1, w2, w3);
            [v31, v32, v33, lambda3] = S_LS_User3_Brutal(H11, H12, H13, H21, H22, H23, H31, H32, H33, g1, g2, g3, v11, v12, v13, v21, v22, v23, n0, w1, w2, w3);
            
          
            %[v11, v12, v13]
            %[v21, v22, v23]
            %[v31, v32, v33]
            Power = [norm(v11)^2+norm(v12)^2+norm(v13)^2 norm(v21)^2+norm(v22)^2+norm(v23)^2 norm(v31)^2+norm(v32)^2+norm(v33)^2]
            Error = [norm(v11-v11_o) norm(v12-v12_o) norm(v13-v13_o)]
            %Lambda = [lambda1 lambda2 lambda3]
            
            %Convergence Detector
            if( (norm(v11-v11_o)+norm(v12-v12_o)+norm(v13-v13_o)+norm(v21-v21_o)+norm(v22-v22_o)+norm(v23-v23_o)+norm(v31-v31_o)+norm(v32-v32_o)+norm(v33-v33_o)) < 10^(-4)  )
            disp('converges!');
            break;%break the for loop if it's true the condition
            end
            
            %{
            subplot(3,1,1);
            drawnow
            plot(W1);
            drawnow
            axis([1 2*Range 0 10])
            drawnow
            
            subplot(3,1,2);
            drawnow
            plot(W2);
            drawnow
            axis([1 2*Range 0 10])
            drawnow
            
            
            subplot(3,1,3);
            drawnow
            plot(W3);
            drawnow
            axis([1 2*Range 0 10])
            drawnow
            %}
            end
            


            %%Forward Training: LS Algorithm
            [g1, g2, g3] = LS(H11, H12, H13, H21, H22, H23, H31, H32, H33, v11, v12, v13, v21, v22, v23, v31, v32, v33, n0);
            %norm(g1)^2
            %norm(g2)^2

            SINR1 = norm(g1'*(H11*v11+H12*v21+H13*v31))^2/(norm(g1'*(H11*v12+H12*v22+H13*v32))^2+norm(g1'*(H11*v13+H12*v23+H13*v33))^2+n0*g1'*g1);
            SINR2 = norm(g2'*(H21*v12+H22*v22+H23*v32))^2/(norm(g2'*(H21*v11+H22*v21+H23*v31))^2+norm(g2'*(H21*v13+H22*v23+H23*v33))^2+n0*g2'*g2);
            SINR3 = norm(g3'*(H31*v13+H32*v23+H33*v33))^2/(norm(g3'*(H31*v12+H32*v22+H33*v32))^2+norm(g3'*(H31*v12+H32*v22+H33*v32))^2+n0*g3'*g3);
            C1(Realization, numiters) = abs(log2(1+SINR1));
            C2(Realization, numiters) = abs(log2(1+SINR2));
            C3(Realization, numiters) = abs(log2(1+SINR3));
            
    end
            
    
end



%% Plot C(bits/channel)
%figure
%hold on

p1=plot(iternums, mean(C1)+mean(C2)+mean(C3),'--');

axis([1 numiters 0 25])

xlabel('Number of iterations')
ylabel('C(bits/channel)')
title('Simple Receivers')
%}

