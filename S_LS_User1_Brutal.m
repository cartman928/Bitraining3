function  [v11, v12, v13, lambda1] = S_LS_User1_Brutal(H11, H12, H13, H21, H22, H23, H31, H32, H33, g1, g2, g3, v21, v22, v23, v31, v32, v33, n0, w1, w2, w3)
%update filters by sudo-LS algorithm 

    %Power Constraint
    P = 1;

    %% brutal-force search for lambda1 
    %First-Round Search:-10~+10(0.1)
    Precision1 = 0.1;
    for n = -100:100

        v11h = ...
          ( g1'*H11*w1 - g1'*H11*(g1'*(H12*v21+H13*v31))'*w1...
            -g2'*H21*(g2'*(H22*v21+H23*v31))'*w2...
            -g3'*H31*(g3'*(H32*v21+H33*v31))'*w3)/(  H11'*g1*g1'*H11*w1 + H21'*g2*g2'*H21*w2 + H31'*g3*g3'*H31*w3 + (n'*Precision1)*eye(2));

        v12h = ...
          ( g2'*H21*w2 - g1'*H11*(g1'*(H12*v22+H13*v32))'*w1...
            -g2'*H21*(g2'*(H22*v22+H23*v32))'*w2...
            -g3'*H31*(g3'*(H32*v22+H33*v32))'*w3)/(  H11'*g1*g1'*H11*w1 + H21'*g2*g2'*H21*w2 + H31'*g3*g3'*H31*w3 + (n'*Precision1)*eye(2));
        
        v13h = ...
          ( g3'*H31*w3 - g1'*H11*(g1'*(H12*v23+H13*v33))'*w1...
            -g2'*H21*(g2'*(H22*v23+H23*v33))'*w2...
            -g3'*H31*(g3'*(H32*v23+H33*v33))'*w3)/(  H11'*g1*g1'*H11*w1 + H21'*g2*g2'*H21*w2 + H31'*g3*g3'*H31*w3 + (n'*Precision1)*eye(2));

        if (P-(norm(v11h)^2+norm(v12h)^2+norm(v13h)^2)) < 0
        W1(n+101) = 10;
        else
        W1(n+101) = P-(norm(v11h)^2+norm(v12h)^2+norm(v13h)^2);
        end
    end
   [M1,I1] = min(W1);
    lambda1 = (I1-101)*Precision1;
    
    %Second-Round Search:-0.1~+0.1(0.001)
    Precision2 = 0.001;
    for n = -100:100

        v11h = ...
          ( g1'*H11*w1 - g1'*H11*(g1'*(H12*v21+H13*v31))'*w1...
            -g2'*H21*(g2'*(H22*v21+H23*v31))'*w2...
            -g3'*H31*(g3'*(H32*v21+H33*v31))'*w3)/(  H11'*g1*g1'*H11*w1 + H21'*g2*g2'*H21*w2 + H31'*g3*g3'*H31*w3 + (n'*Precision2+lambda1)*eye(2));

        v12h = ...
          ( g2'*H21*w2 - g1'*H11*(g1'*(H12*v22+H13*v32))'*w1...
            -g2'*H21*(g2'*(H22*v22+H23*v32))'*w2...
            -g3'*H31*(g3'*(H32*v22+H33*v32))'*w3)/(  H11'*g1*g1'*H11*w1 + H21'*g2*g2'*H21*w2 + H31'*g3*g3'*H31*w3 + (n'*Precision2+lambda1)*eye(2));
        
        v13h = ...
          ( g3'*H31*w3 - g1'*H11*(g1'*(H12*v23+H13*v33))'*w1...
            -g2'*H21*(g2'*(H22*v23+H23*v33))'*w2...
            -g3'*H31*(g3'*(H32*v23+H33*v33))'*w3)/(  H11'*g1*g1'*H11*w1 + H21'*g2*g2'*H21*w2 + H31'*g3*g3'*H31*w3 + (n'*Precision2+lambda1)*eye(2));
        
        W2(n+101) = abs(P-(norm(v11h)^2+norm(v12h)^2+norm(v13h)^2));
    end
    [M2,I2] = min(W2);
    lambda1 = (I2-101)*Precision2+lambda1;

%Second-Round Search:-0.001~+0.001(0.00001)
    Precision3 = 0.00001;
    for n = -100:100

        v11h = ...
          ( g1'*H11*w1 - g1'*H11*(g1'*(H12*v21+H13*v31))'*w1...
            -g2'*H21*(g2'*(H22*v21+H23*v31))'*w2...
            -g3'*H31*(g3'*(H32*v21+H33*v31))'*w3)/(  H11'*g1*g1'*H11*w1 + H21'*g2*g2'*H21*w2 + H31'*g3*g3'*H31*w3 + (n'*Precision3+lambda1)*eye(2));

        v12h = ...
          ( g2'*H21*w2 - g1'*H11*(g1'*(H12*v22+H13*v32))'*w1...
            -g2'*H21*(g2'*(H22*v22+H23*v32))'*w2...
            -g3'*H31*(g3'*(H32*v22+H33*v32))'*w3)/(  H11'*g1*g1'*H11*w1 + H21'*g2*g2'*H21*w2 + H31'*g3*g3'*H31*w3 + (n'*Precision3+lambda1)*eye(2));
        
        v13h = ...
          ( g3'*H31*w3 - g1'*H11*(g1'*(H12*v23+H13*v33))'*w1...
            -g2'*H21*(g2'*(H22*v23+H23*v33))'*w2...
            -g3'*H31*(g3'*(H32*v23+H33*v33))'*w3)/(  H11'*g1*g1'*H11*w1 + H21'*g2*g2'*H21*w2 + H31'*g3*g3'*H31*w3 + (n'*Precision3+lambda1)*eye(2));
        
        W3(n+101) = abs(P-(norm(v11h)^2+norm(v12h)^2+norm(v13h)^2));
    end
    [M3,I3] = min(W3);
    lambda1 = (I3-101)*Precision3+lambda1;
    
    
    
%%  
        v11h = ...
          ( g1'*H11*w1 - g1'*H11*(g1'*(H12*v21+H13*v31))'*w1...
            -g2'*H21*(g2'*(H22*v21+H23*v31))'*w2...
            -g3'*H31*(g3'*(H32*v21+H33*v31))'*w3)/(  H11'*g1*g1'*H11*w1 + H21'*g2*g2'*H21*w2 + H31'*g3*g3'*H31*w3 + lambda1*eye(2));

        v12h = ...
          ( g2'*H21*w2 - g1'*H11*(g1'*(H12*v22+H13*v32))'*w1...
            -g2'*H21*(g2'*(H22*v22+H23*v32))'*w2...
            -g3'*H31*(g3'*(H32*v22+H33*v32))'*w3)/(  H11'*g1*g1'*H11*w1 + H21'*g2*g2'*H21*w2 + H31'*g3*g3'*H31*w3 + lambda1*eye(2));
        
        v13h = ...
          ( g3'*H31*w3 - g1'*H11*(g1'*(H12*v23+H13*v33))'*w1...
            -g2'*H21*(g2'*(H22*v23+H23*v33))'*w2...
            -g3'*H31*(g3'*(H32*v23+H33*v33))'*w3)/(  H11'*g1*g1'*H11*w1 + H21'*g2*g2'*H21*w2 + H31'*g3*g3'*H31*w3 + lambda1*eye(2));
        
        v11 = v11h';
        v12 = v12h';
        v13 = v13h';    
end

    
    
