function  [v21, v22, v23] = S_LS_User2_Brutal(H11, H12, H13, H21, H22, H23, H31, H32, H33, g1, g2, g3, v11, v12, v13, v31, v32, v33, n0, w1, w2, w3)
%update filters by sudo-LS algorithm 

    %Power Constraint
    P = 1;

    %brutal-force search for lambda1 
    Range = 10000;
    Precision = 0.0001;
    %-1/+1
    for n = -Range:Range

        v21h = ...
          ( g1'*H12*w1 - g1'*H12*(g1'*(H11*v11+H13*v31))'*w1...
            -g2'*H22*(g2'*(H21*v11+H23*v31))'*w2...
            -g3'*H32*(g3'*(H31*v11+H33*v31))'*w3)/(  H12'*g1*g1'*H12*w1 + H22'*g2*g2'*H22*w2 + H32'*g3*g3'*H32*w3 + (n'*Precision)*eye(2));

        v22h = ...
          ( g2'*H22*w2 - g1'*H12*(g1'*(H11*v12+H13*v32))'*w1...
            -g2'*H22*(g2'*(H21*v12+H23*v32))'*w2...
            -g3'*H32*(g3'*(H31*v12+H33*v32))'*w3)/(  H12'*g1*g1'*H12*w1 + H22'*g2*g2'*H22*w2 + H32'*g3*g3'*H32*w3 + (n'*Precision)*eye(2));
        
        v23h = ...
          ( g3'*H32*w3 - g1'*H12*(g1'*(H11*v13+H13*v33))'*w1...
            -g2'*H22*(g2'*(H21*v13+H23*v33))'*w2...
            -g3'*H32*(g3'*(H31*v13+H33*v33))'*w3)/(  H12'*g1*g1'*H12*w1 + H22'*g2*g2'*H22*w2 + H32'*g3*g3'*H32*w3 + (n'*Precision)*eye(2));

        W(n+Range+1) = abs(norm(v21h)^2+norm(v22h)^2+norm(v23h)^2-P);
    end
    
    %plot(W);
    %axis([1 2*Range 0 0.1]);
    
    [M,I] = min(W);
    lambda2 = (I-Range-1)*Precision;
    
    
        v21h = ...
          ( g1'*H12*w1 - g1'*H12*(g1'*(H11*v11+H13*v31))'*w1...
            -g2'*H22*(g2'*(H21*v11+H23*v31))'*w2...
            -g3'*H32*(g3'*(H31*v11+H33*v31))'*w3)/(  H12'*g1*g1'*H12*w1 + H22'*g2*g2'*H22*w2 + H32'*g3*g3'*H32*w3 + lambda2*eye(2));

        v22h = ...
          ( g2'*H22*w2 - g1'*H12*(g1'*(H11*v12+H13*v32))'*w1...
            -g2'*H22*(g2'*(H21*v12+H23*v32))'*w2...
            -g3'*H32*(g3'*(H31*v12+H33*v32))'*w3)/(  H12'*g1*g1'*H12*w1 + H22'*g2*g2'*H22*w2 + H32'*g3*g3'*H32*w3 + lambda2*eye(2));
        
        v23h = ...
          ( g3'*H32*w3 - g1'*H12*(g1'*(H11*v13+H13*v33))'*w1...
            -g2'*H22*(g2'*(H21*v13+H23*v33))'*w2...
            -g3'*H32*(g3'*(H31*v13+H33*v33))'*w3)/(  H12'*g1*g1'*H12*w1 + H22'*g2*g2'*H22*w2 + H32'*g3*g3'*H32*w3 + lambda2*eye(2));
        
        v21 = v21h';
        v22 = v22h';
        v23 = v23h';
    
end

    
    
