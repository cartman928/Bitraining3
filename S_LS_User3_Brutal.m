function  [v31, v32, v33, lambda3, W3] = S_LS_User2_Brutal(H11, H12, H13, H21, H22, H23, H31, H32, H33, g1, g2, g3, v11, v12, v13, v21, v22, v23, n0, w1, w2, w3, Range, Precision)
%update filters by sudo-LS algorithm 

    %Power Constraint
    P = 1;

    %brutal-force search for lambda1 

    for n = -Range:Range

        v31h = ...
          ( g1'*H13*w1 - g1'*H13*(g1'*(H11*v11+H12*v21))'*w1...
            -g2'*H23*(g2'*(H21*v11+H22*v21))'*w2...
            -g3'*H33*(g3'*(H31*v11+H32*v21))'*w3)/(  H13'*g1*g1'*H13*w1 + H23'*g2*g2'*H23*w2 + H33'*g3*g3'*H33*w3 + (n'*Precision)*eye(2));

        v32h = ...
          ( g2'*H23*w2 - g1'*H13*(g1'*(H11*v12+H12*v22))'*w1...
            -g2'*H23*(g2'*(H21*v12+H22*v22))'*w2...
            -g3'*H33*(g3'*(H31*v12+H32*v22))'*w3)/(  H13'*g1*g1'*H13*w1 + H23'*g2*g2'*H23*w2 + H33'*g3*g3'*H33*w3 + (n'*Precision)*eye(2));
        
        v33h = ...
          ( g3'*H33*w3 - g1'*H13*(g1'*(H11*v13+H12*v23))'*w1...
            -g2'*H23*(g2'*(H21*v13+H22*v23))'*w2...
            -g3'*H33*(g3'*(H31*v13+H32*v23))'*w3)/(  H13'*g1*g1'*H13*w1 + H23'*g2*g2'*H23*w2 + H33'*g3*g3'*H33*w3 + (n'*Precision)*eye(2));

        W(n+1+Range) = abs(norm(v31h)^2+norm(v32h)^2+norm(v33h)^2-P);
        W3(n+1+Range) = norm(v31h)^2+norm(v32h)^2+norm(v33h)^2-P;
    end
    
    %plot(W);
    %axis([1 2*Range 0 0.001]);
    
    
    [M,I] = min(W);
    lambda3 = (I-Range-1)*Precision;
    
        v31h = ...
          ( g1'*H13*w1 - g1'*H13*(g1'*(H11*v11+H12*v21))'*w1...
            -g2'*H23*(g2'*(H21*v11+H22*v21))'*w2...
            -g3'*H33*(g3'*(H31*v11+H32*v21))'*w3)/(  H13'*g1*g1'*H13*w1 + H23'*g2*g2'*H23*w2 + H33'*g3*g3'*H33*w3 + lambda3*eye(2));

        v32h = ...
          ( g2'*H23*w2 - g1'*H13*(g1'*(H11*v12+H12*v22))'*w1...
            -g2'*H23*(g2'*(H21*v12+H22*v22))'*w2...
            -g3'*H33*(g3'*(H31*v12+H32*v22))'*w3)/(  H13'*g1*g1'*H13*w1 + H23'*g2*g2'*H23*w2 + H33'*g3*g3'*H33*w3 + lambda3*eye(2));
        
        v33h = ...
          ( g3'*H33*w3 - g1'*H13*(g1'*(H11*v13+H12*v23))'*w1...
            -g2'*H23*(g2'*(H21*v13+H22*v23))'*w2...
            -g3'*H33*(g3'*(H31*v13+H32*v23))'*w3)/(  H13'*g1*g1'*H13*w1 + H23'*g2*g2'*H23*w2 + H33'*g3*g3'*H33*w3 + lambda3*eye(2));
        
        v31 = v31h';
        v32 = v32h';
        v33 = v33h';
        
        %norm(v31)^2+norm(v32)^2+norm(v33)^2
        
        a=1;
    
end

    
    
