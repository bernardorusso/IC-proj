// Function Distortion Matrix

function [A,Q] = dist_mat(Q,R)
    M = [-R, 1, 0, 0, 0, Q(1,3), 0;
          0, 1, 0, 0, Q(1,2), 0, 0;
          0, 0,-R, 0, 0, 0, Q(2,4);
          0, 0, 0, 1, Q(2,2), 0, 0;
          0, 0, 0, 0, 1,-1, 1;
          0, 0,-R, 1, 0, Q(2,3), 0;
         -R, 0, 0, 0, 0, 0, Q(1,4)]
    
    b = [0; 0; 0; 0; 1; 0; 0]
    [L U P] = LU_fact(M)
    x = linear_sist(L,U,P,b)
    A = [x(1) x(2) 0;
         x(3) x(4) 0]
    A(3,:) = [(x(7)-1)/R, x(7)-x(6), 1]
    Q(3,:) = [1, 1, 1, 1]
    for i=2:4
        Q(:,i) = x(i+3)*Q(:,i)
    end
endfunction
