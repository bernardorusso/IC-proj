// Function Distortion Matrix

function [A] = dist_mat(Q,R)
    x = Q(1,2:4)-Q(1,1)
    y = Q(2,2:4)-Q(2,1)
    M = [1.,0.,0.,0.,-x(1),0.;
         1.,R ,0.,0.,-x(2),-x(2)*R;
         0.,R ,0.,0.,0.,-x(3)*R;
         0.,0.,1.,0.,-y(1),0.;
         0.,0.,1.,R ,-y(2),-y(2)*R;
         0.,0.,0.,R ,0.,-y(3)*R];
    b = [x(1);x(2);x(3);y(1);y(2);y(3)];
    M=make_permute(make_permute(make_permute(M,5,3),4,6),5,6);
    b=make_permute(make_permute(make_permute(b,5,3),4,6),5,6);
    [L U P] = LU_fact(M)
    x = linear_sist(L,U,P,b)
    A = [x(1),x(2),0;
         x(3),x(4),0;
         x(5),x(6),1]
endfunction
