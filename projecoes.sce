// O problema

// A = [a  b  0
//      c  d  0   -> considerando que não há translação
//      m  n  1]

// Q = [0 x1 x2 x3  => [0  x1/z1  x2/z2  x3/z3
//      0 y1 y2 y3      0  y1/z1  y2/z2  y3/z3]
//      1 z1 z2 z3]

// P = [0  0  B  B  => [0  0  B  B
//      0 -H -H  0      0 -H -H  0]
//      1  1  1  1]

// Desenvolvendo AP = Q, temos:
//    -bH/z1      =   x1/z1
//    -dH/z1      =   y1/z1  
//   (aB-bH)/z2   =   x2/z2  
//   (cB-dH)/z2   =   y2/z2  
//     aB/z3      =   x3/z3  
//     cB/z3      =   y3/z3  
//     1-nH       =    z1
//    1+mB-nH     =    z2    ->  (1-nH)+mB = z1+mB=z2
//     1+mB       =    z3    ->  1+mB = 1+(z2-z1) = z3 -> z1-z2+z3 = 1//
//
// Onde nós conhecemos os termos xi/zi e yi/zi para todo i (que são os vértices
// do quadrilátero da foto)
//
// Denotaremos xi/zi e yi/zi por Xi e Yi respenctivamente 
//
// Podemos transformar o problema no sistema linear
//
// [-B   H   0   0   0  X2   0    [a     [0
//   0   H   0   0  X1   0   0     b      0
//   0   0  -B   0   0   0  Y3     c      0
//   0   0   0   H  Y1   0   0  *  d   =  0
//   0   0   0   0   1  -1   1     z1     1
//   0   0  -B   H   0  Y2   0     z2     0
//  -B   0   0   0   0   0  X3]    z3]    0]
//
// Nós não podemos resolver o sistema, uma vez que não conhecemos os valores
// de B e H, mas, se pedirmos ao usuário uma estimativa da razão R de B/H,
// podemos, sem perda de generalidade, substituir B por R e H por 1 em todas
// as equações, resultando no sistema:
//
// [-R   1   0   0   0  X2   0    [a     [0
//   0   1   0   0  X1   0   0     b      0
//   0   0  -R   0   0   0  Y3     c      0
//   0   0   0   1  Y1   0   0  *  d   =  0
//   0   0   0   0   1  -1   1     z1     1
//   0   0  -R   1   0  Y2   0     z2     0
//  -R   0   0   0   0   0  X3]    z3]    0]
//
// Uma vez resolvido o sistema, teremos os valores necessários para resolver
// 
// [ 0  1  * [m  = [1-z1
//  -R  1     n]    1-z2
//  -R  0]          1-z3]
//
// Então, temos todos os valores da matriz A e a última linha de Q
// Como A é invertível, A*P=Q -> P=(A^-1)*Q


function [P,A] = ajusta(Q,R)
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
    P = inv(A)*Q
endfunction

// Resolvendo os Exemplos
[P1 A1] = ajusta(Q1, 0.5);
[P2 A2] = ajusta(Q2, 2.0);
[P3 A3] = ajusta(Q3, 1.0);
