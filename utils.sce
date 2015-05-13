// Funções auxiliares
function [L,U,P]=LU_fact(A)
    [n n]=size(A)
    P=eye(n,n)
    for i=1:(n-1)
        j=i
        M=A(i,i)
        for l=0:n-i-1
            if abs(A((n-l),i))>M then
                j=n-l
                M=abs(A((n-l),i))
            end
        A=make_permute(A,i,j)
        P=make_permute(P,i,j)
        end
        for k=(i+1):n
            A(k,(i+1):n) = A(k,(i+1):n)-(A(i,(i+1):n)*(A(k,i)/A(i,i)))
            A(k,i) = A(k,i)/A(i,i)
        end
    end
    L=tril(A,-1)+eye(n,n)
    U=triu(A)
endfunction

function [A]=make_permute(A,i,j)
    A([i j],:)=A([j i],:)
endfunction

function x=linear_sist(L,U,P,b)
    n=length(b)
    y=zeros(n,1)
    x=zeros(n,1)
    b=P*b
    for i=1:n
        y(i,:)=b(i,:)-L(i,1:i-1)*y(1:i-1,:)
    end
    for j=1:n
        i=n-(j-1)
        x(i,:)=(y(i,:)-U(i,i+1:n)*x(i+1:n,:))/U(i,i)
    end
endfunction
