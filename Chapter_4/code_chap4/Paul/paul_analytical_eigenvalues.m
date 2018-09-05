function lambda=paul_analytical_eigenvalues(C,N)

    lambda=zeros(N,1);
    for k=0:(N-1)
        lambda(k+1)=(k+1)*(1-2/(C+1))^(k+1)*(2/(C+1)+1/(k+1));
    end

end