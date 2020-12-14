function [] = testGS(a,b,c,d)

 % this function test the GS method for 4 different grid sizes (a,b,c,d).
 % It plots the residual for each case after each GS iteration
 % it returns when the residual for the last system is smaller than 1E-10
 % or after the maxnb of iteration has been reached

mata = construct_laplace_matrix(a)
matb = construct_laplace_matrix(b) ;
matc = construct_laplace_matrix(c) ;
matd = construct_laplace_matrix(d) ;

solutiona = rand(a,1);
solutionb = rand(b,1);
solutionc = rand(c,1);
solutiond = rand(d,1);

%The rhs is set to be 1 everywhere and the boundary condition is
%homogeneous (first and last point)
rhsa = ones(a,1); rhsa(1,1) = 0; rhsa(a,1) =0;
rhsb = ones(b,1); rhsb(1,1) = 0; rhsb(b,1) =0;
rhsc = ones(c,1); rhsc(1,1) = 0; rhsc(c,1) =0;
rhsd = ones(d,1); rhsd(1,1) = 0; rhsd(d,1) =0;

close all
figure(1)
norm_residual =1;
ite = 0;
while norm_residual>1E-10 && ite<10
    
    subplot(2,2,1)
    hold on
    xa = linspace(0,1,a);
    solutiona = GaussSeidel(mata,rhsa,solutiona);
    resa = abs(mata*solutiona-rhsa);
    plot(xa(2:a-2),log10(resa(2:a-2)));
        
    subplot(2,2,2)
    hold on    
    xb = linspace(0,1,b);
    solutionb = GaussSeidel(matb,rhsb,solutionb);
    resb = abs(matb*solutionb-rhsb);
    plot(xb(2:b-2),log10(resb(2:b-2)));
    
    subplot(2,2,3)
    hold on
    xc = linspace(0,1,c);
    solutionc = GaussSeidel(matc,rhsc,solutionc);
    resc = abs(matc*solutionc-rhsc);
    plot(xc(2:c-2),log10(resc(2:c-2)));
    
    subplot(2,2,4)
    hold on
    solutiond = GaussSeidel(matd,rhsd,solutiond);
    xd = linspace(0,1,d);
    resd = abs(matd*solutiond-rhsd);
    plot(xd(2:d-2),log10(resd(2:d-2)));
    
    % lets monitor convergence based off the norm othe residual of this
    % system
    norm_residual = mean(resd);
    ite = ite+1;
    %pause(3)
end

ite_needed_to_converge = ite
end