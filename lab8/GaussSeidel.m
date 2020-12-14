function [sols] = GaussSeidel(matrix,rhs,initial)

N = size(matrix);
N = N(1,1);

sols = initial;

%% this is the gauss seidel function

% it takes a sytems and perform 5 smoothing operations

for ktr = 1:5
    for itr=1:N    
        sols(itr) = rhs(itr);
        for  jtr=1:N       
        %% todo: implement GS smoothing    
            if itr ~= jtr
                sols(itr) = sols(itr) - matrix(itr,jtr)*sols(jtr);
            end
        end
        sols(itr) = sols(itr)/matrix(itr,itr);
    
    end
end
%vec = sols;
end