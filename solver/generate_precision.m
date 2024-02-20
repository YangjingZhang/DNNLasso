%% Generate the ground truth precision matrix
%% Input: dimension d and figure type
%% type = 1, Gamma = BB^T + D, where
%%           Prob(B_ij=-1)=(1-rho)/2,Prob(B_ij=1)=(1-rho)/2,Prob(B_ij=0)=rho
%%           the number of nonzero elements of B is 10*d
%%           D is diagonal with entries uniformly in [0,0.1] + 0.0001
%%           see (Yoon and Kim, 2021) random graph  
%% type = 2, see (Yoon and Kim, 2021) graph with clusters
function Gamma = generate_precision(d,type)
B = zeros(d,d);
switch type
    case 1
        no_of_nz = min(10*d,d^2); %number of nonzero entries
        ind_of_nz = randperm(d^2,no_of_nz);
        B(ind_of_nz) = [ones(no_of_nz/2,1);-ones(no_of_nz/2,1)];
        Gamma = B*B' + diag(rand(1,d)/10) + 0.0001*eye(d);
    case 2
        if d < 100
            no_of_blocks = 1;
        elseif d <= 200
            no_of_blocks = 5;
        else
            no_of_blocks = 10;
        end
        size_of_block = round(d/no_of_blocks);
        no_of_nz = round(d/2)*2;
        for i = 1:no_of_blocks           
            if i == no_of_blocks
                size_of_block = d - size_of_block*(no_of_blocks - 1);
            end
            ind_of_nz = randperm(size_of_block^2,no_of_nz);
            B0 = zeros(size_of_block,size_of_block);
            B0(ind_of_nz) = [ones(no_of_nz/2,1);-ones(no_of_nz/2,1)];
            if i == no_of_blocks
                B((size_of_block*(i - 1) + 1):end,((size_of_block*(i - 1) + 1):end)) = ...
                    B0*B0' + diag(rand(1,size_of_block)/10) + 0.0001*eye(size_of_block);
            else
                B((size_of_block*(i - 1) + 1):size_of_block*i,(size_of_block*(i - 1) + 1):size_of_block*i) = ...
                    B0*B0' + diag(rand(1,size_of_block)/10) + 0.0001*eye(size_of_block);
            end
        end
        Gamma = B;
end
end
