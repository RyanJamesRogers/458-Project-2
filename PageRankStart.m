clear;
clc;


%Defining Incidence Matrix
M = [1 0 0 1 ...
     0 1 1 1 ...
     0 1 1 1 ...
     ];
 
[r,c]= size(M);

col_sums = sum(M,1);


create_A(M,r,c,col_sums);
 
function A = create_A(M,r,c,col_sums)
    for i= 1:r
        for j= 1:c
            if M(i,j)==1
                A(i,j) = 1/col_sums(j);
            else
                A(i,j) = 0;
            end
        end
    end  
end


    

    
        
       

%Suppose you have n webpages that link to each other. Let Nj be
%the number of different webpages to which webpage j links.




