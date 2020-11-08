% clear;clc;
% %Question 2a
% A = [0,0,1/5,0,0,0,0,0,0,0,0,0;
%     1,0,1/5,0,0,0,0,0,0,0,0,0;
%     0,0,0,0,0,0,1/2,0,0,0,0,0;
%     0,0,1/5,0,1/2,0,0,0,0,0,0,0;
%     0,0,1/5,0,0,0,0,0,0,0,0,0;
%     0,0,1/5,0,0,0,0,0,0,0,0,0;
%     0,0,0,1,1/2,0,0,0,0,0,0,0;
%     0,0,0,0,0,0,1/2,0,1,0,0,0;
%     0,0,0,0,0,0,0,1,0,0,0,0;
%     0,0,0,0,0,0,0,0,0,0,1,0;
%     0,0,0,0,0,0,0,0,0,0,0,1;
%     0,0,0,0,0,0,0,0,0,1,0,0];
% 
% % Question 2b
% B = [0,1/12,1/5,0,0,1/12,0,0,0,0,0,0;
%     1,1/12,1/5,0,0,1/12,0,0,0,0,0,0;
%     0,1/12,0,0,0,1/12,1/2,0,0,0,0,0;
%     0,1/12,1/5,0,1/2,1/12,0,0,0,0,0,0;
%     0,1/12,1/5,0,0,1/12,0,0,0,0,0,0;
%     0,1/12,1/5,0,0,1/12,0,0,0,0,0,0;
%     0,1/12,0,1,1/2,1/12,0,0,0,0,0,0;
%     0,1/12,0,0,0,1/12,1/2,0,1,0,0,0;
%     0,1/12,0,0,0,1/12,0,1,0,0,0,0;
%     0,1/12,0,0,0,1/12,0,0,0,0,1,0;
%     0,1/12,0,0,0,1/12,0,0,0,0,0,1;
%     0,1/12,0,0,0,1/12,0,0,0,1,0,0];
% u0 = ones(12,1)*1/12;
% 
% %Part 2d
% different_a=[];
% different_a_iters=[];
% for alpha=.05:.1:.95
%        C_1 = create_C(B,alpha,u0);
%        [result,iter]=page_rank(C_1,u0);
%        different_a=[different_a,result];
%        different_a_iters = [different_a_iters,iter];
% end
% disp("Results for different A values")
% disp(different_a)
% disp(different_a_iters)
% 
% possible_us=[u0];
% u1=[1/6,0,1/6,0,1/6,0,1/6,0,1/6,0,1/6,0]';
% u2=[1/4,0,0,1/4,0,0,1/4,0,0,1/4,0,0]';
% u3=[1/3,0,0,0,1/3,0,0,0,1/3,0,0,0]';
% u4=[1/2,0,0,0,0,0,1/2,0,0,0,0,0,]';
% possible_us=[possible_us,u1,u2,u3,u4];
% different_u=[];
% different_u_iters=[];
% for j=1:size(possible_us,2)
%        alpha0=.8;
%        u_j=possible_us(:,j);
%        C_1 = create_C(B,alpha0,u_j);
%        [result1,iter1]=page_rank(C_1,u_j);
%        different_u=[different_u,result1];
%        different_u_iters = [different_u_iters,iter1];
% end
% disp("Results for differnt U values")
% disp(different_u)
% disp(different_u_iters)
% 
% 
% %As we can see from these results, the page rank algorithm depends on both 
% % U and Alpha. The higher alpha gets, the more concentrated the importance
% % vector becomes, and using a starting vector that is closer and closer to
% % just an [1,0,...,0] makes the importance vector more concentrated as well
% % The rate of convergence also goes down as alpha goes up and as U gets
% % more concentrated
% %Part 2e
% %Does not converge for any starting distribution
% disp("Outcome for non-stochastic matrix") 
% disp(page_rank(A,u0))
% 
% %Part 2f
% %Does not converge for any starting distribution
% disp("Outcome for a matrix with cycles")
% disp(page_rank(B,u0))

%Question 3a
incidence=load("incidencematrix.mat");
C= load_matrix(incidence.M);
u=(1/1000)*ones(1000,1);
[import,iters]=page_rank(C,u);
pages=[[1:length(import)]' , import];
top_10=sortrows(pages,2);
disp(top_10(1:11,[1,2]))

%Part 3b
%My current thought is to buy link from the cheapest 50/100/200 sites and
%see if we can get it into the top 10
bottom_fifty=top_10(end-51:end,1);
c1=incidence.M;
c1(end+1,:)=0;
c1(:,end+1)=0;
for i=1:length(bottom_fifty)
    c1(1001,bottom_fifty(i))=1;
    c1(1,1001)=1;
end
u1=(1/1001)*ones(1001,1);
[import1,iters]=page_rank(c1,u1);
pages1=[[1:length(import1)]' , import1];
top_10=sortrows(pages1,2);
disp(find(top_10(:,1)==1001))
x=1;

%Question 1
function [importance,iters] = page_rank(L,importance)
    stopping=false;
    iters=0;
    while ~stopping && iters<10000
        iters=iters+1;
        prev_importance=importance;
        importance=L*importance;
        if norm(importance-prev_importance)<=.00001
            stopping=true;
        end
    end 
    if iters>9999
        disp("Too many iterations!")
    end
    
end
function [C]=create_C(B,alpha,u)
    [rows,cols]=size(B);
    e=ones(1,cols);
    C=alpha*B + (1-alpha)*u*e;
end
%Question 3a
function [C] = load_matrix(M)     
    
    [rows,cols]=size(M);
    %Create A from M
    A=zeros(rows,cols);
    column_sums=sum(M,1);
    for i=1:rows
        for j=1:cols
            if M(i,j)==1
                A(i,j)=1/column_sums(j);
            end
        end
    end
    %Modify A to create B
    for k=cols:-1:1
        if sum(A(:,k))==0
            A(:,k)=1/rows;
        end
    end
    %Weight B to create C
    alpha=.95;
    u=(1/rows) *ones(rows,1);
    C = create_C(A,alpha,u);
end
