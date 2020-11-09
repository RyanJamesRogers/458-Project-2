clear;
clc;
% M is not a link matrix but instead needs to be turned into one
R = load('incidencematrix.mat');
M = struct2mat(R);
cost = 0;

x = (1/1000)*ones(1000,1);
rnkDes = GetC(.3,x,M);

while getRank(1001,rnkDes)>0.1*length(M)
    M = addpage(M);
    M = addconnection(1001,length(M),M);
    rnkDes = GetC(.3,x,M);
    cost = cost+1000;
    
end
cost








function val = getRank(site, C)
    x = (1/1000)*ones(1000,1);
    order = 1:1:1000;
    x_ordering = [order', x];
    final = sortrows(x_ordering,2,'descend');
    val 
end




function m_new = addPage(M)
    [r,c]=size(M);
    m_new=M;
    m_new(end+1,:)=0;
    m_new(:,end+1)=0;
end

function dubnew = addconnection(site,from,m_new)
    m_new(site,from)=1;
    m_new(from,site)=1;
    dubnew=m_new;
end

function [C] = GetC(alpha,u,A)

% Use equation given to calculate C
C = alpha*A + (1-alpha)*u*ones(1,1000);
end
    
function [M,n]= struct2mat(R)
%STRUCT2MAT	Converts a structure into a matrix.
%	[X,n]= STRUCT2MAT(S) converts a structre S into a numeric matrix X. The 
%	contents of each *numeric* field of S (either a vector or a matrix) will 
%	form 1 column of X. Fieldnames are returned in cell array 'n'. If the 
%	fields of S aren't of the same length, the columns of X will be padded 
%	with NaN.
%	Example:
%	s= struct('a',['string of letters'],'b',[1 2; 3 4],'c',[1 2 3 4 5 6 7 8 9])
%	[x,n]= struct2mat(s)
%	Author: F. de Castro
% Extract field names
fn= fieldnames(R);
% Identify numeric & find out maximum length
len=   zeros(1,numel(fn));
isnum= zeros(1,numel(fn));
maxlen= 0;
for j= 1:numel(fn)
	if isnumeric(R.(fn{j})) 
		isnum(j)= j;
		len(j)= numel(R.(fn{j}));
		if len(j) > maxlen, maxlen= numel(R.(fn{j})); end
	end
end
isnum= isnum(isnum~=0);
ncol= numel(isnum);
% Preallocation
M= NaN(maxlen, ncol);
n= cell(ncol,1);
% Take numeric fields and field names
for j= 1:ncol
	M(1:len(isnum(j)),j)= R.(fn{isnum(j)})(:);
	n(j)= {fn{isnum(j)}};
end
M=reshape(M,1000,1000)';
end
