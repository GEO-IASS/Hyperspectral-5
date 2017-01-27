function [index_set, abundances, rec, err]=AAM(x,L,numit)
% AAM Alternating angle minimization implementation of the MESMA algorithm
%
% This function executes the AAM algorithm on a set of libraries. For every
% pixel in x, an endmember set is constructed containing zero or one
% endmembers from each library. The indices and abundances are returned. 
% The algorithm employed is iterative alternating angle minimization.
%
% Input:  X: the N data points to unmix in d dimensions, (d,N)
%         L: spectral library, cell array of p elements of size (d,~)
%         K: Optional number of iterations. Default 3
% Output: index_set identifying the endmembers from each library, (p,N)
%         abundances with respect to these endmembers, (p,N)
%         reconstruct contains the reconstructed spectral (d,N)
%         err contains the reconstruction error (Euclidean distance) (1,N)
%
%
% Rob Heylen, 2016, University of Antwerp.



% Turn this on to activate sanity checks on the libraries and input values.
% Turn this off if you are sure there are no doubles in the libraries, and
% you do not want to check if pixels are contained in the libraries.
if nargin==2
    numit=3;
end

[d,num]=size(x);
p=numel(L);
for i=1:p
    N(i)=size(L{i},2);
end
mask=zeros(2^p-1,p);
err=zeros(2^p-1,num);
for setcnt=1:2^p-1
    setmask=logical(de2bi(setcnt,p));
    mask(setcnt,:)=setmask;
    q=sum(setmask);
    idx=find(setmask);
    if q==1
        t=zeros(N(idx),num);
        for i=1:N(idx)
            t(i,:)=sqrt(sum((x-L{idx}(:,i)*ones(1,num)).^2));
        end
        [err(setcnt,:),I{setcnt}]=min(t);
        A{setcnt}=ones(1,num);
    else
        clear Lt;
        for i=1:q
            Lt{i}=L{idx(i)};
        end
        [I{setcnt}, A{setcnt}, ~, err(setcnt,:)] = AAM_no_subset( x, Lt, numit); 
    end
end
[err,J]=min(err);
abundances=zeros(p,num);
index_set=zeros(p,num);
rec=zeros(size(x));

for i=1:num
    q=sum(mask(J(i),:));
    idx=find(mask(J(i),:));
    index_set(idx,i)=I{J(i)}(:,i);
    E=zeros(d,q);
    for j=1:q
        E(:,j)=L{idx(j)}(:,I{J(i)}(j,i));
    end
    abundances(idx,i)=A{J(i)}(:,i);
    rec(:,i)=E*A{J(i)}(:,i);
end
end
