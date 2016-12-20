function a=randa(N,p)
    a=rand(1,N,p);
    a=-log(a);
    a=a./(sum(a,2)*ones(1,p));
    a=a';
end