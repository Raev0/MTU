function num=Loc(MAP,ELE)
[L,W]=size(MAP);
temp=0;
if length(ELE)==2 && W>2;
    for n=3:W
        ELE=[ELE,' '];
    end
end


for n=1:L
        obj=sum(MAP(n,1:W)==ELE);
        if obj==W
            temp=n;
        end
end
num=temp;
end