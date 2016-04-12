function DC_value=FindDC(A)
    [m,n]=size(A);
    Ini_value=zeros(1,m-1);
    if m<2 || n<5
    else
        for i=2:m
            for j=5:n
               if isempty(A{i,j})
               else
                   Ini_value(i-1)=A{i,j};
               end
            end
        end
    end
    DC_value=Ini_value;
end
