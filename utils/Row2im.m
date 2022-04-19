function k = Row2im(A,imSize, winSize,flag)

if nargin < 4
    flag = 1;
end

Nx = imSize(1);
Ny = imSize(2);
Nc = imSize(3:end);
Bx = winSize(1);
By = winSize(2);
A = reshape(A,[size(A,1),Bx*By,Nc]);

count = 1;
k = zeros(imSize);
w = k;
    for y = 1 : By 
        for x = 1 : Bx
            k(x:Nx-Bx+x,y:Ny-By+y,:) = k(x:Nx-Bx+x,y:Ny-By+y,:) +...
                reshape(A(:,count,:),[Nx-Bx+1,Ny-By+1,prod(Nc)]);
            w(x:Nx-Bx+x,y:Ny-By+y,:) = w(x:Nx-Bx+x,y:Ny-By+y,:) + 1; %% w is the number of repeating
            count = count + 1;
        end
    end
if flag == 1
    k = k./w;
end

end
