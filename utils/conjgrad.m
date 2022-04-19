function [x] = conjgrad(A,b,x,MaxIter)
    r=b-A*x;
    p=r;
    rsold=r'*r; 
    for i=1:MaxIter
        Ap=A*p;
        alpha=rsold/(p'*Ap);
        x=x+alpha*p;
        r=r-alpha*Ap;
        rsnew=r'*r;
        if sqrt(rsnew)<1e-10
              break;
        end
        p=r+rsnew/rsold*p;
        rsold=rsnew;
    end
end


% function [x] = conjgrad(A,b,x)
% %     r=b-A*x;
% %     p=r;
% %     rsold=r'*r;
% %  
% %     for i=1:1e6
% %         Ap=A*p;
% %         alpha=rsold/(p'*Ap);
% %         x=x+alpha*p;
% %         r=r-alpha*Ap;
% %         rsnew=r'*r;
% %         if sqrt(rsnew)<1e-10
% %               break;
% %         end
% %         p=r+rsnew/rsold*p;
% %         rsold=rsnew;
% %     
%     epsilon = 1.0e-9;
%     n = length(b);
%     r = b - feval(A,x); s=r;
%     for numIter = 1:n
%         u = feval(func,s);
%         alpha = dot(s,r)/dot(s,u);
%         x = x + alpha*s;
%         r = b - feval(A,x);
%         if sqrt(dot(r,r)) < epsilon
%             return
%         else 
%             beta = - dot(r,u)/dot(s,u);
%             s = r + beta*s;
%         end
%     end
%     error('Too many iterations');
% end