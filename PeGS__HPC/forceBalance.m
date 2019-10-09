%%forceBalance.m
%%This imposes force balance on a given set of forces and their alpha and
%%beta angles as defined by fig. 4.11 in James Puckett's thesis.
%%This is James' matlab port of peDiskSolves force balance code.
%%I just wrapped to be compatible with PEGS

%%Last Edit on 2017/07/11 by Jonathan Kollmer (jekollme@ncsu.edu)

function [alpha,f] = forceBalance(f,alpha,beta)
    clear obj
    
    z = length(alpha);
    obj.f = f(:)';
    obj.alpha = alpha(:)';
    obj.beta = beta(:)';
    obj.z = z;
    obj.x = [f;alpha];% 0.1*ones(1,2*obj.z-3);
    
    if obj.z==2
        dbeta = (obj.beta(2)-obj.beta(1))/2;
        obj.f(1)     = obj.x(1);
        obj.alpha(1) = acos(sin(dbeta));
        if (obj.alpha(1)<pi/2)          %why is this backward????
            obj.alpha(1)=acos(sin(-dbeta));
        end
        if (obj.alpha(1)>pi/2) %cant pull
            obj.alpha(1)=obj.alpha(1)-pi;
        end
        if (obj.alpha(1)<-pi/2)%cant pull
            obj.alpha(1)=obj.alpha(1)+pi;
        end
        obj.f(2)     = obj.f(1);
        obj.alpha(2) = -obj.alpha(1);

    else                
        Z = obj.z;                
        F = obj.f;
        B = obj.beta-obj.beta(Z);  %rotate disc so that beta(z)=0
        A = obj.alpha;                
        u  = 2*Z-3+1;%number of unknowns                
        cnt=1;
        for i=1:Z %distribute x to unknowns
            if cnt<u
                A(i)= obj.x(cnt);
                cnt=cnt+1;
            end
            if  cnt<u%what's y?, should be u right?? why y in c++ code?
                F(i)= obj.x(cnt); 
                cnt=cnt+1;
            end
        end
        %solve for f(z-1),f(z),alpha(z)
        za  = Z-1;
        zb  = Z;
        ii = 1:(za-1);
        fsa2 = sum(F(ii).*sin(A(ii)));
        fsab2= sum(F(ii).*sin(A(ii)+B(ii)));
        fcab2= sum(F(ii).*cos(A(ii)+B(ii)));
        F(za)= (-fsa2+fsab2)./(sin(A(za))-sin(A(za)+B(za))); %z-2 %% JO: / = ./ ??
        fsa1 	 = fsa2+ F(za)*sin(A(za));
        fsab1    = fsab2+F(za)*sin(A(za)+B(za));
        fcab1    = fcab2+F(za)*cos(A(za)+B(za));          
        F(zb) 	 = sqrt( fsab1^2 + fcab1^2 );
        %           if (f[zb]>0)
        A(zb)= asin( -fsa1/F(zb) );
        obj.alpha = A;
        obj.f     = F;                
    end
    alpha = obj.alpha(:)';
    f = obj.f(:)';
    
    
    
   %%%BEGIN SANITY CHECKS %%%%%%%  
    
   if obj.z>1
       Fx = sum(f.*cos(obj.alpha+obj.beta));
       Fy = sum(f.*sin(obj.alpha+obj.beta));
       if length(Fx)>1 || length(Fy)>1
           display('uh-oh, something bad has happened');
       end
   
       for k=1:length(f)
           if (isnan(f(k))) %for some reason result sometimes gets to be NAN, not sure why this happens
               f(k) = 0; %temproary fix is to set it zero then
               display('Warning: ForceBalance encountered a NAN value in f(k) where none was expected. Setting to 0 instead.');
           end 
%            if (f(k)<0) %for some reason, not sure why this happens
%                f(k) = 0; %temproary fix is to set it zero then
%                display('Warning: ForceBalance encountered a negative value in f(k) where none was expected. Setting to 0 instead.');
%            end 
           if (isreal(f(k))==0) %for some reason, not sure why this happens
               f(k) = 0; %temproary fix is to set it zero then
               display('Warning: ForceBalance encountered a complex value in f(k) where none was expected. Setting to 0 instead.');
           end   
           if (isnan(alpha(k))) %for some reason result sometimes gets to be NAN, not sure why this happens
               alpha(k) = 0; %temproary fix is to set it zero then
               display('Warning: ForceBalance encountered a NAN value in a(k) where none was expected. Setting to 0 instead.');
           end 
%            if (a(k)<0) %for some reason, not sure why this happens
%                a(k) = 0; %temproary fix is to set it zero then
%                display('Warning: ForceBalance encountered a negative value in a(k) where none was expected. Setting to 0 instead.');
%            end 
           if (isreal(alpha(k))==0) %for some reason, not sure why this happens
               alpha(k) = 0; %temproary fix is to set it zero then
               display('Warning: ForceBalance encountered a complex value in a(k) where none was expected. Setting to 0 instead.');
           end    
       end
   end
   %%%END SANITY CHECKS %%%%%%%    
    
end  