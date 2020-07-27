function   [UPDATE, Dismax] = checklistB(r_cutE, r_list, xpFsaved, ypFsaved, xpFold, ypFold)

% Call to check if update is needed M30, 2012

Dismax = 0.0;

% start the calculation of actual maximum displacement 

Dis_X = abs(xpFold - xpFsaved);   % x displacement compared to saved 
Dis_Y = abs(ypFold - ypFsaved);    % y displacement compared to saved
 
 maxDis_X = max(Dis_X);
 maxDis_Y = max(Dis_Y);
 
 if (maxDis_X > Dismax) 
     Dismax = maxDis_X;
 end

 if(maxDis_Y > Dismax)
     Dismax = maxDis_Y;
 end
 
 DisCheck = 2.0*sqrt(2.0*Dismax*Dismax);  

%  conservative skin crossing 
%  Based ON ELASTIC CUT-OFF criterion as its the longer one and yields 
%  smaller criteria for crossing i.e, the list is updates as many times as
% the elastic interaction neighbours need to be revised
 
 if (DisCheck > (r_list - r_cutE))
     UPDATE = 1;
 else
     UPDATE = 0;
 end
 
end

