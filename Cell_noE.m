   %  Program BD_N with Verlet neighbour lists - no elasticity
   
   %  Bidisperse spheres - Types A and B
   %  A has diameter eq 1 (diameter sigma = 1)  and velocity v0A
   %  VEL(N) AND RANCELL(N) - reflect this difference
   
   % Periodic boundary conditions in 2D,
   % Boxlength of {-L/2, L/2} in simulation units centered at origin 
   
   % (xpU, ypU, TpU) are unfolded coordinates
   % (xpF, ypF, TpF) are  folded coordinates
   
   % Time scaled with (3/D_r)
   
   % SET RANDOM NUMBER before using PARFOR
 
   % LOAD SAVED CONFIGURATION IF NEEDED

   % load ('./*.mat')

   % DEFINE PARAMETERS 
   
   %  cellrunID = 1;
 
   %   fileIDF = fopen('cellF1.dat','w');
   %   fileIDU = fopen('cellU1.dat','w');
     
   N = 40*40;                                                               
   M = 40;                                                                       
   B_T = sqrt(1.0);                                                          
   B_R = sqrt(1.0);                                                              
   v0 = 20.0;                             

   % SET MEMORY

   vel = zeros(1,N);    
   radcell = zeros(1,N);                                          

   % ELASTIC PARAMETERS

   mu = 0.4;                                                               
   alpha = 0.0;                % No ELASTICITY                                
   
   con1 = alpha*(1.+mu)*mu/pi;                  
   con2 = alpha*(1.+mu)*(6.+mu)/pi;   
   
   % LJ parameters - simulations in sigma units

   sig = 1.0;                  % Diameter of cell                                           
   epsLJ  = 1.0;               % Energy scale                                             
   sizedif = 0.04;             % Diameter diff 
   
   % Size diff is diff in width = diameter
   % Radius of A is sigma/2 - sizediff/2
   % Radius of B is sigma/2 + sizediff/2

   rsigA = 0.5*(sig - sizedif);                                          
   rsigB = 0.5*(sig + sizedif);                                           
 
   v0A = v0;                                                                 
   v0B = v0;                                                                 
   
   % Set radius of cells  
        
    Nhalf = N/2;                                                             
     
     for Ncell = 1: Nhalf;
      radcell(Ncell) = rsigA;     % radcell is RADIUS
      vel(Ncell) = v0A;           % vel is SPEED
     end

     for Ncell = Nhalf: N;
       radcell(Ncell) = rsigB;
       vel(Ncell) = v0B;
     end

   % set cut off parameters 

   EHSDAB = (rsigA + rsigB);               
   EHSDAA = (rsigA + rsigA);               
   EHSDBB = (rsigB + rsigB);               

   MAXHSD = sig + sizedif;     % This is equal to EHSDBB

   r_cut =  1.2*MAXHSD;        % Maximum cut-off LJ   
   r_cutE = 1.2*MAXHSD;        % Cut off - Elasticity                                              
   r_list = 3.7*MAXHSD;        % Cut off - List                                              
   
   % MAKE SURE r_list > max(r_cut, r_cutE)

   nselect = 2;                                                              
   MaxD = 0.02;                                                              
   dtmax = 10^-4;                                                                 
   nT = 2*10^6;                                                          

   delta = 0.2;                                                               
   Nsave = 500;                       
   Nprint = 500;
   
%   ALLOCATE MEMORY 

     xpFold = zeros(1,N);
     ypFold = zeros(1,N);
     theFold = zeros(1,N);
 
     xpU = xpFold;
     ypU = ypFold;
     theU = theFold;  

     disX = zeros(1,N);
     disY = zeros(1,N);
     disT=  zeros(1,N);

 % initialize INITIAL CONFIGURATION 

    t = 0.d0;                                                              
    savecount = 0;     

    [xpFold, ypFold, theFold, BoxL, densA, densB] = diluteB(N, delta, sig, rsigA, rsigB);

     Boxhalf = BoxL/2.0; 
     Area = BoxL*BoxL;  
     
%     50 % is CELL A and 50% is CELL B 
%     Thus the formulae:
%     AcellA = (N/2.0)*pi*(rsigA)*(rsigA);
%     AcellB = (N/2.0)*pi*(rsigB)*(rsigB);
%     densA = AcellA/Area;
%     densB = AcellB/Area;

      density = densA + densB;
      
%     xpFold = xpFnew;
%     ypFold = ypFnew;
%     theFold = theFnew;
    
    xpU = xpFold;                                
    ypU = ypFold;                               
    theU = theFold;
    
    xp1 = xpFold;
    yp1 = ypFold;
    thep1 = theFold;
    
    tp1 = 0.0;
    
    figure(1)
    plot(xpFold(1:Nhalf), ypFold(1:Nhalf),'o', 'MarkerFaceColor','y', 'MarkerSize', 10);
    axis([ -Boxhalf  +Boxhalf  -Boxhalf  +Boxhalf]); 
    
    hold on;
    plot(xpFold(Nhalf+1:N), ypFold(Nhalf+1:N),'o', 'MarkerFaceColor','y', 'MarkerSize', 10);
    hold on;
    quiver(xpFold, ypFold, cos(theFold), sin(theFold), 0); 
    hold off; 
    
    
  

 % INITIALIZE TOTAL FORCES AND TORQUES
  
      FTx = zeros(1,N);
      FTy = zeros(1,N);
      TTz = zeros(1,N);
 
 % CREATE LISTS FOR FIRST USE      
   
 [POINT, LIST, xpFsaved, ypFsaved, theFsaved, timesaved] = LISTUB(r_list, N, xpFold, ypFold, theFold, BoxL, t);
        
 % CALCULATE TOTAL FORCE USING LISTS 
     
          % Lennard Jones force
          % LJ subroutine will have the forces in terms of SCALED POSITIONS
          
          
 [LJFTx, LJFTy, LJTTz] = ForceLJB(POINT, LIST, r_cut,  N, xpFold, ypFold, theFold, BoxL, radcell);

        FTx = LJFTx;
        FTy = LJFTy;
        TTz = LJTTz;
    
        
       % GENERATE trial displacement for FIRST TIME STEP
  
 disX = (sqrt(2.0)*B_T*sqrt(dtmax)).*randn(1,N);
 disY = (sqrt(2.0)*B_T*sqrt(dtmax)).*randn(1,N);
 disT = (sqrt(2.0)*B_R*sqrt(dtmax)).*randn(1,N);
 
       % FIRST STEP - MOVE CELLS 

 xpFnew = xpFold + (dtmax.*(FTx + vel.*cos(theFold))) + disX; 
 ypFnew = ypFold + (dtmax.*(FTy + vel.*sin(theFold))) + disY; 
 theFnew = theFold + (dtmax.*TTz) + disT;
 
 xpU = xpFnew;         % Unfolded coordinates
 ypU = ypFnew;         % Unfolded coordinates
 theU = theFnew;       % Unfolded coordinates
 
 % FIRST STEP - ENFORCE PERIODICITY 

     for icell = 1: N;     

  if (xpFnew(icell)  > +Boxhalf)   
      xpFnew(icell) = xpFnew(icell) - BoxL;                
  end 
    
  if (xpFnew(icell)  < -Boxhalf)   
      xpFnew(icell) = xpFnew(icell) + BoxL;                
  end 
              
  if (ypFnew(icell)  > +Boxhalf)   
      ypFnew(icell) = ypFnew(icell) - BoxL;                
  end 
    
  if (ypFnew(icell)  < -Boxhalf)   
     ypFnew(icell) = ypFnew(icell) + BoxL;                
  end      

     end
     
 % UPDATE OLD VALUES FOR USE 

    xpFold = xpFnew;
    ypFold = ypFnew;
    theFold = theFnew;
    
    %    fprintf(fileIDF, '%d %f %f %f\n', N, B_T, B_R, mu);
    %    fprintf(fileIDF, '%f %f %f %f\n',v0,alpha,delta,BoxL);
    
    %    fprintf(fileIDU, '%d %f %f %f\n', N, B_T, B_R, mu);
    %    fprintf(fileIDU, '%f %f %f %f\n',v0,alpha,delta,BoxL);
   
 % START TIME STEPPER - use ADAPTIVE EXPLICIT EULER for time stepper
    
    for nTime = 1: nT;
  
 % UPDATE LIST if needed 
   
     [UPDATE, Dismax] = checklistB(r_cutE, r_list, xpFsaved, ypFsaved, xpFold, ypFold);
 
    if (UPDATE == 1)      
    count_update = nTime;   
    [POINT, LIST, xpFsaved, ypFsaved, theFsaved, timesaved] = LISTUB(r_list, N, xpFold, ypFold, theFold, BoxL, t);
    end
    
 % GENERATE RANDOM NUMBERS for noise 
    
     Xran = (B_T*sqrt(2.0)).*randn(1,N);
     Yran = (B_T*sqrt(2.0)).*randn(1,N);
     Tran = (B_R*sqrt(2.0)).*randn(1,N);

 % DEFINE ARRAYS

  CT = cos(theFold);
  ST = sin(theFold);

 % CALCULATE TOTAL FORCE USING LISTS BASED ON OLD positions
     
     % Lennard Jones force

     [LJFTx, LJFTy, LJTTz] = ForceLJB(POINT, LIST, r_cut, N, xpFold, ypFold, theFold, BoxL, radcell);


        FTx = LJFTx;
        FTy = LJFTy;
        TTz = LJTTz;

        
     % TENTATIVE DISPLACEMENT MADE
    
     % note that max random displacement is (1.4142*0.5)*(sqrt(dtmax)) and
     % that is less than 10^-2 = MaxD since 1.4142/2 = 0.77 approximately
     % and sqrt(dtmax) is less than 10^-2.

     %  FTDx = (dtmax.*FTx) + (dtmax.*vel.*CT);
     %  FTDy = (dtmax.*FTy) + (dtmax.*vel.*ST);

     MTx = max(abs(dtmax.*FTx));
     MTy = max(abs(dtmax.*FTy));
      
     MaxMax = max(MTx,MTy);
      
     % CORRECTED DISPLACEMENT MADE with modified dt - based on deterministic force displacement    
      
     if (MaxMax > MaxD);
          dt = dtmax*(MaxD/MaxMax);
     else
          dt = dtmax;
     end
     
     FDX = dt.*(FTx + vel.*CT); 
     FDY = dt.*(FTy + vel.*ST); 
     TDT = dt.*TTz;     
    
     disX = sqrt(dt).*Xran;
     disY = sqrt(dt).*Yran;   
     disT = sqrt(dt).*Tran;
 
     xpU = xpU + FDX + disX; 
     ypU = ypU + FDY + disY; 
     theU = theU + TDT + disT; 
   
     xpFnew = xpFold + FDX + disX; 
     ypFnew = ypFold + FDY + disY;
     theFnew = theFold + TDT + disT; 

      % ENFORCE PERIODICITY     
      
      for icell = 1: N;     

  if (xpFnew(icell)  > +Boxhalf)   
      xpFnew(icell) = xpFnew(icell) - BoxL;                
  end 
    
  if (xpFnew(icell)  < -Boxhalf)   
      xpFnew(icell) = xpFnew(icell) + BoxL;                
  end 
              
  if (ypFnew(icell)  > +Boxhalf)   
      ypFnew(icell) = ypFnew(icell) - BoxL;                
  end 
    
  if (ypFnew(icell)  < -Boxhalf)   
     ypFnew(icell) = ypFnew(icell) + BoxL;                
  end      

      end
      
  % UPDATE CURRENT TIME

    t = t + dt;

  % save data, fig, mat and also binned velocity data to get averages
  % to do the binning move back to actual bins using BoxL  
                                         
        if (mod(nTime,Nsave) == 0);
   
   savecount = savecount + 1;    
   [xp2, yp2, thep2, tp2] =  savedataB(t, xpU, ypU, theU, xpFnew, ypFnew, theFnew, v0, B_R, B_T, alpha, mu, N, savecount, BoxL, radcell, vel, xp1, yp1, thep1, tp1);

  % reset values for next save and print
 
        xp1 = xp2;
        yp1 = yp2;
        tp1 = tp2;

       end
         
 % UPDATE VARIABLES FOR NEXT RUN 
    
        xpFold = xpFnew;
        ypFold = ypFnew;
        theFold = theFnew;     
    
    end

    % fclose(fileIDF);
    % fclose(fileIDU);
    
 % end of the number of time steps specified 


