classdef caco3_main
    % trying a simple diagenesis
    % irregular grid
    % separate into subroutines
    %====================================
    % cpp options (see also defines.h)
    % sense     :  not doing any signal change experiments, used for lysocline and CaCO3 burial calculations
    % biotest   :  examining different styles of biotubation
    % track2    :  tracking signals with multipe CaCO3 species at different time steps
    % size      :  two types of CaCO3 species with different sizes
    % nonrec    :  not storing the profile files but only CaCO3 conc. and burial flux at the end of simulation
    % nondisp   :  not displaying the results
    % showiter  :  showing each iteration on display
    % sparse    :  use sparse matrix solver for caco3 and co2 system
    % recgrid   :  recording the grid to be used for making transition matrix in LABS
    % allnobio  :  assuming no bioturbation for all caco3 species
    % allturbo2 :  assuming homogeneous bio-mixing for all caco3 species
    % alllabs   :  assuming non-local mixing from LABS for all caco3 species
    % ===================================
    properties (Constant)
        
        nspcc = 4;  % default: 4; if def_sense: 12;  if def_track2: 42; if def_size: 8    % number of CaCO3 species
        ztot=500d0; % if def_sense: 50.0d0;   % cm , total sediment thickness
        nz = 100;       % grid number
        tol = 1d-8; % 1d-6    % tolerance of error
        nrec = 15;      % total recording time of sediment profiles
        
        %% Constants
        cai = 10.3d-3;      % mol kg-1 calcium conc. seawater
        fact = 1d-3;  % w/r factor to facilitate calculation
        % :: rhosed = 2.09d0 % g/cm3 sediment particle density assuming opal
        rhosed = 2.6d0;  % g/cm3 sediment particle density assming kaolinite
        rhoom = 1.2d0;  % g/cm3 organic particle density
        rhocc = 2.71d0;  % g/cm3 organic particle density
        mom = 30d0;  % g/mol OM assuming CH2O
        % :: msed = 87.11d0 % g/mol arbitrary sediment g/mol assuming opal (SiO2•n(H2O) )
        msed = 258.16d0;  % g/mol arbitrary sediment g/mol assuming kaolinite ( 	Al2Si2O5(OH)4 )
        mcc = 100d0;  % g/mol CaCO3
        ox2om = 1.3d0;  % o2/om ratio for om decomposition (Emerson 1985; Archer 1991)
        % om2cc = 0.666d0;   % rain ratio of organic matter to calcite
        om2cc = 0.7d0;   % for signal tracking exp; rain ratio of organic matter to calcite
        % :: om2cc = 0.5d0  % rain ratio of organic matter to calcite
        ncc = 4.5d0;   % (Archer et al. 1989) reaction order for caco3 dissolution
        
        o2th = 0d0;  % threshold oxygen level below which not to calculate
        dev = 1d-6;  % deviation addumed
        zml_ref = 12d0;  % a referece mixed layer depth
        ccx_th = 1d-300;  % threshold caco3 conc. (mol cm-3) below which calculation is not conducted
        omx_th = 1d-300;  % threshold om    conc. (mol cm-3) below which calculation is not conducted
        
        
        %% OM degradation
        % komi = 2d0;       % /yr  % a reference om degradation rate const.; MUDS
        komi = 0.5d0;       % /yr  % arbitrary
        % komi = 0.1d0;     % /yr  % Canfield 1994
        % komi = 0.06d0;    % /yr  % ?? Emerson 1985? who adopted relatively slow decomposition rate
        
        %% options, compare defines.h
        def_test = false;       % using 'test' directory
        def_biotest = false;   	% testing 5kyr signal change event
        def_size = false;       % testting two size caco3 simulation
        
        def_sense = false;       % without signal tracking
        def_track2 = false;   	% using method2 to track signals (default 42 species)
        def_nondisp = true;       % won't showing results on display
        def_nonrec = false;       % define not recording profiles?
        def_showiter = false;      % show co2 iterations & error in calccaco3sys
        def_sparse = false;       % using sparse matrix solve (you need UMFPACK)
        
        def_allnobio = false;       % with/without bioturbation
        def_allturbo2 = false;      % all turbo2 mixing
        def_alllabs = false;        % all labs mixing
        def_allnonlocal = false; 	% ON if assuming non-local mixing (i.e., if labs or turbo2 is ON)
        
        def_oxonly = false;      % enabling only oxic degradation of om
        
        def_nodissolve = false;     % assuming no caco3 dissolution
        
        def_fullclump = false;      % allowing full 2 substitution including 17o
        
        def_mocsy = false;          % using mocsy for caco3 thermodynamics
        MOCSY_USE_PRECISION = 2;    % digit used for mocsy
        
        def_recgrid = false;    % recording the grid to be used for making transition matrix in LABS
        
        %        kcci = 10.0d0*365.25d0      % /yr; a reference caco3 dissolution rate const.
        kcci = 1d0*365.25d0         % /yr; cf., 0.15 to 30 d-1 Emerson and Archer (1990) 0.1 to 10 d-1 in Archer 1991
        %        kcci = 0d0*365.25d0        % /yr; if def_nodissolve 
    end
    
    properties
        
        
        %%
        z; dz;   % depth, individual sediment layer thickness
        
        poroi = 0.8d0   % a reference porosity
        poro;           % porosity
        sporo; sporoi; sporof; ; porof  % solid volume fraction (1 - poro), i and f denote the top and bottom values of variables
        
        zox;
    end
    
    methods(Static)
        
        function obj = caco3_main(arg1)
            % set default values
            %            fprintf('in caco3_main \n');
            %            obj.nz = arg1;
        end
        
        function [dz,z] = makegrid(beta,nz,ztot, def_recgrid)
            
            for iz=1:nz
                z(iz) = iz*ztot/nz;  % regular grid
                if (iz==1)
                    dz(iz) = ztot*log((beta+(z(iz)/ztot)^2d0)/(beta-(z(iz)/ztot)^2d0))/log((beta+1d0)/(beta-1d0));
                end
                if (iz~=1)
                    dz(iz) = ztot*log((beta+(z(iz)/ztot)^2d0)/(beta-(z(iz)/ztot)^2d0))/log((beta+1d0)/(beta-1d0)) - sum(dz(1:iz-1));
                end
            end
            
            for iz=1:nz  % depth is defined at the middle of individual layers
                if (iz==1)
                    z(iz)=dz(iz)*0.5d0;
                end
                if (iz~=1)
                    z(iz) = z(iz-1)+dz(iz-1)*0.5d0 + 0.5d0*dz(iz);
                end
            end
            
            %% ~~~~~~~~~~~~~ saving grid for LABS ~~~~~~~~~~~~~~~~~~~~~~
            if(def_recgrid)
                %            open(unit=100, file='C:/cygwin64/home/YK/LABS/1dgrid.txt',action='write',status='unknown')
                str = sprintf('.LABS/1dgrid.txt');
                file_tmp = fopen(str,'wt');
                
                for iz = 1:nz
                    %                 write(100,*) dz(iz)
                    fprintf(file_tmp,'%17.16e \n',dz(iz));
                end
                fclose(file_tmp);
            end
            
        end
        
        function [poro, porof, sporof, sporo, sporoi] = getporosity(z, poroi, nz)
            %             nz=100;
            %             poroi = 0.8d0;   % a reference porosity
            poro = poroi;  % constant porosity
            % ----------- Archer's parameterization
            calgg = 0.0d0;  % caco3 in g/g (here 0 is assumed )
            pore_max =  1d0 - ( 0.483d0 + 0.45d0 * calgg) / 2.5d0;  % porosity at the bottom
            exp_pore = 0.25d0*calgg + 3.d0 *(1d0-calgg);  % scale depth of e-fold decrease of porosity
            poro = exp(-z/exp_pore) * (1.0d0-pore_max) + pore_max;
            porof = pore_max;  % porosity at the depth
            porof = poro(nz);  % this assumes zero-porosity gradient at the depth; these choices do not affect the calculation
            sporof = 1d0-porof;  %  volume fraction of solids at bottom depth
            % ------------------
            % cf., poro = poro_0*exp(-z(iz)/poro_scale)  % Hydrate modeling parameterization where poro_0 = 0.69 ... poro_scale = 2000 (m)
            sporoi = 1d0-poroi;  % volume fraction of solids at the seawater-sediment interface (SWI)
            sporo = 1d0 - poro;  %  volume fraction of solids
            
            %             for iz=1:global_variable.nz
            %                 fprintf('%17.16e %17.16e %17.16e\n',z(iz),poro(iz),sporo(iz));
            %             end
        end
        
        function [keq1    ,keq2   ,keqcc  ,co3sat, dif_dic    ,dif_alk    ,dif_o2 ,kom, kcc, global_var] = coefs(tmp,sal,dep, global_var)
            
            ff = global_var.poro.*global_var.poro;       % representing tortuosity factor
            
            dif_dic0 = (151.69d0 + 7.93d0*tmp); % cm2/yr at 2 oC (Huelse et al. 2018)
            dif_alk0 = (151.69d0 + 7.93d0*tmp); % cm2/yr  (Huelse et al. 2018)
            dif_o20 =  (348.62d0 + 14.09d0*tmp);
            
            dif_dic = dif_dic0*ff;  % reflecting tortuosity factor
            dif_alk = dif_alk0*ff;
            dif_o2 = dif_o20*ff;
            
            kom = global_var.komi * ones(1, global_var.nz);  % assume reference values for all reaction terms
            kcc = global_var.kcci * ones(global_var.nz, global_var.nspcc);
            
            if(global_var.def_size)
                % assume stronger dissolution for fine species (1-4)
                kcc(:,1) = kcci*10d0;
                kcc(:,2) = kcci*10d0;
                kcc(:,3) = kcci*10d0;
                kcc(:,4) = kcci*10d0;
            end
            
            keq1 = caco3_therm.calceq1(tmp,sal,dep);    % carbonic acid dissociation const. function called from caco3_therm.f90
            keq2 = caco3_therm.calceq2(tmp,sal,dep);    % bicarbonate dissociation const. function called from caco3_therm.f90
            
            keqcc = caco3_therm.calceqcc(tmp,sal,dep);  % calcite solubility function called from caco3_therm.f90
            co3sat = keqcc/global_var.cai;         % co3 conc. at calcite saturation
            
            
        end %coefs
        
        function [omflx, detflx, ccflx] = flxstat(om2cc, ccflxi, mcc, nspcc)
            %% determine steady state flux
            
            omflx = om2cc*ccflxi;  % om rain = rain ratio x caco3 rain
            detflx = (1d0/9d0)*ccflxi*mcc; % 90% of mass flux becomes inorganic C; g cm-2 yr-1
            ccflx = ones(nspcc,1)*ccflxi/nspcc;  %  rains of individual caco3 species (vector of nspcc species) is equivalently distributed as default
            
        end
        
        function [w, wi] = burial_pre(detflx,ccflx, msed,mvsed,mvcc, poroi, nz)
            % initial guess of burial profile, requiring porosity profile
            % w = burial rate, wi = burial rate initial guess
            
            % burial rate w from rain fluxes represented by volumes
            % initial guess assuming a box representation (this guess is accurate when there is no caco3 dissolution occurring)
            % om is not considered as it gets totally depleted
            
            wi = (detflx./msed.*mvsed + sum(ccflx).*mvcc)./(1d0-poroi);
            w = wi * ones(1,nz);
            
        end
        
        function age = dep2age(dz,w,nz)
            % depth -age conversion
            
            dage = dz./w;  % time spans of individual sediment layers
            age = 0.0d0;
            for iz=1:nz  % assigning ages to depth in the same way to assign depths to individual grids
                if (iz==1)
                    age(iz)=dage(iz)*0.5d0;
                end
                if (iz~=1)
                    age(iz) = age(iz-1)+dage(iz-1)*0.5d0 + 0.5d0*dage(iz);
                end
            end
            
        end
        
        function [up, dwn, cnr, adf] = calcupwindscheme(w, nz)
            % ------------ determine variables to realize advection
            %  upwind scheme
            %  up  ---- burial advection at grid i = sporo(i)*w(i)*(some conc. at i) - sporo(i-1)*w(i-1)*(some conc. at i - 1)
            %  dwn ---- burial advection at grid i = sporo(i+1)*w(i+1)*(some conc. at i+1) - sporo(i)*w(i)*(some conc. at i)
            %  cnr ---- burial advection at grid i = sporo(i+1)*w(i+1)*(some conc. at i+1) - sporo(i-1)*w(i-1)*(some conc. at i - 1)
            %  when burial rate is positive, scheme need to choose up, i.e., up = 1.
            %  when burial rate is negative, scheme need to choose dwn, i.e., dwn = 1.
            %  where burial change from positive to negative or vice versa, scheme chooses cnr, i.e., cnr = 1. for the mass balance sake
            
            up = zeros(1, nz);
            dwn= zeros(1, nz);
            cnr = zeros(1, nz);
            adf= ones(1, nz);
            
            for iz=1:nz
                if (iz==1)
                    if (w(iz)>=0d0 && w(iz+1)>=0d0)  % positive burial
                        up(iz) = 1;
                    elseif (w(iz)<=0d0 && w(iz+1)<=0d0) then  % negative burial
                        dwn(iz) = 1;
                    else   %  where burial sign changes
                        if (~(w(iz)*w(iz+1) <=0d0))
                            fprintf('error #1 in calcupwindscheme');
                            return
                        end
                        cnr(iz) = 1;
                    end
                elseif (iz==nz)
                    if (w(iz)>=0d0 && w(iz-1)>=0d0)
                        up(iz) = 1;
                    elseif (w(iz)<=0d0 && w(iz-1)<=0d0)
                        dwn(iz) = 1;
                    else
                        if (~(w(iz)*w(iz-1) <= 0d0))
                            fprintf('error #2 in calcupwindscheme');
                            return
                        end
                        cnr(iz) = 1;
                    end
                else
                    if (w(iz) >=0d0)
                        if (w(iz+1)>=0d0 && w(iz-1)>=0d0)
                            up(iz) = 1;
                        else
                            cnr(iz) = 1;
                        end
                    else
                        if (w(iz+1)<=0d0 && w(iz-1)<=0d0)
                            dwn(iz) = 1;
                        else
                            cnr(iz) = 1;
                        end
                    end
                end
            end
            
        end
        
        
        function [trans,izrec,izrec2,izml,nonlocal] = make_transmx(labs,nspcc,turbo2,nobio,dz,sporo,nz,z, zml_ref, def_size)
            % % make transition matrix
            
            % allocate output arrrays
            trans = zeros(nz, nz, nspcc+2);
            nonlocal = false(1, nspcc + 2);     % % initial assumption
            % allocate local arrrays
            zml = zeros(1, nspcc+2);
            translabs = zeros(nz, nz);
            translabs_tmp = zeros(nz, nz);
            dbio = zeros(1, nz);
            transdbio = zeros(nz, nz);   % transition matrix to realize Fickian mixing with biodiffusion coefficient dbio which is defined just above
            transturbo2 = zeros(nz, nz);    % transition matrix for random mixing
            
            %             Dominik: Todo: test
            %~~~~~~~~~~~~ loading transition matrix from LABS ~~~~~~~~~~~~~~~~~~~~~~~~
            if (any(labs))
                %                translabs = 0d0
                nlabs = 7394;  % number of labs transition matrices to be read
                %    print*, 'loading LABS ' , labs
                for ilabs=1:nlabs
                    % %  Dominik: todo: Run and test with input
                    % %                    translabs_tmp=0d0  % transition matrix to be read
                    % %                   write(dumchr(1),'(i0)') ilabs*2000
                    % %                    open(unit=file_tmp, file='C:/Users/YK/Desktop/biot-res/trans-test-1kyr-frq-20190315/mix/' ...
                    % %                        //'transmtx-'//trim(adjustl(dumchr(1)))//'.txt',action='read',status='unknown');
                    % %                    for iz = 1:nz
                    % %                        read(file_tmp,*) translabs_tmp(iz,:)  % reading
                    % %                    end
                    % %                    close(file_tmp)
                    translabs = translabs + translabs_tmp;  % adding up all transition matrices
                end
                translabs = translabs/real(nlabs); % and averaging all transition matrices
            end
            
            
            if (true)   % devided by the time duration when transition matrices are created in LABS and weakening by a factor
                % if (.false.) then
                translabs = translabs *365.25d0/10d0*1d0/2.3d0;
                % translabs = translabs *365.25d0/10d0*1d0/10d0
            end
            %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            zml=zml + zml_ref; % mixed layer depth assumed to be a reference value at first
            
            zrec = 1.1d0*max(zml);  % depth where recording is made, aimed at slightly below the bottom of mixed layer
            zrec2 = 2.0d0*max(zml);  % depth where recording is made ver. 2, aimed at 2 time bottom depth of mixed layer
            
            
            if(def_size)
                zml(2+1)=20d0;   % fine species have larger mixed layers
                zml(2+2)=20d0;   % note that total number of solid species is 2 + nspcc including om, clay and nspcc of caco3; thus index has '2+'
                zml(2+3)=20d0;
                zml(2+4)=20d0;
                zrec = 1.1d0*min(zml);  % first recording is made below minimum depth of mixed layer
                zrec2 = 1.1d0*max(zml); % second recording is made below maximum depth of mixed layer
            end
            
            for iz=1:nz     % determine grid locations where signal recording is made
                if (z(iz)<=zrec)
                    izrec = iz;
                end
                if (z(iz)<=zrec2)
                    izrec2 = iz;
                end
            end
            
            % above nonlocal = false(1, nspcc + 2);     % % initial assumption
            for isp=1:nspcc+2
                if (turbo2(isp) || labs(isp))
                    nonlocal(isp)= true;    % % if mixing is made by turbo2 or labs, then nonlocal
                end
                dbio = zeros(1, nz);
                for iz = 1:nz
                    if (z(iz) <=zml(isp))
                        dbio(iz) =  0.15d0;   %  within mixed layer 150 cm2/kyr (Emerson, 1985)
                        izml = iz;   % determine grid of bottom of mixed layer
                    else
                        dbio(iz) =  0d0; % no biodiffusion in deeper depths
                    end
                end
                
                % transition matrix to realize Fickian mixing with biodiffusion coefficient dbio which is defined just above
                for iz = 1:izml
                    if (iz==1)
                        transdbio(iz,iz) = 0.5d0*(sporo(iz)*dbio(iz)+sporo(iz+1)*dbio(iz+1))*(-1d0)/(0.5d0*(dz(iz)+dz(iz+1)));
                        transdbio(iz+1,iz) = 0.5d0*(sporo(iz)*dbio(iz)+sporo(iz+1)*dbio(iz+1))*(1d0)/(0.5d0*(dz(iz)+dz(iz+1)));
                    elseif (iz==izml)
                        transdbio(iz,iz) = 0.5d0*(sporo(iz)*dbio(iz)+sporo(iz-1)*dbio(iz-1))*(-1d0)/(0.5d0*(dz(iz)+dz(iz-1)));
                        transdbio(iz-1,iz) = 0.5d0*(sporo(iz)*dbio(iz)+sporo(iz-1)*dbio(iz-1))*(1d0)/(0.5d0*(dz(iz)+dz(iz-1)));
                    else
                        transdbio(iz,iz) = 0.5d0*(sporo(iz)*dbio(iz)+sporo(iz-1)*dbio(iz-1))*(-1d0)/(0.5d0*(dz(iz)+dz(iz-1)))  ...
                            + 0.5d0*(sporo(iz)*dbio(iz)+sporo(iz+1)*dbio(iz+1))*(-1d0)/(0.5d0*(dz(iz)+dz(iz+1)));
                        transdbio(iz-1,iz) = 0.5d0*(sporo(iz)*dbio(iz)+sporo(iz-1)*dbio(iz-1))*(1d0)/(0.5d0*(dz(iz)+dz(iz-1)));
                        transdbio(iz+1,iz) = 0.5d0*(sporo(iz)*dbio(iz)+sporo(iz+1)*dbio(iz+1))*(1d0)/(0.5d0*(dz(iz)+dz(iz+1)));
                    end
                end
                
                %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                % transition matrix for random mixing
                transturbo2 = transturbo2 + 0.0015d0/1d0;  % arbitrary assumed probability
                for iz=1:izml  % when i = j, transition matrix contains probabilities with which particles are moved from other layers of sediment
                    transturbo2(iz,iz)=-0.0015d0*(izml-1)/1d0;
                end
                
                if (turbo2(isp))
                    translabs = transturbo2;   % translabs temporarily used to represents nonlocal mixing
                end
                trans(:,:,isp) = transdbio(:,:);  %  firstly assume local mixing implemented by dbio
                
                if (nonlocal(isp))
                    trans(:,:,isp) = translabs(:,:);  % if nonlocal, replaced by either turbo2 mixing or labs mixing
                end
                if (nobio(isp))
                    trans(:,:,isp) = 0d0;  % if assuming no bioturbation, transition matrix is set at zero
                end
            end % for
            
            % even when all are local Fickian mixing, mixing treatment must be the same as in case of nonlocal
            % if mixing intensity and depths are different between different species
            if (all(~nonlocal))
                for isp=1:nspcc+2-1
                    if (any(trans(:,:,isp+1)~=trans(:,:,isp)))
                        nonlocal=true;
                    end
                end
            end
            
        end
        
        function [omx, izox, kom] = omcalc(oxic,anoxic,o2x,om, komi,nz,sporo,sporoi,sporof, w,wi,dt,up,dwn,cnr,adf,trans,nspcc,labs,turbo2,nonlocal,omflx,poro,dz, o2th)
            
            % allocate output arrays
            omx = zeros(1, nz);         %  om conc. mol cm-3 sld;
            kom = zeros(1, nz);         %  degradation rate consts. for each nz grids
            % izox: integer for grid number of zox
            
            % allocate local array ( down further down)
            
            
            %  amx and ymx correspond to A and B in Ax = B
            %  dgesv subroutine of BLAS returns the solution x in ymx
            %  emx stores errors (not required for dgesv); ipiv is required for dgesv
            %  E.g., difference form of governing equation at grid 2 can be expressed as
            %
            %               sporo(2)*(omx(2)-om(2))/dt + (sporo(2)*w(2)*omx(2)-sporo(1)*w(1)*omx(1))/dz(2) + sporo(2)*kom(2)*omx(2) = 0
            %
            %  In above difference equation, w(2) and w(1) are assumed to be positive for illustration purpose (and dropping adf, up, dwn and cnr terms)
            %  and omx(2) and omx(1) are unknowns to be solved.
            %  x contains omx(1), omx(2), ...., omx(nz), i.e., nz unknowns
            %  and thus the above equation fills the matrix A as
            %
            %               A(2,1) =  (-sporo(1)*w(1)*1)/dz(2)
            %               A(2,2) =  sporo(2)*(1)/dt + (sporo(2)*w(2)*1)/dz(2) + sporo(2)*kom(2)*1
            %
            %  and the matrix B as
            %
            %               - B(2) = sporo(2)*(-om(2))/dt
            %
            %  Matrices A and B are filled in this way. Note again amx and ymx correspond A and B, respectively.
            
            izox_calc_done = false;
            if (oxic)
                for iz=1:nz
                    if (o2x(iz) > o2th)
                        kom(iz) = komi;
                        % izox = iz
                        if (~izox_calc_done)
                            izox = iz;      % set grid number down
                        end
                    else% unless anoxi degradation is allowed, om cannot degradate below zox
                        kom(iz) = 0d0
                        if (anoxic)
                            kom(iz) = komi;
                        end
                        izox_calc_done = true;
                    end
                end
            end
            
            nsp=1;               % number of species considered here; 1, only om
            nmx = nz*nsp;        % # of col (... row) of matrix A to in linear equations Ax = B to be solved, each species has nz (# of grids) unknowns
            
            %  allocate local matrices used to solve linear system Ax = B
            amx = zeros(nmx, nmx);          % amx corresponds to A in Ax = B, but ymx is also x when Ax = B is solved. emx is array of error
            ymx = zeros(1, nmx);            % ymx correspond to B in Ax = B, but ymx is also x when Ax = B is solved
            emx = zeros(1, nmx);            % emx is array of error
            ipiv = zeros(1, nmx);           % matrix used to solve linear system Ax = B
            
            %            amx = 0d0;
            %            ymx = 0d0;
            
            for iz = 1:nz
                row = 1 + (iz-1)*nsp;    % row number is obtained from grid number; here simply gird 1 corresponds to row 1
                if (iz == 1)            % need to reflect upper boundary, rain flux; and be careful that iz - 1 does not exit
                    ymx(row) = sporo(iz)*(-om(iz))/dt ...     % time change term
                        - omflx/dz(1);      % rain flux term
                    amx(row,row) = ( sporo(iz)*(1d0)/dt ...  % time change term
                        + adf(iz)*up(iz)*(sporo(iz)*w(iz)*1d0-sporoi*wi*0d0)/dz(1)   ...                        % advection terms
                        + adf(iz)*dwn(iz)*(sporo(iz+1)*w(iz+1)*0d0-sporo(iz)*w(iz)*1d0)/dz(1)   ...
                        + adf(iz)*cnr(iz)*(sporo(iz+1)*w(iz+1)*0d0-sporoi*wi*0d0)/dz(1)   ...
                        + sporo(iz)*kom(iz));   %  rxn term
                    % matrix filling at grid iz but for unknwon at grid iz + 1 (here only advection terms)
                    amx(row,row+nsp) =  (adf(iz)*dwn(iz)*(sporo(iz+1)*w(iz+1)*1d0-sporo(iz)*w(iz)*0d0)/dz(1)   ...
                        + adf(iz)*cnr(iz)*(sporo(iz+1)*w(iz+1)*1d0-sporoi*wi*0d0)/dz(1)    );
                elseif (iz == nz)    % need to reflect lower boundary; none; but must care that iz + 1 does not exist
                    ymx(row) = 0d0 + sporo(iz)*(-om(iz))/dt ;  % time change term
                    amx(row,row) = ( sporo(iz)*(1d0)/dt ... % time change term
                        + adf(iz)*up(iz)*(sporo(iz)*w(iz)*1d0-sporo(iz-1)*w(iz-1)*0d0)/dz(iz)  ...      % advection terms
                        + adf(iz)*dwn(iz)*(sporof*w(iz)*1d0-sporo(iz)*w(iz)*1d0)/dz(iz)  ...
                        + adf(iz)*cnr(iz)*(sporof*w(iz)*1d0-sporo(iz-1)*w(iz-1)*0d0)/dz(iz)  ...
                        + sporo(iz)*kom(iz) );  % rxn term
                    % filling matrix at grid iz but for unknown at grid iz-1 (only advection terms)
                    amx(row,row-nsp) = ( adf(iz)*up(iz)*(sporof*w(iz)*0d0-sporo(iz-1)*w(iz-1)*1d0)/dz(iz)  ...
                        + adf(iz)*cnr(iz)*(sporof*w(iz)*0d0-sporo(iz-1)*w(iz-1)*1d0)/dz(iz)  );
                else % do not have to care about boundaries
                    ymx(row) = 0d0 + sporo(iz)*(0d0-om(iz))/dt ; % time change term
                    amx(row,row) = (sporo(iz)*(1d0)/dt ...     % time change term
                        + adf(iz)*up(iz)*(sporo(iz)*w(iz)*1d0-sporo(iz-1)*w(iz-1)*0d0)/dz(iz)  ...      % advection terms
                        + adf(iz)*dwn(iz)*(sporo(iz+1)*w(iz+1)*0d0-sporo(iz)*w(iz)*1d0)/dz(iz)  ...
                        + adf(iz)*cnr(iz)*(sporo(iz+1)*w(iz+1)*0d0-sporo(iz-1)*w(iz-1)*0d0)/dz(iz)  ...
                        + sporo(iz)*kom(iz)   );        % rxn term
                    % filling matrix at grid iz but for unknown at grid iz+1 (only advection terms)
                    amx(row,row+nsp) =  (adf(iz)*dwn(iz)*(sporo(iz+1)*w(iz+1)*1d0-sporo(iz)*w(iz)*0d0)/dz(iz)  ...
                        + adf(iz)*cnr(iz)*(sporo(iz+1)*w(iz+1)*1d0-sporo(iz-1)*w(iz-1)*0d0)/dz(iz) );
                    % filling matrix at grid iz but for unknown at grid iz-1 (only advection terms)
                    amx(row,row-nsp) =  ( adf(iz)*up(iz)*(sporo(iz)*w(iz)*0d0-sporo(iz-1)*w(iz-1)*1d0)/dz(iz)  ...
                        + adf(iz)*cnr(iz)*(sporo(iz+1)*w(iz+1)*0d0-sporo(iz-1)*w(iz-1)*1d0)/dz(iz)  );
                end
                % diffusion terms are reflected with transition matrices
                if (turbo2(1) || labs(1))
                    for iiz = 1:nz
                        col = 1 + (iiz-1)*nsp ;
                        if (trans(iiz,iz,1)==0d0)
                            continue        % cycle in fortran
                        end
                        amx(row,col) = amx(row,col) -trans(iiz,iz,1)/dz(iz)*dz(iiz)*(1d0-poro(iiz));
                    end
                else
                    for iiz = 1:nz
                        col = 1 + (iiz-1)*nsp;
                        if (trans(iiz,iz,1)==0d0)
                            continue        % cycle in fortran
                        end
                        amx(row,col) = amx(row,col) -trans(iiz,iz,1)/dz(iz);
                    end
                end
            end % for
            
            ymx = - ymx';  % I have filled matrix B without changing signs; here I change signs at once & make a vector
            
            %            call dgesv(nmx,int(1),amx,nmx,ipiv,ymx,nmx,infobls)
            omx = matrixDivide(amx,ymx);
            %             omx = ymx; % now passing the solution to unknowns omx
            
            
        end
        
        function[omadv,omdec,omdif,omrain,omres,omtflx] = calcflxom(omflx,sporo,om,omx,dt,w,dz,z,nz,turbo2,labs,poro,up,dwn,cnr,adf,rho,mom,trans,kom,sporof) % ,nonlocal,sporoi,wi,nspcc)
            %% calculating the fluxes relevant to om diagenesis (and checking the calculation satisfies the difference equations )
            
            % preallocate output variables
            % fluxes, adv, dec, dif, res, t and rain denote burial, decomposition, dissolution, diffusion, residual, time change and rain fluxes, respectively
            omadv = 0d0;
            omdec = 0d0;
            omdif = 0d0;
            omres = 0d0;
            omtflx = 0d0;
            omrain = 0d0;
            
            nsp = 1;        % independent chemical variables (better specified as global var???)
            
            for iz = 1:nz
                row = 1 + (iz-1)*nsp;
                if (iz == 1)
                    omtflx = omtflx + sporo(iz)*(omx(iz)-om(iz))/dt*dz(iz) ;
                    omadv = omadv + adf(iz)*up(iz)*(sporo(iz)*w(iz)*omx(iz)-0d0)/dz(iz)*dz(iz)  ...
                        + adf(iz)*dwn(iz)*(sporo(iz+1)*w(iz+1)*omx(iz+1)-sporo(iz)*w(iz)*omx(iz))/dz(iz)*dz(iz)  ...
                        + adf(iz)*cnr(iz)*(sporo(iz+1)*w(iz+1)*omx(iz+1)-0d0)/dz(iz)*dz(iz);
                    omdec = omdec + sporo(iz)*kom(iz)*omx(iz)*dz(iz);
                    omrain = omrain - omflx/dz(1)*dz(iz);
                elseif (iz == nz)
                    omtflx = omtflx + sporo(iz)*(omx(iz)-om(iz))/dt*dz(iz);
                    omadv = omadv + adf(iz)*up(iz)*(sporo(iz)*w(iz)*omx(iz)-sporo(iz-1)*w(iz-1)*omx(iz-1))/dz(iz)*dz(iz) ...
                        + adf(iz)*dwn(iz)*(sporof*w(iz)*omx(iz)-sporo(iz)*w(iz)*omx(iz))/dz(iz)*dz(iz) ...
                        + adf(iz)*cnr(iz)*(sporof*w(iz)*omx(iz)-sporo(iz-1)*w(iz-1)*omx(iz-1))/dz(iz)*dz(iz);
                    omdec = omdec + sporo(iz)*kom(iz)*omx(iz)*dz(iz);
                else
                    omtflx = omtflx + sporo(iz)*(omx(iz)-om(iz))/dt*dz(iz);
                    omadv = omadv + adf(iz)*up(iz)*(sporo(iz)*w(iz)*omx(iz)-sporo(iz-1)*w(iz-1)*omx(iz-1))/dz(iz)*dz(iz)  ...
                        + adf(iz)*dwn(iz)*(sporo(iz+1)*w(iz+1)*omx(iz+1)-sporo(iz)*w(iz)*omx(iz))/dz(iz)*dz(iz)  ...
                        + adf(iz)*cnr(iz)*(sporo(iz+1)*w(iz+1)*omx(iz+1)-sporo(iz-1)*w(iz-1)*omx(iz-1))/dz(iz)*dz(iz);
                    omdec = omdec + sporo(iz)*kom(iz)*omx(iz)*dz(iz);
                end
                if (turbo2(1)||labs(1))
                    for iiz = 1:nz
                        if (trans(iiz,iz,1)==0d0)
                            continue        % cycle in fortran
                        end
                        omdif = omdif -trans(iiz,iz,1)/dz(iz)*dz(iiz)*(1d0-poro(iiz))*dz(iz)*omx(iiz);
                    end
                else
                    for iiz = 1:nz
                        if (trans(iiz,iz,1)==0d0)
                            continue        % cycle in fortran
                        end
                        omdif = omdif -trans(iiz,iz,1)/dz(iz)*dz(iz)*omx(iiz);  % check previous versions
                    end
                end
            end % for
            
            omres = omadv + omdec + omdif + omrain + omtflx;    % this is residual flux should be zero equations are exactly satisfied
            
            if (any(omx<0d0))  % if negative om conc. is detected, need to stop
                % fprintf('negative om, stop \n');
                str = sprintf('ERROR_MSG_NEGATIVE_OM.txt');
                file_tmp = fopen(str,'wt');
                for iz=1:nz
                    %        write(file_tmp,*)(trans(iz,iiz,isp),iiz=1,nz)
                    fprintf(file_tmp,'%17.16e \t %17.16e \t %17.16e \t %17.16e \t %17.16e \t %17.16e \t %17.16e \t ', z(iz),omx(iz)*mom/rho(iz)*100d0,w(iz),up(iz),dwn(iz),cnr(iz),adf(iz));
                    fprintf(file_tmp,'\n');
                end
                fclose(file_tmp);
                msg = 'Error: NEGATIVE OM, STOP.';
                error(msg)
            end
            if (any(isnan(omx)))  % if NAN, ... the same ... stop
                fprintf('OM concentr. omx:\n');
                fprintf(omx);
                msg = 'Error: NaN OM, STOP.';
                error(msg)
            end
            
            
        end
        
        function o2x = o2calc_ox(izox,nz,poro,o2,kom,omx,sporo,dif_o2,dz,dt, ox2om, o2i)
            %% o2 calculation when o2 penetration depth (zox) is the same as bottom depth.
            
            nsp=1;               % number of species considered here; 1, only om
            nmx = nz*nsp;        % # of col (... row) of matrix A to in linear equations Ax = B to be solved, each species has nz (# of grids) unknowns
            
            %  allocate local matrices used to solve linear system Ax = B
            amx = zeros(nmx, nmx);          % amx corresponds to A in Ax = B, but ymx is also x when Ax = B is solved. emx is array of error
            ymx = zeros(1, nmx);            % ymx correspond to B in Ax = B, but ymx is also x when Ax = B is solved
            
            
            for iz = 1:nz
                row = 1 + (iz-1)*nsp;
                if (iz == 1) % be careful about upper boundary
                    ymx(row) = ( poro(iz)*(0d0-o2(iz))/dt ...          % time change term
                        - ((poro(iz)*dif_o2(iz)+poro(iz+1)*dif_o2(iz+1))*0.5d0*(0d0)/(0.5d0*(dz(iz)+dz(iz+1))) ...   % diffusion term
                        - poro(iz)*dif_o2(iz)*(0d0-o2i*1d-6/1d3)/dz(iz))/dz(iz) ...     % rxn term
                        + sporo(iz)*ox2om*kom(iz)*omx(iz)  );
                    amx(row,row) = ( poro(iz)*(1d0)/dt ...   % time change term
                        - ((poro(iz)*dif_o2(iz)+poro(iz+1)*dif_o2(iz+1))*0.5d0*(-1d0)/(0.5d0*(dz(iz)+dz(iz+1))) ...     % diffusion term
                        - poro(iz)*dif_o2(iz)*(1d0)/dz(iz))/dz(iz)  );
                    % filling matrix at grid iz but for unknown at grid iz+1 (only diffusion term)
                    amx(row,row+nsp) = (-((poro(iz)*dif_o2(iz)+poro(iz+1)*dif_o2(iz+1))*0.5d0*(1d0)/(0.5d0*(dz(iz)+dz(iz+1))) - 0d0)/dz(iz) );
                elseif (iz == nz)      % be careful about lower boundary
                    ymx(row) = (0d0 ...
                        + poro(iz)*(0d0-o2(iz))/dt ...  % time change term
                        - (0d0 - 0.5d0*(poro(iz)*dif_o2(iz)+poro(iz-1)*dif_o2(iz-1))*(0d0)/(0.5d0*(dz(iz-1)+dz(iz))))/dz(iz) ...  % diffusion term
                        + sporo(iz)*ox2om*kom(iz)*omx(iz)  );    % rxn term
                    amx(row,row) = ( poro(iz)*(1d0)/dt ... % time change term
                        - (0d0 - 0.5d0*(poro(iz)*dif_o2(iz)+poro(iz-1)*dif_o2(iz-1))*(1d0)/(0.5d0*(dz(iz-1)+dz(iz))))/dz(iz) );     % diffusion term
                    % filling matrix at grid iz but for unknown at grid iz-1 (only diffusion term)
                    amx(row,row-nsp) = ( -(0d0 - 0.5d0*(poro(iz)*dif_o2(iz)+poro(iz-1)*dif_o2(iz-1))*(-1d0)/(0.5d0*(dz(iz-1)+dz(iz))))/dz(iz) );
                else
                    ymx(row) = ( 0d0 ...
                        + poro(iz)*(0d0-o2(iz))/dt ... % time change term
                        - (0.5d0*(poro(iz+1)*dif_o2(iz+1)+poro(iz)*dif_o2(iz))*(0d0)/(0.5d0*(dz(iz+1)+dz(iz))) ...    % diffusion term
                        - 0.5d0*(poro(iz)*dif_o2(iz)+poro(iz-1)*dif_o2(iz-1))*(0d0)/(0.5d0*(dz(iz)+dz(iz-1))))/dz(iz)  ...
                        + sporo(iz)*ox2om*kom(iz)*omx(iz)  );   % rxn term
                    amx(row,row) = ( poro(iz)*(1d0)/dt ... % time change term
                        - (0.5d0*(poro(iz+1)*dif_o2(iz+1)+poro(iz)*dif_o2(iz))*(-1d0)/(0.5d0*(dz(iz+1)+dz(iz))) ... % diffusion term
                        - 0.5d0*(poro(iz)*dif_o2(iz)+poro(iz-1)*dif_o2(iz-1))*(1d0)/(0.5d0*(dz(iz)+dz(iz-1))))/dz(iz)  );
                    % filling matrix at grid iz but for unknown at grid iz+1 (only diffusion term)
                    amx(row,row+nsp) = (-(0.5d0*(poro(iz+1)*dif_o2(iz+1)+poro(iz)*dif_o2(iz))*(1d0)/(0.5d0*(dz(iz+1)+dz(iz))) - 0d0)/dz(iz)  );
                    % filling matrix at grid iz but for unknown at grid iz-1 (only diffusion term)
                    amx(row,row-nsp) = (- (0d0 - 0.5d0*(poro(iz)*dif_o2(iz)+poro(iz-1)*dif_o2(iz-1))*(-1d0)/(0.5d0*(dz(iz)+dz(iz-1))))/dz(iz) );
                end
            end
            
            ymx = - ymx';  % sign & transpose change; see above for the case of om
            
            %             fprintf('amx   \n');
            %             % showing parameters relevant to burial on screen
            %             for iz=1:nz
            %                 fprintf('%17.16e \t %17.16e \n',amx(:,iz));
            %             end
            % call dgesv(nmx,int(1),amx,nmx,ipiv,ymx,nmx,infobls) % solving
            % o2x = ymx % passing solutions to unknowns
            o2x = matrixDivide(amx,ymx);
            
            %             fprintf('ymx   , o2x \n');
            %             % showing parameters relevant to burial on screen
            %             for iz=1:nz
            %                 fprintf('%17.16e \t %17.16e \n',ymx(iz), o2x(iz));
            %             end
            
        end
        
        
        function [o2dec,o2dif,o2tflx,o2res] = calcflxo2_ox(nz,sporo,kom,omx,dz,poro,dif_o2,dt,o2,o2x, ox2om,o2i)
            %% fluxes relevant to o2 (at the same time checking the satisfaction of difference equations)
            
            o2dec = 0d0;
            o2dif = 0d0;
            o2tflx = 0d0;
            
            for iz = 1:nz
                if (iz == 1)
                    o2dec = o2dec + sporo(iz)*ox2om*kom(iz)*omx(iz)*dz(iz);
                    o2tflx = o2tflx + (o2x(iz)-o2(iz))/dt*dz(iz)*poro(iz);
                    o2dif = o2dif - ((poro(iz)*dif_o2(iz)+poro(iz+1)*dif_o2(iz+1))*0.5d0*(o2x(iz+1)-o2x(iz))/(0.5d0*(dz(iz)+dz(iz+1))) ...
                        - poro(iz)*dif_o2(iz)*(o2x(iz)-o2i*1d-6/1d3)/dz(iz))/dz(iz) *dz(iz);
                elseif (iz == nz)
                    o2dec = o2dec + (1d0-poro(iz))*ox2om*kom(iz)*omx(iz)/poro(iz)*dz(iz)*poro(iz);
                    o2tflx = o2tflx + (o2x(iz)-o2(iz))/dt*dz(iz)*poro(iz);
                    o2dif = o2dif ...
                        - (0d0 - 0.5d0*(poro(iz)*dif_o2(iz)+poro(iz-1)*dif_o2(iz-1))*(o2x(iz)-o2x(iz-1))/(0.5d0*(dz(iz-1)+dz(iz))))/dz(iz) ...
                        *dz(iz);
                else
                    o2dec = o2dec + (1d0-poro(iz))*ox2om*kom(iz)*omx(iz)/poro(iz)*dz(iz)*poro(iz);
                    o2tflx = o2tflx + (o2x(iz)-o2(iz))/dt*dz(iz)*poro(iz);
                    o2dif = o2dif ...
                        - (0.5d0*(poro(iz+1)*dif_o2(iz+1)+poro(iz)*dif_o2(iz))*(o2x(iz+1)-o2x(iz))/(0.5d0*(dz(iz+1)+dz(iz))) ...
                        - 0.5d0*(poro(iz)*dif_o2(iz)+poro(iz-1)*dif_o2(iz-1))*(o2x(iz)-o2x(iz-1))/(0.5d0*(dz(iz)+dz(iz-1))))/dz(iz)  ...
                        *dz(iz);
                end
            end
            
            o2res = o2dec + o2dif + o2tflx;  % residual flux
            
        end
        
        function o2x = o2calc_sbox(izox,nz,poro,o2,kom,omx,sporo,dif_o2,dz,dt, ox2om, o2i)
            %% o2 calculation when o2 is depleted within the calculation domain.
            
            nsp=1;               % number of species considered here; 1, only om
            nmx = nz*nsp;        % # of col (... row) of matrix A to in linear equations Ax = B to be solved, each species has nz (# of grids) unknowns
            
            %  allocate local matrices used to solve linear system Ax = B
            amx = zeros(nmx, nmx);          % amx corresponds to A in Ax = B, but ymx is also x when Ax = B is solved. emx is array of error
            ymx = zeros(1, nmx);            % ymx correspond to B in Ax = B, but ymx is also x when Ax = B is solved
            
            
            for iz = 1:nz
                row = 1 + (iz-1)*nsp;
                if (iz == 1)
                    ymx(row) = ( poro(iz)*(0d0-o2(iz))/dt ... % time change
                        - ((poro(iz)*dif_o2(iz)+poro(iz+1)*dif_o2(iz+1))*0.5d0*(0d0)/(0.5d0*(dz(iz)+dz(iz+1))) ...   % diffusion
                        - poro(iz)*dif_o2(iz)*(0d0-o2i*1d-6/1d3)/dz(iz))/dz(iz)  ...
                        + sporo(iz)*ox2om*kom(iz)*omx(iz) );    % rxn
                    amx(row,row) = (poro(iz)*(1d0)/dt ...        % time change
                        - ((poro(iz)*dif_o2(iz)+poro(iz+1)*dif_o2(iz+1))*0.5d0*(-1d0)/(0.5d0*(dz(iz)+dz(iz+1))) ...  % diffusion
                        - poro(iz)*dif_o2(iz)*(1d0)/dz(iz))/dz(iz)	);
                    
                    % filling matrix at grid iz but for unknown at grid iz+1 (only diffusion term)
                    amx(row,row+nsp) = ( - ((poro(iz)*dif_o2(iz)+poro(iz+1)*dif_o2(iz+1))*0.5d0*(1d0)/(0.5d0*(dz(iz)+dz(iz+1))) - 0d0)/dz(iz) );
                elseif (iz>1 && iz<= izox)
                    ymx(row) = ( 0d0 ...
                        + poro(iz)*(0d0-o2(iz))/dt ...    % time change
                        - (0.5d0*(poro(iz+1)*dif_o2(iz+1)+poro(iz)*dif_o2(iz))*(0d0)/(0.5d0*(dz(iz+1)+dz(iz))) ...        % diffusion
                        - 0.5d0*(poro(iz)*dif_o2(iz)+poro(iz-1)*dif_o2(iz-1))*(0d0)/(0.5d0*(dz(iz)+dz(iz-1))))/dz(iz)  ...
                        + sporo(iz)*ox2om*kom(iz)*omx(iz)  );   % rxn
                    amx(row,row) = ( poro(iz)*(1d0)/dt ...            % time change
                        - (0.5d0*(poro(iz+1)*dif_o2(iz+1)+poro(iz)*dif_o2(iz))*(-1d0)/(0.5d0*(dz(iz+1)+dz(iz))) ... % diffusion
                        - 0.5d0*(poro(iz)*dif_o2(iz)+poro(iz-1)*dif_o2(iz-1))*(1d0)/(0.5d0*(dz(iz)+dz(iz-1))))/dz(iz) );
                    
                    % filling matrix at grid iz but for unknown at grid iz+1 (only diffusion term)
                    amx(row,row+nsp) = ( - (0.5d0*(poro(iz+1)*dif_o2(iz+1)+poro(iz)*dif_o2(iz))*(1d0)/(0.5d0*(dz(iz+1)+dz(iz))) - 0d0)/dz(iz)  );
                    % filling matrix at grid iz but for unknown at grid iz-1 (only diffusion term)
                    amx(row,row-nsp) = ( - (0d0 - 0.5d0*(poro(iz)*dif_o2(iz)+poro(iz-1)*dif_o2(iz-1))*(-1d0)/(0.5d0*(dz(iz)+dz(iz-1))))/dz(iz) );
                elseif (iz> izox)   % at lower than zox; zero conc. is forced
                    ymx(row) =  0d0;
                    amx(row,row) = 1d0;
                end
            end
            
            ymx = - ymx';  % change signs & transpose
            
            %                         call dgesv(nmx,int(1),amx,nmx,ipiv,ymx,nmx,infobls) % solving
            %                         o2x = ymx % passing solution to variable
            o2x = matrixDivide(amx,ymx);
            
        end
        
        
        function [o2dec,o2dif,o2tflx,o2res] = calcflxo2_sbox(nz,sporo,kom,omx,dz,poro,dif_o2,dt,o2,o2x,izox, ox2om, o2i)
            %% o2 calculation of fluxes relevant to oxygen
            
            o2dec = 0d0;
            o2dif = 0d0;
            o2tflx = 0d0;
            
            for iz = 1:nz
                if (iz == 1)
                    o2dec = o2dec + sporo(iz)*ox2om*kom(iz)*omx(iz)*dz(iz);
                    o2tflx = o2tflx + (o2x(iz)-o2(iz))/dt*dz(iz)*poro(iz);
                    o2dif = o2dif - ((poro(iz)*dif_o2(iz)+poro(iz+1)*dif_o2(iz+1))*0.5d0*(o2x(iz+1)-o2x(iz))/(0.5d0*(dz(iz)+dz(iz+1))) ...
                        - poro(iz)*dif_o2(iz)*(o2x(iz)-o2i*1d-6/1d3)/dz(iz))/dz(iz) *dz(iz);
                elseif (iz>1 && iz<=izox)
                    o2dec = o2dec + (1d0-poro(iz))*ox2om*kom(iz)*omx(iz)/poro(iz)*dz(iz)*poro(iz);
                    o2tflx = o2tflx + (o2x(iz)-o2(iz))/dt*dz(iz)*poro(iz);
                    o2dif = o2dif ...
                        - (0.5d0*(poro(iz+1)*dif_o2(iz+1)+poro(iz)*dif_o2(iz))*(o2x(iz+1)-o2x(iz))/(0.5d0*(dz(iz+1)+dz(iz))) ...
                        - 0.5d0*(poro(iz)*dif_o2(iz)+poro(iz-1)*dif_o2(iz-1))*(o2x(iz)-o2x(iz-1))/(0.5d0*(dz(iz)+dz(iz-1))))/dz(iz) ...
                        *dz(iz);
                end
            end
            
            o2res = o2dec + o2dif + o2tflx;  % residual flux
            
        end
        
        
        function [ccx,dicx,alkx,rcc,dt, flg_500, itr] = calccaco3sys(ccx,dicx,alkx,rcc,dt, nspcc,dic,alk,dep,sal,tmp,labs,turbo2,nonlocal,sporo,sporoi,sporof,poro,dif_alk,dif_dic, ...
                w,up,dwn,cnr,adf,dz,trans,cc,oxco2,anco2,co3sat,kcc,ccflx,ncc, nz, ...
                tol, poroi, flg_500, fact, alki,dici,ccx_th , def_nonrec, def_sparse, def_showiter, def_sense)
            
            %% ~~~~~~~~~~~~~~~~~~~~~~ CaCO3 solid, ALK and DIC  calculation ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            
            % some input variables
            % ncc = 4.5d0   % (Archer et al. 1989) reaction order for caco3 dissolution
            % kcc(nz,nspcc)   caco3 dissolution rate consts.
            
            %  allocate output matrices
            %            rcc = zeros(nz, nspcc);         % dissolution rate of caco3
            %            ccx = zeros(nz, nspcc);        %  concentration of caco3, mol cm-3 sld;
            %            dicx = zeros(1, nz);       % mol cm-3 porewater; dic, x denotes dummy variables
            %            alkx = zeros(1, nz);       % mol cm-3 porewater; alk, x denotes dummy variables
            
            
            %  allocate local matrices
            prox = zeros(1,nz);
            co2x = zeros(1,nz);
            hco3x = zeros(1,nz);
            co3x = zeros(1,nz);
            dco3_ddic = zeros(1,nz);
            dco3_dalk = zeros(1,nz);
            drcc_dcc = zeros(nz, nspcc);
            
            drcc_dco3 = zeros(nz, nspcc);
            drcc_ddic = zeros(nz, nspcc);
            drcc_dalk = zeros(nz, nspcc);
            
            drcc_domega = zeros(nz, nspcc);
            domega_dalk = zeros(1,nz);
            domega_ddic = zeros(1,nz);
            omega = zeros(1,nz);
            
            %       Here the system is non-linear and thus Newton's method is used (e.g., Steefel and Lasaga, 1994).
            %
            %       Problem: f(x) = 0  (note that x is array, containing N unknowns, here N = nmx)
            %       Expanding around x0 (initial/previous guess),
            %           f(x) =  f(x0) + f'(x0)*(x-x0) + O((x-x0)**2)
            %       where f'(x0) represents a Jacobian matrix.
            %       Solution is sought by iteration
            %           x = x0 - f'(x0)^-1*f(x0) or x- x0 = - f'(x0)^-1*f(x0)
            %       More practically by following steps.
            %           (1) solving Ax = B where A = f'(x0) and B = -f(x0), which gives delta (i.e., x - x0) (again these are arrays)
            %           (2) update solution by x = x0 + delta
            %           (3) replace x as x0
            %       Three steps are repeated until relative solution difference (i.e., delta) becomes negligible.
            %
            %       Now matrices A and B are Jacobian and array that contains f(x), respectively, represented by amx and ymx in this code
            %
            %       E.g., if equation at grid iz for caco3 is given by (for simplicity it cuts off several terms)
            %           (sporo(iz)*ccx(iz)-sporo(iz)*cc(iz))/dt + (sporo(iz)*w(iz)*ccx(iz)-sporo(iz-1)*w(iz-1)*ccx(iz-1))/dz(iz) + rcc(iz) = 0
            %       Then,
            %           B(row) = - left-hand side
            %       where row is the row number of caco3 at grid iz, i.e., row = 1+(iz-1)*nsp + isp -1 where isp = 1,..., nspcc,
            %       and
            %           A(row,row) = -dB(row)/dccx(iz) = (sporo(iz))/dt + (sporo(iz)*w(iz)*1)/dz(iz) + drcc_dcc(iz)
            %           A(row,row-nsp) = -dB(row)/dccx(iz-1) = (-sporo(iz-1)*w(iz-1)*1)/dz(iz)
            %           A(row,row+nspcc-1+1) = -dB(row)/ddic(iz) = drcc_ddic(iz)
            %           A(row,row+nspcc-1+2) = -dB(row)/dalk(iz) = drcc_dalk(iz)
            %       Something like this.
            %       Note, however, the present code uses ln (conc.) as x following Steefel and Lasaga (1994). So treatment is correspondingly a bit different.
            %       E.g.,
            %           dB(row)/dln(alk(iz)) = dB(row)/dalk(iz)*dalk(iz)/dln(alk(iz)) =  dB(row)/dalk(iz) * alkx(iz) = drcc_dalk(iz)*alkx(iz)
            %           ln x = ln x0 + delta, or, x = x0*exp(delta)
            %
            %       See e.g., Steefel and Lasaga (1994) for more details.
            
            
            flg_500 = false;
            error = 1d4;
            itr = 0;
            
            nsp = 2 + nspcc;  % now considered species are dic, alk and nspcc of caco3
            nmx = nz*nsp;  % col (and row) of matrix; the same number of unknowns
            
            %  allocate local matrices used to solve linear system Ax = B
            amx = zeros(nmx, nmx);          % amx corresponds to A in Ax = B, but ymx is also x when Ax = B is solved. emx is array of error
            ymx = zeros(1, nmx);            % ymx correspond to B in Ax = B, but ymx is also x when Ax = B is solved
            %            emx = zeros(1, nmx);            % emx is array of error
            %            ipiv = zeros(1, nmx);           % matrix used to solve linear system Ax = B
            dumx = zeros(nmx, nmx);          % amx corresponds to A in Ax = B, but ymx is also x when Ax = B is solved. emx is array of error
            
            
            %%% TODO: translate from here
            while(error > tol)
                amx = zeros(nmx, nmx);          % amx corresponds to A in Ax = B, but ymx is also x when Ax = B is solved. emx is array of error
                ymx = zeros(1, nmx);            % ymx correspond to B in Ax = B, but ymx is also x when Ax = B is solved
                
                % calling subroutine from caco3_therm.f90 to calculate aqueous co2 species
                % call calcspecies(dicx,alkx,temp,sal,dep,prox,co2x,hco3x,co3x,nz,infosbr)
                [prox,co2x,hco3x,co3x,infosbr] = caco3_therm.calcspecies(dicx,alkx,tmp,sal,dep);
                %                co3x = co3x';
                if (infosbr==1) % which means error in calculation
                    dt=dt/10d0;
                    if(def_sense)
                        % go to 500
                        flg_500= true;
                        return
                    else
                        % in fortran stop
                        msg = 'Error:infosbr==1 - error in calcspecies, STOP.';
                        error(msg)
                    end
                end
                % calling subroutine from caco3_therm.f90 to calculate derivatives of co3 wrt alk and dic
                % call calcdevs(dicx,alkx,temp,sal,dep,nz,infosbr,dco3_dalk,dco3_ddic)
                [dco3_dalk,dco3_ddic, infosbr] = caco3_therm.calcdevs(dicx,alkx,tmp,sal,dep);
                if (infosbr==1) % which means error in calculation
                    dt=dt/10d0;
                    if(def_sense)
                        % go to 500
                        flg_500= true;
                        return
                    else
                        % in fortran stop
                        msg = 'Error:infosbr==1 - error in calcdevs, STOP.';
                        error(msg)
                    end
                end
                
                for isp=1:nspcc
                    % calculation of dissolution rate for individual species
                    rcc(:,isp) = kcc(:,isp).*ccx(:,isp).*abs(1d0-co3x(:).*1d3/co3sat).^(ncc).*((1d0-co3x(:)*1d3/co3sat)>0d0);    % ((1d0-co3x(:)*1d3/co3sat)>0d0): 1 if true; 0 if false
                    % calculation of derivatives of dissolution rate wrt conc. of caco3 species, dic and alk
                    drcc_dcc(:,isp) = kcc(:,isp).*abs(1d0-co3x(:).*1d3/co3sat).^ncc.*((1d0-co3x(:)*1d3/co3sat)>0d0);    % ((1d0-co3x(:)*1d3/co3sat)>0d0): 1 if true; 0 if false
                    drcc_dco3(:,isp) = kcc(:,isp).*ccx(:,isp)*ncc.*abs(1d0-co3x(:)*1d3/co3sat).^(ncc-1d0) .*((1d0-co3x(:)*1d3/co3sat)>0d0) *(-1d3/co3sat);    % ((1d0-co3x(:)*1d3/co3sat)>0d0): 1 if true; 0 if false
                    drcc_ddic(:,isp) = drcc_dco3(:,isp).*dco3_ddic(:);
                    drcc_dalk(:,isp) = drcc_dco3(:,isp).*dco3_dalk(:);
                end
                
                % todo from here
                % % TODO: add mocsy option here 
                
                for iz = 1:nz
                    row = 1 + (iz-1)*nsp;
                    if (iz == 1)  % when upper condition must be taken account; *** comments for matrix filling are given only in this case
                        for isp = 1:nspcc  % multiple caco3 species
                            % put f(x) for isp caco3 species
                            ymx(row+isp-1) =  + sporo(iz)*(ccx(iz,isp)-cc(iz,isp))/dt ...
                                - ccflx(isp)/dz(1) ...
                                + adf(iz)*up(iz)*(sporo(iz)*w(iz)*ccx(iz,isp)-0d0)/dz(1)  ...
                                + adf(iz)*dwn(iz)*(sporo(iz+1)*w(iz+1)*ccx(iz+1,isp)-sporo(iz)*w(iz)*ccx(iz,isp))/dz(1)  ...
                                + adf(iz)*cnr(iz)*(sporo(iz+1)*w(iz+1)*ccx(iz+1,isp)-0d0)/dz(1)  ...
                                + sporo(iz)*rcc(iz,isp);
                            % derivative of f(x) wrt isp caco3 conc. at grid iz in ln
                            amx(row+isp-1,row+isp-1) = ( sporo(iz)*(1d0)/dt ...
                                + adf(iz)*up(iz)*(sporo(iz)*w(iz)*1d0-0d0)/dz(1)   ...
                                + adf(iz)*dwn(iz)*(0d0-sporo(iz)*w(iz)*1d0)/dz(1)  ...
                                + sporo(iz)* drcc_dcc(iz,isp)  ...
                                )* ccx(iz,isp);
                            
                            % derivative of f(x) wrt isp caco3 conc. at grid iz+1 in ln
                            amx(row+isp-1,row+isp-1+nsp) =  (adf(iz)*dwn(iz)*(sporo(iz+1)*w(iz+1)*1d0-0d0)/dz(1)  ...
                                + adf(iz)*cnr(iz)*(sporo(iz+1)*w(iz+1)*1d0-0d0)/dz(1)  ...
                                )*ccx(iz+1,isp);
                            % derivative of f(x) wrt dic conc. at grid iz in ln
                            amx(row+isp-1,row+nspcc) = (...
                                + sporo(iz)*drcc_ddic(iz,isp)  ...
                                )*dicx(iz);
                            % derivative of f(x) wrt alk conc. at grid iz in ln
                            amx(row+isp-1,row+nspcc+1) = ( sporo(iz)*drcc_dalk(iz,isp)  ...
                                )*alkx(iz);
                            % DIC
                            % derivative of f(x) for dic at iz wrt isp caco3 conc. at grid iz in ln
                            amx(row+nspcc,row+isp-1) = (...
                                - (1d0-poro(iz))*drcc_dcc(iz,isp)  ...
                                )*ccx(iz,isp)*fact;
                            % ALK
                            % derivative of f(x) for alk at iz wrt isp caco3 conc. at grid iz in ln
                            amx(row+nspcc+1,row+isp-1) = (...
                                - 2d0* (1d0-poro(iz))*drcc_dcc(iz,isp)  ...
                                )*ccx(iz,isp)*fact;
                        end
                        %  DIC
                        % put f(x) for dic at iz
                        ymx(row+nspcc) = ( ...
                            + poro(iz)*(dicx(iz)-dic(iz))/dt ...
                            - ((poro(iz)*dif_dic(iz)+poro(iz+1)*dif_dic(iz+1))*0.5d0*(dicx(iz+1)-dicx(iz))/(0.5d0*(dz(iz)+dz(iz+1))) ...
                            - poro(iz)*dif_dic(iz)*(dicx(iz)-dici*1d-6/1d3)/dz(iz))/dz(iz)  ...
                            - oxco2(iz) ...
                            - anco2(iz) ...
                            - (1d0-poro(iz))*sum(rcc(iz,:))  ...
                            )*fact;
                        % put derivative of f(x) for dic at iz wrt dic at iz in ln
                        amx(row+nspcc,row+nspcc) = (...
                            + poro(iz)*(1d0)/dt ...
                            - ((poro(iz)*dif_dic(iz)+poro(iz+1)*dif_dic(iz+1))*0.5d0*(-1d0)/(0.5d0*(dz(iz)+dz(iz+1))) ...
                            - poro(iz)*dif_dic(iz)*(1d0)/dz(iz))/dz(iz)...
                            - (1d0-poro(iz))*sum(drcc_ddic(iz,:))  ...
                            )*dicx(iz)*fact;
                        % put derivative of f(x) for dic at iz wrt dic at iz+1 in ln
                        amx(row+nspcc,row+nspcc+nsp) = (...
                            - ((poro(iz)*dif_dic(iz)+poro(iz+1)*dif_dic(iz+1))*0.5d0*(1d0)/(0.5d0*(dz(iz)+dz(iz+1))) ...
                            - 0d0)/dz(iz)...
                            )*dicx(iz+1)*fact;
                        % put derivative of f(x) for dic at iz wrt alk at iz in ln
                        amx(row+nspcc,row+nspcc+1) = ( ...
                            - (1d0-poro(iz))*sum(drcc_dalk(iz,:))  ...
                            )*alkx(iz)*fact;
                        % ALK
                        % put f(x) for alk at iz
                        ymx(row+nspcc+1) = (...
                            + poro(iz)*(alkx(iz)-alk(iz))/dt ...
                            - ((poro(iz)*dif_alk(iz)+poro(iz+1)*dif_alk(iz+1))*0.5d0*(alkx(iz+1)-alkx(iz))/(0.5d0*(dz(iz)+dz(iz+1))) ...
                            - poro(iz)*dif_alk(iz)*(alkx(iz)-alki*1d-6/1d3)/dz(iz))/dz(iz) ...
                            - anco2(iz) ...
                            - 2d0* (1d0-poro(iz))*sum(rcc(iz,:))  ...
                            )*fact;
                        % put derivative of f(x) for alk at iz wrt alk at iz in ln
                        amx(row+nspcc+1,row+nspcc+1) = (...
                            + poro(iz)*(1d0)/dt ...
                            - ((poro(iz)*dif_alk(iz)+poro(iz+1)*dif_alk(iz+1))*0.5d0*(-1d0)/(0.5d0*(dz(iz)+dz(iz+1))) ...
                            - poro(iz)*dif_alk(iz)*(1d0)/dz(iz))/dz(iz)  ...
                            - 2d0* (1d0-poro(iz))*sum(drcc_dalk(iz,:))  ...
                            )*alkx(iz)*fact;
                        % put derivative of f(x) for alk at iz wrt alk at iz+1 in ln
                        amx(row+nspcc+1,row+nspcc+1+nsp) = (...
                            - ((poro(iz)*dif_alk(iz)+poro(iz+1)*dif_alk(iz+1))*0.5d0*(1d0)/(0.5d0*(dz(iz)+dz(iz+1))) ...
                            - 0d0)/dz(iz)...
                            )*alkx(iz+1)*fact;
                        % put derivative of f(x) for alk at iz wrt dic at iz in ln
                        amx(row+nspcc+1,row+nspcc) = (...
                            - 2d0* (1d0-poro(iz))*sum(drcc_ddic(iz,:))  ...
                            )*dicx(iz)*fact;
                    elseif (iz == nz)   % need be careful about lower boundary condition; no diffusive flux from the bottom
                        for isp=1:nspcc
                            ymx(row+isp-1) = ...
                                + sporo(iz)*(ccx(iz,isp)-cc(iz,isp))/dt ...
                                + adf(iz)*up(iz)*(sporo(iz)*w(iz)*ccx(iz,isp)-sporo(iz-1)*w(iz-1)*ccx(iz-1,isp))/dz(iz)  ...
                                + adf(iz)*cnr(iz)*(sporof*w(iz)*ccx(iz,isp)-sporo(iz-1)*w(iz-1)*ccx(iz-1,isp))/dz(iz)  ...
                                + adf(iz)*dwn(iz)*(sporof*w(iz)*ccx(iz,isp)-sporo(iz)*w(iz)*ccx(iz,isp))/dz(iz)  ...
                                + sporo(iz)*rcc(iz,isp);
                            amx(row+isp-1,row+isp-1) = (...
                                + sporo(iz)*(1d0)/dt ...
                                + adf(iz)*up(iz)*(sporo(iz)*w(iz)*1d0-0d0)/dz(iz)  ...
                                + adf(iz)*cnr(iz)*(sporof*w(iz)*1d0-0d0)/dz(iz)  ...
                                + adf(iz)*dwn(iz)*(sporof*w(iz)*1d0-sporo(iz)*w(iz)*1d0)/dz(iz)  ...
                                + sporo(iz)*drcc_dcc(iz,isp)   ...
                                )*ccx(iz,isp);
                            amx(row+isp-1,row+isp-1-nsp) = ( ...
                                + adf(iz)*up(iz)*(0d0-sporo(iz-1)*w(iz-1)*1d0)/dz(iz)  ...
                                + adf(iz)*cnr(iz)*(0d0-sporo(iz-1)*w(iz-1)*1d0)/dz(iz)  ...
                                )*ccx(iz-1,isp);
                            amx(row+isp-1,row+nspcc) = (...
                                + sporo(iz)*drcc_ddic(iz,isp) ...
                                )*dicx(iz);
                            amx(row+isp-1,row+nspcc+1) = (...
                                + sporo(iz)*drcc_dalk(iz,isp) ...
                                )*alkx(iz);
                            
                            %DIC
                            amx(row+nspcc,row+isp-1) = (...
                                - sporo(iz)*drcc_dcc(iz,isp)  ...
                                )*ccx(iz,isp)*fact;
                            %ALK
                            amx(row+nspcc+1,row+isp-1) = (...
                                - 2d0*sporo(iz)*drcc_dcc(iz,isp)  ...
                                )*ccx(iz,isp)*fact;
                        end
                        % DIC
                        ymx(row+nspcc) = (...
                            + poro(iz)*(dicx(iz)-dic(iz))/dt ...
                            - (0d0 - 0.5d0*(poro(iz)*dif_dic(iz)+poro(iz-1)*dif_dic(iz-1))*(dicx(iz)-dicx(iz-1))  ...
                            /(0.5d0*(dz(iz-1)+dz(iz))))/dz(iz) ...
                            - oxco2(iz) ...
                            - anco2(iz) ...
                            - sporo(iz)*sum(rcc(iz,:))  ...
                            )*fact;
                        amx(row+nspcc,row+nspcc) = ( ...
                            + poro(iz)*(1d0)/dt ...
                            - (0d0 - 0.5d0*(poro(iz)*dif_dic(iz)+poro(iz-1)*dif_dic(iz-1))*(1d0)/(0.5d0*(dz(iz-1)+dz(iz))))/dz(iz) ...
                            - sporo(iz)*sum(drcc_ddic(iz,:))  ...
                            )*dicx(iz)*fact;
                        amx(row+nspcc,row+nspcc-nsp) = ( ...
                            - (0d0 - 0.5d0*(poro(iz)*dif_dic(iz)+poro(iz-1)*dif_dic(iz-1))*(-1d0)/(0.5d0*(dz(iz-1)+dz(iz))))/dz(iz) ...
                            ) * dicx(iz-1)*fact;
                        amx(row+nspcc,row+nspcc+1) = (...
                            - sporo(iz)*sum(drcc_dalk(iz,:))  ...
                            )*alkx(iz)*fact;
                        % ALK
                        ymx(row+nspcc+1) = ( ...
                            + poro(iz)*(alkx(iz)-alk(iz))/dt ...
                            - (0d0 - 0.5d0*(poro(iz)*dif_alk(iz)+poro(iz-1)*dif_alk(iz-1))*(alkx(iz)-alkx(iz-1))/(0.5d0*(dz(iz-1)+dz(iz))))/dz(iz) ...
                            - anco2(iz) ...
                            - 2d0*sporo(iz)*sum(rcc(iz,:))  ...
                            )*fact;
                        amx(row+nspcc+1,row+nspcc+1) = ( ...
                            + poro(iz)*(1d0)/dt ...
                            - (0d0 - 0.5d0*(poro(iz)*dif_alk(iz)+poro(iz-1)*dif_alk(iz-1))*(1d0)/(0.5d0*(dz(iz-1)+dz(iz))))/dz(iz) ...
                            - 2d0*sporo(iz)*sum(drcc_dalk(iz,:))  ...
                            )*alkx(iz)*fact;
                        amx(row+nspcc+1,row+nspcc+1-nsp) = ( ...
                            - (0d0 - 0.5d0*(poro(iz)*dif_alk(iz)+poro(iz-1)*dif_alk(iz-1))*(-1d0)/(0.5d0*(dz(iz-1)+dz(iz))))/dz(iz) ...
                            ) * alkx(iz-1)*fact;
                        amx(row+nspcc+1,row+nspcc) = (...
                            - 2d0*sporo(iz)*sum(drcc_ddic(iz,:))  ...
                            )*dicx(iz)*fact;
                    else %  do not have to be careful abount boundary conditions
                        for isp=1:nspcc
                            ymx(row+isp-1) = ...
                                + sporo(iz)*(ccx(iz,isp)-cc(iz,isp))/dt ...
                                + adf(iz)*up(iz)*(sporo(iz)*w(iz)*ccx(iz,isp)-sporo(iz-1)*w(iz-1)*ccx(iz-1,isp))/dz(iz)  ...
                                + adf(iz)*dwn(iz)*(sporo(iz+1)*w(iz+1)*ccx(iz+1,isp)-sporo(iz)*w(iz)*ccx(iz,isp))/dz(iz)  ...
                                + adf(iz)*cnr(iz)*(sporo(iz+1)*w(iz+1)*ccx(iz+1,isp)-sporo(iz-1)*w(iz-1)*ccx(iz-1,isp))/dz(iz)  ...
                                + sporo(iz)*rcc(iz,isp);
                            amx(row+isp-1,row+isp-1) = (...
                                + sporo(iz)*(1d0)/dt ...
                                + adf(iz)*up(iz)*(sporo(iz)*w(iz)*1d0-0d0)/dz(iz)  ...
                                + adf(iz)*dwn(iz)*(0d0-sporo(iz)*w(iz)*1d0)/dz(iz)  ...
                                + sporo(iz)*drcc_dcc(iz,isp)  ...
                                )*ccx(iz,isp);
                            amx(row+isp-1,row+isp-1+nsp) =  (...
                                + adf(iz)*dwn(iz)*(sporo(iz+1)*w(iz+1)*1d0-0d0)/dz(iz)  ...
                                + adf(iz)*cnr(iz)*(sporo(iz+1)*w(iz+1)*1d0-0d0)/dz(iz)  ...
                                )*ccx(iz+1,isp);
                            amx(row+isp-1,row+isp-1-nsp) =  (...
                                + adf(iz)*up(iz)*(0d0-sporo(iz-1)*w(iz-1)*1d0)/dz(iz)  ...
                                + adf(iz)*cnr(iz)*(0d0-sporo(iz-1)*w(iz-1)*1d0)/dz(iz)  ...
                                )*ccx(iz-1,isp);
                            amx(row+isp-1,row+nspcc) = (...
                                + sporo(iz)*drcc_ddic(iz,isp)  ...
                                )*dicx(iz);
                            amx(row+isp-1,row+nspcc+1) = (...
                                + sporo(iz)*drcc_dalk(iz,isp) ...
                                )*alkx(iz);
                            % DIC
                            amx(row+nspcc,row+isp-1) = (...
                                - sporo(iz)*drcc_dcc(iz,isp)  ...
                                )*ccx(iz,isp)*fact;
                            % ALK
                            amx(row+nspcc+1,row+isp-1) = (...
                                - 2d0*sporo(iz)*drcc_dcc(iz,isp)  ...
                                )*ccx(iz,isp)*fact;
                        end
                        % DIC
                        ymx(row+nspcc) = ( ...
                            + poro(iz)*(dicx(iz)-dic(iz))/dt ...
                            - (0.5d0*(poro(iz+1)*dif_dic(iz+1)+poro(iz)*dif_dic(iz))*(dicx(iz+1)-dicx(iz))/(0.5d0*(dz(iz+1)+dz(iz))) ...
                            - 0.5d0*(poro(iz)*dif_dic(iz)+poro(iz-1)*dif_dic(iz-1))*(dicx(iz)-dicx(iz-1))/(0.5d0*(dz(iz)+dz(iz-1))))/dz(iz)  ...
                            - oxco2(iz) ...
                            - anco2(iz) ...
                            - sporo(iz)*sum(rcc(iz,:))  ...
                            )*fact;
                        amx(row+nspcc,row+nspcc) = (...
                            + poro(iz)*(1d0)/dt ...
                            - (0.5d0*(poro(iz+1)*dif_dic(iz+1)+poro(iz)*dif_dic(iz))*(-1d0)/(0.5d0*(dz(iz+1)+dz(iz))) ...
                            - 0.5d0*(poro(iz)*dif_dic(iz)+poro(iz-1)*dif_dic(iz-1))*(1d0)/(0.5d0*(dz(iz)+dz(iz-1))))/dz(iz)  ...
                            - sporo(iz)*sum(drcc_ddic(iz,:))  ...
                            )*dicx(iz)*fact;
                        amx(row+nspcc,row+nspcc+nsp) = (...
                            - (0.5d0*(poro(iz+1)*dif_dic(iz+1)+poro(iz)*dif_dic(iz))*(1d0)/(0.5d0*(dz(iz+1)+dz(iz))) ...
                            - 0d0)/dz(iz)  ...
                            )*dicx(iz+1)*fact;
                        amx(row+nspcc,row+nspcc-nsp) = (...
                            - (0d0 ...
                            - 0.5d0*(poro(iz)*dif_dic(iz)+poro(iz-1)*dif_dic(iz-1))*(-1d0)/(0.5d0*(dz(iz)+dz(iz-1))))/dz(iz) ...
                            )*dicx(iz-1)*fact;
                        amx(row+nspcc,row+nspcc+1) = (...
                            - sporo(iz)*sum(drcc_dalk(iz,:))  ...
                            )*alkx(iz)*fact;
                        % ALK
                        ymx(row+nspcc+1) = (...
                            + poro(iz)*(alkx(iz)-alk(iz))/dt ...
                            - (0.5d0*(poro(iz+1)*dif_alk(iz+1)+poro(iz)*dif_alk(iz))*(alkx(iz+1)-alkx(iz))/(0.5d0*(dz(iz+1)+dz(iz))) ...
                            - 0.5d0*(poro(iz)*dif_alk(iz)+poro(iz-1)*dif_alk(iz-1))*(alkx(iz)-alkx(iz-1))/(0.5d0*(dz(iz)+dz(iz-1))))/dz(iz) ...
                            - anco2(iz) ...
                            - 2d0*sporo(iz)*sum(rcc(iz,:))  ...
                            ) *fact;
                        amx(row+nspcc+1,row+nspcc+1) = (...
                            + poro(iz)*(1d0)/dt ...
                            - (0.5d0*(poro(iz+1)*dif_alk(iz+1)+poro(iz)*dif_alk(iz))*(-1d0)/(0.5d0*(dz(iz+1)+dz(iz))) ...
                            - 0.5d0*(poro(iz)*dif_alk(iz)+poro(iz-1)*dif_alk(iz-1))*(1d0)/(0.5d0*(dz(iz)+dz(iz-1))))/dz(iz)  ...
                            - 2d0*sporo(iz)*sum(drcc_dalk(iz,:))  ...
                            )*alkx(iz)*fact;
                        amx(row+nspcc+1,row+nspcc+1+nsp) = ( ...
                            - (0.5d0*(poro(iz+1)*dif_alk(iz+1)+poro(iz)*dif_alk(iz))*(1d0)/(0.5d0*(dz(iz+1)+dz(iz))) ...
                            - 0d0)/dz(iz)  ...
                            )*alkx(iz+1)*fact;
                        amx(row+nspcc+1,row+nspcc+1-nsp) = (...
                            - (0d0 ...
                            - 0.5d0*(poro(iz)*dif_alk(iz)+poro(iz-1)*dif_alk(iz-1))*(-1d0)/(0.5d0*(dz(iz)+dz(iz-1))))/dz(iz) ...
                            )*alkx(iz-1)*fact;
                        amx(row+nspcc+1,row+nspcc) = (...
                            - 2d0*sporo(iz)*sum(drcc_ddic(iz,:))  ...
                            )*dicx(iz)*fact;
                    end
                    % diffusion terms are filled with transition matrices
                    for isp=1:nspcc
                        if (turbo2(isp+2) || labs(isp+2))
                            for iiz = 1:nz
                                col = 1 + (iiz-1)*nsp;
                                if (trans(iiz,iz,isp+2)==0d0)
                                    continue        % cycle in fortran
                                end
                                amx(row+isp-1,col+isp-1) = amx(row+isp-1,col+isp-1) ...
                                    - trans(iiz,iz,isp+2)/dz(iz)*dz(iiz)*(1d0-poro(iiz))*ccx(iiz,isp);
                                ymx(row+isp-1) = ymx(row+isp-1) ...
                                    - trans(iiz,iz,isp+2)/dz(iz)*dz(iiz)*(1d0-poro(iiz))*ccx(iiz,isp);
                            end
                        else
                            for iiz = 1:nz
                                col = 1 + (iiz-1)*nsp;
                                if (trans(iiz,iz,isp+2)==0d0)
                                    continue        % cycle in fortran
                                end
                                amx(row+isp-1,col+isp-1) = amx(row+isp-1,col+isp-1) -trans(iiz,iz,isp+2)/dz(iz)*ccx(iiz,isp);   % here Hauptdiagonale
                                ymx(row+isp-1) = ymx(row+isp-1) - trans(iiz,iz,isp+2)/dz(iz)*ccx(iiz,isp);
                            end
                        end
                    end
                end % for loop
                
                ymx = - ymx';  % because I put f(x) into ymx (=B), minus sign need be added & transpose
                
                if(~def_nonrec)
                    if (any(isnan(ymx)))
                        
                        % fprintf('negative om, stop \n');
                        str = sprintf('chk_ymx_pre.txt');
                        file_tmp = fopen(str,'wt');
                        for iz=1:nmx
                            %        write (file_tmp,*) ymx(iz)      write(file_tmp,*)(trans(iz,iiz,isp),iiz=1,nz)
                            fprintf(file_tmp,'%17.16e', ymx(iz));
                            fprintf(file_tmp,'\n');
                        end
                        fclose(file_tmp);
                        msg = 'Error: NAN in ymx, STOP.';
                        error(msg)
                        %                         print*,'NAN in ymx'
                        %                         open(unit=file_tmp,file=trim(adjustl(workdir))//'chk_ymx_pre.txt',status = 'unknown')
                        %                         do iz = 1, nmx
                        %                         write (file_tmp,*) ymx(iz)
                        %                         enddo
                        %                         close(file_tmp)
                        %                         stop
                    end
                end
                
                if(~def_sparse)
                    % using non-sparse solver
                    %                     fprintf('ymx before \n');
                    %                  	for iz=1:nmx
                    %                         fprintf('%17.16e \n', ymx(iz));
                    %                     end
                    %                 fprintf('end amx(1,1) %17.16e \n', amx(1,1));
                    %                 str = sprintf('matlab_chk_amx.txt');
                    %                 file_tmp = fopen(str,'wt');
                    %
                    %                 for iz=1:nmx
                    %                     %        write(file_tmp,*)(trans(iz,iiz,isp),iiz=1,nz)
                    %                     fprintf(file_tmp,'%17.16e\t',amx(iz,:));
                    %                     fprintf(file_tmp,'\n');
                    %                 end
                    %                 fclose(file_tmp);
                    
                    %                    call dgesv(nmx,int(1),amx,nmx,ipiv,ymx,nmx,infobls)
                    ymx = matrixDivide(amx,ymx);
%                                         fprintf('ymx after \n');
%                                          for iz=1:nmx
%                                              fprintf('%17.16e \n', ymx(iz));
%                                          end
                    
                else
                    %  slowest way of using sparse matrix solver
                    
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    %%%%%% Dominik TODO: implement solution with UMFPACK %%%%%%%%
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                end
                
                for iz = 1:nz
                    row = 1+(iz-1)*nsp;
                    for isp=1:nspcc
                        if (ymx(row+isp-1)>10d0) % this help conversion
                            ccx(iz,isp) = ccx(iz,isp)*1.5d0;
                        elseif (ymx(row+isp-1)<-10d0)  % this help conversion
                            ccx(iz,isp) = ccx(iz,isp)*0.5d0;
                        else
                            ccx(iz,isp) = ccx(iz,isp)*exp(ymx(row+isp-1));
                        end
                        if (ccx(iz,isp)<ccx_th)  % too small trancate value and not be accounted for error
                            ccx(iz,isp)=ccx_th;
                            ymx(row+isp-1) = 0d0;
                        end
                    end
                    if (ymx(row+nspcc)>10d0)
                        dicx(iz)=dicx(iz)*1.5d0;
                    elseif (ymx(row+nspcc)<-10d0)
                        dicx(iz)=dicx(iz)*0.5d0;
                    else
                        dicx(iz) = dicx(iz)*exp(ymx(row+nspcc));
                    end
                    if (ymx(row+nspcc+1)>10d0)
                        alkx(iz) = alkx(iz)*1.5d0;
                    elseif (ymx(row+nspcc+1)<-10d0)
                        alkx(iz) = alkx(iz)*0.5d0;
                    else
                        alkx(iz) = alkx(iz)*exp(ymx(row+nspcc+1));
                    end
                    if (dicx(iz)<1d-100)
                        ymx(row+nspcc) = 0d0;
                    end
                    if (alkx(iz)<1d-100)
                        ymx(row+nspcc+1) = 0d0;
                    end
                end
                
                %              	fprintf('ccx after \n');
                %                  ccx
                
                error = max(exp(abs(ymx))) - 1d0;
                itr = itr + 1;
                if(def_showiter)
                    fprintf('co2 iteration, error %i \t %17.16e \n',itr, error);
                end
                
                
                %  if negative or NAN calculation stops
                if (any(ccx<0d0))
                    ccx
                    msg = 'negative ccx, stop';
                    error(msg)
                    %     print*,'negative ccx, stop'
                    %     print*,ccx
                    %     stop
                end
                if (any(isnan(ccx)))
                    ccx
                    msg = 'nan ccx, stop';
                    error(msg)
                end
                
                if (any(dicx<0d0))
                    dicx
                    msg = 'negative dicx, stop';
                    error(msg)
                end
                if (any(isnan(dicx)))
                    dicx
                    msg = 'nan dicx, stop';
                    error(msg)
                end
                
                if (any(alkx<0d0))
                    alkx
                    msg = 'negative alkx, stop';
                    error(msg)
                end
                if (any(isnan(alkx)))
                    alkx
                    msg = 'nan alkx, stop';
                    error(msg)
                end
                
            end     % while(error > tol)
            
        end     % function calccaco3sys
        
        function[cctflx,ccdis,ccdif,ccadv,ccrain,ccres,alktflx,alkdis,alkdif,alkdec,alkres, dictflx,dicdis,dicdif,dicres,dicdec, dw] = ...
                calcflxcaco3sys(dw, nspcc,ccx,cc, ccflx, dt,dz,rcc,adf,up,dwn,cnr,w,dif_alk,dif_dic,dic,dicx,alk,alkx,oxco2,anco2,trans, turbo2,labs,nonlocal,sporof,it,nz,poro,sporo, ...
                dici,alki,mvcc,tol)
            %% calculation of fluxes relevant to caco3 and co2 system
            
            nsp = nspcc+2;          % independent chemical variables
            
            % allocate output matrices / variables
            cctflx = zeros(1, nspcc);
            ccdis = zeros(1, nspcc);
            ccdif = zeros(1, nspcc);
            ccadv = zeros(1, nspcc);
            ccrain = zeros(1, nspcc);
            ccres = zeros(1, nspcc);
            
            dictflx = 0d0;
            dicdis = 0d0;
            dicdif = 0d0;
            dicdec = 0d0;
            dicres = 0d0;
            
            alktflx = 0d0;
            alkdis = 0d0;
            alkdif = 0d0;
            alkdec = 0d0;
            alkres = 0d0;
            
            for iz = 1:nz
                row = 1 + (iz-1)*nsp;
                if (iz == 1)
                    for isp=1:nspcc
                        cctflx(isp) = cctflx(isp) + (1d0-poro(iz))*(ccx(iz,isp)-cc(iz,isp))/dt *dz(iz);
                        ccdis(isp) = ccdis(isp)  + (1d0-poro(iz))*rcc(iz,isp) *dz(iz);
                        ccrain(isp) = ccrain(isp) - ccflx(isp)/dz(1)*dz(iz);
                        ccadv(isp) = ccadv(isp) + adf(iz)*up(iz)*(sporo(iz)*w(iz)*ccx(iz,isp)-0d0)/dz(1) * dz(iz) ...
                            + adf(iz)*dwn(iz)*(sporo(iz+1)*w(iz+1)*ccx(iz+1,isp)-sporo(iz)*w(iz)*ccx(iz,isp))/dz(1) * dz(iz)  ...
                            + adf(iz)*cnr(iz)*(sporo(iz+1)*w(iz+1)*ccx(iz+1,isp)-0d0)/dz(1) * dz(iz);
                    end
                    %  DIC
                    dictflx = dictflx +(dicx(iz)-dic(iz))/dt*dz(iz)*poro(iz);
                    dicdif = dicdif - ((poro(iz)*dif_dic(iz)+poro(iz+1)*dif_dic(iz+1))*0.5d0*(dicx(iz+1)-dicx(iz))/(0.5d0*(dz(iz)+dz(iz+1))) ...
                        - poro(iz)*dif_dic(iz)*(dicx(iz)-dici*1d-6/1d3)/dz(iz))/dz(iz)*dz(iz);
                    dicdec = dicdec - oxco2(iz)*dz(iz) - anco2(iz)*dz(iz);
                    dicdis = dicdis - sum(rcc(iz,:))*sporo(iz)*dz(iz) ;
                    % ALK
                    alktflx = alktflx + (alkx(iz)-alk(iz))/dt*dz(iz)*poro(iz);
                    alkdif = alkdif - ((poro(iz)*dif_alk(iz)+poro(iz+1)*dif_alk(iz+1))*0.5d0*(alkx(iz+1)-alkx(iz))/(0.5d0*(dz(iz)+dz(iz+1))) ...
                        - poro(iz)*dif_alk(iz)*(alkx(iz)-alki*1d-6/1d3)/dz(iz))/dz(iz)*dz(iz);
                    alkdec = alkdec - anco2(iz)*dz(iz) ;
                    alkdis = alkdis - 2d0* sporo(iz)*sum(rcc(iz,:))*dz(iz) ;
                elseif (iz == nz)
                    for isp=1:nspcc
                        cctflx(isp) = cctflx(isp) + sporo(iz)*(ccx(iz,isp)-cc(iz,isp))/dt *dz(iz);
                        ccdis(isp) = ccdis(isp)  + sporo(iz)*rcc(iz,isp) *dz(iz);
                        ccadv(isp) = ccadv(isp) ...
                            + adf(iz)*up(iz)*(sporo(iz)*w(iz)*ccx(iz,isp)-sporo(iz-1)*w(iz-1)*ccx(iz-1,isp))/dz(iz) * dz(iz)  ...
                            + adf(iz)*cnr(iz)*(sporof*w(iz)*ccx(iz,isp)-sporo(iz-1)*w(iz-1)*ccx(iz-1,isp))/dz(iz) * dz(iz)  ...
                            + adf(iz)*dwn(iz)*(sporof*w(iz)*ccx(iz,isp)-sporo(iz)*w(iz)*ccx(iz,isp))/dz(iz) * dz(iz)  ;
                    end
                    % DIC
                    dictflx = dictflx +(dicx(iz)-dic(iz))/dt*dz(iz)*poro(iz) ;
                    dicdif = dicdif - (0d0 ...
                        - 0.5d0*(poro(iz)*dif_dic(iz)+poro(iz-1)*dif_dic(iz-1))*(dicx(iz)-dicx(iz-1))/(0.5d0*(dz(iz-1)+dz(iz))) ...
                        )/dz(iz)*dz(iz);
                    dicdec = dicdec - oxco2(iz)*dz(iz) - anco2(iz)*dz(iz) ;
                    dicdis = dicdis - sporo(iz)*sum(rcc(iz,:))*dz(iz) ;
                    % ALK
                    alktflx = alktflx + (alkx(iz)-alk(iz))/dt*dz(iz)*poro(iz);
                    alkdif = alkdif - (0d0 ...
                        - 0.5d0*(poro(iz)*dif_alk(iz)+poro(iz-1)*dif_alk(iz-1))*(alkx(iz)-alkx(iz-1))/(0.5d0*(dz(iz-1)+dz(iz))))/dz(iz)*dz(iz);
                    alkdec = alkdec - anco2(iz)*dz(iz);
                    alkdis = alkdis - 2d0* sporo(iz)*sum(rcc(iz,:))*dz(iz);
                else
                    for isp=1:nspcc
                        cctflx(isp) = cctflx(isp) + sporo(iz)*(ccx(iz,isp)-cc(iz,isp))/dt *dz(iz);
                        ccdis(isp) = ccdis(isp)  + sporo(iz)*rcc(iz,isp) *dz(iz);
                        ccadv(isp) = ccadv(isp) ...
                            + adf(iz)*up(iz)*(sporo(iz)*w(iz)*ccx(iz,isp)-sporo(iz-1)*w(iz-1)*ccx(iz-1,isp))/dz(iz) * dz(iz)  ...
                            + adf(iz)*dwn(iz)*(sporo(iz+1)*w(iz+1)*ccx(iz+1,isp)-sporo(iz)*w(iz)*ccx(iz,isp))/dz(iz) * dz(iz)  ...
                            + adf(iz)*cnr(iz)*(sporo(iz+1)*w(iz+1)*ccx(iz+1,isp)-sporo(iz-1)*w(iz-1)*ccx(iz-1,isp))/dz(iz) * dz(iz)  ;
                    end
                    % DIC
                    dictflx = dictflx +(dicx(iz)-dic(iz))/dt*dz(iz)*poro(iz) ;
                    dicdif = dicdif - (0.5d0*(poro(iz+1)*dif_dic(iz+1)+poro(iz)*dif_dic(iz))*(dicx(iz+1)-dicx(iz))/(0.5d0*(dz(iz+1)+dz(iz))) ...
                        - 0.5d0*(poro(iz)*dif_dic(iz)+poro(iz-1)*dif_dic(iz-1))*(dicx(iz)-dicx(iz-1))/(0.5d0*(dz(iz)+dz(iz-1))) ...
                        )/dz(iz)*dz(iz);
                    dicdec = dicdec - oxco2(iz)*dz(iz) - anco2(iz)*dz(iz) ;
                    dicdis = dicdis - sporo(iz)*sum(rcc(iz,:))*dz(iz) ;
                    % ALK
                    alktflx = alktflx + (alkx(iz)-alk(iz))/dt*dz(iz)*poro(iz);
                    alkdif = alkdif - (0.5d0*(poro(iz+1)*dif_alk(iz+1)+poro(iz)*dif_alk(iz))*(alkx(iz+1)-alkx(iz))/(0.5d0*(dz(iz+1)+dz(iz))) ...
                        - 0.5d0*(poro(iz)*dif_alk(iz)+poro(iz-1)*dif_alk(iz-1))*(alkx(iz)-alkx(iz-1))/(0.5d0*(dz(iz)+dz(iz-1))))/dz(iz)*dz(iz);
                    alkdec = alkdec - anco2(iz)*dz(iz);
                    alkdis = alkdis - 2d0* sporo(iz)*sum(rcc(iz,:))*dz(iz) ;
                end %if
                for isp=1:nspcc
                    if (labs(isp+2) || turbo2(isp+2))
                        for iiz = 1:nz
                            if (trans(iiz,iz,isp+2)==0d0)
                                continue        % cycle in fortran
                            end
                            ccdif(isp) = ccdif(isp) -trans(iiz,iz,isp+2)/dz(iz)*dz(iiz)*(1d0-poro(iiz))*dz(iz)*ccx(iiz,isp);
                        end
                    else
                        for iiz = 1:nz
                            if (trans(iiz,iz,isp+2)==0d0)
                                continue        % cycle in fortran
                            end
                            ccdif(isp) = ccdif(isp) -trans(iiz,iz,isp+2)/dz(iz)*dz(iz)*ccx(iiz,isp);
                        end
                    end %if
                    if (labs(isp+2) || turbo2(isp+2))
                        for iiz = 1:nz
                            if (trans(iiz,iz,isp+2)==0d0)
                                continue        % cycle in fortran
                            end
                            dw(iz) = dw(iz) - mvcc*(-trans(iiz,iz,isp+2)/dz(iz)*dz(iiz)*(1d0-poro(iiz))*ccx(iiz,isp));
                        end
                    else
                        if (nonlocal(isp+2))
                            for iiz = 1:nz
                                if (trans(iiz,iz,isp+2)==0d0)
                                    continue        % cycle in fortran
                                end
                                dw(iz) = dw(iz) -mvcc*(-trans(iiz,iz,isp+2)/dz(iz)*ccx(iiz,isp));
                            end
                        end %if
                    end %if
                end %for
                dw(iz) = dw(iz) -(1d0-poro(iz))*mvcc*sum(rcc(iz,:));
            end  %% end first for
            
            % residual fluxes
            ccres = cctflx +  ccdis +  ccdif + ccadv + ccrain;
            dicres = dictflx + dicdis + dicdif + dicdec ;
            alkres = alktflx + alkdis + alkdif + alkdec ;
            
            % if (abs(alkres)/max(alktflx,alkdis ,alkdif , alkdec) > tol*10d0) then   % if residual fluxes are relatively large, record just in case
            % print*,'not enough accuracy in co2 calc:stop',abs(alkres)/max(alktflx,alkdis ,alkdif , alkdec)
            % write(file_err,*)it,'not enough accuracy in co2 calc:stop',abs(alkres)/max(alktflx,alkdis ,alkdif , alkdec)  &
            % ,alkres, alktflx,alkdis , alkdif , alkdec
            
            if (abs(alkres)/max(abs(ccflx)) > tol*10d0)    % if residual fluxes are relatively large, record just in case
                fprintf('not enough accuracy in co2 calc:stop: %17.16e \n',abs(alkres)/max(abs(ccflx)));
                fprintf('not enough accuracy in co2 calc:stop: alkres, alktflx, alkdis , alkdif , alkdec %17.16e \t %17.16e \t %17.16e \t %17.16e \t %17.16e \t \n',alkres, alktflx, alkdis , alkdif , alkdec);
                %                 write(file_err,*)it,'not enough accuracy in co2 calc:stop',abs(alkres)/maxval(abs(ccflx))  &
                %                 ,alkres, alktflx,alkdis , alkdif , alkdec
            end %if
            
            
        end
        
        function ptx = claycalc(nz,sporo,pt,dt,w,dz,detflx,adf,up,dwn,cnr,trans, nspcc,labs,turbo2,nonlocal,poro,sporof, msed, def_nonrec)
            
            nsp = 1;        %  independent chemical variables, only consider clay
            nmx = nz*nsp;   % matrix is linear and solved like om and o2, so see comments there for calculation procedures
            
            % allocate output vector
            ptx = zeros(1, nz);     % mol cm-3 sld; clay conc.
            
            %  allocate local matrices used to solve linear system Ax = B
            amx = zeros(nmx, nmx);          % amx corresponds to A in Ax = B, but ymx is also x when Ax = B is solved. emx is array of error
            ymx = zeros(1, nmx);            % ymx correspond to B in Ax = B, but ymx is also x when Ax = B is solved
            %            emx = zeros(1, nmx);            % emx is array of error
            %            ipiv = zeros(1, nmx);           % matrix used to solve linear system Ax = B
            dumx = zeros(nmx, nmx);          % amx corresponds to A in Ax = B, but ymx is also x when Ax = B is solved. emx is array of error
            
            error = 1d4;
            itr = 0;
            
            
            for iz = 1:nz
                row = 1 + (iz-1)*nsp;
                if (iz == 1)
                    ymx(row) = ...
                        + sporo(iz)*(-pt(iz))/dt ...
                        - detflx/msed/dz(iz);
                    amx(row,row) = (...
                        + sporo(iz)*(1d0)/dt ...
                        + adf(iz)*up(iz)*(sporo(iz)*w(iz)*1d0-0d0)/dz(iz)   ...
                        + adf(iz)*dwn(iz)*(sporo(iz+1)*w(iz+1)*0d0-sporo(iz)*w(iz)*1d0)/dz(iz)   ...
                        + adf(iz)*cnr(iz)*(sporo(iz+1)*w(iz+1)*0d0-0d0)/dz(iz)   ...
                        );
                    amx(row,row+nsp) =  (...
                        + adf(iz)*dwn(iz)*(sporo(iz+1)*w(iz+1)*1d0-sporo(iz)*w(iz)*0d0)/dz(iz)   ...
                        + adf(iz)*cnr(iz)*(sporo(iz+1)*w(iz+1)*1d0-0d0)/dz(iz)   ...
                        );
                elseif (iz == nz)
                    ymx(row) = ...
                        + sporo(iz)*(-pt(iz))/dt;
                    amx(row,row) = (...
                        + sporo(iz)*(1d0)/dt ...
                        + adf(iz)*up(iz)*(sporo(iz)*w(iz)*1d0-0d0)/dz(iz)  ...
                        + adf(iz)*cnr(iz)*(sporof*w(iz)*1d0-0d0)/dz(iz)  ...
                        + adf(iz)*dwn(iz)*(sporof*w(iz)*1d0-sporo(iz)*w(iz)*1d0)/dz(iz)  ...
                        );
                    amx(row,row-nsp) = ( ...
                        + adf(iz)*up(iz)*(0d0-sporo(iz-1)*w(iz-1)*1d0)/dz(iz)  ...
                        + adf(iz)*cnr(iz)*(0d0-sporo(iz-1)*w(iz-1)*1d0)/dz(iz)  ...
                        );
                else
                    ymx(row) = ...
                        + sporo(iz)*(-pt(iz))/dt;
                    amx(row,row) = (...
                        + sporo(iz)*(1d0)/dt ...
                        + adf(iz)*up(iz)*(sporo(iz)*w(iz)*1d0-0d0)/dz(iz)  ...
                        + adf(iz)*dwn(iz)*(sporo(iz+1)*w(iz+1)*0d0-sporo(iz)*w(iz)*1d0)/dz(iz)  ...
                        + adf(iz)*cnr(iz)*(sporo(iz+1)*w(iz+1)*0d0-0d0)/dz(iz)  ...
                        );
                    amx(row,row+nsp) =  (...
                        + adf(iz)*dwn(iz)*(sporo(iz+1)*w(iz+1)*1d0-sporo(iz)*w(iz)*0d0)/dz(iz)  ...
                        + adf(iz)*cnr(iz)*(sporo(iz+1)*w(iz+1)*1d0-0d0)/dz(iz)  ...
                        );
                    amx(row,row-nsp) =  (...
                        + adf(iz)*up(iz)*(0d0-sporo(iz-1)*w(iz-1)*1d0)/dz(iz)  ...
                        + adf(iz)*cnr(iz)*(0d0-sporo(iz-1)*w(iz-1)*1d0)/dz(iz)  ...
                        );
                end
                if (labs(2) || turbo2(2))
                    for iiz = 1:nz
                        col = 1 + (iiz-1)*nsp;
                        if (trans(iiz,iz,2)==0d0)
                            continue        % cycle in fortran
                        end
                        amx(row,col) = amx(row,col) -trans(iiz,iz,2)/dz(iz)*dz(iiz)*(1d0-poro(iiz));
                    end
                else
                    for iiz = 1:nz
                        col = 1 + (iiz-1)*nsp;
                        if (trans(iiz,iz,2)==0d0)
                            continue        % cycle in fortran
                        end
                        amx(row,col) = amx(row,col) -trans(iiz,iz,2)/dz(iz);
                    end
                end
            end     % first fortran do
            
            ymx = - ymx';
            
            if(~def_nonrec)
                if (any(isnan(ymx)))
                    str = sprintf('chk_ymx_pre_pt.txt');
                    file_tmp = fopen(str,'wt');
                    for iz=1:nmx
                        %        write (file_tmp,*) ymx(iz)      write(file_tmp,*)(trans(iz,iiz,isp),iiz=1,nz)
                        fprintf(file_tmp,'%17.16e', ymx(iz));
                        fprintf(file_tmp,'\n');
                    end
                    fclose(file_tmp);
                    msg = 'NAN in pre ymx:pt, stop';
                    error(msg)
                end
            end
            
            %                call dgesv(nmx,int(1),amx,nmx,ipiv,ymx,nmx,infobls)
            ymx = matrixDivide(amx,ymx);
            
            if(~def_nonrec)
                if (any(isnan(amx)))
                    str = sprintf('chk_amx_pt.txt');
                    file_tmp = fopen(str,'wt');
                    for iz=1:nmx
                        fprintf(file_tmp,'%17.16e', amx(iz));
                        fprintf(file_tmp,'\n');
                    end
                    fclose(file_tmp);
                    msg = 'NAN in amx:pt, stop';
                    error(msg)
                end
                
                if (any(isnan(ymx)))
                    str = sprintf('chk_ymx_pt.txt');
                    file_tmp = fopen(str,'wt');
                    for iz=1:nmx
                        fprintf(file_tmp,'%17.16e', ymx(iz));
                        fprintf(file_tmp,'\n');
                    end
                    fclose(file_tmp);
                    msg = 'NAN in ymx:pt, stop';
                    error(msg)
                end
            end
            
            ptx = ymx;    	% mol cm-3 sld; clay conc.
            
        end
        
        
        
        function [pttflx,ptdif,ptadv,ptres,ptrain, dw] = calcflxclay(dw, nz,sporo,ptx,pt,dt,dz,detflx,w,adf,up,dwn,cnr,sporof,trans,turbo2,labs,nonlocal,poro, msed,mvsed)
            
            % allocate output variables (clay fluxes )
            pttflx = 0d0;
            ptdif = 0d0;
            ptadv = 0d0;
            ptres = 0d0;
            ptrain = 0d0;
            
            nsp = 1;        %  independent chemical variables, only consider clay
            
            for iz = 1:nz
                row = 1 + (iz-1)*nsp;
                if (iz == 1)
                    pttflx = pttflx + sporo(iz)*(ptx(iz)-pt(iz))/dt*dz(iz);
                    ptrain = ptrain - detflx/msed;
                    ptadv = ptadv ...
                        + adf(iz)*up(iz)*(sporo(iz)*w(iz)*ptx(iz)-0d0)/dz(iz)*dz(iz)  ...
                        + adf(iz)*dwn(iz)*(sporo(iz+1)*w(iz+1)*ptx(iz+1)-sporo(iz)*w(iz)*ptx(iz))/dz(iz)*dz(iz)  ...
                        + adf(iz)*cnr(iz)*(sporo(iz+1)*w(iz+1)*ptx(iz+1)-0d0)/dz(iz)*dz(iz)  ;
                elseif (iz == nz)
                    pttflx = pttflx + (1d0-poro(iz))*(ptx(iz)-pt(iz))/dt*dz(iz);
                    ptadv = ptadv ...
                        + adf(iz)*up(iz)*(sporo(iz)*w(iz)*ptx(iz)-sporo(iz-1)*w(iz-1)*ptx(iz-1))/dz(iz)*dz(iz)  ...
                        + adf(iz)*cnr(iz)*(sporof*w(iz)*ptx(iz)-sporo(iz-1)*w(iz-1)*ptx(iz-1))/dz(iz)*dz(iz)  ...
                        + adf(iz)*dwn(iz)*(sporof*w(iz)*ptx(iz)-sporo(iz)*w(iz)*ptx(iz))/dz(iz)*dz(iz)  ;
                else
                    pttflx = pttflx + (1d0-poro(iz))*(ptx(iz)-pt(iz))/dt*dz(iz);
                    ptadv = ptadv ...
                        + adf(iz)*up(iz)*(sporo(iz)*w(iz)*ptx(iz)-sporo(iz-1)*w(iz-1)*ptx(iz-1))/dz(iz)*dz(iz)  ...
                        + adf(iz)*dwn(iz)*(sporo(iz+1)*w(iz+1)*ptx(iz+1)-sporo(iz)*w(iz)*ptx(iz))/dz(iz)*dz(iz)  ...
                        + adf(iz)*cnr(iz)*(sporo(iz+1)*w(iz+1)*ptx(iz+1)-sporo(iz-1)*w(iz-1)*ptx(iz-1))/dz(iz)*dz(iz);
                end
                if(turbo2(2) || labs(2))
                    for iiz = 1:nz
                        if (trans(iiz,iz,2)==0d0)
                            continue        % cycle in fortran
                        end
                        ptdif = ptdif -trans(iiz,iz,2)*ptx(iiz)/dz(iz)*dz(iiz)*dz(iz);
                    end
                else
                    for iiz = 1:nz
                        if (trans(iiz,iz,2)==0d0)
                            continue        % cycle in fortran
                        end
                        ptdif = ptdif -trans(iiz,iz,2)*ptx(iiz)/dz(iz)    ...
                            *dz(iz);
                    end
                end
                if(turbo2(2) || labs(2))
                    for iiz = 1:nz
                        if (trans(iiz,iz,2)==0d0)
                            continue        % cycle in fortran
                        end
                        dw(iz) = dw(iz) - mvsed*(-trans(iiz,iz,2)*ptx(iiz)/dz(iz)*dz(iiz)*(1d0-poro(iiz)));
                    end
                else
                    if (nonlocal(2))
                        for iiz = 1:nz
                            if (trans(iiz,iz,2)==0d0)
                                continue        % cycle in fortran
                            end
                            dw(iz) = dw(iz) - mvsed*(-trans(iiz,iz,2)*ptx(iiz)/dz(iz));
                        end
                    end
                end
            end     % first fortran do
            
            ptres = pttflx + ptdif + ptadv + ptrain;
        end
        
        function [rho, frt] = getsldprop(nz, omx, ptx, ccx, nspcc, w, up, dwn, cnr, adf, z, mom, msed, mcc, mvom, mvsed, mvcc)
            %% calculate solid property, rho (density) and frt (total vol.frac)
            
            for iz=1:nz
                rho(iz) = omx(iz)*mom + ptx(iz)*msed +  sum(ccx(iz,:))*mcc;  % calculating bulk density
                frt(iz) = omx(iz)*mvom + ptx(iz)*mvsed + sum(ccx(iz,:))*mvcc;  % calculation of total vol. fraction of solids
            end
            
            % check error for density (rho)
            if (any(rho<0d0))  % if negative density stop ....
                str = sprintf('NEGATIVE_RHO.txt');
                file_tmp = fopen(str,'wt');
                for iz=1:nz
                    %        write (file_tmp,*) ymx(iz)      write(file_tmp,*)(trans(iz,iiz,isp),iiz=1,nz)
                    fprintf(file_tmp,'%17.16e \t %17.16e \t %17.16e \t %17.16e \t %17.16e \t %17.16e \t %17.16e \n ', z(iz),rho(iz),w(iz),up(iz),dwn(iz),cnr(iz),adf(iz));
                    %                    fprintf(file_tmp,'\n');
                end
                fclose(file_tmp);
                msg = 'negative density, STOP.';
                error(msg)
            end
        end
        
        function [w,wi] = burial_calc(detflx,ccflx,nspcc,omflx,dw,dz,poro,nz, msed,mvsed,mvcc,mvom,poroi)
            %% calculate new burial velocity
            
            % w is up dated by solving
            %           d(1 - poro)*w/dz = dw
            % note that dw has recorded volume changes by reactions and non-local mixing (see Eqs. B2 and B6 in ms)
            % finite difference form is
            %           if (iz/=1) {(1-poro(iz))*w(iz)-(1-poro(iz-1))*w(iz-1)}/dz(iz) = dw(iz)
            %           if (iz==1) (1-poro(iz))*w(iz) = total volume flux + dw(iz)*dz(iz)
            % which leads to the following calculations
            
            wi = (detflx/msed*mvsed + sum(ccflx)*mvcc +omflx*mvom)/(1d0-poroi);  % upper value; (1d0-poroi) is almost meaningless, see below
            for iz=1:nz
                if (iz==1)
                    w(iz)=((1d0-poroi)*wi + dw(iz)*dz(iz))/(1d0-poro(iz));
                else
                    w(iz)=((1d0-poro(iz-1))*w(iz-1) + dw(iz)*dz(iz))/(1d0-poro(iz));
                end
            end
            
            
        end
        
        function recordprofile(itrec, nz,z,age,pt,msed,wi,rho,cc,ccx,dic,dicx,alk,alkx,co3,co3x,co3sat ...
                ,rcc,pro,o2x,oxco2,anco2,om,mom,mcc,d13c_ocni,d18o_ocni,up,dwn,cnr,adf,nspcc,ptx,w,frt,prox,omx,d13c_blk,d18o_blk)
            %% write sediment profiles in files
            
            if (itrec==0)
                %    open(unit=file_tmp,file=trim(adjustl(workdir))//'ptx-'//trim(adjustl(dumchr(1)))//'.txt',action='write',status='replace')
                str = sprintf('./01_recprofile_fickian_nosig_default/matlab_ptx-%3.3i.txt',itrec);
                file_tmp = fopen(str,'wt');
                for iz = 1:nz
                    %        write(file_tmp,*) z(iz),age(iz),pt(iz)*msed/2.5d0*100,0d0,1d0  ,wi
                    fprintf(file_tmp,'%17.16e \t %17.16e \t %17.16e \t %17.16e \t %17.16e \t %17.16e \n ',z(iz),age(iz),pt(iz)*msed/2.5d0*100,0d0,1d0  ,wi);
                end
                fclose(file_tmp);
                
                %                open(unit=file_tmp,file=trim(adjustl(workdir))//'ccx-'//trim(adjustl(dumchr(1)))//'.txt',action='write',status='replace')
                str = sprintf('./01_recprofile_fickian_nosig_default/matlab_ccx-%3.3i.txt',itrec);
                file_tmp = fopen(str,'wt');
                for iz = 1:nz
                    %                 write(file_tmp,*) z(iz),age(iz),sum(cc(iz,:))*100d0/2.5d0*100d0, dic(iz)*1d3, alk(iz)*1d3, co3(iz)*1d3-co3sat &
                    %                 , sum(rcc(iz,:)),-log10(pro(iz))
                    fprintf(file_tmp,'%17.16e \t %17.16e \t %17.16e \t %17.16e \t %17.16e \t %17.16e \t %17.16e \t %17.16e \n ',...
                        z(iz),age(iz),sum(cc(iz,:))*100d0/2.5d0*100d0, dic(iz)*1d3, alk(iz)*1d3, co3(iz)*1d3-co3sat, sum(rcc(iz,:)),-log10(pro(iz) ));
                end
                fclose(file_tmp)
                
                %                open(unit=file_tmp,file=trim(adjustl(workdir))//'o2x-'//trim(adjustl(dumchr(1)))//'.txt',action='write',status='replace')
                str = sprintf('./01_recprofile_fickian_nosig_default/matlab_o2x-%3.3i.txt',itrec);
                file_tmp = fopen(str,'wt');
                for iz = 1:nz
                    %                 write(file_tmp,*) z(iz),age(iz),o2x(iz)*1d3, oxco2(iz), anco2(iz)
                    fprintf(file_tmp,'%17.16e \t %17.16e \t %17.16e \t %17.16e \t %17.16e \n ', z(iz),age(iz),o2x(iz)*1d3, oxco2(iz), anco2(iz));
                end
                fclose(file_tmp)
                
                %                open(unit=file_tmp,file=trim(adjustl(workdir))//'omx-'//trim(adjustl(dumchr(1)))//'.txt',action='write',status='replace')
                str = sprintf('./01_recprofile_fickian_nosig_default/matlab_omx-%3.3i.txt',itrec);
                file_tmp = fopen(str,'wt');
                for iz = 1:nz
                    %                 write(file_tmp,*) z(iz),age(iz),om(iz)*mom/2.5d0*100d0
                    fprintf(file_tmp,'%17.16e \t %17.16e \t %17.16e \n ', z(iz),age(iz),om(iz)*mom/2.5d0*100d0);
                end
                fclose(file_tmp)
                
                %                open(unit=file_tmp,file=trim(adjustl(workdir))//'ccx_sp-'//trim(adjustl(dumchr(1)))//'.txt' ,action='write',status='replace')
                str = sprintf('./01_recprofile_fickian_nosig_default/matlab_ccx_sp-%3.3i.txt',itrec);
                fmt=[repmat('%17.16e \t',1,nspcc+2) '\n'];
                file_tmp = fopen(str,'wt');
                for iz = 1:nz
                    %                write(file_tmp,*) z(iz),age(iz),(cc(iz,isp)*mcc/2.5d0*100d0,isp=1,nspcc)
                    fprintf(file_tmp,fmt,z(iz), age(iz),cc(iz,:)*mcc/2.5d0*100d0);
                end
                fclose(file_tmp)
                
                %                open(unit=file_tmp,file=trim(adjustl(workdir))//'sig-'//trim(adjustl(dumchr(1)))//'.txt' ,action='write',status='replace')
                str = sprintf('./01_recprofile_fickian_nosig_default/matlab_sig-%3.3i.txt',itrec);
                file_tmp = fopen(str,'wt');
                for iz = 1:nz
                    %                write(file_tmp,*) z(iz),age(iz),d13c_ocni,d18o_ocni
                    fprintf(file_tmp,'%17.16e \t %17.16e \t %17.16e \t %17.16e \n ', z(iz),age(iz),d13c_ocni,d18o_ocni);
                end
                fclose(file_tmp)
                
                %                open(unit=file_tmp,file=trim(adjustl(workdir))//'bur-'//trim(adjustl(dumchr(1)))//'.txt' ,action='write',status='replace')
                str = sprintf('./01_recprofile_fickian_nosig_default/matlab_bur-%3.3i.txt',itrec);
                file_tmp = fopen(str,'wt');
                for iz = 1:nz
                    %                write(file_tmp,*) z(iz),age(iz),w(iz),up(iz),dwn(iz),cnr(iz),adf(iz)
                    fprintf(file_tmp,'%17.16e \t %17.16e \t %17.16e \t %17.16e \t %17.16e \t %17.16e \t %17.16e \n ', z(iz),age(iz),w(iz),up(iz),dwn(iz),cnr(iz),adf(iz));
                end
                fclose(file_tmp)
            else
                
                %    open(unit=file_tmp,file=trim(adjustl(workdir))//'ptx-'//trim(adjustl(dumchr(1)))//'.txt',action='write',status='replace')
                str = sprintf('./01_recprofile_fickian_nosig_default/matlab_ptx-%3.3i.txt',itrec);
                file_tmp = fopen(str,'wt');
                for iz = 1:nz
                    %        write(file_tmp,*) z(iz),age(iz),pt(iz)*msed/2.5d0*100,0d0,1d0  ,wi
                    fprintf(file_tmp,'%17.16e \t %17.16e \t %17.16e \t %17.16e \t %17.16e \t %17.16e \n ',z(iz),age(iz),ptx(iz)*msed/rho(iz)*100d0,rho(iz),frt(iz) ,w(iz));
                end
                fclose(file_tmp);
                
                %                open(unit=file_tmp,file=trim(adjustl(workdir))//'ccx-'//trim(adjustl(dumchr(1)))//'.txt',action='write',status='replace')
                str = sprintf('./01_recprofile_fickian_nosig_default/matlab_ccx-%3.3i.txt',itrec);
                file_tmp = fopen(str,'wt');
                for iz = 1:nz
                    %                 write(file_tmp,*) z(iz),age(iz),sum(cc(iz,:))*100d0/2.5d0*100d0, dic(iz)*1d3, alk(iz)*1d3, co3(iz)*1d3-co3sat &
                    %                 , sum(rcc(iz,:)),-log10(pro(iz))
                    fprintf(file_tmp,'%17.16e \t %17.16e \t %17.16e \t %17.16e \t %17.16e \t %17.16e \t %17.16e \t %17.16e \n ',...
                        z(iz),age(iz),sum(ccx(iz,:))*mcc/rho(iz)*100d0, dicx(iz)*1d3, alkx(iz)*1d3, co3x(iz)*1d3-co3sat, sum(rcc(iz,:)),-log10(prox(iz)) );
                end
                fclose(file_tmp)
                
                %                open(unit=file_tmp,file=trim(adjustl(workdir))//'o2x-'//trim(adjustl(dumchr(1)))//'.txt',action='write',status='replace')
                str = sprintf('./01_recprofile_fickian_nosig_default/matlab_o2x-%3.3i.txt',itrec);
                file_tmp = fopen(str,'wt');
                for iz = 1:nz
                    %                 write(file_tmp,*) z(iz),age(iz),o2x(iz)*1d3, oxco2(iz), anco2(iz)
                    fprintf(file_tmp,'%17.16e \t %17.16e \t %17.16e \t %17.16e \t %17.16e \n ', z(iz),age(iz),o2x(iz)*1d3, oxco2(iz), anco2(iz));
                end
                fclose(file_tmp)
                
                %                open(unit=file_tmp,file=trim(adjustl(workdir))//'omx-'//trim(adjustl(dumchr(1)))//'.txt',action='write',status='replace')
                str = sprintf('./01_recprofile_fickian_nosig_default/matlab_omx-%3.3i.txt',itrec);
                file_tmp = fopen(str,'wt');
                for iz = 1:nz
                    %                 write(file_tmp,*) z(iz),age(iz),om(iz)*mom/2.5d0*100d0
                    fprintf(file_tmp,'%17.16e \t %17.16e \t %17.16e \n ', z(iz),age(iz),omx(iz)*mom/rho(iz)*100d0);
                end
                fclose(file_tmp)
                
                %                open(unit=file_tmp,file=trim(adjustl(workdir))//'ccx_sp-'//trim(adjustl(dumchr(1)))//'.txt' ,action='write',status='replace')
                str = sprintf('./01_recprofile_fickian_nosig_default/matlab_ccx_sp-%3.3i.txt',itrec);
                fmt=[repmat('%17.16e \t',1,nspcc+2) '\n'];
                file_tmp = fopen(str,'wt');
                for iz = 1:nz
                    %                write(file_tmp,*) z(iz),age(iz),(cc(iz,isp)*mcc/2.5d0*100d0,isp=1,nspcc)
                    fprintf(file_tmp,fmt,z(iz), age(iz),ccx(iz,:)*mcc/rho(iz)*100d0);
                end
                fclose(file_tmp)
                
                %                open(unit=file_tmp,file=trim(adjustl(workdir))//'sig-'//trim(adjustl(dumchr(1)))//'.txt' ,action='write',status='replace')
                str = sprintf('./01_recprofile_fickian_nosig_default/matlab_sig-%3.3i.txt',itrec);
                file_tmp = fopen(str,'wt');
                for iz = 1:nz
                    %                write(file_tmp,*) z(iz),age(iz),d13c_ocni,d18o_ocni
                    fprintf(file_tmp,'%17.16e \t %17.16e \t %17.16e \t %17.16e \n ', z(iz),age(iz),d13c_blk(iz),d18o_blk(iz));
                end
                fclose(file_tmp)
                
                %                open(unit=file_tmp,file=trim(adjustl(workdir))//'bur-'//trim(adjustl(dumchr(1)))//'.txt' ,action='write',status='replace')
                str = sprintf('./01_recprofile_fickian_nosig_default/matlab_bur-%3.3i.txt',itrec);
                file_tmp = fopen(str,'wt');
                for iz = 1:nz
                    %                write(file_tmp,*) z(iz),age(iz),w(iz),up(iz),dwn(iz),cnr(iz),adf(iz)
                    fprintf(file_tmp,'%17.16e \t %17.16e \t %17.16e \t %17.16e \t %17.16e \t %17.16e \t %17.16e \n ', z(iz),age(iz),w(iz),up(iz),dwn(iz),cnr(iz),adf(iz));
                end
                fclose(file_tmp)
                
            end
            
        end
        
        function [rectime, cntrec, time_spn, time_trs, time_aft] = recordtime(nrec,wi, ztot, def_biotest, def_sense, def_nonrec)
            %% set recording time
            
            rectime = zeros(1, nrec);
            
            %%% +++ tracing experiment
            time_spn = ztot / wi *50d0; % yr  % spin-up duration, 50 times the shortest residence time possible (assuming no caco3 dissolution)
            time_trs = 50d3; %  duration of signal change event
            time_aft = time_trs*3d0;  % duration of simulation after the event
            
            if(def_biotest)
                time_trs = 5d3;   %  smaller duration of event assumed
                time_aft = time_trs*10d0;
            end
            % distributing recording time in 3 different periods
            for itrec=1:nrec/3
                rectime(itrec)=itrec*time_spn/real(nrec/3);
            end
            for itrec=nrec/3+1:nrec/3*2
                rectime(itrec)=rectime(nrec/3)+(itrec-nrec/3)*time_trs/real(nrec/3);
            end
            for itrec=nrec/3*2+1:nrec
                rectime(itrec)=rectime(nrec/3*2)+(itrec-nrec/3*2)*time_aft/real(nrec/3);
            end
            if(def_sense)
                time_trs = 0d0;   %  there is no event needed
                time_aft = 0d0;
                % time_spn = 2d0*time_spn  % first run without dissolution
                for itrec=1:nrec
                    rectime(itrec)=itrec*time_spn/real(nrec);
                end
            end
            %%% ++++
            cntrec = 1;  % rec number (increasing with recording done )
            if(~def_nonrec)
                % open(unit=file_tmp,file=trim(adjustl(workdir))//'rectime.txt',action='write',status='unknown')
                str = sprintf('./01_recprofile_fickian_nosig_default/matlab_rectime.txt');
                file_tmp = fopen(str,'wt');
                for itrec=1:nrec
                    %write(file_tmp,*) rectime(itrec)  % recording when records are made
                    fprintf(file_tmp,'%17.16e \n ', rectime(itrec));
                end
                fclose(file_tmp);
            end
            
            
        end
        
        function [d13c_sp,d18o_sp] = sig2sp_pre(d13c_ocni,d13c_ocnf,d18o_ocni,d18o_ocnf, def_sense, def_size, nspcc)
            %% end-member signal assignment
            
            if(~def_sense) % without signal tracking
                %  four end-member caco3 species interpolation
                d13c_sp(1)=d13c_ocni;
                d18o_sp(1)=d18o_ocni;
                d13c_sp(2)=d13c_ocni;
                d18o_sp(2)=d18o_ocnf;
                d13c_sp(3)=d13c_ocnf;
                d18o_sp(3)=d18o_ocni;
                d13c_sp(4)=d13c_ocnf;
                d18o_sp(4)=d18o_ocnf;
            else
                d18o_sp = zeros(1, nspcc);
                d13c_sp = zeros(1, nspcc);
            end
            
            if(def_size)
                % assuming two differently sized caco3 species and giving additional 4 species their end-member isotopic compositions
                d13c_sp(5)=d13c_ocni;
                d18o_sp(5)=d18o_ocni;
                d13c_sp(6)=d13c_ocni;
                d18o_sp(6)=d18o_ocnf;
                d13c_sp(7)=d13c_ocnf;
                d18o_sp(7)=d18o_ocni;
                d13c_sp(8)=d13c_ocnf;
                d18o_sp(8)=d18o_ocnf;
            end
        end
        
        function [dt] = timestep(time,nt_spn,nt_trs,nt_aft, time_spn,time_trs, dt)
            %% set timestep: smaller when closer to & during event
            
            if (time <= time_spn)   % spin-up
                if (time+dt> time_spn)
                    dt=time_trs/nt_trs;     %/5000d0   % when close to 'event', time step needs to get smaller
                else
                    dt = time_spn/nt_spn;   % /800d0 % otherwise larger time step is better to fasten calculation
                end
            elseif (time>time_spn && time<=time_spn+time_trs)  % during event
                dt = time_trs/nt_trs;           % /5000d0
            elseif (time>time_spn+time_trs)
                dt=time_trs/nt_aft; %1000d0 % not too large time step
            end
        end
        
        function [d13c_ocn,d18o_ocn,ccflx,d18o_sp,d13c_sp] = signal_flx(time, time_spn,time_trs,d13c_ocni,d13c_ocnf,d18o_ocni,d18o_ocnf ...
                ,ccflx,ccflxi,d18o_sp,d13c_sp,it,nspcc, flxfini,flxfinf, def_track2, def_size, def_biotest)
            
            % local varaiable
            flxfin = 0.0;       % flux ratio of fine particles:
            cntsp = 0;          % counting caco3 species numbers when signal is tracked by method 2, signal is assigned species as time goes so counting the number of species whose signal has been assigned
            % flux fractions used for assigning flux values to realized isotope input changes
            flxfrc = zeros(1, nspcc);
            flxfrc2 = zeros(1, nspcc);
            
            if (time <= time_spn)   % spin-up
                d13c_ocn = d13c_ocni;  % take initial values
                d18o_ocn = d18o_ocni;
                ccflx = zeros(1, nspcc);
                ccflx(1) = ccflxi;  % raining only species with initial signal values
                if(def_track2)
                    cntsp=1;  % when signal is tracked by method 2, signal is assigned species as time goes so counting the number of species whose signal has been assigned
                    d18o_sp(cntsp)=d18o_ocn;
                    d13c_sp(cntsp)=d13c_ocn;
                    ccflx = zeros(1, nspcc);
                    ccflx(cntsp) = ccflxi;
                end
                if(def_size)
                    ccflx = zeros(1, nspcc);
                    flxfin = flxfini;   % flux ratio of fine particles: i and f denote initial and final values
                    ccflx(1) = (1d0-flxfin)*ccflxi; % fine species with initial signals
                    ccflx(5) = flxfin*ccflxi;       % coarse species with initial signals
                end
            elseif (time>time_spn && time<=time_spn+time_trs)   % during event
                d13c_ocn = d13c_ocni + (time-time_spn)*(d13c_ocnf-d13c_ocni)/time_trs; % single step
                d18o_ocn = d18o_ocni + (time-time_spn)*(d18o_ocnf-d18o_ocni)/time_trs;  % single step
                % creating spike (go from initial to final and come back to initial again taking half time of event duration )
                % this shape is assumed for d18o and fine (& coarse) caco3 flux changes
                if (time-time_spn<=time_trs/2d0)
                    d18o_ocn = d18o_ocni + (time-time_spn)*(d18o_ocnf-d18o_ocni)/time_trs*2d0;
                    flxfin = flxfini + (time-time_spn)*(flxfinf-flxfini)/time_trs*2d0;
                else
                    d18o_ocn = 2d0*d18o_ocnf - d18o_ocni - (time-time_spn)*(d18o_ocnf-d18o_ocni)/time_trs*2d0;
                    flxfin = 2d0*flxfinf - flxfini - (time-time_spn)*(flxfinf-flxfini)/time_trs*2d0;
                end
                if(~def_biotest)
                    % assuming a d13 excursion with shifts occurring with 1/10 times event duration
                    % the same assumption for depth change
                    if (time-time_spn<=time_trs/10d0)
                        d13c_ocn = d13c_ocni + (time-time_spn)*(d13c_ocnf-d13c_ocni)/time_trs*10d0;
                    elseif(time-time_spn>time_trs/10d0 && time-time_spn<=time_trs/10d0*9d0)
                        d13c_ocn = d13c_ocnf;
                    elseif  (time-time_spn>time_trs/10d0*9d0)
                        d13c_ocn = 10d0*d13c_ocnf - 9d0*d13c_ocni  - (time-time_spn)*(d13c_ocnf-d13c_ocni)/time_trs*10d0;
                    end
                end %#if
                if (~(d13c_ocn>=d13c_ocnf && d13c_ocn<=d13c_ocni))  % check if calculated d13c and d18o are within the assumed ranges
                    fprintf('error in d13c: %17.16e',d13c_ocn);
                    msg = 'error in d13c, stop';
                    error(msg)
                end
                % calculating fractions of flux assigned to individual caco3 species in case of interpolation of 2 signal inputs by 4 species
                % NEED to extend to allow tacking of any number of signals
                flxfrc(1) = abs(d13c_ocnf-d13c_ocn)/(abs(d13c_ocnf-d13c_ocn)+abs(d13c_ocni-d13c_ocn));
                flxfrc(2) = abs(d13c_ocni-d13c_ocn)/(abs(d13c_ocnf-d13c_ocn)+abs(d13c_ocni-d13c_ocn));
                flxfrc(3) = abs(d18o_ocnf-d18o_ocn)/(abs(d18o_ocnf-d18o_ocn)+abs(d18o_ocni-d18o_ocn));
                flxfrc(4) = abs(d18o_ocni-d18o_ocn)/(abs(d18o_ocnf-d18o_ocn)+abs(d18o_ocni-d18o_ocn));
                while(2>1)  % do
                    % case flxfrc2(1)=0
                    flxfrc2=zeros(1,nspcc);
                    flxfrc2(2) = flxfrc(1);
                    flxfrc2(3) = flxfrc(3);
                    flxfrc2(4) = flxfrc(2)-flxfrc2(3);
                    if (all(flxfrc2>=0d0))
                        break   % exit in fortran (terminates while-loop)
                    end
                    flxfrc2=zeros(1,nspcc);
                    % case flxfrc2(2)=0
                    flxfrc2(1) = flxfrc(1);
                    flxfrc2(3) = flxfrc(3)-flxfrc2(1);
                    flxfrc2(4) = flxfrc(4);
                    if (all(flxfrc2>=0d0));
                        break   % exit in fortran (terminates while-loop)
                    end
                    flxfrc2=zeros(1,nspcc);
                    % case flxfrc2(3)=0
                    flxfrc2(1) = flxfrc(3);
                    flxfrc2(2) = flxfrc(1)-flxfrc2(1);
                    flxfrc2(4) = flxfrc(2);
                    if (all(flxfrc2>=0d0));
                        break   % exit in fortran (terminates while-loop)
                    end
                    flxfrc2=zeros(1,nspcc);
                    % case flxfrc2(4)=0
                    flxfrc2(2) = flxfrc(4);
                    flxfrc2(1) = flxfrc(1)-flxfrc2(2);
                    flxfrc2(3) = flxfrc(2);
                    if (all(flxfrc2>=0d0))
                        break   % exit in fortran (terminates while-loop)
                    end
                    %                                             print*,'error' % should not come here
                    %                                             stop
                    msg = 'error in signal_flx, stop';  % should not come here
                    error(msg)
                end
                if(def_track2)
                    % tracking as time goes
                    % assignment of new caco3 species is conducted so that nspcc species is used up within 5000 interations
                    % timing of new species assignment can change at different time period during the event
                    if (time-time_spn<=time_trs/10d0)
                        if (mod(it,floor(5000/(nspcc-2)))==0)
                            cntsp=cntsp+1;           % new species assinged
                            d18o_sp(cntsp)=d18o_ocn;
                            d13c_sp(cntsp)=d13c_ocn;
                            ccflx = zeros(1, nspcc);
                            ccflx(cntsp) = ccflxi;
                        end
                    elseif (time-time_spn>time_trs/10d0 && time-time_spn<=time_trs/10d0*9d0)
                        if (mod(it,floor(5000/(nspcc-2)))==0)
                            cntsp=cntsp+1;
                            d18o_sp(cntsp)=d18o_ocn;
                            d13c_sp(cntsp)=d13c_ocn;
                            ccflx = zeros(1, nspcc);
                            ccflx(cntsp) = ccflxi;
                        end
                    elseif  (time-time_spn>time_trs/10d0*9d0)
                        if (mod(it,floor(5000/(nspcc-2)))==0)
                            cntsp=cntsp+1;
                            d18o_sp(cntsp)=d18o_ocn;
                            d13c_sp(cntsp)=d13c_ocn;
                            ccflx = zeros(1, nspcc);
                            ccflx(cntsp) = ccflxi;
                        end
                    end
                else % # def_track2
                    % usually using 4 species interpolation for 2 signals
                    % NEED to be extended so that any number of signals can be tracked with corresponding minimum numbers of caco3 species
                    for isp=1:nspcc
                        ccflx(isp)=flxfrc2(isp)*ccflxi;
                    end
                end %#if def_track2
                
                if(def_size)
                    for isp=1:4 % fine species
                        ccflx(isp)=flxfrc2(isp)*ccflxi*(1d0-flxfin);
                    end
                    for isp=5:8 % coarse species
                        ccflx(isp)=flxfrc2(isp-4)*ccflxi*flxfin;
                    end
                end %#if
                
            elseif (time>time_spn+time_trs)     % after event
                d13c_ocn = d13c_ocni;        % now again initial values
                d18o_ocn = d18o_ocni;
                ccflx = zeros(1, nspcc);
                ccflx(1) = ccflxi;  % allowing only caco3 flux with initial signal values
                if(def_track2)
                    if (cntsp+1~=nspcc)  % checking used caco3 species number is enough and necessary
                        % print*,'fatal error in counting',cntsp
                        % stop
                        fprintf('fatal error in counting  caco3 species numbers : %i',cntsp);
                        msg = 'fatal error in counting  caco3 species numbers, stop';
                        error(msg)
                    end
                    d18o_sp(cntsp+1)=d18o_ocn;
                    d13c_sp(cntsp+1)=d13c_ocn;
                    ccflx = 0d0;
                    ccflx(cntsp+1) = ccflxi;
                end %# if
                if(def_size)
                    ccflx = zeros(1, nspcc);
                    flxfin=flxfini;
                    ccflx(1) = (1d0-flxfin)*ccflxi;  % fine species
                    ccflx(5) = flxfin*ccflxi;  % coarse species
                end %# if
                if(def_biotest)
                    d13c_ocn = d13c_ocnf;   % finish with final value
                    d18o_ocn = d18o_ocni;
                    ccflx = zeros(1, nspcc);
                    ccflx(3) = ccflxi;  % caco3 species with d13c_ocnf and d18o_ocni
                end %# if
            end
            
        end
        
        function [dep] = bdcnd(time, time_spn, time_trs, depi, depf, def_biotest)
            
            if (time <= time_spn)    % spin-up
                dep = depi;  % initial depth
            elseif (time>time_spn && time<=time_spn+time_trs)  % during event
                if(~def_biotest)
                    % assuming a d13 excursion with shifts occurring with 1/10 times event duration
                    % the same assumption for depth change
                    if (time-time_spn<=time_trs/10d0)
                        dep = depi + (depf-depi)*(time-time_spn)/time_trs*10d0;
                    elseif (time-time_spn>time_trs/10d0 && time-time_spn<=time_trs/10d0*9d0)
                        dep = depf;
                    elseif  (time-time_spn>time_trs/10d0*9d0)
                        dep = 10d0*depf-9d0*depi - (depf-depi)*(time-time_spn)/time_trs*10d0;
                    end %if
                else %#
                    dep = depi;
                end %#if
            elseif (time>time_spn+time_trs)
                dep = depi;
            end
        end
        
        
        function sigrec(w,file_sigmlyid,file_sigmlydid,file_sigbtmid,time,age,izrec,d13c_blk,d13c_blkc  ...
                ,d13c_blkf,d18o_blk,d18o_blkc,d18o_blkf,ccx,mcc,rho,ptx,msed,izrec2,nz, def_size)
            %% recording signals at 3 different depths (btm of mixed layer, 2xdepths of btm of mixed layer and btm depth of calculation domain)
            
            
            if(~def_size)
                if (all(w>=0d0))    % not recording when burial is negative
                    
                    % 	write(file_sigmly,*) time-age(izrec),d13c_blk(izrec),d18o_blk(izrec) &
                    %   	,sum(ccx(izrec,:))*mcc/rho(izrec)*100d0,ptx(izrec)*msed/rho(izrec)*100d0
                    fprintf(file_sigmlyid,'%17.16e \t %17.16e \t %17.16e \t %17.16e \t %17.16e \n ',...
                        time-age(izrec),d13c_blk(izrec),d18o_blk(izrec), sum(ccx(izrec,:))*mcc/rho(izrec)*100d0, ptx(izrec)*msed/rho(izrec)*100d0 );
                    % 	write(file_sigmlyd,*) time-age(izrec2),d13c_blk(izrec2),d18o_blk(izrec2) &
                    %    	,sum(ccx(izrec2,:))*mcc/rho(izrec2)*100d0,ptx(izrec2)*msed/rho(izrec2)*100d0
                    fprintf(file_sigmlydid,'%17.16e \t %17.16e \t %17.16e \t %17.16e \t %17.16e \n ',...
                        time-age(izrec2),d13c_blk(izrec2),d18o_blk(izrec2), sum(ccx(izrec2,:))*mcc/rho(izrec2)*100d0,ptx(izrec2)*msed/rho(izrec2)*100d0 );
                    %  	write(file_sigbtm,*) time-age(nz),d13c_blk(nz),d18o_blk(nz) &
                    %       ,sum(ccx(nz,:))*mcc/rho(nz)*100d0,ptx(nz)*msed/rho(nz)*100d0
                    fprintf(file_sigbtmid,'%17.16e \t %17.16e \t %17.16e \t %17.16e \t %17.16e \n ',...
                        time-age(nz),d13c_blk(nz),d18o_blk(nz), sum(ccx(nz,:))*mcc/rho(nz)*100d0,ptx(nz)*msed/rho(nz)*100d0 );
                end
            else
                if (all(w>=0d0))    % not recording when burial is negative
                    % 	    write(file_sigmly,*) time-age(izrec),d13c_blk(izrec),d18o_blk(izrec) &
                    %        ,sum(ccx(izrec,:))*mcc/rho(izrec)*100d0,ptx(izrec)*msed/rho(izrec)*100d0  &
                    %        ,d13c_blkf(izrec),d18o_blkf(izrec),sum(ccx(izrec,1:4))*mcc/rho(izrec)*100d0  &
                    %       ,d13c_blkc(izrec),d18o_blkc(izrec),sum(ccx(izrec,5:8))*mcc/rho(izrec)*100d0
                    fprintf(file_sigmlyid,'%17.16e \t %17.16e \t %17.16e \t %17.16e \t %17.16e \t %17.16e \t %17.16e \t %17.16e \t %17.16e \t %17.16e \t %17.16e \n ',...
                        time-age(izrec),d13c_blk(izrec),d18o_blk(izrec) ...
                        ,sum(ccx(izrec,:))*mcc/rho(izrec)*100d0,ptx(izrec)*msed/rho(izrec)*100d0  ...
                        ,d13c_blkf(izrec),d18o_blkf(izrec),sum(ccx(izrec,1:4))*mcc/rho(izrec)*100d0  ...
                        ,d13c_blkc(izrec),d18o_blkc(izrec),sum(ccx(izrec,5:8))*mcc/rho(izrec)*100d0  );
                    %     write(file_sigmlyd,*) time-age(izrec2),d13c_blk(izrec2),d18o_blk(izrec2) &
                    %         ,sum(ccx(izrec2,:))*mcc/rho(izrec2)*100d0,ptx(izrec2)*msed/rho(izrec2)*100d0  &
                    %         ,d13c_blkf(izrec2),d18o_blkf(izrec2),sum(ccx(izrec2,1:4))*mcc/rho(izrec2)*100d0  &
                    %         ,d13c_blkc(izrec2),d18o_blkc(izrec2),sum(ccx(izrec2,5:8))*mcc/rho(izrec2)*100d0
                    fprintf(file_sigmlydid,'%17.16e \t %17.16e \t %17.16e \t %17.16e \t %17.16e \t %17.16e \t %17.16e \t %17.16e \t %17.16e \t %17.16e \t %17.16e \n ',...
                        time-age(izrec2),d13c_blk(izrec2),d18o_blk(izrec2) ...
                        ,sum(ccx(izrec2,:))*mcc/rho(izrec2)*100d0,ptx(izrec2)*msed/rho(izrec2)*100d0 ...
                        ,d13c_blkf(izrec2),d18o_blkf(izrec2),sum(ccx(izrec2,1:4))*mcc/rho(izrec2)*100d0  ...
                        ,d13c_blkc(izrec2),d18o_blkc(izrec2),sum(ccx(izrec2,5:8))*mcc/rho(izrec2)*100d0 );
                    %     write(file_sigbtm,*) time-age(nz),d13c_blk(nz),d18o_blk(nz) &
                    %         ,sum(ccx(nz,:))*mcc/rho(nz)*100d0,ptx(nz)*msed/rho(nz)*100d0 &
                    %         ,d13c_blkf(nz),d18o_blkf(nz),sum(ccx(nz,1:4))*mcc/rho(nz)*100d0  &
                    %         ,d13c_blkc(nz),d18o_blkc(nz),sum(ccx(nz,5:8))*mcc/rho(nz)*100d0
                    fprintf(file_sigbtmid,'%17.16e \t %17.16e \t %17.16e \t %17.16e \t %17.16e \t %17.16e \t %17.16e \t %17.16e \t %17.16e \t %17.16e \t %17.16e \n ',...
                        time-age(nz),d13c_blk(nz),d18o_blk(nz) ...
                        ,sum(ccx(nz,:))*mcc/rho(nz)*100d0,ptx(nz)*msed/rho(nz)*100d0 ...
                        ,d13c_blkf(nz),d18o_blkf(nz),sum(ccx(nz,1:4))*mcc/rho(nz)*100d0  ...
                        ,d13c_blkc(nz),d18o_blkc(nz),sum(ccx(nz,5:8))*mcc/rho(nz)*100d0 );
                end
                
            end
            
        end
        
    end
    
end