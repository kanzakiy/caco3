classdef caco3_test
    
    % testing the functions and subroutines
    properties
    end
    
    methods(Static)
        
        function bc = caco3_set_boundary_cond(global_var)
            
            bc.ccflxi = 10d-6;      % mol (CaCO3) cm-2 yr-1  % a reference caco3 flux; Emerson and Archer (1990)
            bc.omflx = 12d-6;       % mol cm-2 yr-1       % a reference om flux; Emerson and Archer (1990)
            bc.detflx = 180d-6;     % g cm-2 yr-1  % a reference detrital flux; MUDS input http://forecast.uchicago.edu/Projects/muds.html
            bc.om = 1d-8*ones(1,global_var.nz);           % assume an arbitrary low conc.
            bc.alki =  2285d0;      % uM  % a reference ALK; MUDS
            bc.dici = 2211d0;       % uM   % a reference DIC; MUDS
            bc.o2i = 165d0;         % uM     % a reference O2 input ; MUDS
            bc.o2 = 0.0;
            
            bc.oxic = true;            % oxic only model of OM degradation by Emerson (1985)
            bc.anoxic = true;          % oxic-anoxic model of OM degradation by Archer (1991)
            if(global_var.def_oxonly)
                bc.anoxic = false;         % oxic only model of OM degradation by Emerson (1985)
            end
        end
        
        function mix_type = caco3_set_mixing(global_var)
            %%  define the type of mixing to be used
            
            % biogenic reworking assumed?
            if(global_var.def_allnobio)
                mix_type.nobio =  true(1, global_var.nspcc + 2);
            else
                mix_type.nobio = false(1, global_var.nspcc + 2);
            end
            
            % mixing info from labs?
            if(global_var.def_alllabs)
                mix_type.labs =  true(1, global_var.nspcc + 2);
            else
                mix_type.labs = false(1, global_var.nspcc + 2);
            end
            
            % random mixing?
            if(global_var.def_allturbo2)
                mix_type.turbo2 =  true(1, global_var.nspcc + 2);
            else
                mix_type.turbo2 = false(1, global_var.nspcc + 2);
            end
            
            % ON if assuming non-local mixing (i.e., if labs or turbo2 is ON)
            if(global_var.def_allnonlocal)
                mix_type.nonlocal =  true(1, global_var.nspcc + 2);
            else
                mix_type.nonlocal = false(1, global_var.nspcc + 2);
            end
            
        end
        
        
        function test_caco3_therm()
            
            % test the functions and subroutines to calculate caco3 thermodynamics
            tmp = 2.0;
            sal = 35.0;
            n = 100;
            
            
            %            disp('     depth      keqcc      keqag')
            fprintf('%4s %17s %17s \n','depth','keqcc','keqag');
            for i=1:n
                dep = 6*i/n;
                keqcc = caco3_therm.calceqcc(tmp,sal,dep);
                keqag = caco3_therm.calceqag(tmp,sal,dep);
                fprintf('%4.3f %17.16e %17.16e\n', dep, keqcc, keqag);
            end
        end
        
        function chk_caco3_therm_sbrtns()
            
            % % check the calculation of aqueous co2 species (and its derivatives) using subroutines in caco3_therm.f90
            tmp = 2.0;
            sal = 35.0;
            dep = 4.0d0; % depth in km
            nz = 100;
            dic = ones(nz,1)*2211.0d0*1.0d-6/1.0d3;  % 2211 uM converted to mol/cm3; all the same value at different grids
            alk = ones(nz,1)*2285.0d0*1.0d-6/1.0d3;
            
            
            % calling subroutine to calculate all aqueous co2 species and pH
            [ph,co2,hco3,co3,info] = caco3_therm.calcspecies(dic,alk,tmp,sal,dep);
            if (info~=0) % if ph cannot be calculated in the subroutine, info=1 is returned
                fprintf('error in calcspecies');
                exit()  % stop in fortran
            end
            
            % calling subroutine to calculate derivatives of co3 conc. wrt dic and alk (defined as dco3_ddic and dco3_dalk, respectively)
            [dco3_dalk,dco3_ddic, info] = caco3_therm.calcdevs(dic,alk,tmp,sal,dep);
            if (info~=0) then % if ph cannot be calculated in the subroutine, info=1 is returned
                fprintf('error in calcdevs')
                exit()  % stop in fortran
            end
            
            % printing results on screen; if you want check with MATLAB version by copy and paste
            for iz=1:nz
                %%%% print*,dic(iz),alk(iz),co2(iz),hco3(iz),co3(iz),ph(iz),dco3_ddic(iz),dco3_dalk(iz)
                fprintf('%17.16e %17.16e %17.16e %17.16e %17.16e %17.16e %17.16e %17.16e\n', dic(iz),alk(iz),co2(iz),hco3(iz),co3(iz),ph(iz),dco3_ddic(iz),dco3_dalk(iz));
                
                % note that the above concentrations (alk,dic,co2,hco3,co3) are all in mol/cm3
                % ph is actually H+ concentration in mol/L
                % dco3_dalk and dco3_ddic should be dimensionless
            end
            
        end
        
        function chk_grid()
            % % check the calculation of sediment grid cells caco3_main.f90
            
            nz = 100;
            
            ztot=500.0d0;
            beta = 1.00000000005d0;  % a parameter to make a grid; closer to 1, grid space is more concentrated around the sediment-water interface (SWI)
            [dz,z] = caco3_main.makegrid(beta,global_var.nz, global_var.ztot, global_var.def_recgrid);
            
            for iz=1:nz
                fprintf('%2.1f %17.16e %17.16e\n', iz,dz(iz),z(iz));
            end
            
        end
        
        function chk_poro()
            % % check the calculation of transition matrix using subroutines in caco3_main.f90
            
            % initialize/define global properties/variables
            global_var = caco3_main;
            % global_var = caco3_main(arg1, arg2);
            
            beta = 1.00000000005d0;  % a parameter to make a grid; closer to 1, grid space is more concentrated around the sediment-water interface (SWI)
            [global_var.dz,global_var.z] = caco3_main.makegrid(beta,global_var.nz, global_var.ztot, global_var.def_recgrid);
            
            
            %            global_var = caco3_main.getporosity(global_var.z, global_var); % assume porosity profile
            [global_var.poro, global_var.porof, global_var.sporof, global_var.sporo, global_var.sporoi] = caco3_main.getporosity(global_var.z, global_var.poroi, global_var.nz);      % assume porosity profile
            
            % showing porosity profile as function of depth
            % porosity              ---> poro
            % volume solid fraction ---> sporo = 1 - poro
            % z                     ---> depth
            fprintf('\n \n');
            for iz=1:global_var.nz
                fprintf('%17.16e %17.16e %17.16e\n',global_var.z(iz),global_var.poro(iz),global_var.sporo(iz));
            end
            
        end
        
        function chk_coefs(tmp_in, sal_in, dep_in)
            % % check the calculation of kinetic and thermodynamic coefficients using subroutines in caco3_test_mod_v5_6.f90
            
            % initialize/define global properties/variables
            global_var = caco3_main;
            tmp = tmp_in;   % 2.0;
            sal = sal_in;   % 35.0;
            dep = dep_in;   % 0.0d0; % depth in km
            
            beta = 1.00000000005d0;  % a parameter to make a grid; closer to 1, grid space is more concentrated around the sediment-water interface (SWI)
            [global_var.dz,global_var.z] = caco3_main.makegrid(beta,global_var.nz, global_var.ztot, global_var.def_recgrid);
            
            
            %            global_var = caco3_main.getporosity(global_var.z, global_var); % assume porosity profile
            [global_var.poro, global_var.porof, global_var.sporof, global_var.sporo, global_var.sporoi] = caco3_main.getporosity(global_var.z, global_var.poroi, global_var.nz);      % assume porosity profile
            
            % determining diffusion coeffs, rate constants for om decomposition and cc dissolution
            % , dissociation constants of co2 and hco3, and co3 conc. at calcite saturation
            % Note that porosity is needed for diffusion coefficient to reflect tortuosity
            [keq1    ,keq2   ,keqcc  ,co3sat, dif_dic    ,dif_alk    ,dif_o2 ,kom, kcc, global_var] = caco3_main.coefs(tmp,sal,dep, global_var);
            
            % showing thermodynamic consts which are not functions of depth
            fprintf('keq1    ,keq2   ,keqcc  ,co3sat, kom  \n');
            fprintf('%17.16e %17.16e %17.16e %17.16e %17.16e\n\n',keq1    ,keq2   ,keqcc  ,co3sat, kom);
            
            
            % showing coeffs against depth
            fprintf('z   ,dif_dic    ,dif_alk    ,dif_o2 \n');
            for iz=1:global_var.nz
                fprintf('%17.16e %17.16e %17.16e %17.16e %17.16e',global_var.z(iz), dif_dic(iz),dif_alk(iz),dif_o2(iz));
                fprintf(' \n');
                % depth [cm], dic, alk and o2 diffusion coefficients [cm2 yr-1], om decompposition rate const. [ yr-1]
            end
            
            % dissolution rate constant are species specific so recording on txt files
            fprintf('z ,    Species 1,   Species 2,    Species 3,    Species 4 \n');
            for iz=1:global_var.nz
                fprintf('%17.16e %17.16e %17.16e %17.16e %17.16e\n\n\n',global_var.z(iz), kcc);
                fprintf(' \n');
            end
            
        end
        
        function chk_bur(tmp_in, sal_in, dep_in)
            % % check the calculation of coefficients for sediment burial using subroutines in caco3_test_mod_v5_6.f90
            % <<<< Maybe first only consider Fickian mixing so switch off all macros in defines.h (you can switch on 'test') >>>>>
            
            % you can check coefficients for sediment burial
            % please copy and paste results to whatever file to be compared with results with MATLAB version
            
            % initialize/define global properties/variables
            global_var = caco3_main;
            
            % initialize the boundary conditions
            bc = caco3_test.caco3_set_boundary_cond(global_var);
            
            tmp = tmp_in;   % 2.0;
            sal = sal_in;   % 35.0;
            dep = dep_in;   % 0.0d0; % depth in km
            
            beta = 1.00000000005d0;  % a parameter to make a grid; closer to 1, grid space is more concentrated around the sediment-water interface (SWI)
            [global_var.dz,global_var.z] = caco3_main.makegrid(beta,global_var.nz, global_var.ztot, global_var.def_recgrid);
            
            
            %            global_var = caco3_main.getporosity(global_var.z, global_var); % assume porosity profile
            [global_var.poro, global_var.porof, global_var.sporof, global_var.sporo, global_var.sporoi] = caco3_main.getporosity(global_var.z, global_var.poroi, global_var.nz);      % assume porosity profile
            %%%%%%%%%%%%% flx assignement and initial guess for burial rate %%%%%%%%%%%%%%%%%%%%%%
            % assume fluxes of om, cc and clay, required to calculate burial velocity
            [omflx, detflx, ccflx] = caco3_main.flxstat(global_var.om2cc, bc.ccflxi, global_var.mcc, global_var.nspcc);
            fprintf('ccflx %17.16e \n', ccflx);
            fprintf('om2cc, ccflxi, detflx, omflx, sum(ccflx) \n');
            fprintf('%17.16e %17.16e %17.16e %17.16e %17.16e \n', global_var.om2cc, bc.ccflxi, detflx, omflx, sum(ccflx));
            fprintf(' \n');
            
            % molar volume (cm3 mol-1) needed for burial rate calculation
            mvom = global_var.mom/global_var.rhoom;  % om
            mvsed = global_var.msed/global_var.rhosed; % clay
            mvcc = global_var.mcc/global_var.rhocc; % caco3
            
            % initial guess of burial profile, requiring porosity profile
            % w = burial rate, wi = burial rate initial guess
            [w, wi] = caco3_main.burial_pre(detflx,ccflx,global_var.msed,mvsed,mvcc,global_var.poroi, global_var.nz);
            
            % % depth -age conversion
            age = caco3_main.dep2age(global_var.dz, w, global_var.nz);
            
            % % determine factors for upwind scheme to represent burial advection
            [up, dwn, cnr, adf] = caco3_main.calcupwindscheme(w, global_var.nz);
            
            fprintf('z   ,w  ,age    ,up    ,dwn ,cnr    ,adf     \n');
            % showing parameters relevant to burial on screen
            for iz=1:global_var.nz
                fprintf('%17.16e %17.16e %17.16e %17.16e %17.16e %17.16e %17.16e \n', global_var.z(iz),w(iz),age(iz),up(iz),dwn(iz),cnr(iz) ,adf(iz));
                % depth [cm], burial rate [cm yr-1], grid age [yr], and factors for burial advection;
                % up, dwn, cnr uses (i-1,i),(i,i+1) and (i-1,i+1), respectively for burial at i.
                % adf is the factor to account for mass balance when cnr is not zero.
            end
        end
        
        function chk_trans(tmp_in, sal_in, dep_in)
            % % check the calculation of transition matrix using subroutines in caco3_test_mod_v5_6.f90
            % you can check transition matrix for individual species
            % please copy and paste results to whatever file to be compared with results with MATLAB version
            
            % initialize/define global properties/variables
            global_var = caco3_main;
            
            % initialize mixing type
            mix_type = caco3_test.caco3_set_mixing(global_var);
            
            % initialize the boundary conditions
            bc = caco3_test.caco3_set_boundary_cond(global_var);
            
            tmp = tmp_in;   % 2.0;
            sal = sal_in;   % 35.0;
            dep = dep_in;   % 0.0d0; % depth in km
            
            beta = 1.00000000005d0;  % a parameter to make a grid; closer to 1, grid space is more concentrated around the sediment-water interface (SWI)
            [global_var.dz,global_var.z] = caco3_main.makegrid(beta,global_var.nz, global_var.ztot, global_var.def_recgrid);
            
            
            %           global_var = caco3_main.getporosity(global_var.z, global_var); % assume porosity profile
            [global_var.poro, global_var.porof, global_var.sporof, global_var.sporo, global_var.sporoi] = caco3_main.getporosity(global_var.z, global_var.poroi, global_var.nz);      % assume porosity profile
            
            % % make transition matrix
            [trans,izrec,izrec2,izml,mix_type.nonlocal] = ...
                caco3_main.make_transmx(mix_type.labs,global_var.nspcc,mix_type.turbo2,mix_type.nobio,global_var.dz,global_var.sporo,global_var.nz,global_var.z, global_var.zml_ref, global_var.def_size);
            % recording transition matrices for individual species (isp=1 --> om, 2-->clay,3~2+nspcc --> caco3 species)
            % in default, transition matrices are the same for all the species
            % file will be created your working directory of a name 'chk_trans_sp-xxx.txt' where xxx denotes species number
            % sizes of matrices are all (nz, nz )
            for isp=1:global_var.nspcc+2
                %     write(dumchr(1),'(i3.3)') isp
                %     open(unit=file_tmp,file='./chk_trans_sp-'//trim(adjustl(dumchr(1)))//'.txt',action='write',status='replace')
                str = sprintf('matlab_chk_trans_sp-%3.3i.txt',isp);
                file_tmp = fopen(str,'wt');
                
                for iz=1:global_var.nz
                    %        write(file_tmp,*)(trans(iz,iiz,isp),iiz=1,nz)
                    fprintf(file_tmp,'%17.16e\t',trans(iz,:,isp));
                    fprintf(file_tmp,'\n');
                end
                fclose(file_tmp);
            end
        end
        
        function chk_om(tmp_in, sal_in, dep_in)
            % % check the calculation of om conc. and flx using subroutines in caco3_test_mod_v5_6.f90
            % you can check om calculation including om conc. and flx
            % please copy and paste results to whatever file to be compared with results with MATLAB version
            % NOTE: here oxygen conc. has to be assume. The model calculate om and o2 iteratively, but not here.
            
            interval =10; % choose a value between 1 to nz; om depth profile is shown with this interval; e.g., if inteval = nz, om conc. at all depths are shown
            % e.g., if interval = 5, om conc. at 5 depths are shown
            
            % initialize/define global properties/variables
            global_var = caco3_main;
            
            % initialize mixing type
            mix_type = caco3_test.caco3_set_mixing(global_var);
            
            % initialize the boundary conditions
            bc = caco3_test.caco3_set_boundary_cond(global_var);
            
            tmp = tmp_in;   % 2.0;
            sal = sal_in;   % 35.0;
            dep = dep_in;   % 0.0d0; % depth in km
            
            beta = 1.00000000005d0;  % a parameter to make a grid; closer to 1, grid space is more concentrated around the sediment-water interface (SWI)
            [global_var.dz,global_var.z] = caco3_main.makegrid(beta,global_var.nz, global_var.ztot, global_var.def_recgrid);
            
            
            %            global_var = caco3_main.getporosity(global_var.z, global_var); % assume porosity profile
            [global_var.poro, global_var.porof, global_var.sporof, global_var.sporo, global_var.sporoi] = caco3_main.getporosity(global_var.z, global_var.poroi, global_var.nz);      % assume porosity profile
            %%%%%%%%%%%%% flx assignement and initial guess for burial rate %%%%%%%%%%%%%%%%%%%%%%
            % assume fluxes of om, cc and clay, required to calculate burial velocity
            [omflx, detflx, ccflx] = caco3_main.flxstat(global_var.om2cc, bc.ccflxi, global_var.mcc, global_var.nspcc);
            fprintf('ccflx %17.16e \n', ccflx);
            fprintf('om2cc, ccflxi, detflx, omflx, sum(ccflx) \n');
            fprintf('%17.16e %17.16e %17.16e %17.16e %17.16e \n', global_var.om2cc, bc.ccflxi, detflx, omflx, sum(ccflx));
            fprintf(' \n');
            
            % molar volume (cm3 mol-1) needed for burial rate calculation
            mvom = global_var.mom/global_var.rhoom;  % om
            mvsed = global_var.msed/global_var.rhosed; % clay
            mvcc = global_var.mcc/global_var.rhocc; % caco3
            
            % initial guess of burial profile, requiring porosity profile
            % w = burial rate, wi = burial rate initial guess
            [w, wi] = caco3_main.burial_pre(detflx,ccflx,global_var.msed,mvsed,mvcc,global_var.poroi, global_var.nz);
            
            % % depth -age conversion
            age = caco3_main.dep2age(global_var.dz, w, global_var.nz);
            
            % % determine factors for upwind scheme to represent burial advection
            [up, dwn, cnr, adf] = caco3_main.calcupwindscheme(w, global_var.nz);
            
            
            % % make transition matrix
            [trans,izrec,izrec2,izml,mix_type.nonlocal] = ...
                caco3_main.make_transmx(mix_type.labs,global_var.nspcc,mix_type.turbo2,mix_type.nobio,global_var.dz,global_var.sporo,global_var.nz,global_var.z, global_var.zml_ref, global_var.def_size);
            
            
            o2 = bc.o2i*1d-6/1d3 * ones(1, global_var.nz);     % o2 conc. in uM converted to mol/cm3
            
            omx = bc.om;
            o2x = o2;
            
            time = 0d0; % model time [yr]
            it = 1; % integration count
            nt = 10; % total integration
            dt = 100d0; % time step [yr]
            
            rho = 2.5d0; % assume here
            
            for it=1:nt
                %    print'(A,i0,A,E11.3,A,E11.3,A)','(it,dt,time)  (',it,',',dt,',',time,')'
                fprintf('(it,dt,time)  %2.2i %11.3e %11.3e\n', it, dt, time);
                
                % %%%%%%%%%%%%%%%  om conc. calculation  % %%%%%%%%%%%%%%%
                % omx: mol cm-3 sld; om conc.
                % izox: integer for grid number of zox
                % kom: degradation rate consts. for each nz grids
                [omx, izox, kom] = ...
                    caco3_main.omcalc(bc.oxic, bc.anoxic,o2x,bc.om, global_var.komi,global_var.nz,global_var.sporo,global_var.sporoi,global_var.sporof, w, wi, dt, up, dwn, cnr, adf,trans, ...
                    global_var.nspcc, mix_type.labs,mix_type.turbo2, mix_type.nonlocal, omflx, global_var.poro, global_var.dz, global_var.o2th);
                % calculating the fluxes relevant to om diagenesis (and checking the calculation satisfies the difference equations )
                [omadv,omdec,omdif,omrain,omres,omtflx] = ...
                    caco3_main.calcflxom(omflx,global_var.sporo,bc.om,omx,dt,w,global_var.dz,global_var.z,global_var.nz,mix_type.turbo2,mix_type.labs, global_var.poro,up,dwn,cnr,adf,rho, global_var.mom,trans,kom,global_var.sporof);
                % % Dominik: what are the next 2 lines for? structure the output?
                % %                     write(dumchr(2),'(i0)') interval
                % %                     dumchr(1)="(A,"//trim(adjustl(dumchr(2)))//"E11.3"//")"
                
                omx_wtpc = omx.*global_var.mom./rho*100d0;
                
                %                     % showing results on screen
                fprintf('~~~~ conc ~~~~ \n');
                % %                      fprintf('z  : \n');
                % %                     fprintf('%17.16e\n', global_var.z(1:interval:100));
                % %                      fprintf('OM  : \n');
                % %                     fprintf('%17.16e\n', omx_wtpc(1:interval:100));
                fprintf('z   ,OM   \n');
                % showing parameters relevant to burial on screen
                for iz=1:interval:global_var.nz
                    fprintf('%17.16e \t %17.16e\n', global_var.z(iz), omx_wtpc(iz));
                end
                fprintf('++++ flx ++++ \n');
                fprintf(' \t tflx \t adv \t dif \t omrxn \t ccrxn \t rain \t res \n')
                fprintf('OM :  \t %17.16e \t %17.16e \t %17.16e \t  %17.16e \t %17.16e \t %17.16e \t %17.16e', omtflx, omadv,  omdif, omdec,0d0,omrain, omres);
                
                fprintf('\n');
                fprintf('\n');
                fprintf('\n');
                %
                bc.om = omx;
                time = time +dt;
                %
            end
            
        end
        
        
        
        function chk_om_o2(tmp_in, sal_in, dep_in)
            % % check the calculation of conc. and flx of om & o2 using subroutines in caco3_test_mod_v5_6.f90
            % you can check calculation of concs. and flxes of om & o2
            % please copy and paste results to whatever file to be compared with results with MATLAB version
            % NOTE: This tests o2 calculation, but o2 calculation should be boring without om, so the code here does iterations at individual time steps as in the whole code.
            
            interval =10; % choose a value between 1 to nz; om depth profile is shown with this interval; e.g., if inteval = nz, om conc. at all depths are shown
            % e.g., if interval = 5, om conc. at 5 depths are shown
            
            % initialize/define global properties/variables
            global_var = caco3_main;
            
            
            % initialize mixing type
            mix_type = caco3_test.caco3_set_mixing(global_var);
            
            % initialize the boundary conditions
            bc = caco3_test.caco3_set_boundary_cond(global_var);
            
            tmp = tmp_in;   % 2.0;
            sal = sal_in;   % 35.0;
            dep = dep_in;   % 0.0d0; % depth in km
            
            beta = 1.00000000005d0;  % a parameter to make a grid; closer to 1, grid space is more concentrated around the sediment-water interface (SWI)
            [global_var.dz,global_var.z] = caco3_main.makegrid(beta,global_var.nz, global_var.ztot, global_var.def_recgrid);
            
            
            %            global_var = caco3_main.getporosity(global_var.z, global_var); % assume porosity profile
            [global_var.poro, global_var.porof, global_var.sporof, global_var.sporo, global_var.sporoi] = caco3_main.getporosity(global_var.z, global_var.poroi, global_var.nz);      % assume porosity profile
            %%%%%%%%%%%%% flx assignement and initial guess for burial rate %%%%%%%%%%%%%%%%%%%%%%
            % assume fluxes of om, cc and clay, required to calculate burial velocity
            [omflx, detflx, ccflx] = caco3_main.flxstat(global_var.om2cc, bc.ccflxi, global_var.mcc, global_var.nspcc);
            %            fprintf('ccflx %17.16e \n', ccflx);
            %            fprintf('om2cc, ccflxi, detflx, omflx, sum(ccflx) \n');
            %            fprintf('%17.16e %17.16e %17.16e %17.16e %17.16e \n', global_var.om2cc, bc.ccflxi, detflx, omflx, sum(ccflx));
            %            fprintf(' \n');
            
            % molar volume (cm3 mol-1) needed for burial rate calculation
            mvom = global_var.mom/global_var.rhoom;  % om
            mvsed = global_var.msed/global_var.rhosed; % clay
            mvcc = global_var.mcc/global_var.rhocc; % caco3
            
            % initial guess of burial profile, requiring porosity profile
            % w = burial rate, wi = burial rate initial guess
            [w, wi] = caco3_main.burial_pre(detflx,ccflx,global_var.msed,mvsed,mvcc,global_var.poroi, global_var.nz);
            
            % % depth -age conversion
            age = caco3_main.dep2age(global_var.dz, w, global_var.nz);
            
            % % determine factors for upwind scheme to represent burial advection
            [up, dwn, cnr, adf] = caco3_main.calcupwindscheme(w, global_var.nz);
            
            
            % % make transition matrix
            [trans,izrec,izrec2,izml,mix_type.nonlocal] = ...
                caco3_main.make_transmx(mix_type.labs,global_var.nspcc,mix_type.turbo2,mix_type.nobio,global_var.dz,global_var.sporo,global_var.nz,global_var.z, global_var.zml_ref, global_var.def_size);
            
            [keq1    ,keq2   ,keqcc  ,co3sat, dif_dic    ,dif_alk    ,dif_o2 ,kom, kcc, global_var] = caco3_main.coefs(tmp,sal,dep, global_var);
            
            bc.o2 = bc.o2i*1d-6/1d3 * ones(1, global_var.nz);     % o2 conc. in uM converted to mol/cm3
            
            omx = bc.om;
            o2x = bc.o2;
            
            time = 0d0; % model time [yr]
            it = 1; % integration count
            nt = 10; % total integration
            dt = 100d0; % time step [yr]
            
            rho = 2.5d0; % assume here
            
            %%%  addition to chk_om.f90 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            global_var.zox = 10d0;  % initial assumption on oxygen penetaration depth [cm]
            %%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            
            for it=1:nt
                %    print'(A,i0,A,E11.3,A,E11.3,A)','(it,dt,time)  (',it,',',dt,',',time,')'
                fprintf('(it,dt,time)  %2.2i %11.3e %11.3e\n', it, dt, time);
                
                itr = 0;        % iteration number for om and o2 calcuation
                error = 1d4;  	% error in ieration for zox
                minerr= 1d4;    % recording minimum relative difference in zox from previously considered zox
                
                while(error > global_var.tol)
                    % %%%%%%%%%%%%%%%  om conc. calculation  % %%%%%%%%%%%%%%%
                    % omx: mol cm-3 sld; om conc.
                    % izox: integer for grid number of zox
                    % kom: degradation rate consts. for each nz grids
                    [omx, izox, kom] = ...
                        caco3_main.omcalc(bc.oxic, bc.anoxic,o2x,bc.om, global_var.komi,global_var.nz,global_var.sporo,global_var.sporoi,global_var.sporof, w, wi, dt, up, dwn, cnr, adf,trans, ...
                        global_var.nspcc, mix_type.labs,mix_type.turbo2, mix_type.nonlocal, omflx, global_var.poro, global_var.dz, global_var.o2th);
                    % calculating the fluxes relevant to om diagenesis (and checking the calculation satisfies the difference equations )
                    [omadv,omdec,omdif,omrain,omres,omtflx] = ...
                        caco3_main.calcflxom(omflx,global_var.sporo,bc.om,omx,dt,w,global_var.dz,global_var.z,global_var.nz,mix_type.turbo2,mix_type.labs, global_var.poro,up,dwn,cnr,adf,rho, global_var.mom,trans,kom,global_var.sporof);
                    
                    fprintf('izox = %i \n',izox);      % sb omcalc calculates izox, which is the deepest grid where o2 >=0.
                    
                    % %                     omx_wtpc = omx.*global_var.mom./rho*100d0;
                    % %                     fprintf('~~~~ conc ~~~~ \n');
                    % %                     fprintf('z   ,OM   \n');
                    % %                     % showing parameters relevant to burial on screen
                    % %                     for iz=1:interval:global_var.nz
                    % %                         fprintf('%17.16e \t %17.16e\n', global_var.z(iz), omx_wtpc(iz));
                    % %                     end
                    % %                     fprintf('++++ flx ++++ \n');
                    % %                     fprintf(' \t tflx \t adv \t dif \t omrxn \t ccrxn \t rain \t res \n')
                    % %                     fprintf('OM :  \t %17.16e \t %17.16e \t %17.16e \t  %17.16e \t %17.16e \t %17.16e \t %17.16e', omtflx, omadv,  omdif, omdec,0d0,omrain, omres);
                    % %     %                fprintf('O2 :  \t %17.16e \t %17.16e \t %17.16e \t  %17.16e \t %17.16e \t %17.16e \t %17.16e', o2tflx, 0d0, o2dif, o2dec, 0d0, 0d0, o2res);
                    
                    
                    if (izox == global_var.nz)     % fully oxic; lower boundary condition ---> no diffusive out flow
                        % o2 calculation when o2 penetration depth (zox) is the same as bottom depth.
                        o2x = caco3_main.o2calc_ox(izox,global_var.nz,global_var.poro,bc.o2,kom,omx,global_var.sporo,dif_o2,global_var.dz,dt, global_var.ox2om, bc.o2i);
                        %  fluxes relevant to o2 (at the same time checking the satisfaction of difference equations)
                        [o2dec,o2dif,o2tflx,o2res] = caco3_main.calcflxo2_ox(global_var.nz,global_var.sporo,kom,omx,global_var.dz,global_var.poro,dif_o2,dt,bc.o2,o2x, global_var.ox2om, bc.o2i);
                    else        %% if oxygen is depleted within calculation domain, lower boundary changes to zero concs.
                        % o2 calculation when o2 is depleted within the calculation domain.
                        o2x = caco3_main.o2calc_sbox(izox,global_var.nz,global_var.poro,bc.o2,kom,omx,global_var.sporo,dif_o2,global_var.dz,dt, global_var.ox2om, bc.o2i);
                        % fluxes relevant to oxygen
                        [o2dec,o2dif,o2tflx,o2res] = caco3_main.calcflxo2_sbox(global_var.nz,global_var.sporo,kom,omx,global_var.dz,global_var.poro,dif_o2,dt,bc.o2,o2x,izox, global_var.ox2om, bc.o2i);
                    end
                    
                    %                 % showing intermediate results on screen
                    %                 omx_wtpc = omx.*global_var.mom./rho*100d0;
                    %                 fprintf('~~~~ conc ~~~~ itr = %i \n', itr +1);
                    %                 fprintf('z   ,OM   , o2\n');
                    %                 % showing parameters relevant to burial on screen
                    %                 for iz=1:interval:global_var.nz
                    %                     fprintf('%17.16e \t %17.16e \t %17.16e \n', global_var.z(iz), omx_wtpc(iz), o2x(iz)*1d3);
                    %                 end
                    %                 fprintf('++++ flx ++++ \n');
                    %                 fprintf(' \t tflx \t adv \t dif \t omrxn \t ccrxn \t rain \t res \n')
                    %                 fprintf('OM :  \t %17.16e \t %17.16e \t %17.16e \t  %17.16e \t %17.16e \t %17.16e \t %17.16e \n', omtflx, omadv,  omdif, omdec,0d0,omrain, omres);
                    %                 fprintf('O2 :  \t %17.16e \t %17.16e \t %17.16e \t  %17.16e \t %17.16e \t %17.16e \t %17.16e', o2tflx, 0d0, o2dif, o2dec, 0d0, 0d0, o2res);%
                    %                 fprintf('\n');
                    %                 fprintf('\n');
                    
                    
                    % update of zox
                    zoxx = 0d0;      % zox dummy variable
                    iz=1;
                    %                    for iz=1:global_var.nz
                    while(iz<global_var.nz+1)
                        if (o2x(iz)<=0d0)
                            break   % exit in fortran (terminates for-loop)
                        end
                        iz = iz+1;
                    end
                    
                    if (iz==global_var.nz+1) % oxygen never gets less than 0
                        zoxx = global_var.ztot; % zox is the bottom depth
                    elseif (iz==1)      % calculating zox interpolating at z=0 with SWI conc. and at z=z(iz) with conc. o2x(iz)
                        zoxx = (global_var.z(iz)*bc.o2i*1d-6/1d3 + 0d0*abs(o2x(iz)))/(bc.o2i*1d-6/1d3+abs(o2x(iz)));
                    else     % calculating zox interpolating at z=z(iz-1) with o2x(iz-1) and at z=z(iz) with conc. o2x(iz)
                        zoxx = (global_var.z(iz)*o2x(iz-1) + global_var.z(iz-1)*abs(o2x(iz)))/(o2x(iz-1)+abs(o2x(iz)));
                    end
                    
                    % error evaluation as relative difference of zox
                    error = abs((global_var.zox -zoxx)/global_var.zox);
                    
                    %                    fprintf( 'itr,zox, zoxx, error %i \t %17.16e \t %17.16e \t %17.16e \n',itr, global_var.zox, zoxx, error);
                    %                    fprintf('~~~~~~~~~~~////~~~~~~~~~~~~~ \n');
                    
                    if (global_var.zox==zoxx)
                        break   % exit in fortran (terminates for-loop)
                    end
                    global_var.zox = 0.5d0*(global_var.zox + zoxx);  % new zox
                    
                    % if iteration reaches 100, error in zox is tested assuming individual grid depths as zox and find where error gets minimized
                    if (itr>=100 && itr <= global_var.nz+99)
                        global_var.zox = global_var.z(itr-99); % zox value in next test
                        if (minerr >=error )	 % if this time error is less than last adopt as optimum
                            if (itr~=100)
                                izox_minerr = itr -100;
                                minerr = error;
                            end
                        end
                    elseif (itr == (global_var.nz+100))    % check last test z(nz)
                        if (minerr >=error )
                            izox_minerr = itr -100;
                            minerr = error;
                        end
                        global_var.zox = z(izox_minerr);  % determine next test which should be most optimum
                    elseif (itr == (global_var.nz+101))  % results should be optimum and thus exit
                        break   % exit in fortran (terminates for-loop)
                    end
                    
                    if (itr > (global_var.nz+101))
                        % in fortran stop
                        msg = 'Error: (itr > (nz+101)), STOP.';
                        error(msg)
                    end
                    
                    itr = itr + 1;
                end
                
                %                     % showing results on screen
                omx_wtpc = omx.*global_var.mom./rho*100d0;
                fprintf('~~~~ conc ~~~~ itr = %i \n', itr);
                fprintf('z   ,OM   , o2\n');
                % showing parameters relevant to burial on screen
                for iz=1:interval:global_var.nz
                    fprintf('%17.16e \t %17.16e \t %17.16e \n', global_var.z(iz), omx_wtpc(iz), o2x(iz)*1d3);
                end
                fprintf('++++ flx ++++ \n');
                fprintf(' \t tflx \t adv \t dif \t omrxn \t ccrxn \t rain \t res \n')
                fprintf('OM :  \t %17.16e \t %17.16e \t %17.16e \t  %17.16e \t %17.16e \t %17.16e \t %17.16e \n', omtflx, omadv,  omdif, omdec,0d0,omrain, omres);
                fprintf('O2 :  \t %17.16e \t %17.16e \t %17.16e \t  %17.16e \t %17.16e \t %17.16e \t %17.16e', o2tflx, 0d0, o2dif, o2dec, 0d0, 0d0, o2res);
                
                fprintf('\n');
                fprintf('\n');
                fprintf('\n');
                %
                bc.om = omx;
                bc.o2 = o2x;
                time = time +dt;
                %
                
            end
            
        end
        
        
        function chk_cc(tmp_in, sal_in, dep_in)
            % just check the calculation of conc. and flx of caco3, dic and alk using subroutines in caco3_test_mod_v5_6.f90
            
            % you can check conc. and flx of caco3, alk, dic
            % please copy and paste results to whatever file to be compared with results with MATLAB version
            % NOTE: here decomposition of om is ignored. The test including decomposition of om is made in another file
            % because burial is not adjusted, cc conc. can goes > 100 wt% in this test
            
            tmp = tmp_in;   % 2.0;
            sal = sal_in;   % 35.0;
            dep = dep_in;   % 0.0d0; % depth in km
            %            dep = 4.0d0;    % km water depth; note that temperature and salinitiy has initially assumed values in globalvariables.mod
            
            flg_500 = false;    % error in calculation?
            
            interval =10; % choose a value between 1 to nz; om depth profile is shown with this interval; e.g., if inteval = nz, om conc. at all depths are shown
            % e.g., if interval = 5, om conc. at 5 depths are shown
            
            % initialize/define global properties/variables
            global_var = caco3_main;
            
            % initialize
            dw = zeros(1, global_var.nz);       % burial rate dummy
            
            % initialize mixing type
            mix_type = caco3_test.caco3_set_mixing(global_var);
            
            % initialize the boundary conditions
            bc = caco3_test.caco3_set_boundary_cond(global_var);
            
            
            beta = 1.00000000005d0;  % a parameter to make a grid; closer to 1, grid space is more concentrated around the sediment-water interface (SWI)
            [global_var.dz,global_var.z] = caco3_main.makegrid(beta,global_var.nz, global_var.ztot, global_var.def_recgrid);
            
            
            %            global_var = caco3_main.getporosity(global_var.z, global_var); % assume porosity profile
            [global_var.poro, global_var.porof, global_var.sporof, global_var.sporo, global_var.sporoi] = caco3_main.getporosity(global_var.z, global_var.poroi, global_var.nz);      % assume porosity profile
            %%%%%%%%%%%%% flx assignement and initial guess for burial rate %%%%%%%%%%%%%%%%%%%%%%
            % assume fluxes of om, cc and clay, required to calculate burial velocity
            [omflx, detflx, ccflx] = caco3_main.flxstat(global_var.om2cc, bc.ccflxi, global_var.mcc, global_var.nspcc);
            %            fprintf('ccflx %17.16e \n', ccflx);
            %            fprintf('om2cc, ccflxi, detflx, omflx, sum(ccflx) \n');
            %            fprintf('%17.16e %17.16e %17.16e %17.16e %17.16e \n', global_var.om2cc, bc.ccflxi, detflx, omflx, sum(ccflx));
            %            fprintf(' \n');
            
            % molar volume (cm3 mol-1) needed for burial rate calculation
            mvom = global_var.mom/global_var.rhoom;  % om
            mvsed = global_var.msed/global_var.rhosed; % clay
            mvcc = global_var.mcc/global_var.rhocc; % caco3
            
            % initial guess of burial profile, requiring porosity profile
            % w = burial rate, wi = burial rate initial guess
            [w, wi] = caco3_main.burial_pre(detflx,ccflx,global_var.msed,mvsed,mvcc,global_var.poroi, global_var.nz);
            
            % % depth -age conversion
            age = caco3_main.dep2age(global_var.dz, w, global_var.nz);
            
            % % determine factors for upwind scheme to represent burial advection
            [up, dwn, cnr, adf] = caco3_main.calcupwindscheme(w, global_var.nz);
            
            
            % % make transition matrix
            [trans,izrec,izrec2,izml,mix_type.nonlocal] = ...
                caco3_main.make_transmx(mix_type.labs,global_var.nspcc,mix_type.turbo2,mix_type.nobio,global_var.dz,global_var.sporo,global_var.nz,global_var.z, global_var.zml_ref, global_var.def_size);
            
            [keq1    ,keq2   ,keqcc  ,co3sat, dif_dic    ,dif_alk    ,dif_o2 ,kom, kcc, global_var] = caco3_main.coefs(tmp,sal,dep, global_var);
            
            %%   INITIAL CONDITIONS %%%%%%%%%%%%%%%%%%%
            cc = 1d-8*ones(global_var.nz, global_var.nspcc);     	% assume an arbitrary low conc.
            dic = bc.dici*1d-6/1d3*ones(1,global_var.nz);           % mol/cm3; factor is added to change uM to mol/cm3
            alk = bc.alki*1d-6/1d3*ones(1,global_var.nz);     % mol/cm3
            
            % calling subroutine to calculate all aqueous co2 species and pH
            [ph,co2,hco3,co3,info] = caco3_therm.calcspecies(dic,alk,tmp,sal,dep);
            
            %  passing to transient variables, subscript x denotes dummy variable used during iteration
            ccx = cc;       %  concentration of caco3, mol cm-3 sld;
            dicx = dic;     % mol cm-3 porewater; dic, x denotes dummy variables
            alkx = alk;     % mol cm-3 porewater; alk, x denotes dummy variables
            
            % this may not be necessary as these individual species assume equilibrium
            co2x = co2;
            hco3x = hco3;
            co3x = co3;
            
            time = 0d0; % model time [yr]
            it = 1; % integration count
            nt = 10; % total integration
            dt = 1d2; % time step [yr]
            
            rho = 2.5d0*ones(1, global_var.nz); % assume here density (this is going to be calculated based on solid phase composition )
            
            oxco2 = zeros(1, global_var.nz);  % oxic degradation of om; here assumed 0 at all time and depth
            anco2 = zeros(1, global_var.nz);   % anoxic degradation of om; here assumed 0 at all time and depth
            
            
            rcc = zeros(global_var.nz, global_var.nspcc);         % dissolution rate of caco3
            
            for it=1:nt
                %    print'(A,i0,A,E11.3,A,E11.3,A)','(it,dt,time)  (',it,',',dt,',',time,')'
                fprintf('(it,dt,time)  %2.2i %11.3e %11.3e\n', it, dt, time);
                
                % call calccaco3sys()
                [ccx,dicx,alkx,rcc,dt, flg_500, itr] = ...
                    caco3_main.calccaco3sys(ccx,dicx,alkx,rcc, dt, global_var.nspcc,dic,alk,dep,sal,tmp,mix_type.labs,mix_type.turbo2,mix_type.nonlocal, ...
                    global_var.sporo, global_var.sporoi, global_var.sporof, global_var.poro, dif_alk, dif_dic, ...
                    w, up, dwn, cnr, adf, global_var.dz, trans, cc, oxco2, anco2, co3sat, kcc, ccflx, global_var.ncc, global_var.nz, ...
                    global_var.tol, global_var.poroi, flg_500, global_var.fact, bc.alki,bc.dici, global_var.ccx_th, global_var.def_nonrec, global_var.def_sparse, global_var.def_showiter, global_var.def_sense);
                % % %% Dom test
                % %                 if(it ==nt)
                % %                     fprintf('~~~~ cccx ~~~~ itr = %i \n', itr);
                % %                 	for isp=1:global_var.nz
                % %                         fprintf('%17.16e \t %17.16e \t %17.16e \t  %17.16e  \n', ccx(isp, 1), ccx(isp, 2), ccx(isp, 3), ccx(isp, 4) );
                % %                     end
                % %                 end
                
                if(flg_500)
                    msg = 'error after calccaco3sys, STOP.';
                    error(msg)
                end
                % calling subroutine to calculate all aqueous co2 species and pH
                % call calcspecies(dicx,alkx,temp,sal,dep,prox,co2x,hco3x,co3x,nz,infosbr)
                [prox,co2x,hco3x,co3x,infosbr] = caco3_therm.calcspecies(dicx,alkx,tmp,sal,dep);
                
                if (infosbr==1)
                    msg = 'error after calcspecies after calccaco3sys, STOP.';
                    error(msg)
                end
                
                %          % calculation of fluxes relevant to caco3 and co2 system
                %             % call calcflxcaco3sys()
                [cctflx,ccdis,ccdif,ccadv,ccrain,ccres,alktflx,alkdis,alkdif,alkdec,alkres, dictflx,dicdis,dicdif,dicres,dicdec, dw] = ...
                    caco3_main.calcflxcaco3sys(dw, global_var.nspcc, ccx, cc, ccflx,dt, global_var.dz, rcc, adf, up, dwn, cnr, w, dif_alk, dif_dic, dic, dicx, alk, alkx, oxco2, anco2, trans, ...
                    mix_type.turbo2, mix_type.labs,mix_type.nonlocal, global_var.sporof, it, global_var.nz, global_var.poro, global_var.sporo, ...
                    bc.dici,bc.alki, mvcc, global_var.tol);
                
                fprintf('~~~~ conc ~~~~ itr = %i \n', itr);
                fprintf('z, \t cc, \t dic, \t alk\n');
                % showing parameters relevant to burial on screen
                for iz=1:interval:global_var.nz
                    fprintf('%17.16e \t %17.16e \t %17.16e \t %17.16e \n', global_var.z(iz), sum(ccx(iz,:)*global_var.mcc)/rho(iz)*100d0, dicx(iz)*1d3, alkx(iz)*1d3);
                end
                fprintf('   ..... multiple cc species ..... \n');
                %    write(dumchr(2),'(i0)') interval
                %    dumchr(1)="(i0.3,':',"//trim(adjustl(dumchr(2)))//"E11.3"//")"
                for isp=1:global_var.nspcc
                    fprintf('\n');
                    fprintf('cc species: %i \n', isp);
                    %        print dumchr(1),isp,(ccx(iz,isp)*mcc/rho(iz)*100d0,iz=1,nz,nz/interval)
                    for iz=1:interval:global_var.nz
                        fprintf('%17.16e \n',ccx(iz,isp)*global_var.mcc/rho(iz)*100d0);
                    end
                end
                fprintf('++++ flx ++++ \n');
                fprintf(' \t tflx \t adv \t dif \t omrxn \t ccrxn \t rain \t res \n')
                fprintf('cc :  \t %17.16e \t %17.16e \t %17.16e \t  %17.16e \t %17.16e \t %17.16e \t %17.16e \n', sum(cctflx),  sum(ccadv), sum(ccdif),0d0,sum(ccdis), sum(ccrain), sum(ccres) );
                fprintf('dic :  \t %17.16e \t %17.16e \t %17.16e \t  %17.16e \t %17.16e \t %17.16e \t %17.16e \n', dictflx, 0d0,dicdif, dicdec,  dicdis, 0d0,dicres );
                fprintf('alk :  \t %17.16e \t %17.16e \t %17.16e \t  %17.16e \t %17.16e \t %17.16e \t %17.16e \n', alktflx, 0d0, alkdif, alkdec, alkdis, 0d0, alkres );
                fprintf('   ..... multiple cc species ..... \n');
                for isp=1:global_var.nspcc
                    fprintf('%i \t %17.16e \t %17.16e \t %17.16e \t  %17.16e \t %17.16e \t %17.16e \t %17.16e \n', isp,cctflx(isp), ccadv(isp), ccdif(isp),0d0,ccdis(isp), ccrain(isp), ccres(isp) );
                end
                
                fprintf('\n');
                fprintf('\n');
                fprintf('\n');
                
                
                
                cc = ccx;
                dic = dicx;
                alk = alkx;
                time = time +dt;
                
                
            end
            
        end
        
        function chk_clay(tmp_in, sal_in, dep_in)
            % check the calculation of conc. and flx of clay using subroutines in caco3_test_mod_v5_6.f90
            
            % you can check conc. and flx of caco3, alk, dic
            % please copy and paste results to whatever file to be compared with results with MATLAB version
            % NOTE: here decomposition of om is ignored. The test including decomposition of om is made in another file
            % because burial is not adjusted, cc conc. can goes > 100 wt% in this test
            
            tmp = tmp_in;   % 2.0;
            sal = sal_in;   % 35.0;
            dep = dep_in;   % 0.0d0; % depth in km
            %            dep = 4.0d0;    % km water depth; note that temperature and salinitiy has initially assumed values in globalvariables.mod
            
            flg_500 = false;    % error in calculation?
            
            interval =10; % choose a value between 1 to nz; om depth profile is shown with this interval; e.g., if inteval = nz, om conc. at all depths are shown
            % e.g., if interval = 5, om conc. at 5 depths are shown
            
            % initialize/define global properties/variables
            global_var = caco3_main;
            
            % initialize
            dw = zeros(1, global_var.nz);       % burial rate dummy
            
            % initialize mixing type
            mix_type = caco3_test.caco3_set_mixing(global_var);
            
            % initialize the boundary conditions
            bc = caco3_test.caco3_set_boundary_cond(global_var);
            
            
            beta = 1.00000000005d0;  % a parameter to make a grid; closer to 1, grid space is more concentrated around the sediment-water interface (SWI)
            [global_var.dz,global_var.z] = caco3_main.makegrid(beta,global_var.nz, global_var.ztot, global_var.def_recgrid);
            
            
            %            global_var = caco3_main.getporosity(global_var.z, global_var); % assume porosity profile
            [global_var.poro, global_var.porof, global_var.sporof, global_var.sporo, global_var.sporoi] = caco3_main.getporosity(global_var.z, global_var.poroi, global_var.nz);      % assume porosity profile
            %%%%%%%%%%%%% flx assignement and initial guess for burial rate %%%%%%%%%%%%%%%%%%%%%%
            % assume fluxes of om, cc and clay, required to calculate burial velocity
            [omflx, detflx, ccflx] = caco3_main.flxstat(global_var.om2cc, bc.ccflxi, global_var.mcc, global_var.nspcc);
            %            fprintf('ccflx %17.16e \n', ccflx);
            %            fprintf('om2cc, ccflxi, detflx, omflx, sum(ccflx) \n');
            %            fprintf('%17.16e %17.16e %17.16e %17.16e %17.16e \n', global_var.om2cc, bc.ccflxi, detflx, omflx, sum(ccflx));
            %            fprintf(' \n');
            
            % molar volume (cm3 mol-1) needed for burial rate calculation
            mvom = global_var.mom/global_var.rhoom;  % om
            mvsed = global_var.msed/global_var.rhosed; % clay
            mvcc = global_var.mcc/global_var.rhocc; % caco3
            
            % initial guess of burial profile, requiring porosity profile
            % w = burial rate, wi = burial rate initial guess
            [w, wi] = caco3_main.burial_pre(detflx,ccflx,global_var.msed,mvsed,mvcc,global_var.poroi, global_var.nz);
            
            % % depth -age conversion
            age = caco3_main.dep2age(global_var.dz, w, global_var.nz);
            
            % % determine factors for upwind scheme to represent burial advection
            [up, dwn, cnr, adf] = caco3_main.calcupwindscheme(w, global_var.nz);
            
            
            % % make transition matrix
            [trans,izrec,izrec2,izml,mix_type.nonlocal] = ...
                caco3_main.make_transmx(mix_type.labs,global_var.nspcc,mix_type.turbo2,mix_type.nobio,global_var.dz,global_var.sporo,global_var.nz,global_var.z, global_var.zml_ref, global_var.def_size);
            
            pt = 1d-8*ones(1, global_var.nz);  % mol cm-3 sld; clay conc., assume an arbitrary low conc.
            
            % passing to transient variables
            ptx = pt;
            
            time = 0d0; % model time [yr]
            it = 1; % integration count
            nt = 10; % total integration
            dt = 1d6; % time step [yr]
            
            rho = 2.5d0*ones(1, global_var.nz); % assume here density (this is going to be calculated based on solid phase composition )
            
            for it=1:nt
                %    print'(A,i0,A,E11.3,A,E11.3,A)','(it,dt,time)  (',it,',',dt,',',time,')'
                fprintf('(it,dt,time)  %2.2i %11.3e %11.3e\n', it, dt, time);
                
                % call claycalc()
                ptx = caco3_main.claycalc(global_var.nz,global_var.sporo,pt,dt,w,global_var.dz,detflx,adf,up,dwn,cnr,trans, ...
                    global_var.nspcc,mix_type.labs,mix_type.turbo2,mix_type.nonlocal,global_var.poro,global_var.sporof, global_var.msed, global_var.def_nonrec );
                
                %               call calcflxclay()
                [pttflx,ptdif,ptadv,ptres,ptrain, dw] = caco3_main.calcflxclay(dw, global_var.nz,global_var.sporo, ptx, pt, dt, global_var.dz, detflx, w, adf, up, dwn, cnr, global_var.sporof, ...
                    trans, mix_type.turbo2, mix_type.labs, mix_type.nonlocal, global_var.poro, global_var.msed, mvsed);
                
                fprintf('~~~~ conc ~~~~\n');
                fprintf('z, \t sed \n');
                % showing parameters relevant to burial on screen
                for iz=1:interval:global_var.nz
                    fprintf('%17.16e \t %17.16e \n', global_var.z(iz), ptx(iz)*global_var.msed/rho(iz)*100d0);
                end
                
                fprintf('++++ flx ++++ \n');
                fprintf(' \t tflx \t adv \t dif \t omrxn \t ccrxn \t rain \t res \n')
                fprintf('sed :  \t %17.16e \t %17.16e \t %17.16e \t  %17.16e \t %17.16e \t %17.16e \t %17.16e \n', pttflx, ptadv,ptdif,  0d0, 0d0, ptrain, ptres );
                
                %                 write(dumchr(2),'(i0)') interval
                %                 dumchr(1)="(A,"//trim(adjustl(dumchr(2)))//"E11.3"//")"
                %                 % showing results on screen
                %                 print*,'~~~~ conc ~~~~'
                %                 print dumchr(1), 'z  :',(z(iz),iz=1,nz,nz/interval)
                %                 print dumchr(1), 'sed:',(ptx(iz)*msed/rho(iz)*100d0,iz=1,nz,nz/interval)
                %                 print*,'++++ flx ++++'
                %                 print'(7A11)', 'tflx','adv','dif','omrxn','ccrxn','rain','res'
                %                 print'(A,7E11.3)', 'sed:',pttflx, ptadv,ptdif,  0d0, 0d0, ptrain, ptres
                %
                %                 print*,''
                %                 print*,''
                %                 print*,''
                
                pt = ptx;
                time = time +dt;
                
            end %do
        end
        
        function chk_om_o2_cc_clay(tmp_in, sal_in, dep_in)
            % % check the calculation of conc. and flx of om, o2, cc & clay using subroutines in caco3_test_mod_v5_6.f90
            
            
            % <<<< Maybe first only consider Fickian mixing so switch off all macros in defines.h (you can switch on 'test') >>>>>
            
            % you can check calculation of concs. and flxes of om, o2, cc & clay
            % please copy and paste results to whatever file to be compared with results with MATLAB version
            % NOTE: Burial is not modified according to reactions and non-local mixing so that solid conc. may go above 100 wt %.
            
            interval =10; % choose a value between 1 to nz; om depth profile is shown with this interval; e.g., if inteval = nz, om conc. at all depths are shown
            % e.g., if interval = 5, om conc. at 5 depths are shown
            
            % initialize/define global properties/variables
            global_var = caco3_main;
            
            
            % initialize mixing type
            mix_type = caco3_test.caco3_set_mixing(global_var);
            
            % initialize the boundary conditions
            bc = caco3_test.caco3_set_boundary_cond(global_var);
            
            dw = zeros(1, global_var.nz);                       % burial rate change
            rcc = zeros(global_var.nz, global_var.nspcc);       % dissolution rate of caco3
            
            tmp = tmp_in;   % 2.0;
            sal = sal_in;   % 35.0;
            dep = dep_in;   % 0.0d0; % depth in km
            
            flg_500 = false;    % error in calculation?
            
            
            beta = 1.00000000005d0;  % a parameter to make a grid; closer to 1, grid space is more concentrated around the sediment-water interface (SWI)
            [global_var.dz,global_var.z] = caco3_main.makegrid(beta,global_var.nz, global_var.ztot, global_var.def_recgrid);
            
            
            %            global_var = caco3_main.getporosity(global_var.z, global_var); % assume porosity profile
            [global_var.poro, global_var.porof, global_var.sporof, global_var.sporo, global_var.sporoi] = caco3_main.getporosity(global_var.z, global_var.poroi, global_var.nz);      % assume porosity profile
            %%%%%%%%%%%%% flx assignement and initial guess for burial rate %%%%%%%%%%%%%%%%%%%%%%
            % assume fluxes of om, cc and clay, required to calculate burial velocity
            [omflx, detflx, ccflx] = caco3_main.flxstat(global_var.om2cc, bc.ccflxi, global_var.mcc, global_var.nspcc);
            %            fprintf('ccflx %17.16e \n', ccflx);
            %            fprintf('om2cc, ccflxi, detflx, omflx, sum(ccflx) \n');
            %            fprintf('%17.16e %17.16e %17.16e %17.16e %17.16e \n', global_var.om2cc, bc.ccflxi, detflx, omflx, sum(ccflx));
            %            fprintf(' \n');
            
            % molar volume (cm3 mol-1) needed for burial rate calculation
            mvom = global_var.mom/global_var.rhoom;  % om
            mvsed = global_var.msed/global_var.rhosed; % clay
            mvcc = global_var.mcc/global_var.rhocc; % caco3
            
            % initial guess of burial profile, requiring porosity profile
            % w = burial rate, wi = burial rate initial guess
            [w, wi] = caco3_main.burial_pre(detflx,ccflx,global_var.msed,mvsed,mvcc,global_var.poroi, global_var.nz);
            
            % % depth -age conversion
            age = caco3_main.dep2age(global_var.dz, w, global_var.nz);
            
            % % determine factors for upwind scheme to represent burial advection
            [up, dwn, cnr, adf] = caco3_main.calcupwindscheme(w, global_var.nz);
            
            
            % % make transition matrix
            [trans,izrec,izrec2,izml,mix_type.nonlocal] = ...
                caco3_main.make_transmx(mix_type.labs,global_var.nspcc,mix_type.turbo2,mix_type.nobio,global_var.dz,global_var.sporo,global_var.nz,global_var.z, global_var.zml_ref, global_var.def_size);
            
            [keq1    ,keq2   ,keqcc  ,co3sat, dif_dic    ,dif_alk    ,dif_o2 ,kom, kcc, global_var] = caco3_main.coefs(tmp,sal,dep, global_var);
            
            %   INITIAL CONDITIONS %
            bc.o2 = bc.o2i*1d-6/1d3 * ones(1, global_var.nz);	% o2 conc. in uM converted to mol/cm3
            cc = 1d-8 * ones(global_var.nz, global_var.nspcc);	% mol cm-3 sld; concentration of caco3, assume an arbitrary low conc.
            dic = bc.dici*1d-6/1d3 * ones(1, global_var.nz);  	% mol/cm3; factor is added to change uM to mol/cm3
            alk = bc.alki*1d-6/1d3 * ones(1, global_var.nz);   	% mol/cm3
            pt = 1d-8*ones(1, global_var.nz);                   % mol cm-3 sld; clay conc., assume an arbitrary low conc.
            
            % calling subroutine to calculate all aqueous co2 species and pH
            [pro,co2,hco3,co3,infosbr] = caco3_therm.calcspecies(dic,alk,tmp,sal,dep);
            
            omx = bc.om;
            o2x = bc.o2;
            ccx = cc;
            dicx = dic;
            alkx = alk ;
            % this may not be necessary as these individual species assume equilibrium
            co2x = co2;
            hco3x = hco3;
            co3x = co3;
            ptx = pt;
            
            time = 0d0; % model time [yr]
            it = 1; % integration count
            nt = 10; % total integration
            dt = 1d6; % time step [yr]
            
            rho = 2.5d0*ones(global_var.nz, 1); % assume here density (this is going to be calculated based on solid phase composition )
            
            oxco2 = zeros(1, global_var.nz);  % oxic degradation of om; here assumed 0 at all time and depth
            anco2 = zeros(1, global_var.nz);   % anoxic degradation of om; here assumed 0 at all time and depth
            
            
            %%%  addition to chk_om.f90 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            global_var.zox = 10d0;  % initial assumption on oxygen penetaration depth [cm]
            %%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            
            for it=1:nt
                %    print'(A,i0,A,E11.3,A,E11.3,A)','(it,dt,time)  (',it,',',dt,',',time,')'
                fprintf('(it,dt,time)  %2.2i %11.3e %11.3e\n', it, dt, time);
                
                itr = 0;        % iteration number for om and o2 calcuation
                error = 1d4;  	% error in ieration for zox
                minerr= 1d4;    % recording minimum relative difference in zox from previously considered zox
                
                while(error > global_var.tol)
                    % %%%%%%%%%%%%%%%  om conc. calculation  % %%%%%%%%%%%%%%%
                    % omx: mol cm-3 sld; om conc.
                    % izox: integer for grid number of zox
                    % kom: degradation rate consts. for each nz grids
                    [omx, izox, kom] = ...
                        caco3_main.omcalc(bc.oxic, bc.anoxic,o2x,bc.om, global_var.komi,global_var.nz,global_var.sporo,global_var.sporoi,global_var.sporof, w, wi, dt, up, dwn, cnr, adf,trans, ...
                        global_var.nspcc, mix_type.labs,mix_type.turbo2, mix_type.nonlocal, omflx, global_var.poro, global_var.dz, global_var.o2th);
                    % calculating the fluxes relevant to om diagenesis (and checking the calculation satisfies the difference equations )
                    [omadv,omdec,omdif,omrain,omres,omtflx] = ...
                        caco3_main.calcflxom(omflx,global_var.sporo,bc.om,omx,dt,w,global_var.dz,global_var.z,global_var.nz,mix_type.turbo2,mix_type.labs, global_var.poro,up,dwn,cnr,adf,rho, global_var.mom,trans,kom,global_var.sporof);
                    
                    fprintf('izox = %i \n',izox);      % sb omcalc calculates izox, which is the deepest grid where o2 >=0.
                    
                    % %                     omx_wtpc = omx.*global_var.mom./rho*100d0;
                    % %                     fprintf('~~~~ conc ~~~~ \n');
                    % %                     fprintf('z   ,OM   \n');
                    % %                     % showing parameters relevant to burial on screen
                    % %                     for iz=1:interval:global_var.nz
                    % %                         fprintf('%17.16e \t %17.16e\n', global_var.z(iz), omx_wtpc(iz));
                    % %                     end
                    % %                     fprintf('++++ flx ++++ \n');
                    % %                     fprintf(' \t tflx \t adv \t dif \t omrxn \t ccrxn \t rain \t res \n')
                    % %                     fprintf('OM :  \t %17.16e \t %17.16e \t %17.16e \t  %17.16e \t %17.16e \t %17.16e \t %17.16e', omtflx, omadv,  omdif, omdec,0d0,omrain, omres);
                    % %     %                fprintf('O2 :  \t %17.16e \t %17.16e \t %17.16e \t  %17.16e \t %17.16e \t %17.16e \t %17.16e', o2tflx, 0d0, o2dif, o2dec, 0d0, 0d0, o2res);
                    
                    
                    if (izox == global_var.nz)     % fully oxic; lower boundary condition ---> no diffusive out flow
                        % o2 calculation when o2 penetration depth (zox) is the same as bottom depth.
                        o2x = caco3_main.o2calc_ox(izox,global_var.nz,global_var.poro,bc.o2,kom,omx,global_var.sporo,dif_o2,global_var.dz,dt, global_var.ox2om, bc.o2i);
                        %  fluxes relevant to o2 (at the same time checking the satisfaction of difference equations)
                        [o2dec,o2dif,o2tflx,o2res] = caco3_main.calcflxo2_ox(global_var.nz,global_var.sporo,kom,omx,global_var.dz,global_var.poro,dif_o2,dt,bc.o2,o2x, global_var.ox2om, bc.o2i);
                    else        %% if oxygen is depleted within calculation domain, lower boundary changes to zero concs.
                        % o2 calculation when o2 is depleted within the calculation domain.
                        o2x = caco3_main.o2calc_sbox(izox,global_var.nz,global_var.poro,bc.o2,kom,omx,global_var.sporo,dif_o2,global_var.dz,dt, global_var.ox2om, bc.o2i);
                        % fluxes relevant to oxygen
                        [o2dec,o2dif,o2tflx,o2res] = caco3_main.calcflxo2_sbox(global_var.nz,global_var.sporo,kom,omx,global_var.dz,global_var.poro,dif_o2,dt,bc.o2,o2x,izox, global_var.ox2om, bc.o2i);
                    end
                    
                    %                 % showing intermediate results on screen
                    %                 omx_wtpc = omx.*global_var.mom./rho*100d0;
                    %                 fprintf('~~~~ conc ~~~~ itr = %i \n', itr +1);
                    %                 fprintf('z   ,OM   , o2\n');
                    %                 % showing parameters relevant to burial on screen
                    %                 for iz=1:interval:global_var.nz
                    %                     fprintf('%17.16e \t %17.16e \t %17.16e \n', global_var.z(iz), omx_wtpc(iz), o2x(iz)*1d3);
                    %                 end
                    %                 fprintf('++++ flx ++++ \n');
                    %                 fprintf(' \t tflx \t adv \t dif \t omrxn \t ccrxn \t rain \t res \n')
                    %                 fprintf('OM :  \t %17.16e \t %17.16e \t %17.16e \t  %17.16e \t %17.16e \t %17.16e \t %17.16e \n', omtflx, omadv,  omdif, omdec,0d0,omrain, omres);
                    %                 fprintf('O2 :  \t %17.16e \t %17.16e \t %17.16e \t  %17.16e \t %17.16e \t %17.16e \t %17.16e', o2tflx, 0d0, o2dif, o2dec, 0d0, 0d0, o2res);%
                    %                 fprintf('\n');
                    %                 fprintf('\n');
                    
                    
                    % update of zox
                    zoxx = 0d0;      % zox dummy variable
                    iz=1;
                    %                    for iz=1:global_var.nz
                    while(iz<global_var.nz+1)
                        if (o2x(iz)<=0d0)
                            break   % exit in fortran (terminates for-loop)
                        end
                        iz = iz+1;
                    end
                    
                    if (iz==global_var.nz+1) % oxygen never gets less than 0
                        zoxx = global_var.ztot; % zox is the bottom depth
                    elseif (iz==1)      % calculating zox interpolating at z=0 with SWI conc. and at z=z(iz) with conc. o2x(iz)
                        zoxx = (global_var.z(iz)*bc.o2i*1d-6/1d3 + 0d0*abs(o2x(iz)))/(bc.o2i*1d-6/1d3+abs(o2x(iz)));
                    else     % calculating zox interpolating at z=z(iz-1) with o2x(iz-1) and at z=z(iz) with conc. o2x(iz)
                        zoxx = (global_var.z(iz)*o2x(iz-1) + global_var.z(iz-1)*abs(o2x(iz)))/(o2x(iz-1)+abs(o2x(iz)));
                    end
                    
                    % error evaluation as relative difference of zox
                    error = abs((global_var.zox -zoxx)/global_var.zox);
                    
                    %                    fprintf( 'itr,zox, zoxx, error %i \t %17.16e \t %17.16e \t %17.16e \n',itr, global_var.zox, zoxx, error);
                    %                    fprintf('~~~~~~~~~~~////~~~~~~~~~~~~~ \n');
                    
                    if (global_var.zox==zoxx)
                        break   % exit in fortran (terminates for-loop)
                    end
                    global_var.zox = 0.5d0*(global_var.zox + zoxx);  % new zox
                    
                    % if iteration reaches 100, error in zox is tested assuming individual grid depths as zox and find where error gets minimized
                    if (itr>=100 && itr <= global_var.nz+99)
                        global_var.zox = global_var.z(itr-99); % zox value in next test
                        if (minerr >=error )	 % if this time error is less than last adopt as optimum
                            if (itr~=100)
                                izox_minerr = itr -100;
                                minerr = error;
                            end
                        end
                    elseif (itr == (global_var.nz+100))    % check last test z(nz)
                        if (minerr >=error )
                            izox_minerr = itr -100;
                            minerr = error;
                        end
                        global_var.zox = z(izox_minerr);  % determine next test which should be most optimum
                    elseif (itr == (global_var.nz+101))  % results should be optimum and thus exit
                        break   % exit in fortran (terminates for-loop)
                    end
                    
                    if (itr > (global_var.nz+101))
                        % in fortran stop
                        msg = 'Error: (itr > (nz+101)), STOP.';
                        error(msg)
                    end
                    
                    itr = itr + 1;
                end
                %~~  OM & O2 calculation END ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                
                % calculation of oxic and anoxic degradation of om (oxco2 and anco2, respectively)
                
                for iz = 1:global_var.nz
                    if (o2x(iz) > global_var.o2th)
                        oxco2(iz) = (1d0-global_var.poro(iz))*kom(iz)*omx(iz);  % aerobic respiration
                    else
                        % o2x(iz) = o2th
                        if (bc.anoxic)
                            anco2(iz) = (1d0-global_var.poro(iz))*kom(iz)*omx(iz);  % anaerobic respiration
                        end
                    end
                end
                
                for iz=1:global_var.nz
                    dw(iz) = dw(iz) -(1d0-global_var.poro(iz))*mvom*kom(iz)*omx(iz);  %% burial rate change need reflect volume change caused by chemical reactions
                    % as well as non-local mixing
                    if (mix_type.turbo2(1) || mix_type.labs(1))
                        for iiz = 1:global_var.nz
                            if (trans(iiz,iz,1)==0d0)
                                continue        % cycle in fortran
                            end
                            dw(iz) = dw(iz) - mvom*(-trans(iiz,iz,1)/global_var.dz(iz)*global_var.dz(iiz)*(1d0-global_var.poro(iiz))*omx(iiz));
                        end
                    else
                        if (mix_type.nonlocal(1))
                            for iiz = 1:global_var.nz
                                if (trans(iiz,iz,1)==0d0)
                                    continue        % cycle in fortran
                                end
                                dw(iz) = dw(iz) - mvom*(-trans(iiz,iz,1)/global_var.dz(iz)*omx(iiz));
                            end
                        end
                    end
                end
                
                
                for iz=1:global_var.nz
                    if (omx(iz)<global_var.omx_th)
                        omx(iz)=global_var.omx_th;  %% truncated at minimum value
                    end
                end
                
                %% calculation of caco3 system
                % call calccaco3sys()
                [ccx,dicx,alkx,rcc,dt, flg_500, itr] = ...
                    caco3_main.calccaco3sys(ccx,dicx,alkx,rcc, dt, global_var.nspcc,dic,alk,dep,sal,tmp,mix_type.labs,mix_type.turbo2,mix_type.nonlocal, ...
                    global_var.sporo, global_var.sporoi, global_var.sporof, global_var.poro, dif_alk, dif_dic, ...
                    w, up, dwn, cnr, adf, global_var.dz, trans, cc, oxco2, anco2, co3sat, kcc, ccflx, global_var.ncc, global_var.nz, ...
                    global_var.tol, global_var.poroi, flg_500, global_var.fact, bc.alki,bc.dici, global_var.ccx_th, global_var.def_nonrec, global_var.def_sparse, global_var.def_showiter, global_var.def_sense);
                
                if(flg_500)
                    msg = 'error after calccaco3sys, STOP.';
                    error(msg)
                end
                % calling subroutine to calculate all aqueous co2 species and pH
                % call calcspecies(dicx,alkx,temp,sal,dep,prox,co2x,hco3x,co3x,nz,infosbr)
                [prox,co2x,hco3x,co3x,infosbr] = caco3_therm.calcspecies(dicx,alkx,tmp,sal,dep);
                
                if (infosbr==1)
                    msg = 'error after calcspecies after calccaco3sys, STOP.';
                    error(msg)
                end
                
                %          % calculation of fluxes relevant to caco3 and co2 system
                %             % call calcflxcaco3sys()
                [cctflx,ccdis,ccdif,ccadv,ccrain,ccres,alktflx,alkdis,alkdif,alkdec,alkres, dictflx,dicdis,dicdif,dicres,dicdec, dw] = ...
                    caco3_main.calcflxcaco3sys(dw, global_var.nspcc, ccx, cc, ccflx,dt, global_var.dz, rcc, adf, up, dwn, cnr, w, dif_alk, dif_dic, dic, dicx, alk, alkx, oxco2, anco2, trans, ...
                    mix_type.turbo2, mix_type.labs,mix_type.nonlocal, global_var.sporof, it, global_var.nz, global_var.poro, global_var.sporo, ...
                    bc.dici,bc.alki, mvcc, global_var.tol);
                
                %~~  caco3 calculations END ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                
                %% clay calculation
                % call claycalc()
                ptx = caco3_main.claycalc(global_var.nz,global_var.sporo,pt,dt,w,global_var.dz,detflx,adf,up,dwn,cnr,trans, ...
                    global_var.nspcc,mix_type.labs,mix_type.turbo2,mix_type.nonlocal,global_var.poro,global_var.sporof, global_var.msed, global_var.def_nonrec );
                
                %               call calcflxclay()
                [pttflx,ptdif,ptadv,ptres,ptrain, dw] = caco3_main.calcflxclay(dw, global_var.nz,global_var.sporo, ptx, pt, dt, global_var.dz, detflx, w, adf, up, dwn, cnr, global_var.sporof, ...
                    trans, mix_type.turbo2, mix_type.labs, mix_type.nonlocal, global_var.poro, global_var.msed, mvsed);
                
                % end of clay calculation
                
                %% showing results on screen
                
                omx_wtpc = omx.*global_var.mom./rho*100d0;
                fprintf('~~~~ conc ~~~~ itr = %i \n', itr);
                fprintf('z, \t   OM, \t o2, \t cc, \t dic, \t alk, \t sed\n');
                % showing parameters relevant to burial on screen
                for iz=1:interval:global_var.nz
                    fprintf('%17.16e \t %17.16e \t %17.16e \t %17.16e \t %17.16e \t %17.16e \t %17.16e \n', ...
                        global_var.z(iz), omx_wtpc(iz), o2x(iz)*1d3, sum(ccx(iz,:)*global_var.mcc)/rho(iz)*100d0, dicx(iz)*1d3, alkx(iz)*1d3, ptx(iz)*global_var.msed/rho(iz)*100d0);
                end
                
                fprintf('   ..... multiple cc species ..... \n');
                %    write(dumchr(2),'(i0)') interval
                %    dumchr(1)="(i0.3,':',"//trim(adjustl(dumchr(2)))//"E11.3"//")"
                for isp=1:global_var.nspcc
                    fprintf('\n');
                    fprintf('cc species: %i \n', isp);
                    %        print dumchr(1),isp,(ccx(iz,isp)*mcc/rho(iz)*100d0,iz=1,nz,nz/interval)
                    for iz=1:interval:global_var.nz
                        fprintf('%17.16e \n',ccx(iz,isp)*global_var.mcc/rho(iz)*100d0);
                    end
                end
                
                fprintf('++++ flx ++++ \n');
                fprintf(' \t tflx \t adv \t dif \t omrxn \t ccrxn \t rain \t res \n')
                fprintf('om :  \t %17.16e \t %17.16e \t %17.16e \t  %17.16e \t %17.16e \t %17.16e \t %17.16e \n', omtflx, omadv,  omdif, omdec,0d0,omrain, omres );
                fprintf('o2 :  \t %17.16e \t %17.16e \t %17.16e \t  %17.16e \t %17.16e \t %17.16e \t %17.16e \n', o2tflx,0d0, o2dif,o2dec, 0d0,0d0,o2res );
                fprintf('cc :  \t %17.16e \t %17.16e \t %17.16e \t  %17.16e \t %17.16e \t %17.16e \t %17.16e \n', sum(cctflx),  sum(ccadv), sum(ccdif),0d0,sum(ccdis), sum(ccrain), sum(ccres) );
                fprintf('dic :  \t %17.16e \t %17.16e \t %17.16e \t  %17.16e \t %17.16e \t %17.16e \t %17.16e \n', dictflx, 0d0,dicdif, dicdec,  dicdis, 0d0,dicres );
                fprintf('alk :  \t %17.16e \t %17.16e \t %17.16e \t  %17.16e \t %17.16e \t %17.16e \t %17.16e \n', alktflx, 0d0, alkdif, alkdec, alkdis, 0d0, alkres );
                fprintf('sed :  \t %17.16e \t %17.16e \t %17.16e \t  %17.16e \t %17.16e \t %17.16e \t %17.16e \n', pttflx, ptadv,ptdif,  0d0, 0d0, ptrain, ptres );
                fprintf('   ..... multiple cc species ..... \n');
                for isp=1:global_var.nspcc
                    fprintf('%i \t %17.16e \t %17.16e \t %17.16e \t  %17.16e \t %17.16e \t %17.16e \t %17.16e \n', isp,cctflx(isp), ccadv(isp), ccdif(isp),0d0,ccdis(isp), ccrain(isp), ccres(isp) );
                end
                
                fprintf('\n');
                fprintf('\n');
                fprintf('\n');
                
                om = omx;
                o2 = o2x;
                cc = ccx;
                dic = dicx;
                alk = alkx;
                pt = ptx;
                time = time +dt;
                
            end
            
        end
        
        
        
        function chk_through(tmp_in, sal_in, dep_in)
            %% check the calculation of conc. and flx of om, o2, cc & clay using subroutines in caco3_test_mod_v5_6.f90, also save profiles in .txt files
            % burial rate is updated at a give time step until it converges well
            
            
            % <<<< Maybe first only consider Fickian mixing so switch off all macros in defines.h (you can switch on 'test') >>>>>
            
            % you can check calculation of concs. and flxes of om, o2, cc & clay, with burial modified at individual time steps
            % please copy and paste results to whatever file to be compared with results with MATLAB version
            % NOTE: Now this code is almost the same as the whole code. The difference is only that this code does not track any signals
            
            interval =10; % choose a value between 1 to nz; om depth profile is shown with this interval; e.g., if inteval = nz, om conc. at all depths are shown
            % e.g., if interval = 5, om conc. at 5 depths are shown
            
            % initialize/define global properties/variables
            global_var = caco3_main;
            
            
            % initialize mixing type
            mix_type = caco3_test.caco3_set_mixing(global_var);
            
            % initialize the boundary conditions
            bc = caco3_test.caco3_set_boundary_cond(global_var);
            
            dw = zeros(1, global_var.nz);                       % burial rate change
            rcc = zeros(global_var.nz, global_var.nspcc);       % dissolution rate of caco3
            
            tmp = tmp_in;   % 2.0;
            sal = sal_in;   % 35.0;
            dep = dep_in;   % 0.0d0; % depth in km
            
            flg_500 = false;    % error in calculation?
            
            
            beta = 1.00000000005d0;  % a parameter to make a grid; closer to 1, grid space is more concentrated around the sediment-water interface (SWI)
            [global_var.dz,global_var.z] = caco3_main.makegrid(beta,global_var.nz, global_var.ztot, global_var.def_recgrid);
            
            
            %            global_var = caco3_main.getporosity(global_var.z, global_var); % assume porosity profile
            [global_var.poro, global_var.porof, global_var.sporof, global_var.sporo, global_var.sporoi] = caco3_main.getporosity(global_var.z, global_var.poroi, global_var.nz);      % assume porosity profile
            %%%%%%%%%%%%% flx assignement and initial guess for burial rate %%%%%%%%%%%%%%%%%%%%%%
            % assume fluxes of om, cc and clay, required to calculate burial velocity
            [omflx, detflx, ccflx] = caco3_main.flxstat(global_var.om2cc, bc.ccflxi, global_var.mcc, global_var.nspcc);
            %            fprintf('ccflx %17.16e \n', ccflx);
            %            fprintf('om2cc, ccflxi, detflx, omflx, sum(ccflx) \n');
            %            fprintf('%17.16e %17.16e %17.16e %17.16e %17.16e \n', global_var.om2cc, bc.ccflxi, detflx, omflx, sum(ccflx));
            %            fprintf(' \n');
            
            % molar volume (cm3 mol-1) needed for burial rate calculation
            mvom = global_var.mom/global_var.rhoom;  % om
            mvsed = global_var.msed/global_var.rhosed; % clay
            mvcc = global_var.mcc/global_var.rhocc; % caco3
            
            % initial guess of burial profile, requiring porosity profile
            % w = burial rate, wi = burial rate initial guess
            [w, wi] = caco3_main.burial_pre(detflx,ccflx,global_var.msed,mvsed,mvcc,global_var.poroi, global_var.nz);
            
            % % depth -age conversion
            age = caco3_main.dep2age(global_var.dz, w, global_var.nz);
            
            % % determine factors for upwind scheme to represent burial advection
            [up, dwn, cnr, adf] = caco3_main.calcupwindscheme(w, global_var.nz);
            
            
            % % make transition matrix
            [trans,izrec,izrec2,izml,mix_type.nonlocal] = ...
                caco3_main.make_transmx(mix_type.labs,global_var.nspcc,mix_type.turbo2,mix_type.nobio,global_var.dz,global_var.sporo,global_var.nz,global_var.z, global_var.zml_ref, global_var.def_size);
            
            [keq1    ,keq2   ,keqcc  ,co3sat, dif_dic    ,dif_alk    ,dif_o2 ,kom, kcc, global_var] = caco3_main.coefs(tmp,sal,dep, global_var);
            
            %   INITIAL CONDITIONS %
            bc.o2 = bc.o2i*1d-6/1d3 * ones(1, global_var.nz);	% o2 conc. in uM converted to mol/cm3
            cc = 1d-8 * ones(global_var.nz, global_var.nspcc);	% mol cm-3 sld; concentration of caco3, assume an arbitrary low conc.
            dic = bc.dici*1d-6/1d3 * ones(1, global_var.nz);  	% mol/cm3; factor is added to change uM to mol/cm3
            alk = bc.alki*1d-6/1d3 * ones(1, global_var.nz);   	% mol/cm3
            pt = 1d-8*ones(1, global_var.nz);                   % mol cm-3 sld; clay conc., assume an arbitrary low conc.
            
            % calling subroutine to calculate all aqueous co2 species and pH
            [pro,co2,hco3,co3,infosbr] = caco3_therm.calcspecies(dic,alk,tmp,sal,dep);
            
            omx = bc.om;
            o2x = bc.o2;
            ccx = cc;
            dicx = dic;
            alkx = alk ;
            ptx = pt;
            % this may not be necessary as these individual species assume equilibrium
            co2x = co2;
            hco3x = hco3;
            co3x = co3;
            
            time = 0d0; % model time [yr]
            it = 1; % integration count
            nt = 20; % total integration
            dt = 1d4; % time step [yr]
            
            rho = 2.5d0*ones(global_var.nz, 1); % assume here density (this is going to be calculated based on solid phase composition )
            
            oxco2 = zeros(1, global_var.nz);  % oxic degradation of om; here initially assumed 0
            anco2 = zeros(1, global_var.nz);   % anoxic degradation of om; here initially assumed 0
            
            
            %%%  addition to chk_om.f90 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            global_var.zox = 10d0;  % initial assumption on oxygen penetaration depth [cm]
            %%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        	help= 1;

            
            for it=1:nt
                
                itr_w = 0;          % # of iteration made to converge w
                err_w_min = 1d4;    % minimum relative difference of w compared to the previous w
                err_f = 0d0;        % relative different of total vol. fraction of solids wrt the previous value
                err_w = 10d0;           % err in burial rate
                %                300 continue % point of restart when burial velocity does not converge
                while(err_w > global_var.tol)
                    
                    fprintf('(it,dt,time)  %2.2i %11.3e %11.3e\n', it, dt, time);
                    
                    dw = zeros(1, global_var.nz);                       % burial rate change
                    
                    itr = 0;        % iteration number for om and o2 calcuation
                    error = 1d4;  	% error in ieration for zox
                    minerr= 1d4;    % recording minimum relative difference in zox from previously considered zox
                    
                    while(error > global_var.tol)
                        % %%%%%%%%%%%%%%%  om conc. calculation  % %%%%%%%%%%%%%%%
                        % omx: mol cm-3 sld; om conc.
                        % izox: integer for grid number of zox
                        % kom: degradation rate consts. for each nz grids
                        [omx, izox, kom] = ...
                            caco3_main.omcalc(bc.oxic, bc.anoxic,o2x,bc.om, global_var.komi,global_var.nz,global_var.sporo,global_var.sporoi,global_var.sporof, w, wi, dt, up, dwn, cnr, adf,trans, ...
                            global_var.nspcc, mix_type.labs,mix_type.turbo2, mix_type.nonlocal, omflx, global_var.poro, global_var.dz, global_var.o2th);
                        % calculating the fluxes relevant to om diagenesis (and checking the calculation satisfies the difference equations )
                        [omadv,omdec,omdif,omrain,omres,omtflx] = ...
                            caco3_main.calcflxom(omflx,global_var.sporo,bc.om,omx,dt,w,global_var.dz,global_var.z,global_var.nz,mix_type.turbo2,mix_type.labs, global_var.poro,up,dwn,cnr,adf,rho, global_var.mom,trans,kom,global_var.sporof);
                        
                        fprintf('izox = %i \n',izox);      % sb omcalc calculates izox, which is the deepest grid where o2 >=0.
                        
                        % %                     omx_wtpc = omx.*global_var.mom./rho*100d0;
                        % %                     fprintf('~~~~ conc ~~~~ \n');
                        % %                     fprintf('z   ,OM   \n');
                        % %                     % showing parameters relevant to burial on screen
                        % %                     for iz=1:interval:global_var.nz
                        % %                         fprintf('%17.16e \t %17.16e\n', global_var.z(iz), omx_wtpc(iz));
                        % %                     end
                        % %                     fprintf('++++ flx ++++ \n');
                        % %                     fprintf(' \t tflx \t adv \t dif \t omrxn \t ccrxn \t rain \t res \n')
                        % %                     fprintf('OM :  \t %17.16e \t %17.16e \t %17.16e \t  %17.16e \t %17.16e \t %17.16e \t %17.16e', omtflx, omadv,  omdif, omdec,0d0,omrain, omres);
                        % %     %                fprintf('O2 :  \t %17.16e \t %17.16e \t %17.16e \t  %17.16e \t %17.16e \t %17.16e \t %17.16e', o2tflx, 0d0, o2dif, o2dec, 0d0, 0d0, o2res);
                        
                        
                        if (izox == global_var.nz)     % fully oxic; lower boundary condition ---> no diffusive out flow
                            % o2 calculation when o2 penetration depth (zox) is the same as bottom depth.
                            o2x = caco3_main.o2calc_ox(izox,global_var.nz,global_var.poro,bc.o2,kom,omx,global_var.sporo,dif_o2,global_var.dz,dt, global_var.ox2om, bc.o2i);
                            %  fluxes relevant to o2 (at the same time checking the satisfaction of difference equations)
                            [o2dec,o2dif,o2tflx,o2res] = caco3_main.calcflxo2_ox(global_var.nz,global_var.sporo,kom,omx,global_var.dz,global_var.poro,dif_o2,dt,bc.o2,o2x, global_var.ox2om, bc.o2i);
                        else        %% if oxygen is depleted within calculation domain, lower boundary changes to zero concs.
                            % o2 calculation when o2 is depleted within the calculation domain.
                            o2x = caco3_main.o2calc_sbox(izox,global_var.nz,global_var.poro,bc.o2,kom,omx,global_var.sporo,dif_o2,global_var.dz,dt, global_var.ox2om, bc.o2i);
                            % fluxes relevant to oxygen
                            [o2dec,o2dif,o2tflx,o2res] = caco3_main.calcflxo2_sbox(global_var.nz,global_var.sporo,kom,omx,global_var.dz,global_var.poro,dif_o2,dt,bc.o2,o2x,izox, global_var.ox2om, bc.o2i);
                        end
                        
                        %                 % showing intermediate results on screen
                        %                 omx_wtpc = omx.*global_var.mom./rho*100d0;
                        %                 fprintf('~~~~ conc ~~~~ itr = %i \n', itr +1);
                        %                 fprintf('z   ,OM   , o2\n');
                        %                 % showing parameters relevant to burial on screen
                        %                 for iz=1:interval:global_var.nz
                        %                     fprintf('%17.16e \t %17.16e \t %17.16e \n', global_var.z(iz), omx_wtpc(iz), o2x(iz)*1d3);
                        %                 end
                        %                 fprintf('++++ flx ++++ \n');
                        %                 fprintf(' \t tflx \t adv \t dif \t omrxn \t ccrxn \t rain \t res \n')
                        %                 fprintf('OM :  \t %17.16e \t %17.16e \t %17.16e \t  %17.16e \t %17.16e \t %17.16e \t %17.16e \n', omtflx, omadv,  omdif, omdec,0d0,omrain, omres);
                        %                 fprintf('O2 :  \t %17.16e \t %17.16e \t %17.16e \t  %17.16e \t %17.16e \t %17.16e \t %17.16e', o2tflx, 0d0, o2dif, o2dec, 0d0, 0d0, o2res);%
                        %                 fprintf('\n');
                        %                 fprintf('\n');
                        
                        
                        % update of zox
                        zoxx = 0d0;      % zox dummy variable
                        iz=1;
                        %                    for iz=1:global_var.nz
                        while(iz<global_var.nz+1)
                            if (o2x(iz)<=0d0)
                                break   % exit in fortran (terminates for-loop)
                            end
                            iz = iz+1;
                        end
                        
                        if (iz==global_var.nz+1) % oxygen never gets less than 0
                            zoxx = global_var.ztot; % zox is the bottom depth
                        elseif (iz==1)      % calculating zox interpolating at z=0 with SWI conc. and at z=z(iz) with conc. o2x(iz)
                            zoxx = (global_var.z(iz)*bc.o2i*1d-6/1d3 + 0d0*abs(o2x(iz)))/(bc.o2i*1d-6/1d3+abs(o2x(iz)));
                        else     % calculating zox interpolating at z=z(iz-1) with o2x(iz-1) and at z=z(iz) with conc. o2x(iz)
                            zoxx = (global_var.z(iz)*o2x(iz-1) + global_var.z(iz-1)*abs(o2x(iz)))/(o2x(iz-1)+abs(o2x(iz)));
                        end
                        
                        % error evaluation as relative difference of zox
                        error = abs((global_var.zox -zoxx)/global_var.zox);
                        
                        %                    fprintf( 'itr,zox, zoxx, error %i \t %17.16e \t %17.16e \t %17.16e \n',itr, global_var.zox, zoxx, error);
                        %                    fprintf('~~~~~~~~~~~////~~~~~~~~~~~~~ \n');
                        
                        if (global_var.zox==zoxx)
                            break   % exit in fortran (terminates for-loop)
                        end
                        global_var.zox = 0.5d0*(global_var.zox + zoxx);  % new zox
                        
                        % if iteration reaches 100, error in zox is tested assuming individual grid depths as zox and find where error gets minimized
                        if (itr>=100 && itr <= global_var.nz+99)
                            global_var.zox = global_var.z(itr-99); % zox value in next test
                            if (minerr >=error )	 % if this time error is less than last adopt as optimum
                                if (itr~=100)
                                    izox_minerr = itr -100;
                                    minerr = error;
                                end
                            end
                        elseif (itr == (global_var.nz+100))    % check last test z(nz)
                            if (minerr >=error )
                                izox_minerr = itr -100;
                                minerr = error;
                            end
                            global_var.zox = z(izox_minerr);  % determine next test which should be most optimum
                        elseif (itr == (global_var.nz+101))  % results should be optimum and thus exit
                            break   % exit in fortran (terminates for-loop)
                        end
                        
                        if (itr > (global_var.nz+101))
                            % in fortran stop
                            msg = 'Error: (itr > (nz+101)), STOP.';
                            error(msg)
                        end
                        
                        itr = itr + 1;
                    end
                    %~~  OM & O2 calculation END ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                    
                    % calculation of oxic and anoxic degradation of om (oxco2 and anco2, respectively)
                    
                    for iz = 1:global_var.nz
                        if (o2x(iz) > global_var.o2th)
                            oxco2(iz) = (1d0-global_var.poro(iz))*kom(iz)*omx(iz);  % aerobic respiration
                        else
                            % o2x(iz) = o2th
                            if (bc.anoxic)
                                anco2(iz) = (1d0-global_var.poro(iz))*kom(iz)*omx(iz);  % anaerobic respiration
                            end
                        end
                    end
                    
                    for iz=1:global_var.nz
                        dw(iz) = dw(iz) -(1d0-global_var.poro(iz))*mvom*kom(iz)*omx(iz);  %% burial rate change need reflect volume change caused by chemical reactions
                        % as well as non-local mixing
                        if (mix_type.turbo2(1) || mix_type.labs(1))
                            for iiz = 1:global_var.nz
                                if (trans(iiz,iz,1)==0d0)
                                    continue        % cycle in fortran
                                end
                                dw(iz) = dw(iz) - mvom*(-trans(iiz,iz,1)/global_var.dz(iz)*global_var.dz(iiz)*(1d0-global_var.poro(iiz))*omx(iiz));
                            end
                        else
                            if (mix_type.nonlocal(1))
                                for iiz = 1:global_var.nz
                                    if (trans(iiz,iz,1)==0d0)
                                        continue        % cycle in fortran
                                    end
                                    dw(iz) = dw(iz) - mvom*(-trans(iiz,iz,1)/global_var.dz(iz)*omx(iiz));
                                end
                            end
                        end
                    end
                    
                    
                    for iz=1:global_var.nz
                        if (omx(iz)<global_var.omx_th)
                            omx(iz)=global_var.omx_th;  %% truncated at minimum value
                        end
                    end
                    
                    %% calculation of caco3 system
                    % call calccaco3sys()
                    [ccx,dicx,alkx,rcc,dt, flg_500, itr] = ...
                        caco3_main.calccaco3sys(ccx,dicx,alkx,rcc, dt, global_var.nspcc,dic,alk,dep,sal,tmp,mix_type.labs,mix_type.turbo2,mix_type.nonlocal, ...
                        global_var.sporo, global_var.sporoi, global_var.sporof, global_var.poro, dif_alk, dif_dic, ...
                        w, up, dwn, cnr, adf, global_var.dz, trans, cc, oxco2, anco2, co3sat, kcc, ccflx, global_var.ncc, global_var.nz, ...
                        global_var.tol, global_var.poroi, flg_500, global_var.fact, bc.alki,bc.dici, global_var.ccx_th, global_var.def_nonrec, global_var.def_sparse, global_var.def_showiter, global_var.def_sense);
                    
                    if(flg_500)
                        msg = 'error after calccaco3sys, STOP.';
                        error(msg)
                    end
                    % calling subroutine to calculate all aqueous co2 species and pH
                    % call calcspecies(dicx,alkx,temp,sal,dep,prox,co2x,hco3x,co3x,nz,infosbr)
                    [prox,co2x,hco3x,co3x,infosbr] = caco3_therm.calcspecies(dicx,alkx,tmp,sal,dep);
                    
                    if (infosbr==1)
                        msg = 'error after calcspecies after calccaco3sys, STOP.';
                        error(msg)
                    end
                    
                    %          % calculation of fluxes relevant to caco3 and co2 system
                    %             % call calcflxcaco3sys()
                    [cctflx,ccdis,ccdif,ccadv,ccrain,ccres,alktflx,alkdis,alkdif,alkdec,alkres, dictflx,dicdis,dicdif,dicres,dicdec, dw] = ...
                        caco3_main.calcflxcaco3sys(dw, global_var.nspcc, ccx, cc, ccflx,dt, global_var.dz, rcc, adf, up, dwn, cnr, w, dif_alk, dif_dic, dic, dicx, alk, alkx, oxco2, anco2, trans, ...
                        mix_type.turbo2, mix_type.labs,mix_type.nonlocal, global_var.sporof, it, global_var.nz, global_var.poro, global_var.sporo, ...
                        bc.dici,bc.alki, mvcc, global_var.tol);
                    
                    %~~  caco3 calculations END ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                    
                    %% clay calculation
                    % call claycalc()
                    ptx = caco3_main.claycalc(global_var.nz,global_var.sporo,pt,dt,w,global_var.dz,detflx,adf,up,dwn,cnr,trans, ...
                        global_var.nspcc,mix_type.labs,mix_type.turbo2,mix_type.nonlocal,global_var.poro,global_var.sporof, global_var.msed, global_var.def_nonrec );
                    
                    %               call calcflxclay()
                    [pttflx,ptdif,ptadv,ptres,ptrain, dw] = caco3_main.calcflxclay(dw, global_var.nz,global_var.sporo, ptx, pt, dt, global_var.dz, detflx, w, adf, up, dwn, cnr, global_var.sporof, ...
                        trans, mix_type.turbo2, mix_type.labs, mix_type.nonlocal, global_var.poro, global_var.msed, mvsed);
                    
                    % end of clay calculation
                    
                    %%%%%%%%%%%%%%
                    
                    % checking for total volume of solids, density and burial velocity
                    
                    % call getsldprop() % get solid property, rho (density) and frt (total vol.frac)
                    [rho, frt] = caco3_main.getsldprop(global_var.nz, omx, ptx, ccx, global_var.nspcc, w, up, dwn, cnr, adf, global_var.z, global_var.mom, global_var.msed, global_var.mcc, mvom, mvsed, mvcc);
                    
                    err_f = max(abs(frt - 1d0));  % new error in total vol. fraction (must be 1 in theory)
                    %% ========= calculation of burial velocity =============================
                    
                    wx = w;  % recording previous burial velocity
                    
                    % call burialcalc() % get new burial velocity
                    [w,wi] = caco3_main.burial_calc(detflx,ccflx,global_var.nspcc,omflx,dw,global_var.dz,global_var.poro,global_var.nz, global_var.msed,mvsed,mvcc,mvom,global_var.poroi);
                    
                    %  determine calculation scheme for advection / factors for upwind scheme to represent burial advection
                    % call calcupwindscheme()
                    [up, dwn, cnr, adf] = caco3_main.calcupwindscheme(w, global_var.nz);
                    
                    
                    % error and iteration evaluation
                    itr_w = itr_w + 1;              % counting iteration for w
                    err_w = max(abs((w-wx)/wx));    % relative difference of w
                    if (err_w<err_w_min)
                        err_w_min= err_w;  % recording minimum relative difference of  w
                        wxx = wx;  % recording w which minimizes deviation of total sld fraction from 1
                    end
                    
                    if (itr_w>100)   % if iteration gets too many (100), force to end with optimum w where error is minimum
                        if (itr_w==101)
                            w = wxx;
                            % %                             go to 300  % Dominik: do just go back to 300 if err_w > tol
                        elseif (itr_w==102)
                            w = wxx;
                            fprintf('not converging w %i \t %17.16e \t %17.16e \n',time, err_w, err_w_min);
                            pause;
                            break;
                            % %                             go to 400
                        end
                    end
                    % %                 if (err_w > tol) go to 300
                end
                % %
                % %                     400 continue
                % %
                fprintf('error in frt: %17.16e \n', max(abs(frt - 1d0)));
                
                
                
                %%%%%%%%%%%%%%
                
                %% showing results on screen
                
                omx_wtpc = omx.*global_var.mom./rho'*100d0;
                fprintf('~~~~ conc ~~~~ itr = %i \n', itr);
                fprintf('z, \t   OM, \t o2, \t cc, \t dic, \t alk, \t sed\n');
                % showing parameters relevant to burial on screen
                for iz=1:interval:global_var.nz
                    fprintf('%17.16e \t %17.16e \t %17.16e \t %17.16e \t %17.16e \t %17.16e \t %17.16e \n', ...
                        global_var.z(iz), omx_wtpc(iz), o2x(iz)*1d3, sum(ccx(iz,:)*global_var.mcc)/rho(iz)*100d0, dicx(iz)*1d3, alkx(iz)*1d3, ptx(iz)*global_var.msed/rho(iz)*100d0);
                end
                
                fprintf('   ..... multiple cc species ..... \n');
                %    write(dumchr(2),'(i0)') interval
                %    dumchr(1)="(i0.3,':',"//trim(adjustl(dumchr(2)))//"E11.3"//")"
                for isp=1:global_var.nspcc
                    fprintf('\n');
                    fprintf('cc species: %i \n', isp);
                    %        print dumchr(1),isp,(ccx(iz,isp)*mcc/rho(iz)*100d0,iz=1,nz,nz/interval)
                    for iz=1:interval:global_var.nz
                        fprintf('%17.16e \n',ccx(iz,isp)*global_var.mcc/rho(iz)*100d0);
                    end
                end
                
                fprintf('++++ flx ++++ \n');
                fprintf(' \t tflx \t adv \t dif \t omrxn \t ccrxn \t rain \t res \n')
                fprintf('om :  \t %17.16e \t %17.16e \t %17.16e \t  %17.16e \t %17.16e \t %17.16e \t %17.16e \n', omtflx, omadv,  omdif, omdec,0d0,omrain, omres );
                fprintf('o2 :  \t %17.16e \t %17.16e \t %17.16e \t  %17.16e \t %17.16e \t %17.16e \t %17.16e \n', o2tflx,0d0, o2dif,o2dec, 0d0,0d0,o2res );
                fprintf('cc :  \t %17.16e \t %17.16e \t %17.16e \t  %17.16e \t %17.16e \t %17.16e \t %17.16e \n', sum(cctflx),  sum(ccadv), sum(ccdif),0d0,sum(ccdis), sum(ccrain), sum(ccres) );
                fprintf('dic :  \t %17.16e \t %17.16e \t %17.16e \t  %17.16e \t %17.16e \t %17.16e \t %17.16e \n', dictflx, 0d0,dicdif, dicdec,  dicdis, 0d0,dicres );
                fprintf('alk :  \t %17.16e \t %17.16e \t %17.16e \t  %17.16e \t %17.16e \t %17.16e \t %17.16e \n', alktflx, 0d0, alkdif, alkdec, alkdis, 0d0, alkres );
                fprintf('sed :  \t %17.16e \t %17.16e \t %17.16e \t  %17.16e \t %17.16e \t %17.16e \t %17.16e \n', pttflx, ptadv,ptdif,  0d0, 0d0, ptrain, ptres );
                fprintf('   ..... multiple cc species ..... \n');
                for isp=1:global_var.nspcc
                    fprintf('%i \t %17.16e \t %17.16e \t %17.16e \t  %17.16e \t %17.16e \t %17.16e \t %17.16e \n', isp,cctflx(isp), ccadv(isp), ccdif(isp),0d0,ccdis(isp), ccrain(isp), ccres(isp) );
                end
                fprintf('==== burial etc ==== \n');
                fprintf('z, \t w, \t rho, \t frc \n');
                % showing parameters relevant to burial on screen
                for iz=1:interval:global_var.nz
                    fprintf('%17.16e \t %17.16e \t %17.16e \t %17.16e \n', ...
                        global_var.z(iz),w(iz), rho(iz), frt(iz));
                end
                
                fprintf('\n');
                fprintf('\n');
                fprintf('\n');
                
                om = omx;
                o2 = o2x;
                cc = ccx;
                dic = dicx;
                alk = alkx;
                pt = ptx;
                time = time +dt;
                
                d13c_ocni = 0d0;
                d18o_ocni = 0d0;
                d13c_blk = zeros(1, global_var.nz);
                d18o_blk = zeros(1, global_var.nz);
%                if(mod(it,800)==0)                    
                caco3_main.recordprofile(it, global_var.nz, global_var.z, age, pt, global_var.msed, wi, rho, cc, ccx, dic, dicx, alk, alkx, co3, co3x, co3sat ...
                    , rcc, pro, o2x, oxco2, anco2, om, global_var.mom, global_var.mcc, d13c_ocni, d18o_ocni, up,dwn, cnr, adf, global_var.nspcc, ptx, w, frt, prox, omx, d13c_blk, d18o_blk)
%                help = help+1;
%                end
            end
            
        end
        
        
        function chk_through_signal(tmp_in, sal_in, dep_in, dep_max)
            %% check the calculation of conc. and flx of om, o2, cc & clay using subroutines in caco3_test_mod_v5_6.f90, also save profiles in .txt files
            % burial rate is updated at a give time step until it converges well
            % dep_max: max depth to be changed to during the experiment
            
            % <<<< Maybe first only consider Fickian mixing so switch off all macros in defines.h (you can switch on 'test') >>>>>
            
            % you can check calculation of concs. and flxes of om, o2, cc & clay, with burial modified at individual time steps
            % please copy and paste results to whatever file to be compared with results with MATLAB version
            % NOTE: Now this code is almost the same as the whole code. The difference is only that this code does not track any signals
            
            interval =10; % choose a value between 1 to nz; om depth profile is shown with this interval; e.g., if inteval = nz, om conc. at all depths are shown
            % e.g., if interval = 5, om conc. at 5 depths are shown
            
            % initialize/define global properties/variables
            global_var = caco3_main;
            
            
            % initialize mixing type
            mix_type = caco3_test.caco3_set_mixing(global_var);
            
            % initialize the boundary conditions
            bc = caco3_test.caco3_set_boundary_cond(global_var);
            
            dw = zeros(1, global_var.nz);                       % burial rate change
            rcc = zeros(global_var.nz, global_var.nspcc);       % dissolution rate of caco3
            
            tmp = tmp_in;   % 2.0;
            sal = sal_in;   % 35.0;
            dep = dep_in;   % 0.0d0; % depth in km
            
            bc.ccflxi = 12d-6;      % mol (CaCO3) cm-2 yr-1 - caco3 flux
%            global_var.om2cc = 0.7d0;   % rain ratio of organic matter to calcite

            % open files for output signal at 3 different depths
            file_sigmly = sprintf('./recprofile_signaltrack_fickian/matlab_sigmly.txt');
            file_sigmlyid = fopen(file_sigmly,'wt');
            file_sigmlyd = sprintf('./recprofile_signaltrack_fickian/matlab_sigmlyd.txt');
            file_sigmlydid = fopen(file_sigmlyd,'wt');
            file_sigbtm = sprintf('./recprofile_signaltrack_fickian/matlab_sigbtm.txt');
            file_sigbtmid = fopen(file_sigbtm,'wt');
            file_bound = sprintf('./recprofile_signaltrack_fickian/matlab_bound.txt');
            file_boundid = fopen(file_bound,'wt');
            
            
            flg_500 = false;    % error in calculation?
            
            
            beta = 1.00000000005d0;     % a parameter to make a grid; closer to 1, grid space is more concentrated around the sediment-water interface (SWI)
            [global_var.dz,global_var.z] = caco3_main.makegrid(beta,global_var.nz, global_var.ztot, global_var.def_recgrid);
            
            
            %            global_var = caco3_main.getporosity(global_var.z, global_var); % assume porosity profile
            [global_var.poro, global_var.porof, global_var.sporof, global_var.sporo, global_var.sporoi] = caco3_main.getporosity(global_var.z, global_var.poroi, global_var.nz);      % assume porosity profile
            %%%%%%%%%%%%% flx assignement and initial guess for burial rate %%%%%%%%%%%%%%%%%%%%%%
            % assume fluxes of om, cc and clay, required to calculate burial velocity
            [omflx, detflx, ccflx] = caco3_main.flxstat(global_var.om2cc, bc.ccflxi, global_var.mcc, global_var.nspcc);
            %            fprintf('ccflx %17.16e \n', ccflx);
                        fprintf('om2cc, ccflxi, detflx, omflx, sum(ccflx) \n');
                        fprintf('%17.16e %17.16e %17.16e %17.16e %17.16e \n', global_var.om2cc, bc.ccflxi, detflx, omflx, sum(ccflx));
                        fprintf(' \n');
            
            % molar volume (cm3 mol-1) needed for burial rate calculation
            mvom = global_var.mom/global_var.rhoom;  % om
            mvsed = global_var.msed/global_var.rhosed; % clay
            mvcc = global_var.mcc/global_var.rhocc; % caco3
            
            % initial guess of burial profile, requiring porosity profile
            % w = burial rate, wi = burial rate initial guess
            [w, wi] = caco3_main.burial_pre(detflx,ccflx,global_var.msed,mvsed,mvcc,global_var.poroi, global_var.nz);
            
            % % depth -age conversion
            age = caco3_main.dep2age(global_var.dz, w, global_var.nz);
            
            % % determine factors for upwind scheme to represent burial advection
            [up, dwn, cnr, adf] = caco3_main.calcupwindscheme(w, global_var.nz);
            
            %%% ~~~~~~~~~~~~~~ set recording time
            % call recordtime()
            [rectime, cntrec, time_spn, time_trs, time_aft] = caco3_main.recordtime(global_var.nrec, wi, global_var.ztot, global_var.def_biotest, global_var.def_sense, global_var.def_nonrec);
            
            % water depth, i and f denote initial and final values
            depi = dep_in;  % depth before event
            depf = dep_max;   % max depth to be changed to
            
            % %  flux ratio of fine particles: i and f denote initial and final values
            flxfini = 0.5d0;  %  total caco3 rain flux for fine species assumed before event
            flxfinf = 0.9d0; %  maximum changed value
            
            % ///////////// isotopes  ////////////////
            % initial value of ocean d13c, final value of ocean d13c / d18o
            d13c_ocni = 2d0;  % initial ocean d13c value
            d13c_ocnf = -1d0; % ocean d13c value with maximum change
            d18o_ocni = 1d0; % initial ocean d18o value
            d18o_ocnf = -1d0; % ocean d18o value with maximum change
            % Dominik - initialize ocean d13c, d18o
            d13c_ocn = 0d0;
            d18o_ocn = 0d0;
            % end-member signal assignment
            % call sig2sp_pre()
            [d13c_sp,d18o_sp] = caco3_main.sig2sp_pre(d13c_ocni,d13c_ocnf,d18o_ocni,d18o_ocnf, global_var.def_sense, global_var.def_size, global_var.nspcc);
            
            
            
            % % make transition matrix
            [trans,izrec,izrec2,izml,mix_type.nonlocal] = ...
                caco3_main.make_transmx(mix_type.labs,global_var.nspcc,mix_type.turbo2,mix_type.nobio,global_var.dz,global_var.sporo,global_var.nz,global_var.z, global_var.zml_ref, global_var.def_size);
            
            [keq1    ,keq2   ,keqcc  ,co3sat, dif_dic    ,dif_alk    ,dif_o2 ,kom, kcc, global_var] = caco3_main.coefs(tmp,sal,dep, global_var);
            
            %   INITIAL CONDITIONS %
            bc.o2 = bc.o2i*1d-6/1d3 * ones(1, global_var.nz);	% o2 conc. in uM converted to mol/cm3
            cc = 1d-8 * ones(global_var.nz, global_var.nspcc);	% mol cm-3 sld; concentration of caco3, assume an arbitrary low conc.
            dic = bc.dici*1d-6/1d3 * ones(1, global_var.nz);  	% mol/cm3; factor is added to change uM to mol/cm3
            alk = bc.alki*1d-6/1d3 * ones(1, global_var.nz);   	% mol/cm3
            pt = 1d-8*ones(1, global_var.nz);                   % mol cm-3 sld; clay conc., assume an arbitrary low conc.
            
            % calling subroutine to calculate all aqueous co2 species and pH
            [pro,co2,hco3,co3,infosbr] = caco3_therm.calcspecies(dic,alk,tmp,sal,dep);
            
            omx = bc.om;
            o2x = bc.o2;
            ccx = cc;
            dicx = dic;
            alkx = alk ;
            ptx = pt;
            % this may not be necessary as these individual species assume equilibrium
            co2x = co2;
            hco3x = hco3;
            co3x = co3;
            
            time = 0d0; % model time [yr]
            int_count = 1; % integration count
            nt = 10; % total integration
            dt = 1d2; % time step [yr]
            
            rho = 2.5d0*ones(global_var.nz, 1); % assume here density (this is going to be calculated based on solid phase composition )
            
            oxco2 = zeros(1, global_var.nz);  % oxic degradation of om; here initially assumed 0
            anco2 = zeros(1, global_var.nz);   % anoxic degradation of om; here initially assumed 0
            
            
            %%%  addition to chk_om.f90 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            global_var.zox = 10d0;  % initial assumption on oxygen penetaration depth [cm]
            %%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            %            for int_count=1:nt
            while(2>1)              %   Here change to unknwon iteration
                
                %% ///////// isotopes & fluxes settings //////////////
                if(~global_var.def_sense)   % do signal tracking
                    
                    % determine time step dt by calling timestep(time,nt_spn,nt_trs,nt_aft,dt) where nt_xx denotes total iteration number
                    nt_spn = 800.0;     % timesteps for spin-up
                    nt_trs = 5000.0;    % timesteps close to & during event (signal transition)
                    nt_aft = 1000.0;    % timesteps after event
                    
                    [dt] = caco3_main.timestep(time, nt_spn, nt_trs, nt_aft, time_spn, time_trs, dt);
                    
                    [d13c_ocn, d18o_ocn, ccflx, d18o_sp, d13c_sp] = ...
                        caco3_main.signal_flx(time, time_spn,time_trs,d13c_ocni,d13c_ocnf,d18o_ocni,d18o_ocnf ...
                        ,ccflx,bc.ccflxi,d18o_sp,d13c_sp,int_count,global_var.nspcc,flxfini,flxfinf, global_var.def_track2, global_var.def_size, global_var.def_biotest);
                    
                    [dep] = caco3_main.bdcnd(time, time_spn, time_trs, depi, depf, global_var.def_biotest);
                end
                
                % isotope signals represented by caco3 rain fluxes
                d18o_flx = sum(d18o_sp(:).*ccflx(:))/bc.ccflxi;
                d13c_flx = sum(d13c_sp(:).*ccflx(:))/bc.ccflxi;
                
                if(~global_var.def_track2)
                    if (abs(d13c_flx - d13c_ocn)>global_var.tol || abs(d18o_flx - d18o_ocn)>global_var.tol)  % check comparability with input signals
                        % in fortran stop
                        fprintf('error in assignment of proxy:\n');
                        fprintf('d18o_ocn, d13c_ocn, d18o_flx, d13c_flx \t %17.16E \t %17.16E \t %17.16E \t %17.16E \t', d18o_ocn, d13c_ocn, d18o_flx, d13c_flx);
                        msg = 'error in assignment of proxy, STOP.';
                        error(msg)
                    end
                end
                
                % <<<<<<<<<<<<<<<<<<<<<  NEW: 05/13/2019  <<<<<<<<<<<<<<<<<<<<< <<<<<<<<<<<<<<<<<<<<< <<<<<<<<<<<<<<<<<<<<<
                if(~global_var.def_size)
                    %  recording fluxes of two types of caco3 separately 
                    % write(file_bound,*) time, d13c_ocn, d18o_ocn, (ccflx(isp),isp=1,nspcc),temp, dep, sal,dici,alki, o2i
                    fmt=[repmat('%17.16e \t',1,global_var.nspcc+9) '\n'];
                    fprintf(file_boundid,fmt, time, d13c_ocn, d18o_ocn, ccflx(:), tmp, dep, sal, bc.dici,bc.alki, bc.o2i);
                else 
                    %  do not record separately 
                    % write(file_bound,*) time, d13c_ocn, d18o_ocn, sum(ccflx(1:4)),sum(ccflx(5:8)),(ccflx(isp),isp=1,nspcc),temp, dep, sal,dici,alki, o2i
                    fmt=[repmat('%17.16e \t',1,global_var.nspcc+11) '\n'];
                    fprintf(file_boundid,fmt, time, d13c_ocn, d18o_ocn, sum(ccflx(1:4)),sum(ccflx(5:8)), ccflx(:), tmp, dep, sal,bc.dici,bc.alki, bc.o2i);
                    
                end  
                % <<<<<<<<<<<<<<<<<<<<<  NEW: 05/13/2019  <<<<<<<<<<<<<<<<<<<<< <<<<<<<<<<<<<<<<<<<<< <<<<<<<<<<<<<<<<<<<<<
                
                %% === temperature & pressure and associated boundary changes ====
                % if temperature is changed during signal change event this affect diffusion coeff etc.
                %    call coefs(  &
                %        dif_dic,dif_alk,dif_o2,kom,kcc,co3sat & % output
                %        ,temp,sal,dep,nz,nspcc,poro,cai  & %  input
                %        )
                [keq1    ,keq2   ,keqcc  ,co3sat, dif_dic    ,dif_alk    ,dif_o2 ,kom, kcc, global_var] = caco3_main.coefs(tmp,sal,dep, global_var);
                %% /////////////////////
                
                
                itr_w = 0;          % # of iteration made to converge w
                err_w_min = 1d4;    % minimum relative difference of w compared to the previous w
                err_f = 0d0;        % relative different of total vol. fraction of solids wrt the previous value
                err_w = 10d0;           % err in burial rate
                
                while(err_w > global_var.tol)       %   	300 continue % point of restart when burial velocity does not converge
                    
                    fprintf('(int_count,dt,time)  %2.2i %11.3e %11.3e\n', int_count, dt, time);
%                     if(int_count == 800)
%                         fprintf('w, \t frc \n');
%                         for iz=1:interval:global_var.nz
%                             fprintf('%17.16e \t %17.16e \n', w(iz), frt(iz));
%                         end
%                     end
                    dw = zeros(1, global_var.nz);                       % burial rate change
                    
                    itr_om_o2 = 0;        % iteration number for om and o2 calcuation
                    zox_error = 1d4;  	% error in ieration for zox
                    minerr= 1d4;    % recording minimum relative difference in zox from previously considered zox
                    
                    while(zox_error > global_var.tol)
                        % %%%%%%%%%%%%%%%  om conc. calculation  % %%%%%%%%%%%%%%%
                        % omx: mol cm-3 sld; om conc.
                        % izox: integer for grid number of zox
                        % kom: degradation rate consts. for each nz grids
                        [omx, izox, kom] = ...
                            caco3_main.omcalc(bc.oxic, bc.anoxic,o2x,bc.om, global_var.komi,global_var.nz,global_var.sporo,global_var.sporoi,global_var.sporof, w, wi, dt, up, dwn, cnr, adf,trans, ...
                            global_var.nspcc, mix_type.labs,mix_type.turbo2, mix_type.nonlocal, omflx, global_var.poro, global_var.dz, global_var.o2th);
                        % calculating the fluxes relevant to om diagenesis (and checking the calculation satisfies the difference equations )
                        [omadv,omdec,omdif,omrain,omres,omtflx] = ...
                            caco3_main.calcflxom(omflx,global_var.sporo,bc.om,omx,dt,w,global_var.dz,global_var.z,global_var.nz,mix_type.turbo2,mix_type.labs, global_var.poro,up,dwn,cnr,adf,rho, global_var.mom,trans,kom,global_var.sporof);
                        
%                        fprintf('izox = %i \n',izox);      % sb omcalc calculates izox, which is the deepest grid where o2 >=0.
                        
                        % %                     omx_wtpc = omx.*global_var.mom./rho*100d0;
                        % %                     fprintf('~~~~ conc ~~~~ \n');
                        % %                     fprintf('z   ,OM   \n');
                        % %                     % showing parameters relevant to burial on screen
                        % %                     for iz=1:interval:global_var.nz
                        % %                         fprintf('%17.16e \t %17.16e\n', global_var.z(iz), omx_wtpc(iz));
                        % %                     end
                        % %                     fprintf('++++ flx ++++ \n');
                        % %                     fprintf(' \t tflx \t adv \t dif \t omrxn \t ccrxn \t rain \t res \n')
                        % %                     fprintf('OM :  \t %17.16e \t %17.16e \t %17.16e \t  %17.16e \t %17.16e \t %17.16e \t %17.16e', omtflx, omadv,  omdif, omdec,0d0,omrain, omres);
                        % %     %                fprintf('O2 :  \t %17.16e \t %17.16e \t %17.16e \t  %17.16e \t %17.16e \t %17.16e \t %17.16e', o2tflx, 0d0, o2dif, o2dec, 0d0, 0d0, o2res);
                        
                        
                        if (izox == global_var.nz)     % fully oxic; lower boundary condition ---> no diffusive out flow
                            % o2 calculation when o2 penetration depth (zox) is the same as bottom depth.
                            o2x = caco3_main.o2calc_ox(izox,global_var.nz,global_var.poro,bc.o2,kom,omx,global_var.sporo,dif_o2,global_var.dz,dt, global_var.ox2om, bc.o2i);
                            %  fluxes relevant to o2 (at the same time checking the satisfaction of difference equations)
                            [o2dec,o2dif,o2tflx,o2res] = caco3_main.calcflxo2_ox(global_var.nz,global_var.sporo,kom,omx,global_var.dz,global_var.poro,dif_o2,dt,bc.o2,o2x, global_var.ox2om, bc.o2i);
                        else        %% if oxygen is depleted within calculation domain, lower boundary changes to zero concs.
                            % o2 calculation when o2 is depleted within the calculation domain.
                            o2x = caco3_main.o2calc_sbox(izox,global_var.nz,global_var.poro,bc.o2,kom,omx,global_var.sporo,dif_o2,global_var.dz,dt, global_var.ox2om, bc.o2i);
                            % fluxes relevant to oxygen
                            [o2dec,o2dif,o2tflx,o2res] = caco3_main.calcflxo2_sbox(global_var.nz,global_var.sporo,kom,omx,global_var.dz,global_var.poro,dif_o2,dt,bc.o2,o2x,izox, global_var.ox2om, bc.o2i);
                        end
                        
                        %                 % showing intermediate results on screen
                        %                 omx_wtpc = omx.*global_var.mom./rho*100d0;
                        %                 fprintf('~~~~ conc ~~~~ itr_om_o2 = %i \n', itr_om_o2 +1);
                        %                 fprintf('z   ,OM   , o2\n');
                        %                 % showing parameters relevant to burial on screen
                        %                 for iz=1:interval:global_var.nz
                        %                     fprintf('%17.16e \t %17.16e \t %17.16e \n', global_var.z(iz), omx_wtpc(iz), o2x(iz)*1d3);
                        %                 end
                        %                 fprintf('++++ flx ++++ \n');
                        %                 fprintf(' \t tflx \t adv \t dif \t omrxn \t ccrxn \t rain \t res \n')
                        %                 fprintf('OM :  \t %17.16e \t %17.16e \t %17.16e \t  %17.16e \t %17.16e \t %17.16e \t %17.16e \n', omtflx, omadv,  omdif, omdec,0d0,omrain, omres);
                        %                 fprintf('O2 :  \t %17.16e \t %17.16e \t %17.16e \t  %17.16e \t %17.16e \t %17.16e \t %17.16e', o2tflx, 0d0, o2dif, o2dec, 0d0, 0d0, o2res);%
                        %                 fprintf('\n');
                        %                 fprintf('\n');
                        
                        
                        % update of zox
                        zoxx = 0d0;      % zox dummy variable
                        iz=1;
                        %                    for iz=1:global_var.nz
                        while(iz<global_var.nz+1)
                            if (o2x(iz)<=0d0)
                                break   % exit in fortran (terminates for-loop)
                            end
                            iz = iz+1;
                        end
                        
                        if (iz==global_var.nz+1) % oxygen never gets less than 0
                            zoxx = global_var.ztot; % zox is the bottom depth
                        elseif (iz==1)      % calculating zox interpolating at z=0 with SWI conc. and at z=z(iz) with conc. o2x(iz)
                            zoxx = (global_var.z(iz)*bc.o2i*1d-6/1d3 + 0d0*abs(o2x(iz)))/(bc.o2i*1d-6/1d3+abs(o2x(iz)));
                        else     % calculating zox interpolating at z=z(iz-1) with o2x(iz-1) and at z=z(iz) with conc. o2x(iz)
                            zoxx = (global_var.z(iz)*o2x(iz-1) + global_var.z(iz-1)*abs(o2x(iz)))/(o2x(iz-1)+abs(o2x(iz)));
                        end
                        
                        % error evaluation as relative difference of zox
                        zox_error = abs((global_var.zox -zoxx)/global_var.zox);
                        
                        %                    fprintf( 'itr_om_o2,zox, zoxx, error %i \t %17.16e \t %17.16e \t %17.16e \n',itr_om_o2, global_var.zox, zoxx, error);
                        %                    fprintf('~~~~~~~~~~~////~~~~~~~~~~~~~ \n');
                        
                        if (global_var.zox==zoxx)
                            break   % exit in fortran (terminates for-loop)
                        end
                        global_var.zox = 0.5d0*(global_var.zox + zoxx);  % new zox
                        
                        % if iteration reaches 100, error in zox is tested assuming individual grid depths as zox and find where error gets minimized
                        if (itr_om_o2>=100 && itr_om_o2 <= global_var.nz+99)
                            global_var.zox = global_var.z(itr_om_o2-99); % zox value in next test
                            if (minerr >=zox_error )	 % if this time error is less than last adopt as optimum
                                if (itr_om_o2~=100)
                                    izox_minerr = itr_om_o2 -100;
                                    minerr = zox_error;
                                end
                            end
                        elseif (itr_om_o2 == (global_var.nz+100))    % check last test z(nz)
                            if (minerr >=zox_error )
                                izox_minerr = itr_om_o2 -100;
                                minerr = zox_error;
                            end
                            global_var.zox = z(izox_minerr);  % determine next test which should be most optimum
                        elseif (itr_om_o2 == (global_var.nz+101))  % results should be optimum and thus exit
                            break   % exit in fortran (terminates for-loop)
                        end
                        
                        if (itr_om_o2 > (global_var.nz+101))
                            % in fortran stop
                            msg = 'Error: (itr_om_o2 > (nz+101)), STOP.';
                            error(msg)
                        end
                        
                        itr_om_o2 = itr_om_o2 + 1;
                    end
                    %~~  OM & O2 calculation END ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                    
                    % calculation of oxic and anoxic degradation of om (oxco2 and anco2, respectively)
                    
                    for iz = 1:global_var.nz
                        if (o2x(iz) > global_var.o2th)
                            oxco2(iz) = (1d0-global_var.poro(iz))*kom(iz)*omx(iz);  % aerobic respiration
                        else
                            % o2x(iz) = o2th
                            if (bc.anoxic)
                                anco2(iz) = (1d0-global_var.poro(iz))*kom(iz)*omx(iz);  % anaerobic respiration
                            end
                        end
                    end
                    
                    for iz=1:global_var.nz
                        dw(iz) = dw(iz) -(1d0-global_var.poro(iz))*mvom*kom(iz)*omx(iz);  %% burial rate change need reflect volume change caused by chemical reactions
                        % as well as non-local mixing
                        if (mix_type.turbo2(1) || mix_type.labs(1))
                            for iiz = 1:global_var.nz
                                if (trans(iiz,iz,1)==0d0)
                                    continue        % cycle in fortran
                                end
                                dw(iz) = dw(iz) - mvom*(-trans(iiz,iz,1)/global_var.dz(iz)*global_var.dz(iiz)*(1d0-global_var.poro(iiz))*omx(iiz));
                            end
                        else
                            if (mix_type.nonlocal(1))
                                for iiz = 1:global_var.nz
                                    if (trans(iiz,iz,1)==0d0)
                                        continue        % cycle in fortran
                                    end
                                    dw(iz) = dw(iz) - mvom*(-trans(iiz,iz,1)/global_var.dz(iz)*omx(iiz));
                                end
                            end
                        end
                    end
                    
                    
                    for iz=1:global_var.nz
                        if (omx(iz)<global_var.omx_th)
                            omx(iz)=global_var.omx_th;  %% truncated at minimum value
                        end
                    end
                    
                    %% calculation of caco3 system
                    % call calccaco3sys()
                    [ccx,dicx,alkx,rcc,dt, flg_500, itr_om_o2] = ...
                        caco3_main.calccaco3sys(ccx,dicx,alkx,rcc, dt, global_var.nspcc,dic,alk,dep,sal,tmp,mix_type.labs,mix_type.turbo2,mix_type.nonlocal, ...
                        global_var.sporo, global_var.sporoi, global_var.sporof, global_var.poro, dif_alk, dif_dic, ...
                        w, up, dwn, cnr, adf, global_var.dz, trans, cc, oxco2, anco2, co3sat, kcc, ccflx, global_var.ncc, global_var.nz, ...
                        global_var.tol, global_var.poroi, flg_500, global_var.fact, bc.alki,bc.dici, global_var.ccx_th, global_var.def_nonrec, global_var.def_sparse, global_var.def_showiter, global_var.def_sense);
                    
                    if(flg_500)
                        msg = 'error after calccaco3sys, STOP.';
                        error(msg)
                    end
                    % calling subroutine to calculate all aqueous co2 species and pH
                    % call calcspecies(dicx,alkx,temp,sal,dep,prox,co2x,hco3x,co3x,nz,infosbr)
                    [prox,co2x,hco3x,co3x,infosbr] = caco3_therm.calcspecies(dicx,alkx,tmp,sal,dep);
                    
                    if (infosbr==1)
                        msg = 'error after calcspecies after calccaco3sys, STOP.';
                        error(msg)
                    end
                    
                    %          % calculation of fluxes relevant to caco3 and co2 system
                    %             % call calcflxcaco3sys()
                    [cctflx,ccdis,ccdif,ccadv,ccrain,ccres,alktflx,alkdis,alkdif,alkdec,alkres, dictflx,dicdis,dicdif,dicres,dicdec, dw] = ...
                        caco3_main.calcflxcaco3sys(dw, global_var.nspcc, ccx, cc, ccflx,dt, global_var.dz, rcc, adf, up, dwn, cnr, w, dif_alk, dif_dic, dic, dicx, alk, alkx, oxco2, anco2, trans, ...
                        mix_type.turbo2, mix_type.labs,mix_type.nonlocal, global_var.sporof, int_count, global_var.nz, global_var.poro, global_var.sporo, ...
                        bc.dici,bc.alki, mvcc, global_var.tol);
                    
                    %~~  caco3 calculations END ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                    
                    %% clay calculation
                    % call claycalc()
                    ptx = caco3_main.claycalc(global_var.nz,global_var.sporo,pt,dt,w,global_var.dz,detflx,adf,up,dwn,cnr,trans, ...
                        global_var.nspcc,mix_type.labs,mix_type.turbo2,mix_type.nonlocal,global_var.poro,global_var.sporof, global_var.msed, global_var.def_nonrec );
                    
                    %               call calcflxclay()
                    [pttflx,ptdif,ptadv,ptres,ptrain, dw] = caco3_main.calcflxclay(dw, global_var.nz,global_var.sporo, ptx, pt, dt, global_var.dz, detflx, w, adf, up, dwn, cnr, global_var.sporof, ...
                        trans, mix_type.turbo2, mix_type.labs, mix_type.nonlocal, global_var.poro, global_var.msed, mvsed);
                    
                    % end of clay calculation
                    
                    %%%%%%%%%%%%%%
                    
                    % checking for total volume of solids, density and burial velocity
                    
                    % call getsldprop() % get solid property, rho (density) and frt (total vol.frac)
                    [rho, frt] = caco3_main.getsldprop(global_var.nz, omx, ptx, ccx, global_var.nspcc, w, up, dwn, cnr, adf, global_var.z, global_var.mom, global_var.msed, global_var.mcc, mvom, mvsed, mvcc);
                    
                    err_f = max(abs(frt - 1d0));  % new error in total vol. fraction (must be 1 in theory)
                    %% ========= calculation of burial velocity =============================
                    
                    wx = w;  % recording previous burial velocity
                    
                    % call burialcalc() % get new burial velocity
                    [w,wi] = caco3_main.burial_calc(detflx,ccflx,global_var.nspcc,omflx,dw,global_var.dz,global_var.poro,global_var.nz, global_var.msed,mvsed,mvcc,mvom,global_var.poroi);
                    
                    %  determine calculation scheme for advection / factors for upwind scheme to represent burial advection
                    % call calcupwindscheme()
                    [up, dwn, cnr, adf] = caco3_main.calcupwindscheme(w, global_var.nz);
                    
                    
                    % error and iteration evaluation
                    itr_w = itr_w + 1;              % counting iteration for w
                    err_w = max(abs((w-wx)./wx));    % relative difference of w
                    if (err_w<err_w_min)
                        err_w_min= err_w;  % recording minimum relative difference of  w
                        wxx = wx;  % recording w which minimizes deviation of total sld fraction from 1
                    end
                    
                    if (itr_w>100)   % if iteration gets too many (100), force to end with optimum w where error is minimum
                        if (itr_w==101)
                            w = wxx;
                            % %                             go to 300  % Dominik: do just go back to 300 if err_w > tol
                        elseif (itr_w==102)
                            w = wxx;
                            fprintf('not converging w %i \t %17.16e \t %17.16e \n',time, err_w, err_w_min);
                            %       pause;
                            break;  % %     go to 400
                        end
                    end
                end     % %     if (err_w > tol) go to 300
                % %
                % %                     400 continue
                % %
                
                %********************************************************************************************************************************  ADDED-START
                
                % % depth -age conversion
                age = caco3_main.dep2age(global_var.dz, w, global_var.nz);
                
                % ---------------------
                %/////// ISOTOPES /////
                %  calculating bulk isotopic composition
                % allocate isotope signal for fine and coarse caco3 species
                d18o_blkf = zeros(1,global_var.nz);
                d13c_blkf = zeros(1,global_var.nz);
                d18o_blkc = zeros(1,global_var.nz);
                d13c_blkc = zeros(1,global_var.nz);
                for iz=1:global_var.nz
                    d18o_blk(iz) = sum(d18o_sp(:)'.*ccx(iz,:))/sum(ccx(iz,:));
                    d13c_blk(iz) = sum(d13c_sp(:)'.*ccx(iz,:))/sum(ccx(iz,:));
                    if(global_var.def_size)
                        d18o_blkf(iz) = sum(d18o_sp(1:4)*ccx(iz,1:4))/sum(ccx(iz,1:4));
                        d13c_blkf(iz) = sum(d13c_sp(1:4)*ccx(iz,1:4))/sum(ccx(iz,1:4));
                        d18o_blkc(iz) = sum(d18o_sp(5:8)*ccx(iz,5:8))/sum(ccx(iz,5:8));
                        d13c_blkc(iz) = sum(d13c_sp(5:8)*ccx(iz,5:8))/sum(ccx(iz,5:8));
                    end
                end
                
                % recording
                
                if (time>=rectime(cntrec))
                    %        call recordprofile(cntrec )
                    caco3_main.recordprofile(cntrec, global_var.nz, global_var.z, age, pt, global_var.msed, wi, rho, cc, ccx, dic, dicx, alk, alkx, co3, co3x, co3sat ...
                        , rcc, pro, o2x, oxco2, anco2, bc.om, global_var.mom, global_var.mcc, d13c_ocni, d18o_ocni, up,dwn, cnr, adf, global_var.nspcc, ptx, w, frt, prox, omx, d13c_blk, d18o_blk)
                    
                    cntrec = cntrec + 1;
                    if (cntrec == global_var.nrec+1)
                        break   % exit in fortran (terminates while(2>1)-loop 
                    end
                end
                
                % recording signals at 3 different depths (btm of mixed layer, 2xdepths of btm of mixed layer and btm depth of calculation domain)
                caco3_main.sigrec(w,file_sigmlyid,file_sigmlydid,file_sigbtmid,time,age,izrec,d13c_blk,d13c_blkc  ...
                    ,d13c_blkf,d18o_blk,d18o_blkc,d18o_blkf,ccx,global_var.mcc,rho,ptx,global_var.msed,izrec2,global_var.nz, global_var.def_size)
                %********************************************************************************************************************************  ADDED-END
                
                fprintf('error in frt: %17.16e \n', max(abs(frt - 1d0)));                                
                
                %%%%%%%%%%%%%%
                
                %% showing results on screen
                
                omx_wtpc = omx.*global_var.mom./rho'*100d0;
                fprintf('~~~~ conc ~~~~ itr_om_o2 = %i \n', itr_om_o2);
                fprintf('z, \t   OM, \t o2, \t cc, \t dic, \t alk, \t sed\n');
                % showing parameters relevant to burial on screen
                for iz=1:interval:global_var.nz
                    fprintf('%17.16e \t %17.16e \t %17.16e \t %17.16e \t %17.16e \t %17.16e \t %17.16e \n', ...
                        global_var.z(iz), omx_wtpc(iz), o2x(iz)*1d3, sum(ccx(iz,:)*global_var.mcc)/rho(iz)*100d0, dicx(iz)*1d3, alkx(iz)*1d3, ptx(iz)*global_var.msed/rho(iz)*100d0);
                end
                
                fprintf('   ..... multiple cc species ..... \n');
                %    write(dumchr(2),'(i0)') interval
                %    dumchr(1)="(i0.3,':',"//trim(adjustl(dumchr(2)))//"E11.3"//")"
                for isp=1:global_var.nspcc
                    fprintf('\n');
                    fprintf('cc species: %i \n', isp);
                    %        print dumchr(1),isp,(ccx(iz,isp)*mcc/rho(iz)*100d0,iz=1,nz,nz/interval)
                    for iz=1:interval:global_var.nz
                        fprintf('%17.16e \n',ccx(iz,isp)*global_var.mcc/rho(iz)*100d0);
                    end
                end
                
                fprintf('++++ flx ++++ \n');
                fprintf(' \t tflx \t adv \t dif \t omrxn \t ccrxn \t rain \t res \n')
                fprintf('om :  \t %17.16e \t %17.16e \t %17.16e \t  %17.16e \t %17.16e \t %17.16e \t %17.16e \n', omtflx, omadv,  omdif, omdec,0d0,omrain, omres );
                fprintf('o2 :  \t %17.16e \t %17.16e \t %17.16e \t  %17.16e \t %17.16e \t %17.16e \t %17.16e \n', o2tflx,0d0, o2dif,o2dec, 0d0,0d0,o2res );
                fprintf('cc :  \t %17.16e \t %17.16e \t %17.16e \t  %17.16e \t %17.16e \t %17.16e \t %17.16e \n', sum(cctflx),  sum(ccadv), sum(ccdif),0d0,sum(ccdis), sum(ccrain), sum(ccres) );
                fprintf('dic :  \t %17.16e \t %17.16e \t %17.16e \t  %17.16e \t %17.16e \t %17.16e \t %17.16e \n', dictflx, 0d0,dicdif, dicdec,  dicdis, 0d0,dicres );
                fprintf('alk :  \t %17.16e \t %17.16e \t %17.16e \t  %17.16e \t %17.16e \t %17.16e \t %17.16e \n', alktflx, 0d0, alkdif, alkdec, alkdis, 0d0, alkres );
                fprintf('sed :  \t %17.16e \t %17.16e \t %17.16e \t  %17.16e \t %17.16e \t %17.16e \t %17.16e \n', pttflx, ptadv,ptdif,  0d0, 0d0, ptrain, ptres );
                fprintf('   ..... multiple cc species ..... \n');
                for isp=1:global_var.nspcc
                    fprintf('%i \t %17.16e \t %17.16e \t %17.16e \t  %17.16e \t %17.16e \t %17.16e \t %17.16e \n', isp,cctflx(isp), ccadv(isp), ccdif(isp),0d0,ccdis(isp), ccrain(isp), ccres(isp) );
                end
                fprintf('==== burial etc ==== \n');
                fprintf('z, \t w, \t rho, \t frc \n');
                % showing parameters relevant to burial on screen
                for iz=1:interval:global_var.nz
                    fprintf('%17.16e \t %17.16e \t %17.16e \t %17.16e \n', ...
                        global_var.z(iz),w(iz), rho(iz), frt(iz));
                end
                
                fprintf('\n');
                fprintf('\n');
                fprintf('\n');
                
                om = omx;
                o2 = o2x;
                cc = ccx;
                dic = dicx;
                alk = alkx;
                pt = ptx;
                time = time +dt;
                
%                 d13c_ocni = 0d0;
%                 d18o_ocni = 0d0;
%                 d13c_blk = zeros(1, global_var.nz);
%                 d18o_blk = zeros(1, global_var.nz);

%                 caco3_main.recordprofile(int_count, global_var.nz, global_var.z, age, pt, global_var.msed, wi, rho, cc, ccx, dic, dicx, alk, alkx, co3, co3x, co3sat ...
%                     , rcc, pro, o2x, oxco2, anco2, om, global_var.mom, global_var.mcc, d13c_ocni, d18o_ocni, up,dwn, cnr, adf, global_var.nspcc, ptx, w, frt, prox, omx, d13c_blk, d18o_blk);
                int_count = 1+int_count;
            end
            
            %********************************************************************************************************************************  ADDED-START
        fclose(file_sigmlyid);      % recording signals etc at just below mixed layer 
        fclose(file_sigmlydid);     % recording signals etc at depths of 2x mixed layer thickness 
        fclose(file_sigbtmid);      % % recording signals etc at bottom of sediment  
%********************************************************************************************************************************  ADDED-END

            
        end
        
        
        
    end
end

