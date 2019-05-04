        function chk_cc_test(tmp_in, sal_in, dep_in)
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
            [global_var.dz,global_var.z] = caco3_main.makegrid(beta,global_var.nz, global_var.ztot);
            
            
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
            dt = 1d6; % time step [yr]
            
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
                    global_var.tol, global_var.poroi, flg_500, global_var.fact, bc.alki,bc.dici, global_var.ccx_th, global_var.def_nonrec, global_var.def_sparse, global_var.def_showiter);
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