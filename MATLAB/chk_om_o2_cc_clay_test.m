        function chk_om_o2_cc_clay_test(tmp_in, sal_in, dep_in)
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
            dt = 100d0; % time step [yr]
            
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
                    global_var.tol, global_var.poroi, flg_500, global_var.fact, bc.alki,bc.dici, global_var.ccx_th, global_var.def_nonrec, global_var.def_sparse, global_var.def_showiter);
                
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