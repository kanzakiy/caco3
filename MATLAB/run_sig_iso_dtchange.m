function run_sig_iso_dtchange(cc_rain_flx_in, rainratio_in, dep_in, dt_in, oxonly_in, folder)
%
%   *******************************************************************   %
%   ***  Running IMP for signal tracking and lysocline exps ***********   %
%   *******************************************************************   %
    
%   cc_rain_flx_in      :   influx of CaCO3  [mol cm-2 yr-1]
%   rainratio_in        :   assumed OM/CaCO3-ratio [-]
%   dep_in              :   water depth [km] (set dep_max, below, to different depth for dissolution experiment)
%   dt_in               :   time step to start simulation with [yr]
%   oxonly_in           :   select oxic only model for OM degradation [true/false]
%   folder              :   path to output folder, will be created 

    % Example run: 
    % run_sig_iso_dtchange(6.0e-5, 1.5, 0.24, 1d8, true, './1207_test')
    
    

    %% function to run the model for signal tracking (def_sense = true; in main)
    %% and single lysocline experiments (def_sense = false; in main)
    %% changes dt automatically if pH can't be calculated but uses w and other results from previous calculation
    %% saves profiles in .txt files
    % burial rate is updated at a give time step until it converges well
    
    dep_max = 5.0d0;    %   max depth to be changed to during the experiment
%    folder = './0_test';
    mkdir(folder);      % create output folder
    interval =10;       % if def_nondisp = false: choose a value between 1 to nz; om depth profile is shown with this interval; e.g., if inteval = 1, conc. at all depths are shown
    % e.g., if interval = 10, om conc. at nz/10 depths are shown
    flag_steadystate = false;
    % initialize/define global properties/variables
    global_var = caco3_main;

    % initialize the boundary conditions
    [bc, global_var] = caco3_main.caco3_set_boundary_cond(global_var);

    % initialize mixing type
    mix_type = caco3_main.caco3_set_mixing(global_var);

    % set to input values to boundary conditions:
    bc.ccflxi = cc_rain_flx_in;      % mol (CaCO3) cm-2 yr-1
    global_var.om2cc = rainratio_in;   % rain ratio of organic matter to calcite
    dep = dep_in;   % 0.0d0; % depth in km
    dt = dt_in; % time step [yr]
    if(oxonly_in)
        bc.anoxic = false;         % oxic only model of OM degradation by Emerson (1985)
    end
    global_var.def_oxonly = oxonly_in;


    dw = zeros(1, global_var.nz);                       % burial rate change
    rcc = zeros(global_var.nz, global_var.nspcc);       % dissolution rate of caco3

    tmp = bc.tmp;   % 2.0;
    sal = bc.sal;   % 35.0;

    %            bc.ccflxi = 60d-6;      % mol (CaCO3) cm-2 yr-1 - caco3 flux
    %            global_var.om2cc = 0.7d0;   % rain ratio of organic matter to calcite

    % open files for output signal at 3 different depths

    file_sigmly = sprintf('%s/matlab_sigmly.txt', folder);
    file_sigmlyid = fopen(file_sigmly,'wt');
    file_sigmlyd = sprintf('%s/matlab_sigmlyd.txt', folder);
    file_sigmlydid = fopen(file_sigmlyd,'wt');
    file_sigbtm = sprintf('%s/matlab_sigbtm.txt', folder);
    file_sigbtmid = fopen(file_sigbtm,'wt');
    file_bound = sprintf('%s/matlab_bound.txt', folder);
    file_boundid = fopen(file_bound,'wt');
    file_frac = sprintf('%s/matlab_frac.txt', folder);
    file_fracid = fopen(file_frac,'wt');


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
    [rectime, cntrec, time_spn, time_trs, time_aft] = caco3_main.recordtime(global_var.nrec, wi, global_var.ztot, global_var.def_biotest, global_var.def_sense, global_var.def_nonrec, folder);

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
    % bc.o2 = bc.o2i*1d-6/1d3 * ones(1, global_var.nz);	% o2 conc. in uM converted to mol/cm3   % YK commented out
    om = 1d-8*ones(1, global_var.nz);                   % mol cm-3 sld; clay conc., assume an arbitrary low conc.  % YK added
    o2 = bc.o2i*1d-6/1d3 * ones(1, global_var.nz);	% o2 conc. in uM converted to mol/cm3   % YK added
    cc = 1d-8 * ones(global_var.nz, global_var.nspcc);	% mol cm-3 sld; concentration of caco3, assume an arbitrary low conc.
    dic = bc.dici*1d-6/1d3 * ones(1, global_var.nz);  	% mol/cm3; factor is added to change uM to mol/cm3
    alk = bc.alki*1d-6/1d3 * ones(1, global_var.nz);   	% mol/cm3
    pt = 1d-8*ones(1, global_var.nz);                   % mol cm-3 sld; clay conc., assume an arbitrary low conc.


    % calling subroutine to calculate all aqueous co2 species and pH
    [pro,co2,hco3,co3,infosbr] = caco3_therm.calcspecies(dic,alk,tmp,sal,dep);

    % omx = bc.om;  % YK commented out
    % o2x = bc.o2;  % YK commented out
    omx = om;  % YK added
    o2x = o2;  % YK added
    ccx = cc;
    dicx = dic;
    alkx = alk ;
    ptx = pt;
    % this may not be necessary as these individual species assume equilibrium
    co2x = co2;
    hco3x = hco3;
    co3x = co3;

    co3i=co3(1);  % recording seawater conc. of co3

    time = 0d0;     % model time [yr]
    int_count = 1;  % integration count
    nt = 10;        % total integration
    dt = dt_in;   %1d2;       % time step [yr]
    warmup_done = false;

    rho = 2.5d0*ones(global_var.nz, 1); % assume here density (this is going to be calculated based on solid phase composition )

    oxco2 = zeros(1, global_var.nz);  % oxic degradation of om; here initially assumed 0
    anco2 = zeros(1, global_var.nz);   % anoxic degradation of om; here initially assumed 0


    %%%  addition to chk_om.f90 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    global_var.zox = 10d0;  % initial assumption on oxygen penetaration depth [cm]
    %%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %%  START OF TIME INTEGLATION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %            for int_count=1:nt
    flag_smaller_dt = false;
    while(2>1)              %   Here change to unknwon iteration

        %% ///////// isotopes & fluxes settings //////////////
        if(~global_var.def_sense)   % do signal tracking

            % determine time step dt by calling timestep(time,nt_spn,nt_trs,nt_aft,dt) where nt_xx denotes total iteration number
            nt_spn = 800;     % timesteps for spin-up % YK modified
            if (~warmup_done)
                if (int_count<11)
                    nt_spn = 80000;
                elseif (int_count<21)
                    nt_spn = 8000;
                else
                    warmup_done = true;
                    time = 0d0;
                    int_count = 1;
                    continue        % cycle in fortran
                end
            end
            nt_trs = 5000;    % timesteps close to & during event (signal transition) % YK modified
            nt_aft = 1000;    % timesteps after event  % YK modified

            [dt] = caco3_main.timestep(time, nt_spn, nt_trs, nt_aft, time_spn, time_trs, dt);

            [d13c_ocn, d18o_ocn, ccflx, d18o_sp, d13c_sp] = ...
                caco3_main.signal_flx(time, time_spn,time_trs,d13c_ocni,d13c_ocnf,d18o_ocni,d18o_ocnf ...
                ,ccflx,bc.ccflxi,d18o_sp,d13c_sp,int_count,global_var.nspcc,flxfini,flxfinf, global_var.def_track2, global_var.def_size, global_var.def_biotest);

            [dep] = caco3_main.bdcnd(time, time_spn, time_trs, depi, depf, global_var.def_biotest);
        else
            if(~flag_smaller_dt)    % set to initial dt, if normal iteration
                dt = dt_in;                % set dt to input time-step
            else
                flag_smaller_dt = false;
            end
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
        w_save = w;          % saving previous w

        while(err_w > global_var.tol)       %   	300 continue % point of restart when burial velocity does not converge
            if(int_count == 16)
                34
            end
            fprintf('(int_count,dt,time)  %2.2i %11.3e %11.3e\n', int_count, dt, time);
            dw = zeros(1, global_var.nz);                       % burial rate change caused by reaction and non-local mixing

            oxco2 = zeros(1, global_var.nz);  % oxic degradation of om; here initially assumed 0  % YK added
            anco2 = zeros(1, global_var.nz);   % anoxic degradation of om; here initially assumed 0 % YK added

            % ~~~~~~~~~ OM & O2 iteration wrt zox ~~~~~~~~~~~~~~~

            itr_om_o2 = 0;        % iteration number for om and o2 calcuation
            zox_error = 1d4;  	% error in ieration for zox
            minerr= 1d4;    % recording minimum relative difference in zox from previously considered zox

            while(2>1)       % (zox_error > global_var.tol) % new 18.06.2019
                % %%%%%%%%%%%%%%%  om conc. calculation  % %%%%%%%%%%%%%%%
                % omx: mol cm-3 sld; om conc.
                % izox: integer for grid number of zox
                % kom: degradation rate consts. for each nz grids

                dt_om_o2 = 1d8;
                dt_om_o2 = dt;

                % calculating zox from assumed/previous o2 profiles
                [izox,kom,zox,kom_ox,kom_an] = caco3_main.calc_zox(bc.oxic, bc.anoxic, global_var.nz, o2x, global_var.o2th, global_var.komi, global_var.ztot, global_var.z, bc.o2i, global_var.dz);

                [omx] = ...
                    caco3_main.omcalc(bc.oxic, bc.anoxic,o2x,om, global_var.komi,global_var.nz ...  % YK changed bc.om to om
                    ,global_var.sporo,global_var.sporoi,global_var.sporof, w, wi, dt_om_o2, up, dwn, cnr, adf,trans, ...
                    global_var.nspcc, mix_type.labs,mix_type.turbo2, mix_type.nonlocal, omflx, global_var.poro, global_var.dz, global_var.o2th, kom);
                % calculating the fluxes relevant to om diagenesis (and checking the calculation satisfies the difference equations )
                [omadv,omdec,omdif,omrain,omres,omtflx] = ...
                    caco3_main.calcflxom(omflx,global_var.sporo,om,omx,dt_om_o2,w,global_var.dz ...  % YK changed bc.om to om
                    ,global_var.z,global_var.nz,mix_type.turbo2,mix_type.labs, global_var.poro ...
                    ,up,dwn,cnr,adf,rho, global_var.mom,trans,kom,global_var.sporof);


                %	if (izox == global_var.nz)     % fully oxic; lower boundary condition ---> no diffusive out flow
                % o2 calculation when o2 penetration depth (zox) is the same as bottom depth.
                o2x = caco3_main.o2calc_ox(izox,global_var.nz,global_var.poro,o2 ...  % YK changed bc.o2 to o2
                    ,kom_ox,omx,global_var.sporo,dif_o2,global_var.dz,dt_om_o2, global_var.ox2om, bc.o2i);
                %  fluxes relevant to o2 (at the same time checking the satisfaction of difference equations)
                [o2dec,o2dif,o2tflx,o2res] = caco3_main.calcflxo2_ox(global_var.nz,global_var.sporo ...
                    ,kom_ox,omx,global_var.dz,global_var.poro,dif_o2,dt_om_o2,o2,o2x, global_var.ox2om, bc.o2i); % YK changed bc.o2 to o2
                %                        else        %% if oxygen is depleted within calculation domain, lower boundary changes to zero concs.
                % o2 calculation when o2 is depleted within the calculation domain.

                %%%% new 18.06.2019
                if (all(o2x>=0d0) && izox == global_var.nz)
                    iizox_errmin = global_var.nz;
                    all_oxic = true;       %  XXXXXXXXXXXXXXXXXXXXXXX  NEW ---- June 19 2019 YK added XXXXXXXXXXXXXXXXXXXXXXX
                    % print *,'all oxic',iizox_errmin
                elseif (any(o2x<0d0))
                    all_oxic = false;      %  XXXXXXXXXXXXXXXXXXXXXXX  NEW ---- June 19 2019 YK added XXXXXXXXXXXXXXXXXXXXXXX
                    error_o2min = 1d4;
                    iizox_errmin = izox;
                    for iizox = 1:global_var.nz
                        % fluxes relevant to oxygen
                        o2x = caco3_main.o2calc_sbox(iizox,global_var.nz,global_var.poro,o2 ...   % YK changed bc.o2 to o2
                            ,kom_ox,omx,global_var.sporo,dif_o2,global_var.dz,dt_om_o2, global_var.ox2om, bc.o2i);
                        %                               call o2calc_sbox(  &
                        %                                                               o2x  & % output
                        %                             ,iizox,nz,poro,o2,kom_ox,omx,sporo,dif_o2,dz,dt_om_o2,ox2om,o2i & % input
                        %                             )
                        if (all(o2x>=0d0))
                            if (abs(o2x(max(iizox-1,1)))<error_o2min)
                                error_o2min = abs(o2x(max(iizox-1,1)));
                                iizox_errmin = iizox;
                                % print*,'find smaller difference',iizox_errmin
                            end
                        end
                    end
                    %%%% new 18.06.2019 until here

                    o2x = caco3_main.o2calc_sbox(iizox_errmin,global_var.nz,global_var.poro,o2 ...   % YK changed bc.o2 to o2
                        ,kom_ox,omx,global_var.sporo,dif_o2,global_var.dz,dt_om_o2, global_var.ox2om, bc.o2i);
                    % fluxes relevant to oxygen
                    [o2dec,o2dif,o2tflx,o2res] = caco3_main.calcflxo2_sbox(global_var.nz,global_var.sporo,kom_ox,omx ...
                        ,global_var.dz,global_var.poro,dif_o2,dt_om_o2,o2,o2x,iizox_errmin, global_var.ox2om, bc.o2i);   % YK changed bc.o2 to o2

                    [iizox_errmin,kom_dum(:,1),zox,kom_dum(:,2),kom_dum(:,3)] = caco3_main.calc_zox(bc.oxic, bc.anoxic, global_var.nz, o2x, global_var.o2th, global_var.komi, global_var.ztot, global_var.z, bc.o2i, global_var.dz);

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
                %                         zoxx = 0d0;      % zox dummy variable
                %                         iz=1;
                %                         %                    for iz=1:global_var.nz
                %                         while(iz<global_var.nz+1)
                %                             if (o2x(iz)<=0d0)
                %                                 break   % exit in fortran (terminates for-loop)
                %                             end
                %                             iz = iz+1;
                %                         end
                %
                %                         if (iz==global_var.nz+1) % oxygen never gets less than 0
                %                             zoxx = global_var.ztot; % zox is the bottom depth
                %                         elseif (iz==1)      % calculating zox interpolating at z=0 with SWI conc. and at z=z(iz) with conc. o2x(iz)
                %                             zoxx = (global_var.z(iz)*bc.o2i*1d-6/1d3 + 0d0*abs(o2x(iz)))/(bc.o2i*1d-6/1d3+abs(o2x(iz)));
                %                         else     % calculating zox interpolating at z=z(iz-1) with o2x(iz-1) and at z=z(iz) with conc. o2x(iz)
                %                             zoxx = (global_var.z(iz)*o2x(iz-1) + global_var.z(iz-1)*abs(o2x(iz)))/(o2x(iz-1)+abs(o2x(iz)));
                %                         end

                % error evaluation as relative difference of zox
                %                        zox_error = abs((global_var.zox -zoxx)/global_var.zox);
                zox_error = abs(izox-iizox_errmin);     % relative difference
                %                    fprintf( 'itr_om_o2,zox, izox, error %i \t %17.16e \t %17.16e \t %17.16e \n',itr_om_o2, global_var.zox, izox, zox_error);
                %                    fprintf('~~~~~~~~~~~////~~~~~~~~~~~~~ \n');

                %                         if (global_var.zox==zoxx)
                %                             break   % exit in fortran (terminates for-loop)
                %                         end
                %                         global_var.zox = 0.5d0*(global_var.zox + zoxx);  % new zox
                %
                %                         % if iteration reaches 100, error in zox is tested assuming individual grid depths as zox and find where error gets minimized
                %                         if (itr_om_o2>=100 && itr_om_o2 <= global_var.nz+99)
                %                             global_var.zox = global_var.z(itr_om_o2-99); % zox value in next test
                %                             if (minerr >=zox_error )	 % if this time error is less than last adopt as optimum
                %                                 if (itr_om_o2~=100)
                %                                     izox_minerr = itr_om_o2 -100;
                %                                     minerr = zox_error;
                %                                 end
                %                             end
                %                         elseif (itr_om_o2 == (global_var.nz+100))    % check last test z(nz)
                %                             if (minerr >=zox_error )
                %                                 izox_minerr = itr_om_o2 -100;
                %                                 minerr = zox_error;
                %                             end
                %                             global_var.zox = z(izox_minerr);  % determine next test which should be most optimum
                %                         elseif (itr_om_o2 == (global_var.nz+101))  % results should be optimum and thus exit
                %                             break   % exit in fortran (terminates for-loop)
                %                         end
                %
                %                         if (itr_om_o2 > (global_var.nz+101))
                %                             % in fortran stop
                %                             msg = 'Error: (itr_om_o2 > (nz+101)), STOP.';
                %                             error(msg)
                %                         end

                %                         if (izox==iizox_errmin)
                %                             break   % exit in fortran (terminates for/while-loop)
                %                         end

                %  VVVVVVVVVVVVVVVVVVVVVVVVVVV  NEW start ---- June 19 2019 YK added VVVVVVVVVVVVVVVVVVVVVVVVVVV
                if (izox==iizox_errmin)
                    if (all_oxic)
                        break   % exit in fortran (terminates for/while-loop)
                    else
                        o2x = caco3_main.o2calc_sbox(izox,global_var.nz,global_var.poro,o2 ...   % YK changed bc.o2 to o2
                            ,kom_ox,omx,global_var.sporo,dif_o2,global_var.dz,dt_om_o2, global_var.ox2om, bc.o2i);
                        % print'(A,I0,10E11.3)', 'o2 :',izox,(o2x(iz)*1d3,iz=1,nz,nz/10)
                        % fluxes relevant to oxygen
                        [o2dec,o2dif,o2tflx,o2res] = caco3_main.calcflxo2_sbox(global_var.nz,global_var.sporo,kom_ox,omx ...
                            ,global_var.dz,global_var.poro,dif_o2,dt_om_o2,o2,o2x,izox, global_var.ox2om, bc.o2i);   % YK changed bc.o2 to o2
                        break   % exit in fortran (terminates for/while-loop)
                    end
                end
                %  VVVVVVVVVVVVVVVVVVVVVVVVVVV  NEW end ---- June 19 2019 YK added VVVVVVVVVVVVVVVVVVVVVVVVVVV

                if (zox_error < minerr )
                    minerr = zox_error;
                else
                    if (izox < global_var.nz && iizox_errmin ==global_var.nz)
                        o2x = caco3_main.o2calc_sbox(izox,global_var.nz,global_var.poro,o2 ...   % YK changed bc.o2 to o2
                            ,kom_ox,omx,global_var.sporo,dif_o2,global_var.dz,dt_om_o2, global_var.ox2om, bc.o2i);
                        % fluxes relevant to oxygen
                        [o2dec,o2dif,o2tflx,o2res] = caco3_main.calcflxo2_sbox(global_var.nz,global_var.sporo,kom_ox,omx ...
                            ,global_var.dz,global_var.poro,dif_o2,dt_om_o2,o2,o2x,izox, global_var.ox2om, bc.o2i);   % YK changed bc.o2 to o2
                        break   % exit in fortran (terminates for/while-loop)
                    end
                end

                itr_om_o2 = itr_om_o2 + 1;
            end
            %~~  OM & O2 calculation END ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            % calculation of oxic and anoxic degradation of om (oxco2 and anco2, respectively)
            fprintf('finising om & o2 \n');

            %                     for iz = 1:global_var.nz
            %                         if (o2x(iz) > global_var.o2th)
            %                             oxco2(iz) = (1d0-global_var.poro(iz))*kom(iz)*omx(iz);  % aerobic respiration
            %                         else
            %                             % o2x(iz) = o2th
            %                             if (bc.anoxic)
            %                                 anco2(iz) = (1d0-global_var.poro(iz))*kom(iz)*omx(iz);  % anaerobic respiration
            %                             end
            %                         end
            %                     end

            oxco2(:) = (1d0-global_var.poro(:)).*kom_ox(:).*omx(:);
            anco2(:) = (1d0-global_var.poro(:)).*kom_an(:).*omx(:);

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
                dt = dt/10;
                if(global_var.def_sense)
                    % was in fortran go to 500
                    ccx = cc;
                    alkx = alk;
                    dicx = dic;
                    w = w_save ;
                    % % determine factors for upwind scheme to represent burial advection
                    [up, dwn, cnr, adf] = caco3_main.calcupwindscheme(w, global_var.nz);
                    % in fortran  go to 600
                    flag_smaller_dt = true;
                    break       % jump out of loop and start with smaller dt (in fortran go to 600)
                end
                %                         msg = 'error after calccaco3sys, STOP.';
                %                         error(msg)
            end
            % calling subroutine to calculate all aqueous co2 species and pH
            % call calcspecies(dicx,alkx,temp,sal,dep,prox,co2x,hco3x,co3x,nz,infosbr)
            [prox,co2x,hco3x,co3x,infosbr] = caco3_therm.calcspecies(dicx,alkx,tmp,sal,dep);

            if (infosbr==1)
                dt = dt/10;
                if(global_var.def_sense)
                    % was in fortran go to 500
                    ccx = cc;
                    alkx = alk;
                    dicx = dic;
                    w = w_save ;
                    % % determine factors for upwind scheme to represent burial advection
                    [up, dwn, cnr, adf] = caco3_main.calcupwindscheme(w, global_var.nz);
                    % in fortran  go to 600
                    flag_smaller_dt = true;
                    break       % jump out of loop and start with smaller dt (in fortran go to 600)
                else
                    msg = 'error after calcspecies after calccaco3sys - calcspecies, STOP.';
                    error(msg)
                end
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
            if(global_var.def_sense)
                if (err_f < global_var.tol)    %  ! if total vol. fraction is near enough to 1, steady-state solution is obtained
                    flag_steadystate = true;
                    break   % in fortran exit  (terminates loop: while(err_w > global_var.tol))
                end
            end %#endif                    %% ========= calculation of burial velocity =============================

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

        if(flag_smaller_dt) % if here bc need to use lower dt, start while from beginning
            continue
        end
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
                , rcc, pro, o2x, oxco2, anco2, bc.om, global_var.mom, global_var.mcc, d13c_ocni, d18o_ocni, up,dwn, cnr, adf, global_var.nspcc, ptx, w, frt, prox, omx, d13c_blk, d18o_blk, folder)

            cntrec = cntrec + 1;
            if (cntrec == global_var.nrec+1)
                break   % exit in fortran (terminates while(2>1)-loop
            end
        end

        if(global_var.def_sense)    % if no signal tracking and steady-state is reached exit calculation
            if(flag_steadystate)
                break
            end
        end

        % recording signals at 3 different depths (btm of mixed layer, 2xdepths of btm of mixed layer and btm depth of calculation domain)
        caco3_main.sigrec(w,file_sigmlyid,file_sigmlydid,file_sigbtmid,time,age,izrec,d13c_blk,d13c_blkc  ...
            ,d13c_blkf,d18o_blk,d18o_blkc,d18o_blkf,ccx,global_var.mcc,rho,ptx,global_var.msed,izrec2,global_var.nz, global_var.def_size)
        %********************************************************************************************************************************  ADDED-END

        fprintf('error in frt: %17.16e \n', max(abs(frt - 1d0)));
        fmt=[repmat('%17.16e \t',1,2) '\n'];   % YK added
        fprintf(file_fracid,fmt, time, max(abs(frt - 1d0))); % YK added
        %%%%%%%%%%%%%%

        %% showing results on screen?
        if(~global_var.def_nondisp)

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
        end

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
    fclose(file_boundid);      % % boundary condition file closing  % YK added
    fclose(file_fracid);      % %   file recording total sld fraction closing  % YK added
    %********************************************************************************************************************************  ADDED-END

    % recording end results for lysoclines and caco3 burial fluxes
    if(oxonly_in)
        sprintf('%s/matlab_frac.txt', folder);
%         str_lys = sprintf('./1207_test/lys_sense_cc-%2.1e_rr-%.2f_oxonly.txt',cc_rain_flx_in, rainratio_in);
%         str_ccbur = sprintf('./1207_test/ccbur_sense_cc-%2.1e_rr-%.2f_oxonly.txt',cc_rain_flx_in, rainratio_in);
        str_lys = sprintf('%s/lys_sense_cc-%2.1e_rr-%.2f_oxonly.txt',folder,cc_rain_flx_in, rainratio_in);
        str_ccbur = sprintf('%s/ccbur_sense_cc-%2.1e_rr-%.2f_oxonly.txt',folder,cc_rain_flx_in, rainratio_in);
    else
%         str_lys = sprintf('./1207_test/lys_sense_cc-%2.1e_rr-%.2f_oxanox.txt',cc_rain_flx_in, rainratio_in);
%         str_ccbur = sprintf('./1207_test/ccbur_sense_cc-%2.1e_rr-%.2f_oxanox.txt',cc_rain_flx_in, rainratio_in);
        str_lys = sprintf('%s/lys_sense_cc-%2.1e_rr-%.2f_oxanox.txt',folder, cc_rain_flx_in, rainratio_in);
        str_ccbur = sprintf('%s/ccbur_sense_cc-%2.1e_rr-%.2f_oxanox.txt',folder, cc_rain_flx_in, rainratio_in);
    end

    file_tmp = fopen(str_lys,'at+');    % the 4th column is the plotted CaCO3 wt%!
    fprintf(file_tmp,'%17.16e \t %17.16e \t %17.16e \t %17.16e \t %17.16e \t %17.16e \t %17.16e \n ', 1d6*(co3i*1d3-co3sat), sum(ccx(1,:))*global_var.mcc/rho(1)*100d0, frt(1), ...
        sum(ccx(global_var.nz,:))*global_var.mcc/rho(global_var.nz)*100d0, frt(global_var.nz),sum(ccx(izml,:))*global_var.mcc/rho(izml)*100d0, frt(izml));
    fclose(file_tmp);

    file_tmp = fopen(str_ccbur,'at+');
    fprintf(file_tmp,'%17.16e \t %17.16e \n ', 1d6*(co3i*1d3-co3sat), 1d6*sum(ccadv));
    fclose(file_tmp);
end
