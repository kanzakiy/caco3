classdef caco3_global_var
    % class to define variables
    
    properties
        
        nz = 100;  % grid number
        % % #ifdef sense
        nspcc = 1;  % number of CaCO3 species
        % % #elif defined track2
        nspcc = 42;
        % % #elif defined size
        nspcc = 8;
        % % #else
        nspcc = 4;
        % % #endif
        nsig = 2;
        cc; ccx;  % mol cm-3 sld; concentration of caco3, subscript x denotes dummy variable used during iteration
        om; omx;  % mol cm-3 sld; om conc.
        pt; ptx; % mol cm-3 sld; clay conc.
        ccflx; d13c_sp;d18o_sp; % flux of caco3, d13c signal of caco3, d18o signal of caco3
        d13c_blk; d18o_blk;  % d13c signal of bulk caco3, d18o signal of bulk caco3
        d13c_blkf; d18o_blkf;  d13c_blkc; d18o_blkc; % subscripts f and c denotes variables of fine and coarse caco3 species, respectively
        d13c_flx; d18o_flx;  % d13c signal averaged over flux values, d18o counterpart
        d13c_ocni; d13c_ocnf; d13c_ocn;  % initial value of ocean d13c, final value of ocean d13c, ocean d13c
        d18o_ocni; d18o_ocnf; d18o_ocn;  % the same as above expect 18o insted of 13c
        sigocn; sigocni; sigocnf; sigflx; sigblk; sigblkf; sigblkc; 
        flxfrc; flxfrc2;   % flux fractions used for assigning flux values to realized isotope input changes
        ccflxi = 10d-6;  % mol (CaCO3) cm-2 yr-1  % a reference caco3 flux; Emerson and Archer (1990)
        omflx = 12d-6;  % mol cm-2 yr-1       % a reference om flux; Emerson and Archer (1990)
        detflx = 180d-6;  % g cm-2 yr-1  % a reference detrital flux; MUDS input http://forecast.uchicago.edu/Projects/muds.html
        alki =  2285d0;  % uM  % a reference ALK; MUDS
        dici = 2211d0;  % uM   % a reference DIC; MUDS
        o2i = 165d0;  % uM     % a reference O2 ; MUDS
        % :: komi = 2d0  % /yr  % a reference om degradation rate const.; MUDS
        komi = 0.5d0;   % /yr  % arbitrary
        % :: komi = 0.1d0  % /yr  % Canfield 1994
        % :: komi = 0.06d0  % /yr  % ?? Emerson 1985? who adopted relatively slow decomposition rate
        % :: kcci = 10.0d0*365.25d0  % /yr; a reference caco3 dissolution rate const.
        % #ifndef nodissolve
        kcci = 1d0*365.25d0;   % /yr  ;cf., 0.15 to 30 d-1 Emerson and Archer (1990) 0.1 to 10 d-1 in Archer 1991
        % #else
        kcci = 0d0*365.25d0;   % /yr
        % #endif
        poroi = 0.8d0;   % a reference porosity
        keqcc = 4.4d-7;    % mol2 kg-2 caco3 solutiblity (Mucci 1983 cited by Emerson and Archer 1990)
        ncc = 4.5d0;    % (Archer et al. 1989) reaction order for caco3 dissolution
        temp = 2d0;   % C a refernce temperature
        sal = 35d0;   % wt o/oo  salinity
        cai = 10.3d-3;  % mol kg-1 calcium conc. seawater
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
        om2cc = 0.666d0;   % rain ratio of organic matter to calcite
        % :: om2cc = 0.5d0  % rain ratio of organic matter to calcite
        o2th = 0d0;  % threshold oxygen level below which not to calculate
        dev = 1d-6;  % deviation addumed
        zml_ref = 12d0;  % a referece mixed layer depth
        ccx_th = 1d-300;  % threshold caco3 conc. (mol cm-3) below which calculation is not conducted
        omx_th = 1d-300;  % threshold om    conc. (mol cm-3) below which calculation is not conducted
        dif_alk0;  dif_dic0;  dif_o20;   % diffusion coefficient of alk, dik and o2 in seawater
        zml;  zrec;  zrec2;   % mixed layer depth, sediment depth where proxy signal is read 1 & 2
        chgf;   % variable to check change in total fraction of solid materials
        flxfin;  flxfini;  flxfinf;   %  flux ratio of fine particles: i and f denote initial and final values
        pore_max;  exp_pore;  calgg;   % parameters to determine porosity in Archer (1991)
        mvom;  mvsed;  mvcc;   % molar volumes (cm3 mol-1) mv_i = m_i/rho_i where i = om, sed and cc for organic matter, clay and caco3, respectively
        keq1;  keq2;   % equilibrium const. for h2co3 and hco3 dissociations and functions to calculate them
        co3sat;  keqag;   % co3 conc. at caco3 saturation and solubility product of aragonite
        zox;  zoxx;   % oxygen penetration depth (cm) and its dummy variable
        dep, depi,depf % water depth, i and f denote initial and final values
        corrf; df, err_f;  err_fx;   %  variables to help total fraction of solid materials to converge 1
        dic;  dicx;  alk;  alkx;  o2;  o2x;  %  mol cm-3 porewater; dic, alk and o2 concs., x denotes dummy variables
        dif_dic;  dif_alk;  dif_o2;  % dic, alk and o2 diffusion coeffs inclueing effect of tortuosity
        dbio;  ff;   % biodiffusion coeffs, and formation factor
        co2;  hco3;  co3;  pro;  % co2, hco3, co3 and h+ concs.
        co2x;  hco3x;  co3x;  prox; co3i;  % dummy variables and initial co3 conc.
        poro;  rho;  frt;   % porositiy, bulk density and total volume fraction of solid materials
        sporo;  sporoi;  porof;  sporof;   % solid volume fraction (1 - poro), i and f denote the top and bottom values of variables
        rcc;   % dissolution rate of caco3
        drcc_dcc;  drcc_ddic; % derivatives of caco3 dissolution rate wrt caco3 and dic concs.
        drcc_dalk;  drcc_dco3;  % derivatives of caco3 dissolution rate wrt alk and co3 concs.
        dco3_ddic;  dco3_dalk;   % derivatives of co3 conc. wrt dic and alk concs.
        ddum;   % dummy variable
        dpro_dalk;  dpro_ddic;   % derivatives of h+ conc. wrt alk and dic concs.
        kcc;   % caco3 dissolution rate consts.
        kom;   % degradation rate consts.
        oxco2; anco2;   % oxic and anoxic om degradation rate
        w;  , wi, dw;  wx;  % burial rate, burial rate initial guess, burial rate change, burial rate dummy
        err_w, wxx;  % err in burial rate, dummy dummy burial rate
        err_f_min;  dfrt_df;  d2frt_df2;  dfrt_dfx;  err_w_min;   % variables to minimize errors in burial rate and total fractions of solid phases
        z;  dz; % depth, individual sediment layer thickness
        eta;  beta;   % parameters to make a grid
        dage;  age;   % individual time span and age of sediment grids
        % #ifdef sense
        ztot = 50d0;  % cm , total sediment thickness
        % #else
        ztot = 500d0;  % cm
        % #endif
        nsp = 3;   % independent chemical variables, this does not have to be decided here
        nmx;       % row (and col) number of matrix created to solve linear difference equations
        amx; ymx; emx;  % amx and ymx correspond to A and B in Ax = B, but ymx is also x when Ax = B is solved. emx is array of error
        ipiv; dumx; % matrix used to solve linear system Ax = B
        infobls;  infosbr;   % variables used to tell errors when calling a subroutine to solve matrix
        error;  error2;  minerr;   %  errors in iterations and minimum error produced
        tol = 1d-6;   % tolerance of error
        iz;  row;  col;  itr;  it;  iiz;  itr_w; itr_f;  % integers for sediment grid, matrix row and col and iteration numbers
        cntsp;   % counting caco3 species numbers
        izrec; izrec2;   % grid number where signal is read
        nt = 1000000;   % maximum interation for time (do not have to be defined )
        nrec = 15;   % total recording time of sediment profiles
        cntrec;  itrec;   % integers used for counting recording time
        izox_minerr =0;   % grid number where error in zox is minimum
        up;  dwn;  cnr;  % advection calc. schemes; up or down wind, or central schemes if 1
        adf;   % factor to make sure mass conversion
        time;  dt = 1d2;   % time and time step
        rectime;  time_max;  % recording time and maximum time
        dumreal;   %  dummy variable
        time_spn;  time_trs;  time_aft;   % time durations of spin-up, signal transition and after transition
        % fluxes, adv, dec, dis, dif, res, t and rain denote burial, decomposition, dissoution, diffusion, residual, time change and rain fluxes, respectively
        omadv, omdec, omdif, omrain, omres, omtflx % fluxes of om
        o2tflx, o2dec, o2dif, o2res  % o2 fluxes
        cctflx(nspcc), ccdis(nspcc), ccdif(nspcc), ccadv(nspcc), ccres(nspcc),ccrain(nspcc) % caco3 fluxes
        dictflx, dicdis, dicdif, dicdec, dicres  % dic fluxes
        alktflx, alkdis, alkdif, alkdec, alkres  % alk fluxes
        pttflx, ptdif, ptadv, ptres, ptrain  % clay fluxes
        trans(nz,nz,nspcc+2)  % transition matrix
        transdbio(nz,nz), translabs(nz,nz) % transition matrices created assuming Fickian mixing and LABS simulation
        transturbo2(nz,nz), translabs_tmp(nz,nz) % transition matrices assuming random mixing and LABS simulation
        character*512 workdir, filechr  % work directory and created file names
        character*25 dumchr(3)  % character dummy variables
        character*25 arg, chr(3,4)  % used for reading variables and dummy variables
        integer(kind=4) dumint(8)  % dummy integer
        integer(kind=4) idp, izox  % integer for depth and grid number of zox
        integer(kind=4) narg, ia  % integers for getting input variables
        integer(kind=4) izml, isp, ilabs, nlabs % grid % # of bottom of mixed layer, % # of caco3 species, % # of labs simulation and total % # of labs simulations
        file_tmp=100,file_ccflx=101,file_omflx=102,file_o2flx=103,file_dicflx=104,file_alkflx=105,file_ptflx=106  %  file % #
        file_err=107,file_ccflxes(nspcc), file_bound=108, file_totfrac=109, file_sigmly=110,file_sigmlyd=111  % file % #
        file_sigbtm=112 % file % #
        external dgesv  % subroutine in BALS library
        oxic = .true.  % oxic only model of OM degradation by Emerson (1985)
        anoxic = .true.  % oxic-anoxic model of OM degradation by Archer (1991)
        nobio(nspcc+2) = .false.  % no biogenic reworking assumed
        turbo2(nspcc+2) = .false.  % random mixing
        labs(nspcc+2) = .false.  % mixing info from LABS
        nonlocal(nspcc+2)  % ON if assuming non-local mixing (i.e., if labs or turbo2 is ON)
        flg_500;
        %%%% pi=4.0d0*atan(1.0d0) %
% %         % when using sparse matrix solver
% %         n; nnz;  % number of row (and col) and total number of component not zero
% %         ai;  ap;  % storing row number and total number of non-zero components by a given col
% %         ax;  bx;  % values of non-zero component and B in Ax = B
% %         control;
% %         i,j,cnt,cnt2
% %         info(90)
% %         integer(kind=8) numeric
% %         integer(kind=4) status
% %         integer(kind=8) symbolic
% %         integer(kind=4) sys
% %         real(kind=8), allocatable :: kai;   % x in Ax = B      
        
    end
    
    
end