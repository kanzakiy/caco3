classdef lysocline_exp
%
%   *******************************************************************   %
%   ***  Running the lysocline experiments ****************************   %
%   *******************************************************************   %
%
%    scripts to run the lysocline experiments with oxic-only and
%    oxic-anoxic OM degradation model as in the manuscript

    properties
    end
    
    methods(Static)
        
        function run_all_lysocline_exp_ox(folder)
            %% do oxic-only lysocline experiment as in manuscript
            
            % specify folder where results should be saved (same for oxic-only and oxic-anoxic for plotting reasons, e.g.:
            % run_all_lysocline_exp_ox('./1207_lysocline')
            
            rainratio = [0.0,0.5,0.6666,1.0,1.5];
            dt = 1d8;
            switch_oxonly = [true,false];
            
            nz = 25;
                        
            for i = 1:nz
                for j = 1:10
                    for k = 1:length(rainratio)
                        for l = 1:1 %length(switch_oxonly)
                            cc_rain_flx = j*6e-6
                            dep = i*6.0/25
                            rr = rainratio(k)
                            run_sig_iso_dtchange(cc_rain_flx, rainratio(k), dep, dt, switch_oxonly(l), folder);
%                            caco3_test.chk_through_sig_iso_dtchange(cc_rain_flx, rainratio(k), dep, dt, switch_oxonly(l), folder);
                        end
                    end
                end
            end
            % set
            
            
        end
        
        function run_all_lysocline_exp_oxanox(folder)
            %% do oxic-anoxic lysocline experiment as in manuscript
            
            % specify folder where results should be saved (same for oxic-only and oxic-anoxic for plotting reasons, e.g.:
            % run_all_lysocline_exp_ox('./1207_lysocline')
            
            rainratio = [0.0,0.5,0.6666,1.0,1.5];
            dt = 1d8;
            switch_oxonly = [true,false];
            
            nz = 25;
                        

            for i = 1:nz
                for j = 1:10
                    for k = 1:length(rainratio)
                        for l = 2:2 %length(switch_oxonly)
                            cc_rain_flx = j*6e-6
                            dep = i*6.0/25
                            rr = rainratio(k)
                            run_sig_iso_dtchange(cc_rain_flx, rainratio(k), dep, dt, switch_oxonly(l), folder);
%                            caco3_test.chk_through_sig_iso_dtchange(cc_rain_flx, rainratio(k), dep, dt, switch_oxonly(l), folder);
                        end
                    end
                end
            end
            % set
            
            
        end        
        
        
    end
    
end