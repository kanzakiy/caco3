classdef caco3_therm
    % functions and subroutines to calculate caco3 thermodynamics
    
    properties
    end
    
    methods(Static)
        
        function calceq1_val = calceq1(tmp,sal,dep) % Millero et al. (2006, MC); Pressure effect from Zeebe and Wolf-Gladrow (2001)
            
            tmp_k=tmp+273.15d0;
            pres = dep*100.0d0; % bar
            
            coeff(1)=13.4191D0;
            coeff(2)=0.0331D0;
            coeff(3)=-5.33D-5;
            coeff(4)=-530.1228D0;
            coeff(5)=-6.103d0;
            coeff(6)=-2.06950d0;
            
            calceq1_val=-126.34048d0+6320.813d0/tmp_k+19.568224*log(tmp_k);
            calceq1_val=calceq1_val+coeff(1)*sqrt(sal)+coeff(2)*sal+coeff(3)*sal^2.0d0 ...
                +(coeff(4)*sqrt(sal)+coeff(5)*sal)/tmp_k   ...
                +coeff(6)*sqrt(sal)*log(tmp_k);
            calceq1_val=10.0d0^(-calceq1_val);
            
            % see Zeebe and Wolf-Gladrow (2001) Appendix A.11)
%           calceq1_val=calceq1_val*exp((-(-15.82d0-0.0219d0*tmp)*pres+0.5d0*(1.13d-3-0.1475d-3*tmp)*pres*pres)/83.131d0/tmp_k);     % wrong
            calceq1_val=calceq1_val*exp((-(-25.50d0+0.1271d0*tmp)*pres+0.5d0*(-3.08d-3+0.0877d-3*tmp)*pres*pres)/83.131d0/tmp_k);     % correct
            
        end %function calceq1
        
        function calceq2_val = calceq2(tmp,sal,dep) % Millero et al. (2006, MC); Pressure effect from Zeebe and Wolf-Gladrow (2001)
            
            tmp_k=tmp+273.15d0;
            pres = dep*100.0d0; % bar
            
            coeff(1)=21.0894D0;
            coeff(2)=0.1248D0;
            coeff(3)=-0.0003687D0;
            coeff(4)=-772.483d0;
            coeff(5)=-20.051D0;
            coeff(6)=-3.32254d0;
            
            calceq2_val=-90.18333d0+5143.692d0/tmp_k+14.613358*log(tmp_k);
            calceq2_val=calceq2_val+coeff(1)*sqrt(sal)+coeff(2)*sal+coeff(3)*sal^2.0d0 ...
                +(coeff(4)*sqrt(sal)+coeff(5)*sal)/tmp_k   ...
                +coeff(6)*sqrt(sal)*log(tmp_k);
            calceq2_val=10.0d0^(-calceq2_val);
            
           % see Zeebe and Wolf-Gladrow (2001) Appendix A.11)
%           calceq2_val=calceq2_val*exp((-(-25.50d0+0.1271d0*tmp)*pres+0.5d0*(-3.08d-3+0.0877d-3*tmp)*pres*pres)/83.131d0/tmp_k);      % wrong
            calceq2_val=calceq2_val*exp((-(-15.82d0-0.0219d0*tmp)*pres+0.5d0*(1.13d-3-0.1475d-3*tmp)*pres*pres)/83.131d0/tmp_k);      % correct
            
            
        end %function calceq2
        
        
        function keqcc = calceqcc(tmp,sal,dep) % Mucci (1983) cited by Zeebe and Wolf-Gladrow (2001)
            %% calcite solubility
            %% See e.g. Zeebe and Wolf-Gladrow (2001) Appendix A.10
            tmp_k=tmp+273.15;
            pres = dep*100.0; % bar
            
            keqcc = -171.9065 - 0.077993*tmp_k + 2839.319/tmp_k ...
                +71.595*log10(tmp_k) ...
                +(-0.77712+0.0028426*tmp_k+178.34/tmp_k)*sal^0.5 ...
                -0.07711*sal+0.0041249*sal^1.5;
            keqcc = 10^keqcc;
            
            % effect of pressure on equil. constants (values see e.g. Zeebe and Wolf-Gladrow (2001) Table A.11.1)
            keqcc=keqcc*exp((-(-48.76+0.5304*tmp)*pres+0.5*(-11.76d-3+0.3692d-3*tmp)*pres*pres)/83.131/tmp_k);
            
        end
        
        function keqag = calceqag(tmp,sal,dep) % Mucci (1983) cited by Zeebe and Wolf-Gladrow (2001)
            %% aragonite solubility
            %% See e.g. Zeebe and Wolf-Gladrow (2001) Appendix A.10
            
            tmp_k=tmp+273.15;
            pres = dep*100.0; % bar
            
            keqag = -171.945 - 0.077993*tmp_k + 2903.293/tmp_k ...
                +71.595*log10(tmp_k) ...
                +(-0.068393+0.0017276*tmp_k+88.135/tmp_k)*sqrt(sal) ...
                -0.10018*sal+0.0059415*sal^1.5;
            keqag = 10^keqag;
            
            % effect of pressure on equil. constants (values see e.g. Zeebe and Wolf-Gladrow (2001) Table A.11.1)
            keqag=keqag*exp((-(-46.00+0.5304*tmp)*pres+0.5*(-11.76d-3+0.3692d-3*tmp)*pres*pres)/83.131/tmp_k);
            
        end
        
        function [ph,co2,hco3,co3,info] = calcspecies(dic,alk,tmp,sal,dep)
            % subroutine to calculate all aqueous co2 species and pH
            
            info=0;
            co2 = zeros(1, length(alk));
            hco3 = zeros(1, length(alk));
            co3 = zeros(1, length(alk));
            
            k1=caco3_therm.calceq1(tmp,sal,dep);
            k2=caco3_therm.calceq2(tmp,sal,dep);
            a=1.0d0;
            b=(1.0d0-dic./alk)*k1;
            c=(1.0d0-2.0d0*dic./alk)*k1*k2;
            ph = (-b+sqrt(b.*b-4.0d0.*a.*c))./2.0d0./a;
            if (any(ph<0d0))
                fprintf('... unable to calculate ph - calcspecies \n');
                % print*,dic
                % print*,alk
                % print*,ph
                % stop
                info=1;
                return
            end
            co2 = alk./(k1./ph+2.0d0.*k1.*k2./ph./ph);
            % hco3=co2*k1/ph
            hco3=alk./(1.0d0+2d0.*k2./ph);
            % co3=hco3*k2/ph
            co3=alk./(ph./k2+2.0d0);
            % ph=-log10(ph)
            
            
        end
        
        
        function [dco3_dalk,dco3_ddic,info] = calcdevs(dic,alk,tmp,sal,dep)
            % subroutine to calculate derivatives of co3 conc. wrt dic and alk (defined as dco3_ddic and dco3_dalk, respectively)
            
            info=0;
            k1=caco3_therm.calceq1(tmp,sal,dep);
            k2=caco3_therm.calceq2(tmp,sal,dep);
            %             fprintf('%17.16e %17.16e\n', k1, k2);
            a=1.0d0;
            b=(1.0d0-dic./alk)*k1;
            c=(1.0d0-2.0d0*dic./alk)*k1*k2;
            ph = (-b+sqrt(b.*b-4.0d0.*a.*c))./2.0d0./a;
            if (any(ph<0d0))
                fprintf('... unable to calculate ph - calcdevs \n');
                % print*,dic
                % print*,alk
                % print*,ph
                % stop
                info=1;
                return
            end
            co2 = alk./(k1./ph+2.0d0.*k1.*k2./ph./ph);
            % hco3=co2*k1/ph
            hco3=alk./(1.0d0+2d0.*k2./ph);
            % co3=hco3*k2/ph
            co3=alk./(ph./k2+2.0d0);
            % ph=-log10(ph)
            
            db_dalk = k1.*(-1d0).*dic.*(-1d0)./alk./alk;
            dc_dalk = k1.*k2.*(-2d0).*dic.*(-1d0)./alk./alk;
            db_ddic = k1.*(-1d0./alk);
            dc_ddic = k1.*k2.*(-2d0./alk);
            dph_dalk = -0.5d0.*db_dalk + 0.5d0.*0.5d0.*(b.*b-4d0.*c).^(-0.5d0).*(2d0.*b.*db_dalk - 4d0.*dc_dalk);
            dph_ddic = -0.5d0.*db_ddic + 0.5d0.*0.5d0.*(b.*b-4d0.*c).^(-0.5d0).*(2d0.*b.*db_ddic - 4d0.*dc_ddic);
            dco3_dalk = 1d0./(ph./k2+2d0) + alk.*(-1d0)./((ph./k2+2d0).^2d0).*(1d0./k2).*dph_dalk;
            dco3_ddic = 0d0./(ph./k2+2d0) + alk.*(-1d0)./((ph./k2+2d0).^2d0).*(1d0./k2).*dph_ddic;
        end
        
    end
    
end