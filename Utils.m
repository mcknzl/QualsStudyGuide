%%% Mattingly's Master Equation GUI
%%%     Useful equations
%%%     All units are in imperial (lb, ft, s) and speed is in knots
classdef Utils
    properties (SetAccess = private, GetAccess = public)
        fps_per_kt = 1.68781;
        g0 = 32.2; % ft/s^2
        Alt = [0	1000	2000	3000	4000	5000	6000	7000	8000	9000	10000	11000	12000	13000	14000	15000	16000	17000	18000	19000	20000	21000	22000	23000	24000	25000	26000	27000	28000	29000	30000	31000	32000	33000	34000	35000	36000	37000	38000	39000	40000	41000	42000	43000	44000	45000	46000	47000	48000	49000	50000	51000	52000	53000	54000	55000];
        Tf = [59	55.4	51.8	48.38	44.78	41.18	37.58	33.98	30.56	26.96	23.36	19.76	16.16	12.56	9.14	5.54	1.94	-1.66	-5.26	-8.68	-12.28	-15.88	-19.48	-23.08	-26.5	-30.1	-33.7	-37.3	-40.9	-44.5	-47.92	-51.52	-55.12	-58.72	-62.32	-65.74	-69.34	-69.7	-69.7	-69.7	-69.7	-69.7	-69.7	-69.7	-69.7	-69.7	-69.7	-69.7	-69.7	-69.7	-69.7	-69.7	-69.7	-69.7	-69.7	-69.7];
        Tc = [15	13	11	9.1	7.1	5.1	3.1	1.1	-0.8	-2.8	-4.8	-6.8	-8.8	-10.8	-12.7	-14.7	-16.7	-18.7	-20.7	-22.6	-24.6	-26.6	-28.6	-30.6	-32.5	-34.5	-36.5	-38.5	-40.5	-42.5	-44.4	-46.4	-48.4	-50.4	-52.4	-54.3	-56.3	-56.5	-56.5	-56.5	-56.5	-56.5	-56.5	-56.5	-56.5	-56.5	-56.5	-56.5	-56.5	-56.5	-56.5	-56.5	-56.5	-56.5	-56.5	-56.5];
        Ppsi = [14.7	14.17	13.67	13.17	12.69	12.23	11.78	11.34	10.92	10.51	10.1	9.72	9.35	8.99	8.63	8.29	7.97	7.65	7.34	7.04	6.75	6.47	6.21	5.95	5.7	5.45	5.22	4.99	4.78	4.57	4.36	4.17	3.98	3.8	3.63	3.46	3.3	3.14	2.99	2.58	2.72	2.86	2.609	2.526	2.443	2.36	2.277	2.194	2.111	2.028	1.945	1.862	1.779	1.696	1.613	1.53];
        PP0 = [1	0.9644	0.9298	0.8962	0.8637	0.832	0.8014	0.7716	0.7428	0.7148	0.6877	0.6614	0.636	0.6113	0.5875	0.5643	0.542	0.5203	0.4994	0.4791	0.4595	0.4406	0.4223	0.4046	0.3876	0.3711	0.3552	0.3398	0.325	0.3107	0.297	0.2837	0.2709	0.2586	0.2467	0.2353	0.2243	0.2138	0.2038	0.1942	0.1851	0.176	0.16629	0.15686	0.14743	0.138	0.12857	0.11914	0.10971	0.10028	0.09085	0.08142	0.07199	0.06256	0.05313	0.0437];
        RR0 = [1	0.9711	0.94728	0.9151	0.8881	0.8617	0.8359	0.8106	0.786	0.762	0.7385	0.7156	0.6932	0.6713	0.65	0.6292	0.609	0.5892	0.5699	0.5511	0.5328	0.515	0.4976	0.4806	0.4642	0.4481	0.4325	0.4173	0.4025	0.3881	0.3741	0.3605	0.3473	0.3345	0.322	0.3099	0.2981	0.2844	0.271	0.2583	0.2462	0.2341	0.22118	0.20864	0.1961	0.18356	0.17102	0.15848	0.14594	0.1334	0.12086	0.10832	0.09578	0.08324	0.0707	0.05816];
        A = [661	659	656	654	652	650	647	645	643	640	638	636	633	631	628	626	624	621	619	616	614	611	609	608	604	602	599	597	594	591	589	586	584	581	579	576	573	573	573	573	573	573	573	573	573	573	573	573	573	573	573	573	573	573	573	573];
        rho0 = 0.002739;
        % Engine Parameters
        TSFC0 = 1.014; % slug/ft^3
        Tnorm = 5744; % lbf
        Tmil = 6318.4; % lbf
        Twet = 8616; % lbf
        % Aerodynamics
        CLmax = 1.4;
        CD0 = 0.02;
        AR = 4;
        e = 0.99;
        K;
        K1 = 0.001;
        CD;
        CL;
        LbyD;
        % Wing Loading
        WS = linspace(1,100,1000); % lb/ft^2
        % Curve fits
        fTf;
        fTc;
        fPpsi;
        fPP0;
        fRR0;
        fA;
    end
    
    methods
        function obj = Utils(varargin)
            obj.fTf = fit(obj.Alt',obj.Tf','linearinterp');
            obj.fTc = fit(obj.Alt',obj.Tc','linearinterp');
            obj.fPpsi = fit(obj.Alt',obj.Ppsi','linearinterp');
            obj.fPP0 = fit(obj.Alt',obj.PP0','linearinterp');
            obj.fRR0 = fit(obj.Alt',obj.RR0','linearinterp');
            obj.fA = fit(obj.Alt',obj.A','linearinterp');
            if nargin>0
                inStruct = varargin{1};
                fieldNames = fieldnames(inStruct);
                for i = 1:numel(fieldNames)
                    try
                        field = fieldNames{i};
                        obj.(field) = inStruct.(field);
                    catch
                        errMsg = ['Invalid field name [' field ']: Please fix input structure and try again'];
                        error(errMsg);
                    end
                end
            end
            % Recalculate things
            obj.CD = linspace(obj.CD0,1);
            obj.K = 1/(pi*obj.e*obj.AR);
            obj.CL = sqrt((obj.CD-obj.CD0)./obj.K);
            obj.LbyD = obj.CL./obj.CD;
        end
        
        function TW = MEQ(obj,b,h,V,n,k1,k2,dhdt,dVdt,type)
            rho = obj.rho0 .* obj.fRR0(h);
            a = obj.getAlpha(V,h,type);
            q = 0.5*rho.*V.^2;
            TW = (b./a).*( (q./(b.*obj.WS)).*( obj.CD0 + (n.*k1.*b./q).*obj.WS + ...
                (n.^2.*k2.*b.^2./q.^2).*(obj.WS).^2) + (1./V).*dhdt + dVdt/obj.g0);
        end
        
        function alpha = getAlpha(obj,V,h,type)
            M = V/(obj.fA(h)*obj.fps_per_kt);
            rr0 = obj.fRR0(h);
            switch type
                case 'wet'
                    alpha = 0.76*(0.907+0.262*abs(M-0.5).^(1.5)).*(rr0^0.7);
                case 'dry'
                    alpha = (0.952+0.3*(M-0.4).^2)*(1^0.7);
            end
        end
        
        function W = weightEst(obj,Winit,payload,pilot)
            We = 2.34*Winit^(-0.13); % empty weight regression
            
        end
        
        function [TW,betaOut] = calcTakeoff(obj,betaIn,h,hClear,kTO,TOFL,T,ksi,mu)
            % T/W
            sigma = obj.fRR0(h);
            rho = obj.rho0 * sigma;
            vTO = sqrt(2*obj.WS./(rho*obj.CLmax));
            rLO = vTO.^2/(obj.g0*0.19);
            thetaRot = acos(1-hClear./rLO);
            clearDist = rLO.*sin(thetaRot);
            TW = (betaIn^2*kTO^2)./(sigma*(TOFL-clearDist)*rho*obj.g0*obj.CLmax).*obj.WS;
            % Beta
            vTOdes = vTO(obj.WS == obj.WSdes);
            alphaDes = obj.getAlpha(vTOdes,h,'dry');
            theta = (T+459)/(obj.fTf(h)+459);
            TSFC = (1.5+0.23*vTOdes/obj.fA(h))*sqrt(theta);
            u = (betaIn/(alphaDes*obj.TWdes))*(0.5*ksi*rho*vTOdes^2/(betaIn*obj.WSdes)+mu);
            betaOut = exp(-(TSFC/3600)*obj.WSdes/(obj.g0*(1-u)));
        end
                
        function TW = calcClimb(obj,betaIn,h,V,ROC,type)
            TW = obj.MEQ(betaIn,h,V,1,obj.K1,obj.K,ROC,0,type);
        end
        
        function TW = calcCruise(obj,betaIn,h,V)
            TW = obj.MEQ(betaIn,h,V,1,obj.K1,obj.K,0,0,'dry');
        end
        
        function TW = calcStdyLvlTurn_Gee(obj,betaIn,h,V,n,type)
            TW = obj.MEQ(betaIn,h,V,n,obj.K1,obj.K,0,0,type);
        end
        function TW = calcStdyLvlTurn_Rate(obj,betaIn,h,V,omega,type)
            n = sqrt(1+(omega*pi*V/(180*obj.g0))^2);
            TW = obj.MEQ(betaIn,h,V,n,obj.K1,obj.K,0,0,type);
        end
        function TW = calcStdyLvlTurn_Radius(obj,betaIn,h,V,r,type)
            n = sqrt(1+(V^2/(obj.g0*r))^2);
            TW = obj.MEQ(betaIn,h,V,n,obj.K1,obj.K,0,0,type);
        end
        
        function TW = calcLvlAccel(obj,betaIn,h,V0,V1,t,type)
            avgAccel = (V1-V0)/t;
            avgV = 0.5*(V0+V1);
            n = 1;
            TW = obj.MEQ(betaIn,h,avgV,n,obj.K1,obj.K,0,avgAccel,type);
        end
        
        function WScrit = calcLanding(obj,betaIn,h,Vapp,kTD)
            Vcrit = Vapp/kTD;
            rho = obj.rho0*obj.fRR0(h);
            WScrit = Vcrit^2*rho*obj.CLmax/(2*betaIn);
        end
                
        
    end
end