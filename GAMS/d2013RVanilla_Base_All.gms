$ontext
This is the DICE-2013R model, version DICE2013R_100413_vanilla.gms, revised from April version.
The vanilla version includes only the optimal and baseline scenarios.
These are determined by setting the "ifopt" control at 1 (optimal) or 0 (baseline).
This version has write ("put") output but does not have subroutines ("include").
A full discussion is included in the "DICE 2013R Manual" on the web at dicemodel.net.
$offtext

$title        DICE-2013R October 2013

set        t  Time periods (5 years per period)                    /1*60/ ;

parameters

**Time Step
        tstep    Years per Period                                    /5/

** If optimal control
        ifopt    If optimized 1 and if base is 0                     /0/

** Preferences
        elasmu   Elasticity of marginal utility of consumption     /  1.45 /
        prstp    Initial rate of social time preference per year   / .015  /

** Population and technology
        gama     Capital elasticity in production function        /.300    /
        pop0     Initial world population (millions)              /6838    /
        popadj   Growth rate to calibrate to 2050 pop projection  /0.134   /
        popasym  Asymptotic population (millions)                 /10500   /
        dk       Depreciation rate on capital (per year)          /.100    /
        q0       Initial world gross output (trill 2005 USD)      /63.69   /
        k0       Initial capital value (trill 2005 USD)           /135     /
        a0       Initial level of total factor productivity       /3.80    /
        ga0      Initial growth rate for TFP per 5 years          /0.079   /
        dela     Decline rate of TFP per 5 years                  /0.006   /

** Emissions parameters
        gsigma1  Initial growth of sigma (continuous per year)        /-0.01   /
        dsig     Decline rate of decarbonization per period           /-0.001  /
        eland0   Carbon emissions from land 2010 (GtCO2 per year)     / 3.3    /
        deland   Decline rate of land emissions (per period)          / .2     /
        e0       Industrial emissions 2010 (GtCO2 per year)           /33.61   /
        miu0     Initial emissions control rate for base case 2010    /.039    /

** Carbon cycle
* Initial Conditions
        mat0   Initial Concentration in atmosphere 2010 (GtC)        /830.4   /
        mu0    Initial Concentration in upper strata 2010 (GtC)      /1527.   /
        ml0    Initial Concentration in lower strata 2010 (GtC)      /10010.  /
        mateq  Equilibrium concentration atmosphere  (GtC)           /588     /
        mueq   Equilibrium concentration in upper strata (GtC)       /1350    /
        mleq   Equilibrium concentration in lower strata (GtC)       /10000   /

* Flow paramaters
        b12      Carbon cycle transition matrix                      /.088/
        b23      Carbon cycle transition matrix                      /0.00250/

* These are for declaration and are defined later
        b11      Carbon cycle transition matrix
        b21      Carbon cycle transition matrix
        b22      Carbon cycle transition matrix
        b32      Carbon cycle transition matrix
        b33      Carbon cycle transition matrix
        sig0     Carbon intensity 2010 (kgCO2 per output 2005 USD 2010)

** Climate model parameters
        t2xco2   Equilibrium temp impact (oC per doubling CO2)    / 2.9   /
        fex0     2010 forcings of non-CO2 GHG (Wm-2)              / 0.25   /
        fex1     2100 forcings of non-CO2 GHG (Wm-2)              / 0.70   /
        tocean0  Initial lower stratum temp change (C from 1900)  /.0068  /
        tatm0    Initial atmospheric temp change (C from 1900)    /0.80   /
        c1       Climate equation coefficient for upper level     /0.098  /
        c3       Transfer coefficient upper to lower stratum      /0.088  /
        c4       Transfer coefficient for lower level             /0.025  /
        fco22x   Forcings of equilibrium CO2 doubling (Wm-2)      /3.8    /

** Climate damage parameters
        a1        Damage intercept                                 /0       /
        a2        Damage quadratic term                            /0.00267 /
        a3        Damage exponent                                  /2.00    /

** Abatement cost
        expcost2  Exponent of control cost function               / 2.8  /
        pback     Cost of backstop 2005$ per tCO2 2010            / 344  /
        gback     Initial cost decline backstop cost per period   / .025 /
        limmiu    Upper limit on control rate after 2150          / 1.2  /
        tnopol    Period before which no emissions controls base  / 45   /
        cprice0   Initial base carbon price (2005$ per tCO2)      / 1.0  /
        gcprice   Growth rate of base carbon price per year       /.02   /

** Participation parameters
        periodfullpart Period at which have full participation           /21  /
        partfract2010  Fraction of emissions under control in 2010       / 1  /
        partfractfull  Fraction of emissions under control at full time  / 1  /

** Availability of fossil fuels
        fosslim        Maximum cumulative extraction fossil fuels (GtC) /6000/

** Scaling and inessential parameters
* Note that these are unnecessary for the calculations but are for convenience
        scale1      Multiplicative scaling coefficient              /0.016408662 /
        scale2      Additive scaling coefficient                    /-3855.106895/ ;

* Program control variables
sets     tfirst(t), tlast(t), tearly(t), tlate(t);

PARAMETERS
        L(t)          Level of population and labor
        al(t)         Level of total factor productivity
        sigma(t)      CO2-equivalent-emissions output ratio
        rr(t)         Average utility social discount rate
        ga(t)         Growth rate of productivity from
        forcoth(t)    Exogenous forcing for other greenhouse gases
        gl(t)         Growth rate of labor
        gcost1        Growth of cost factor
        gsig(t)       Change in sigma (cumulative improvement of energy efficiency)
        etree(t)      Emissions from deforestation
        cost1(t)      Adjusted cost for backstop
        partfract(t)  Fraction of emissions in control regime
        lam           Climate model parameter
        gfacpop(t)    Growth factor population
        pbacktime(t)  Backstop price
        optlrsav      Optimal long-run savings rate used for transversality
        scc(t)        Social cost of carbon
        cpricebase(t) Carbon price in base case  ;

* Program control definitions
        tfirst(t) = yes$(t.val eq 1);
        tlast(t)  = yes$(t.val eq card(t));

* Parameters for long-run consistency of carbon cycle
        b11 = 1 - b12;
        b21 = b12*MATEQ/MUEQ;
        b22 = 1 - b21 - b23;
        b32 = b23*mueq/mleq;
        b33 = 1 - b32 ;

* Further definitions of parameters
        sig0 = e0/(q0*(1-miu0));
        lam = fco22x/ t2xco2;
        L("1") = pop0;
        loop(t, L(t+1)=L(t););
        loop(t, L(t+1)=L(t)*(popasym/L(t))**popadj ;);

        ga(t)=ga0*exp(-dela*5*((t.val-1)));
        al("1") = a0; loop(t, al(t+1)=al(t)/((1-ga(t))););

        gsig("1")=gsigma1; loop(t,gsig(t+1)=gsig(t)*((1+dsig)**tstep) ;);
        sigma("1")=sig0;   loop(t,sigma(t+1)=(sigma(t)*exp(gsig(t)*tstep)););

        pbacktime(t)=pback*(1-gback)**(t.val-1);
        cost1(t) = pbacktime(t)*sigma(t)/expcost2/1000;

        etree(t) = eland0*(1-deland)**(t.val-1);
        rr(t) = 1/((1+prstp)**(tstep*(t.val-1)));
        forcoth(t) = fex0+ (1/18)*(fex1-fex0)*(t.val-1)$(t.val lt 19)+ (fex1-fex0)$(t.val ge 19);
        optlrsav = (dk + .004)/(dk + .004*elasmu + prstp)*gama;

        partfract(t)$(ord(T)>periodfullpart) = partfractfull;
        partfract(t)$(ord(T)<periodfullpart+1) = partfract2010+(partfractfull-partfract2010)*(ord(t)-1)/periodfullpart;

        partfract("1")= partfract2010;

        cpricebase(t)= cprice0*(1+gcprice)**(5*(t.val-1));


VARIABLES
        MIU(t)          Emission control rate GHGs
        FORC(t)         Increase in radiative forcing (watts per m2 from 1900)
        TATM(t)         Increase temperature of atmosphere (degrees C from 1900)
        TOCEAN(t)       Increase temperatureof lower oceans (degrees C from 1900)
        MAT(t)          Carbon concentration increase in atmosphere (GtC from 1750)
        MU(t)           Carbon concentration increase in shallow oceans (GtC from 1750)
        ML(t)           Carbon concentration increase in lower oceans (GtC from 1750)
        E(t)            Total CO2 emissions (GtCO2 per year)
        EIND(t)         Industrial emissions (GtCO2 per year)
        C(t)            Consumption (trillions 2005 US dollars per year)
        K(t)            Capital stock (trillions 2005 US dollars)
        CPC(t)          Per capita consumption (thousands 2005 USD per year)
        I(t)            Investment (trillions 2005 USD per year)
        S(t)            Gross savings rate as fraction of gross world product
        RI(t)           Real interest rate (per annum)
        Y(t)            Gross world product net of abatement and damages (trillions 2005 USD per year)
        YGROSS(t)       Gross world product GROSS of abatement and damages (trillions 2005 USD per year)
        YNET(t)         Output net of damages equation (trillions 2005 USD per year)
        DAMAGES(t)      Damages (trillions 2005 USD per year)
        DAMFRAC(t)      Damages as fraction of gross output
        ABATECOST(t)    Cost of emissions reductions  (trillions 2005 USD per year)
        MCABATE(t)      Marginal cost of abatement (2005$ per ton CO2)
        CCA(t)          Cumulative industrial carbon emissions (GTC)
        PERIODU(t)      One period utility function
        CPRICE(t)       Carbon price (2005$ per ton of CO2)
        CEMUTOTPER(t)   Period utility
        UTILITY         Welfare function
;

NONNEGATIVE VARIABLES  MIU, TATM, MAT, MU, ML, Y, YGROSS, C, K, I;

EQUATIONS
*Emissions and Damages
        EEQ(t)           Emissions equation
        EINDEQ(t)        Industrial emissions
        CCACCA(t)        Cumulative carbon emissions

        FORCE(t)         Radiative forcing equation
        DAMFRACEQ(t)     Equation for damage fraction
        DAMEQ(t)         Damage equation

        ABATEEQ(t)       Cost of emissions reductions equation
        MCABATEEQ(t)     Equation for MC abatement
        CARBPRICEEQ(t)   Carbon price equation from abatement

*Climate and carbon cycle
        MMAT(t)          Atmospheric concentration equation
        MMU(t)           Shallow ocean concentration
        MML(t)           Lower ocean concentration
        TATMEQ(t)        Temperature-climate equation for atmosphere
        TOCEANEQ(t)      Temperature-climate equation for lower oceans

*Economic variables
        YGROSSEQ(t)      Output gross equation
        YNETEQ(t)        Output net of damages equation
        YY(t)            Output net equation

        CC(t)            Consumption equation
        CPCE(t)          Per capita consumption definition

        SEQ(t)           Savings rate equation
        KK(t)            Capital balance equation
        RIEQ(t)          Interest rate equation

* Utility
        CEMUTOTPEREQ(t)  Period utility
        PERIODUEQ(t)     Instantaneous utility function equation
        UTIL             Objective function      ;

** Equations of the model
*Emissions and Damages
 eeq(t)..             E(t)           =E= EIND(t) + etree(t);
 eindeq(t)..          EIND(t)        =E= sigma(t) * YGROSS(t) * (1-(MIU(t)));
 ccacca(t+1)..        CCA(t+1)       =E= CCA(t)+ EIND(t)*5/3.666;

 force(t)..           FORC(t)        =E= fco22x * ((log((MAT(t)/588.000))/log(2))) + forcoth(t);
 damfraceq(t) ..      DAMFRAC(t)     =E= (a1*TATM(t))+(a2*TATM(t)**a3) ;
 dameq(t)..           DAMAGES(t)     =E= YGROSS(t) * DAMFRAC(t);

 abateeq(t)..         ABATECOST(t)   =E= YGROSS(t) * cost1(t) * (MIU(t)**expcost2) * (partfract(t)**(1-expcost2));
 mcabateeq(t)..       MCABATE(t)     =E= pbacktime(t) * MIU(t)**(expcost2-1);
 carbpriceeq(t)..     CPRICE(t)      =E= pbacktime(t) * (MIU(t)/partfract(t))**(expcost2-1);

*Climate and carbon cycle
 mmat(t+1)..          MAT(t+1)       =E= MAT(t)*b11 + MU(t)*b21 + (E(t)*(5/3.666));
 mml(t+1)..           ML(t+1)        =E= ML(t)*b33  + MU(t)*b23;
 mmu(t+1)..           MU(t+1)        =E= MAT(t)*b12 + MU(t)*b22 + ML(t)*b32;
 tatmeq(t+1)..        TATM(t+1)      =E= TATM(t) + c1 * ((FORC(t+1)-(fco22x/t2xco2)*TATM(t))-(c3*(TATM(t)-TOCEAN(t))));
 toceaneq(t+1)..      TOCEAN(t+1)    =E= TOCEAN(t) + c4*(TATM(t)-TOCEAN(t));

*Economic variables
 ygrosseq(t)..        YGROSS(t)      =E= (al(t)*(L(t)/1000)**(1-GAMA))*(K(t)**GAMA);
 yneteq(t)..          YNET(t)        =E= YGROSS(t)*(1-damfrac(t));
 yy(t)..              Y(t)           =E= YNET(t) - ABATECOST(t);

 cc(t)..              C(t)           =E= Y(t) - I(t);
 cpce(t)..            CPC(t)         =E= 1000 * C(t) / L(t);

 seq(t)..             I(t)           =E= S(t) * Y(t);
 kk(t+1)..            K(t+1)         =L= (1-dk)**tstep * K(t) + tstep * I(t);
 rieq(t+1)..          RI(t)          =E= (1+prstp) * (CPC(t+1)/CPC(t))**(elasmu/tstep) - 1;

*Utility
 cemutotpereq(t)..    CEMUTOTPER(t)  =E= PERIODU(t) * L(t) * rr(t);
 periodueq(t)..       PERIODU(t)     =E= ((C(T)*1000/L(T))**(1-elasmu)-1)/(1-elasmu)-1;
 util..               UTILITY        =E= tstep * scale1 * sum(t,  CEMUTOTPER(t)) + scale2 ;

*Resource limit
CCA.up(t)       = fosslim;

* Control rate limits
MIU.up(t)            = limmiu*partfract(t);
MIU.up(t)$(t.val<30) = 1;

**  Upper and lower bounds for stability

K.LO(t)         = 1;
MAT.LO(t)       = 10;
MU.LO(t)        = 100;
ML.LO(t)        = 1000;
C.LO(t)         = 2;
TOCEAN.UP(t)    = 20;
TOCEAN.LO(t)    = -1;
TATM.UP(t)      = 40;
CPC.LO(t)      = .01;

* Control variables
* Savings rate for asympotic equilibrium
S.FX(t)$(t.val>50) = optlrsav;

* Base carbon price if base, otherwise optimized
* Warning: If parameters are changed, the next equation might make base case infeasible.
* If so, reduce tnopol so that don't run out of resources.
cprice.up(t)$(ifopt=0) = cpricebase(t);
cprice.up(t)$(t.val>tnopol) = 1000;
cprice.up('1')=cpricebase('1');

* Initial conditions
CCA.FX(tfirst)    = 90;
K.FX(tfirst)      = k0;
MAT.FX(tfirst)    = mat0;
MU.FX(tfirst)     = mu0;
ML.FX(tfirst)     = ml0;
TATM.FX(tfirst)   = tatm0;
TOCEAN.FX(tfirst) = tocean0;

** Solution options
option iterlim = 99900;
option reslim = 99999;
option solprint = on;
option limrow = 0;
option limcol = 0;
model  CO2 /all/;
solve co2 maximizing utility using nlp;
solve co2 maximizing utility using nlp;
solve co2 maximizing utility using nlp;

** POST-SOLVE
* Calculate social cost of carbon
scc(t) = -1000*eeq.m(t)/cc.m(t);

** Display at bottom of output for visual inspection
option decimals=2;
display tatm.l,scc,utility.l,cprice.l,y.l, cpc.l;
option decimals=6;
display ri.l,utility.l;

*Describes a file labeled 'results' with the filename "result.txt" in the current directory
file results /results.txt/;     results.nd = 10 ; results.nw = 0 ; results.pw=1200; results.pc=5;
put results;

put "Year" ;
Loop (T, put (2005+(TSTEP*T.val)) );
put / "Period";
Loop (T, put T.val);
put / "Model Status";
put CO2.modelstat;
put / "** World Economic Parameters" ;
put / "* OUTPUT AND CAPTIAL ACCUMULATION" ;
put / "Capital Share" ;
put GAMA;
put / "Depreciation" ;
put DK;
put / "Initial Output (2005 US International $, trillions)" ;
put YGROSS.l("1");
put / "Initial K (Capital Stock)" ;
put K.l("1");
put / "Initial A (Total Factor Productivity)" ;
put AL("1");

put / "* WELFARE" ;
put / "Time Preference plus Risk (per year)" ;
Loop (T, put PRSTP);
put / "Time Premium" ;
put PRSTP;
put / "Social Time Preference Factor" ;
Loop (T, put power((1/(1+PRSTP)),(TSTEP*(1+(ord(T)-2)))));
put / "Elasticity of MU of Consumption" ;
put ELASMU;

put / "* PRODUCTIVITY" ;
put / "Total Factor Productivity" ;
Loop (T, put AL(T));
put / "Initial Total Factor Productivity" ;
put AL("1");
put / "Rate of Growth of Productivity" ;
Loop (T, put GA(T));

put / "* DAMAGE FUNCTION" ;
put / "Damage Coefficient on Temperature" ;
put a1;
put / "Damage Coefficient on Temperature ^2" ;
put a2;
put / "Exponent on Damages" ;
put A3;

put / "* ABATEMENT COST" ;
put / "Abatement cost function coefficient" ;
Loop (T, put (((pbacktime(T)*SIGMA(T)/EXPCOST2)/1000)));
put / " Price of Backstop Technology (1000 USD per ton C)" ;
put pbacktime("1");
put / " Backstop Price (1000 USD per ton CO2)" ;
Loop (T, put pbacktime(T));
put / "Upper limit of Control Rate" ;
put LIMMIU;
put / "Exponent of control cost function" ;
put EXPCOST2;

put / "* EMISSIONS" ;
put / "Sigma (Industrial MTCO2 per 1000 USD [2000]" ;
Loop (T, put SIGMA(T));
put / "Initial Sigma" ;
put SIGMA("1");
put / "Growth Rate of Sigma" ;
Loop (T, put "-0.01000");
put / "Carbon Emissions from Land Use Change" ;
Loop (T, put ETREE(T));
put / "Initial Carbon Emissions from Land" ;
put ELAND0;

put / "* PARTICIPATION" ;
put / "Proportion of World in Participation" ;
Loop (T, put partfract(T));

put / "* SCALING" ;
put / "Multiplicative Scaling Coefficient in Utility Function" ;
put scale1;
put / "Additive Scaling Coefficient in Utility Function" ;
put scale2;

put / "* POPULATION" ;
put / "Population (millions)" ;
Loop (T, put L(T));

put / "** Global Parameters" ;
put / "* CARBON LIMITS" ;
put / "Maximum Carbon Resources (GtC)" ;
put FOSSLIM;

put / "* CARBON CYCLE INITIAL CONDITIONS" ;
put / "Initial Concentration of CO2 (GtC).." ;
put / "...Atmospheric, Year 2007" ;
put MAT0;
put / "...Atmospheric, Year 2012" ;
put "Not Supported ?";
put / "...Biosphere and Shallow Oceans, Year 2008" ;
put MU0;
put / "...Deep Oceans, Year 2008" ;
put ML0;

put / "* CARBON CYCLE PARAMETERS" ;
put / "AtA b11" ;
put b11;
put / "BtA b21" ;
put b21;
put / "AtB b12" ;
put b12;
put / "BtB b22" ;
put b22;
put / "DtB b32" ;
put b32;
put / "BtD b23" ;
put b23;
put / "DtD b33" ;
put b33;

put / "* TEMPERATURE DATA" ;
put / "Exogenous Forcing (Watts per Square Meter)" ;
Loop (T, put FORCOTH(T));
put / "2000 Forcings, non-CO2 GHG" ;
put FEX0;
put / "2100 Forcings, non-CO2 GHG" ;
put FEX1;
put / "Initial Temperature (2008-2011)" ;
put TATM0;
put / "Initial Temperature of Deep Oceans (deg C above 1900 level)" ;
put TOCEAN0;

put / "* CLIMATE MODEL PARAMETERS" ;
put / "Temperature Sensitivity Coefficient (temp increase per doubling CO2)" ;
put T2XCO2;
put / "Forcings at CO2 doubling (Watts per Meter 2)" ;
put FCO22X;

put / "* CLIMATE MODULE TRANSITION PARAMETERS" ;
put / "Speed of Adjustment Parameter for Atmospheric Temperature" ;
put C1;
put / "Coefficient of Heat Loss from Atmosphere to Oceans" ;
put C3;
put / "Coefficient of Heat Gain by Deep Oceans" ;
put C4;

put / "** Global Environmental Variables" ;
put / "* CARBON CYCLE" ;
put / "Atmospheric concentration of carbon (GTC)" ;
Loop (T, put MAT.l(T));
put / "Atmospheric concentration of carbon (ppm)" ;
Loop (T, put (MAT.l(T)/2.13));
put / "Concentration in biosphere and upper oceans (GTC)" ;
Loop (T, put MU.l(T));
put / "Concentration in deep oceans (GTC)" ;
Loop (T, put ML.l(T));

put / "* CUMULATIVE EMISSIONS" ;
put / " Cumulative Emissions to date" ;
Loop (T, put CCA.l(T));
put / "Ratio To Max" ;
Loop (T, put (CCA.l(T)/FOSSLIM));

put / "* CLIMATE MODULE" ;
put / "Atmospheric Temperature (deg C above preindustrial)" ;
Loop (T, put TATM.l(T));
put / "Total Increase in Forcing (Watts per Meter2, preindustrial)" ;
Loop (T, put FORC.l(T));
put / "Lower Ocean Temperature (deg C above preindustrial)" ;
Loop (T, put TOCEAN.l(T));

put / "** Economic Endogenous Variables" ;
put / "* OUTPUT" ;
put / "Gross Output (trillion USD)" ;
Loop (T, put YGROSS.l(T));
put / "Climate Damages (fraction of gross output)" ;
Loop (T, put DAMFRAC.l(T));
put / "Climate damages (trillion USD)" ;
Loop (T, put DAMAGES.l(T));
put / "Output post-damages yet pre-abatement" ;
Loop (T, put YNET.l(T));
put / "Abatement cost (fraction of gross output)" ;
Loop (T, put "0");
put / "Abatement cost (trillion USD)" ;
Loop (T, put ABATECOST.l(T));
put / "Output (Net of Damages and Abatement, trillion USD pa) " ;
Loop (T, put Y.l(T));

put / "* EMISSIONS" ;
put / "Total Carbon Emissions (GTCO2 per year)" ;
Loop (T, put E.l(T));
put / "Industrial Emissions (GTCO2 per year)" ;
Loop (T, put EIND.l(T));
put / "World Emissions Intensity (sigma)" ;
Loop (T, put (EIND.l(T) / YGROSS.l(T)));

put / "* CAPITAL ACCUMULATION" ;
put / "Gross Investment (trillion 2005USD per year)" ;
Loop (T, put I.l(T));
put / "Capital (trillion 2005USD per year)" ;
Loop (T, put K.l(T));

put / "* DISCOUNT RATES" ;
put / "Utility Discount Rate (per year)" ;
Loop (T, put PRSTP);
put / "First Period gross MPK" ;
Loop (T, put ( GAMA * YGROSS.l(T) / K.l(T) ));

put / "* WELFARE" ;
put / "Consumption (trillion USD per year)" ;
Loop (T, put C.l(T));
put / "Consumption Per Capita (thousand USD per year)" ;
Loop (T, put CPC.l(T));
put / "Utility of Consumption" ;
Loop (T, put PERIODU.l(T));
put / "Total Discounted Utility (trillions 2000USD)" ;
put UTILITY.l;

put / "** CONTROL VARIABLES" ;
put / "Savings Rate (proportion of gross output)" ;
Loop (T, put S.l(T));
put / "Carbon Price (per t CO2)" ;
Loop (T, put cprice.l(T));
put / "Carbon Price (per t CO2) [repeat]" ;
Loop (T, put cprice.l(T));
put / "Carbon Price (per t C)" ;
Loop (T, put (cprice.l(T)*3.666) );
put / "Emissions Control Rate (total)" ;
Loop (T, put MIU.l(T));
put / "Emissions Control Rate (participants)" ;
Loop (T, put ( (cprice.l(T)/pbacktime(T))**(1/(expcost2-1)) )  );
put / "Carbon Price (Global Average)" ;
Loop (T, put (cprice.l(T)*partfract(T)) );
put / "Interest Rate (Real Rate of Return)" ;
Loop (T, put RI.l(T));
put / "Social Cost of Carbon" ;
Loop (T, put scc(T));

putclose;
