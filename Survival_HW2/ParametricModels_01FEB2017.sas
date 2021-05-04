/* SAS Code for Parametric Survival Models lecture, February 1, 2017 */
/* Sandra Gardner, Ph.D. */

options nodate nonumber;

/* modify the following libname statement as necessary */
libname sda 'C:\Users\Faizan\OneDrive\OneDrive\Survival Analysis\Survival_HW2';

ods rtf file='C:\Users\Faizan\OneDrive\OneDrive\Survival Analysis\Survival_HW2';


/* choose some hypothetical values for gamma, scale parameters to illustrate different shapes for hazard functions */

/* hypothetical exponential */
data hexp;
  do t=0 to 5 by .01;
    lambdah=0.6;
    chexph = lambdah*t; 
	sexph = exp(-chexph);
    lambdal=0.3;
    chexpl = lambdal*t;
	sexpl = exp(-chexpl);
    output;
  end;
run;

symbol1 c=red v=point i=join r=1 l=1;
symbol2 c=blue v=point i=join r=1 l=3;
axis1 order=(0 to 2 by .1) minor=none label=('h(t)');
axis2 order=(0 to 5 by 1) minor=none label=('years');
axis3 order=(0 to 2 by .1) minor=none label=('H(t)');
axis4 order=(0 to 1 by .2) minor=none label=('S(t)');

proc gplot data=hexp;
  plot (lambdah lambdal)*t/legend overlay vaxis=axis1 haxis=axis2;
  label lambdah='lambda=0.6'
        lambdal='lambda=0.3';
  title 'Exponential hazard plots';
run;
quit;

proc gplot data=hexp;
  plot (chexph chexpl)*t/legend overlay vaxis=axis3 haxis=axis2;
  label chexph='lambda=0.6'
        chexpl='lambda=0.3';
  title 'Exponential cumulative hazard plots';
run;
quit;

proc gplot data=hexp;
  plot (sexph sexpl)*t/legend overlay vaxis=axis4 haxis=axis2;
  label sexph='lambda=0.6'
        sexpl='lambda=0.3';
  title 'Exponential survival plots';
run;
quit;


/* hypothetical weibull */
data hweibull;
  lambda=.3;
  do t=0 to 5 by .01;
    gamma=1.5;
    hweibullh = gamma*lambda*(t**(gamma-1)); 
	chweibullh = lambda*(t**gamma); 
	sweibullh = exp(-chweibullh);
    gamma=1.0;
    if t>0 then hweibull1 = gamma*lambda*(t**(gamma-1)); 
	chweibull1 = lambda*(t**gamma); 
	sweibull1 = exp(-chweibull1);
    gamma=0.5;
    if t>0 then hweibulll = gamma*lambda*(t**(gamma-1)); 
	chweibulll = lambda*(t**gamma); 
	sweibulll = exp(-chweibulll);
    output;
  end;
run;

symbol1 c=red v=point i=join r=1 l=1;
symbol2 c=blue v=point i=join r=1 l=3;
symbol3 c=black v=point i=join r=1 l=5;
axis1 order=(0 to 2 by .1) minor=none label=('h(t)');
axis2 order=(0 to 5 by 1) minor=none label=('years');
axis3 order=(0 to 2 by .1) minor=none label=('H(t)');
axis4 order=(0 to 1 by .2) minor=none label=('S(t)');

proc gplot data=hweibull;
  plot (hweibullh hweibull1 hweibulll)*t/legend overlay vaxis=axis1 haxis=axis2;
  label hweibullh='gamma=1.5'
        hweibull1='gamma=1.0'
        hweibulll='gamma=0.5';
  title 'Weibull hazard plots - lambda=.3';
run;
quit;

proc gplot data=hweibull;
  plot (chweibullh chweibull1 chweibulll)*t/legend overlay vaxis=axis3 haxis=axis2;
  label chweibullh='gamma=1.5'
        chweibull1='gamma=1.0'
        chweibulll='gamma=0.5';
  title 'Weibull cumulative hazard plots - lambda=.3';
run;
quit;

proc gplot data=hweibull;
  plot (sweibullh sweibull1 sweibulll)*t/legend overlay vaxis=axis4 haxis=axis2;
  label sweibullh='gamma=1.5'
        sweibull1='gamma=1.0'
        sweibulll='gamma=0.5';
  title 'Weibull survival plots - lambda=.3';
run;
quit;

data hweibull;
  lambda=.6;
  do t=0 to 5 by .01;
    gamma=1.5;
    hweibullh = gamma*lambda*(t**(gamma-1)); 
	chweibullh = lambda*(t**gamma); 
	sweibullh = exp(-chweibullh);
    gamma=1.0;
    if t>0 then hweibull1 = gamma*lambda*(t**(gamma-1)); 
	chweibull1 = lambda*(t**gamma); 
	sweibull1 = exp(-chweibull1);
    gamma=0.5;
    if t>0 then hweibulll = gamma*lambda*(t**(gamma-1)); 
	chweibulll = lambda*(t**gamma); 
	sweibulll = exp(-chweibulll);
    output;
  end;
run;

proc gplot data=hweibull;
  plot (hweibullh hweibull1 hweibulll)*t/legend overlay vaxis=axis1 haxis=axis2;
  label hweibullh='gamma=1.5'
        hweibull1='gamma=1.0'
        hweibulll='gamma=0.5';
  title 'Weibull hazard plots - lambda=.6';
run;
quit;

proc gplot data=hweibull;
  plot (chweibullh chweibull1 chweibulll)*t/legend overlay vaxis=axis3 haxis=axis2;
  label chweibullh='gamma=1.5'
        chweibull1='gamma=1.0'
        chweibulll='gamma=0.5';
  title 'Weibull cumulative hazard plots - lambda=.6';
run;
quit;

proc gplot data=hweibull;
  plot (sweibullh sweibull1 sweibulll)*t/legend overlay vaxis=axis4 haxis=axis2;
  label sweibullh='gamma=1.5'
        sweibull1='gamma=1.0'
        sweibulll='gamma=0.5';
  title 'Weibull survival plots - lambda=.6';
run;
quit;

/* hypothetical log normal */
data hlnorm;
  u=0;
  pi=3.14159; 
  do t=0 to 5 by .01;
    scale=1.5;
    if t>0 then slnormh=1-probnorm((log(t)-u)/scale);
    if t>0 then flnormh=1/(sqrt(2*pi)*t*scale)*exp(-1/2*((log(t)-u)/scale)**2);
    hlnormh=flnormh/slnormh;
	chlnormh=-log(slnormh);
    scale=1.0;
    if t>0 then slnorm1=1-probnorm((log(t)-u)/scale);
    if t>0 then flnorm1=1/(sqrt(2*pi)*t*scale)*exp(-1/2*((log(t)-u)/scale)**2);
    hlnorm1=flnorm1/slnorm1;
	chlnorm1=-log(slnorm1);
    scale=0.5;
    if t>0 then slnorml=1-probnorm((log(t)-u)/scale);
    if t>0 then flnorml=1/(sqrt(2*pi)*t*scale)*exp(-1/2*((log(t)-u)/scale)**2);
    hlnorml=flnorml/slnorml;
 	chlnorml=-log(slnorml);
   output;
  end;
run;

symbol1 c=red v=point i=join r=1 l=1;
symbol2 c=blue v=point i=join r=1 l=3;
symbol3 c=black v=point i=join r=1 l=5;
axis1 order=(0 to 2 by .1) minor=none label=('h(t)');
axis2 order=(0 to 5 by 1) minor=none label=('years');
axis3 order=(0 to 2 by .1) minor=none label=('H(t)');
axis4 order=(0 to 1 by .2) minor=none label=('S(t)');

proc gplot data=hlnorm;
  plot (hlnormh hlnorm1 hlnorml)*t/legend overlay vaxis=axis1 haxis=axis2;
  label hlnormh='scale=1.5'
        hlnorm1='scale=1.0'
        hlnorml='scale=0.5';
  title 'Log normal hazard plots - u=0';
run;
quit;

proc gplot data=hlnorm;
  plot (chlnormh chlnorm1 chlnorml)*t/legend overlay vaxis=axis3 haxis=axis2;
  label chlnormh='scale=1.5'
        chlnorm1='scale=1.0'
        chlnorml='scale=0.5';
  title 'Log normal cumulative hazard plots - u=0';
run;
quit;

proc gplot data=hlnorm;
  plot (slnormh slnorm1 slnorml)*t/legend overlay vaxis=axis4 haxis=axis2;
  label slnormh='scale=1.5'
        slnorm1='scale=1.0'
        slnorml='scale=0.5';
  title 'Log normal survival plots - u=0';
run;
quit;

data hlnorm;
  u=0.5;
  pi=3.14159; 
  do t=0 to 5 by .01;
    scale=1.5;
    if t>0 then slnormh=1-probnorm((log(t)-u)/scale);
    if t>0 then flnormh=1/(sqrt(2*pi)*t*scale)*exp(-1/2*((log(t)-u)/scale)**2);
    hlnormh=flnormh/slnormh;
	chlnormh=-log(slnormh);
    scale=1.0;
    if t>0 then slnorm1=1-probnorm((log(t)-u)/scale);
    if t>0 then flnorm1=1/(sqrt(2*pi)*t*scale)*exp(-1/2*((log(t)-u)/scale)**2);
    hlnorm1=flnorm1/slnorm1;
	chlnorm1=-log(slnorm1);
    scale=0.5;
    if t>0 then slnorml=1-probnorm((log(t)-u)/scale);
    if t>0 then flnorml=1/(sqrt(2*pi)*t*scale)*exp(-1/2*((log(t)-u)/scale)**2);
    hlnorml=flnorml/slnorml;
 	chlnorml=-log(slnorml);
   output;
  end;
run;


symbol1 c=red v=point i=join r=1 l=1;
symbol2 c=blue v=point i=join r=1 l=3;
symbol3 c=black v=point i=join r=1 l=5;
axis1 order=(0 to 2 by .1) minor=none label=('h(t)');
axis2 order=(0 to 5 by 1) minor=none label=('years');
axis3 order=(0 to 2 by .1) minor=none label=('H(t)');
axis4 order=(0 to 1 by .2) minor=none label=('S(t)');

proc gplot data=hlnorm;
  plot (hlnormh hlnorm1 hlnorml)*t/legend overlay vaxis=axis1 haxis=axis2;
  label hlnormh='scale=1.5'
        hlnorm1='scale=1.0'
        hlnorml='scale=0.5';
  title 'Log normal hazard plots - u=.5';
run;
quit;

proc gplot data=hlnorm;
  plot (chlnormh chlnorm1 chlnorml)*t/legend overlay vaxis=axis3 haxis=axis2;
  label chlnormh='scale=1.5'
        chlnorm1='scale=1.0'
        chlnorml='scale=0.5';
  title 'Log normal cumulative hazard plots - u=.5';
run;
quit;

proc gplot data=hlnorm;
  plot (slnormh slnorm1 slnorml)*t/legend overlay vaxis=axis4 haxis=axis2;
  label slnormh='scale=1.5'
        slnorm1='scale=1.0'
        slnorml='scale=0.5';
  title 'Log normal survival plots - u=.5';
run;
quit;

/* hypothetical log logistic */
data hllog;
  alpha=1;
  do t=0 to 5 by .01;
    gamma=1/0.5;
    sllogh=1/(1+alpha*t**gamma);
    if t>0 then fllogh=alpha*gamma*t**(gamma-1)/(1+alpha*t**gamma)**2;
    hllogh=fllogh/sllogh;
 	chllogh=-log(sllogh);
    gamma=1;
    sllog1=1/(1+alpha*t**gamma);
    if t>0 then fllog1=alpha*gamma*t**(gamma-1)/(1+alpha*t**gamma)**2;
    hllog1=fllog1/sllog1;
 	chllog1=-log(sllog1);
    gamma=1/1.5;
    sllogl=1/(1+alpha*t**gamma);
    if t>0 then fllogl=alpha*gamma*t**(gamma-1)/(1+alpha*t**gamma)**2;
    hllogl=fllogl/sllogl;
 	chllogl=-log(sllogl);
    output;
  end;
run;

symbol1 c=red v=point i=join r=1 l=1;
symbol2 c=blue v=point i=join r=1 l=3;
symbol3 c=black v=point i=join r=1 l=5;
axis1 order=(0 to 2 by .1) minor=none label=('h(t)');
axis2 order=(0 to 5 by 1) minor=none label=('years');
axis3 order=(0 to 2 by .1) minor=none label=('H(t)');
axis4 order=(0 to 1 by .2) minor=none label=('S(t)');

proc gplot data=hllog;
  plot (hllogh hllog1 hllogl)*t/legend overlay vaxis=axis1 haxis=axis2;
  label hllogh='gamma=1.5'
        hllog1='gamma=1.0'
        hllogl='gamma=0.5';
  title 'Log logistic hazard plots - alpha=1';
run;
quit;

proc gplot data=hllog;
  plot (chllogh chllog1 chllogl)*t/legend overlay vaxis=axis3 haxis=axis2;
  label chllogh='gamma=1.5'
        chllog1='gamma=1.0'
        chllogl='gamma=0.5';
  title 'Log logistic cumulative hazard plots - alpha=1';
run;
quit;

proc gplot data=hllog;
  plot (sllogh sllog1 sllogl)*t/legend overlay vaxis=axis4 haxis=axis2;
  label sllogh='gamma=1.5'
        sllog1='gamma=1.0'
        sllogl='gamma=0.5';
  title 'Log logistic survival plots - alpha=1';
run;
quit;

data hllog;
  alpha=2;
  do t=0 to 5 by .01;
    gamma=1/0.5;
    sllogh=1/(1+alpha*t**gamma);
    if t>0 then fllogh=alpha*gamma*t**(gamma-1)/(1+alpha*t**gamma)**2;
    hllogh=fllogh/sllogh;
 	chllogh=-log(sllogh);
    gamma=1;
    sllog1=1/(1+alpha*t**gamma);
    if t>0 then fllog1=alpha*gamma*t**(gamma-1)/(1+alpha*t**gamma)**2;
    hllog1=fllog1/sllog1;
 	chllog1=-log(sllog1);
    gamma=1/1.5;
    sllogl=1/(1+alpha*t**gamma);
    if t>0 then fllogl=alpha*gamma*t**(gamma-1)/(1+alpha*t**gamma)**2;
    hllogl=fllogl/sllogl;
 	chllogl=-log(sllogl);
    output;
  end;
run;

proc gplot data=hllog;
  plot (hllogh hllog1 hllogl)*t/legend overlay vaxis=axis1 haxis=axis2;
  label hllogh='gamma=1.5'
        hllog1='gamma=1.0'
        hllogl='gamma=0.5';
  title 'Log logistic hazard plots - alpha=2';
run;
quit;

proc gplot data=hllog;
  plot (chllogh chllog1 chllogl)*t/legend overlay vaxis=axis3 haxis=axis2;
  label chllogh='gamma=1.5'
        chllog1='gamma=1.0'
        chllogl='gamma=0.5';
  title 'Log logistic cumulative hazard plots - alpha=2';
run;
quit;

proc gplot data=hllog;
  plot (sllogh sllog1 sllogl)*t/legend overlay vaxis=axis4 haxis=axis2;
  label sllogh='gamma=1.5'
        sllog1='gamma=1.0'
        sllogl='gamma=0.5';
  title 'Log logistic survival plots - alpha=2';
run;
quit;

/* example data set: brain */

/* temporary formats */
proc format;
  value yn
        0='No' 
        1='Yes';
  value agef
        0='<50' 
        1='>=50';
run;

proc univariate data=sda.brain normal;
  var age;
  histogram age/normal;
  title 'Distribution age';
run;

proc freq data=sda.brain;
  tables treat age50;
  title 'Covariate summary';
run;

symbol1 c=red l=1;
* note adding ls and lls plots and output data set with K-M estimates;
* output data set: osurv;
/* overall median      27.430       23.140      31.430	

Mean    Standard Error
44.528             3.285
*/
proc lifetest data=sda.brain plot=(s, ls, lls) outsurv=osurv notable;
  time weeks*event(0);
  title 'LifeTest: Overall Survival';
run;


proc lifetest data=sda.brain plot=(s(atrisk outside)) outsurv=osurv notable;
  time weeks*event(0);
  title 'LifeTest: Overall Survival (Newer options)';
run;

/* suggested plots from Allison: Survival Analysis Using the SAS System
   if logit vs log(t) or probit(Prob(failed at t)) vs log(t) is linear, loglogistic or lognormal model may be good fit */
data osurv2;
  set osurv; /* output data set from Proc Lifetest from overall survival */
  where survival^=1 and survival^=0 and weeks^=0;
  logit=log((1-survival)/survival);	/* - odds of survival, i.e.  -log(survival/(1-survival)) */
  lnorm=probit(1-survival);
  lweeks=log(weeks);
run;
 
symbol1 c=red value=point i=join l=1;
axis1 minor=none label=('probit');
axis2 minor=none label=('log(t)');
proc gplot data=osurv2;
  plot lnorm*lweeks/vaxis=axis1 haxis=axis2;
  title 'Probit(CDF) Plot';
run;
quit; 

/* log weeks (events only) */
/*
              Fitted Normal Distribution for lweeks

                Parameters for Normal Distribution

                  Parameter   Symbol   Estimate

                  Mean        Mu       3.218779
                  Std Dev     Sigma     0.80768


         Goodness-of-Fit Tests for Normal Distribution

 Test                  ----Statistic-----   ------p Value------

 Kolmogorov-Smirnov    D       0.05447400   Pr > D        0.137
 Cramer-von Mises      W-Sq    0.06230297   Pr > W-Sq    >0.250
 Anderson-Darling      A-Sq    0.42332028   Pr > A-Sq    >0.250

*/
proc univariate data=sda.brain normal;
  where event=1;
  var lweeks;
  histogram lweeks/normal;
  cdfplot lweeks/normal;
  title 'Distribution log(weeks) - events only';
run;

proc univariate data=sda.brain normal;
  where event=1;
  var weeks;
  histogram weeks/lognormal;
  histogram weeks/gamma(theta=est);
  cdfplot weeks/lognormal;
  cdfplot weeks/gamma(theta=est);
  title 'Distribution weeks - events only';
run;

data wdist;
  pi=3.14159; 
  do w=-10 to 10 by .1;
     e=exp(w-exp(w));
	 n=exp(-w**2/2)/(sqrt(2*pi));
	 l=exp(w)/(1+exp(w))**2;
     output;
  end;
run;

symbol1 c=red v=point i=join r=1 l=1;
symbol2 c=blue v=point i=join r=1 l=3;
symbol3 c=black v=point i=join r=1 l=5;
axis1 order=(0 to .5 by .1) minor=none label=('f(w)');
axis2 order=(-10 to 10 by 2) minor=none label=('w');
proc gplot data=wdist;
  plot (e n l)*w/overlay legend vaxis=axis1 haxis=axis2;
  title 'Error distributions';
  label e='Extreme Value' n='Normal' l='Logistic';
run;
quit;

symbol1 c=red value=point i=join l=1;
axis1 minor=none label=('logit');
axis2 minor=none label=('log(t)');
proc gplot data=osurv2;
  plot logit*lweeks/vaxis=axis1 haxis=axis2;
  title 'Logit(CDF) Plot';
run;
quit;

/* ODS plots with hazard plot */
proc lifetest data=sda.brain plot=(survival hazard logsurv loglogs) notable;
  time weeks*event(0);
  title 'LifeTest: Overall Survival: hazard';
run;

/* ODS plots with hazard plot - reduced range */
proc lifetest data=sda.brain plot=(hazard(gridu=125)) notable;
  time weeks*event(0);
  title 'LifeTest: Overall Survival: hazard';
run;


/* introducing the lifereg procedure */

/* 
Variable             Sum
event                207
weeks               9426

  note closed form solution for event rate can be calculated from this summary data
  i.e. 207/9426=0.02196 and ln(0.02196)=-3.8185 
*/
proc means data=sda.brain sum maxdec=0;
  var event weeks;
  title 'Data Summary';
run;

proc lifereg data=sda.brain;
  model weeks*event(0)=/d=exponential;
  title 'LifeReg: Overall Survival - Exponential';
run;

/* probplots */
ods output FitStatistics=expfit;
proc lifereg data=sda.brain;
  model weeks*event(0)=/d=exponential;
  probplot;
  inset;
  title 'LifeReg: Overall Survival - Probability Plot (Exponential)';
run;

ods output FitStatistics=weibfit;
proc lifereg data=sda.brain;
  model weeks*event(0)=/d=weibull;
  probplot;
  inset;
  title 'LifeReg: Overall Survival - Probability Plot (Weibull)';
run;

ods output FitStatistics=lnormfit;
proc lifereg data=sda.brain;
  model weeks*event(0)=/d=lnormal;
  probplot;
  inset;
  title 'LifeReg: Overall Survival - Probability Plot (Log normal)';
run;

ods output FitStatistics=llogfit;
proc lifereg data=sda.brain;
  model weeks*event(0)=/d=llogistic;
  probplot;
  inset;
  title 'LifeReg: Overall Survival - Probability Plot (Log logistic)';
run;

ods output FitStatistics=gammafit;
proc lifereg data=sda.brain;
  model weeks*event(0)=/d=gamma;
  probplot;
  inset;
  title 'LifeReg: Overall Survival - Probability Plot (Gamma)';
run;

/* compare models */
proc transpose data=expfit out=texpfit;
  id criterion;
  var value;
run;

proc transpose data=weibfit out=tweibfit;
  id criterion;
  var value;
run;

proc transpose data=lnormfit out=tlnormfit;
  id criterion;
  var value;
run;

proc transpose data=llogfit out=tllogfit;
  id criterion;
  var value;
run;

proc transpose data=gammafit out=tgammafit;
  id criterion;
  var value;
run;

data allfit;
  label model='Model';
  set texpfit(drop=_name_ in=exp)
      tweibfit(drop=_name_ in=weib)
      tlnormfit(drop=_name_ in=lnorm)
      tllogfit(drop=_name_ in=llog)
      tgammafit(drop=_name_ in=gamma);
  if exp then model='Exponential';
     else if weib then model='Weibull';
     else if lnorm then model='LogNormal';
     else if llog then model='LogLogistic';
     else if gamma then model='Gamma';
run;

proc print data=allfit noobs;
  title 'Model Comparison';
run;


/* simple models with no covariates, i.e. intercept model */

/*

                 Fit Statistics

-2 Log Likelihood                        662.275
AIC (smaller is better)                  664.275

                          Standard   95% Confidence     Chi-
Parameter     DF Estimate    Error       Limits       Square Pr > ChiSq

Intercept      1   3.8185   0.0695   3.6823   3.9547 3018.22     <.0001
Scale          0   1.0000   0.0000   1.0000   1.0000
Weibull Scale  1  45.5348   3.1649  39.7357  52.1803
Weibull Shape  0   1.0000   0.0000   1.0000   1.0000


                    Lagrange Multiplier Statistics

                Parameter     Chi-Square    Pr > ChiSq

                Scale             0.6216        0.4304
*/
proc lifereg data=sda.brain;
  model weeks*event(0)=/d=exponential;
  title 'LifeReg: Overall Survival - Exponential';
run;

/* can use likelihood ratio test to compare Exponential to Weibull*/
/* does confidence interval for shape contain 1? */
/* 

Weibull

-2 Log Likelihood                        661.693
AIC (smaller is better)                  665.693

          Analysis of Maximum Likelihood Parameter Estimates

                          Standard   95% Confidence     Chi-
Parameter     DF Estimate    Error       Limits       Square Pr > ChiSq

Intercept      1   3.8321   0.0692   3.6965   3.9676 3069.48     <.0001
Scale          1   0.9608   0.0498   0.8679   1.0636   * 1/shape
Weibull Scale  1  46.1571   3.1925  40.3054  52.8583
Weibull Shape  1   1.0408   0.0540   0.9402   1.1522
*/
proc lifereg data=sda.brain;
  model weeks*event(0)=/d=weibull;
  title 'LifeReg: Overall Survival - Weibull';
run;

/*

-2 Log Likelihood                        608.002
AIC (smaller is better)                  612.002

                          Standard   95% Confidence     Chi-
Parameter     DF Estimate    Error       Limits       Square Pr > ChiSq

Intercept      1   3.3653   0.0644   3.2389   3.4916 2726.77     <.0001
Scale          1   0.9553   0.0479   0.8659   1.0540

*/
proc lifereg data=sda.brain;
  model weeks*event(0)=/d=lnormal;
  title 'LifeReg: Overall Survival - LogNormal';
run;

/*

-2 Log Likelihood                        604.338
AIC (smaller is better)                  608.338

                          Standard   95% Confidence     Chi-
Parameter     DF Estimate    Error       Limits       Square Pr > ChiSq

Intercept      1   3.3233   0.0625   3.2008   3.4458 2828.98     <.0001
Scale          1   0.5398   0.0315   0.4815   0.6052
*/
proc lifereg data=sda.brain;
  model weeks*event(0)=/d=llogistic;
  title 'LifeReg: Overall Survival - LogLogistic';
run;

/*
-2 Log Likelihood                        600.334
AIC (smaller is better)                  606.334

                          Standard   95% Confidence     Chi-
Parameter     DF Estimate    Error       Limits       Square Pr > ChiSq

Intercept      1   3.1407   0.1006   2.9434   3.3380  973.70     <.0001
Scale          1   0.9272   0.0479   0.8379   1.0259
Shape          1  -0.4929   0.1733  -0.8326  -0.1533
*/
proc lifereg data=sda.brain;
  model weeks*event(0)=/d=gamma;
  title 'LifeReg: Overall Survival - Gamma';
run;

/* Likelihood ratio test - Gamma & LogNormal */
data probchi;
 llm1=608.002;
 llm2=600.334;
 lrt=llm1-llm2;
 p=1-probchi(lrt,1);
proc print noobs;
 title 'Likelihood rato test';
run;

/* superimpose estimated parametric model over K-M curve */
/* need to understand how model is parameterized (no covariates, 'average' curve) */
data osurv2;
  set osurv;
  where _censor_ in (0,.);  /* SAS created variable = 0 when an event weeks */

  /* calculate estimate S(ti) for each model */

  /* exponential - note change in sign for beta estimate */
  sexp=exp(-exp(-3.8185)*weeks);

  /* weibull  - note change in sign for beta estimate */
  sweibull=exp(-exp(-1.0408*3.8321)*weeks**1.0408);

  /* alternate weibull using extreme value distribution */
  sev=exp( -exp( (lt-3.8321) /0.9608) );

  /* lognormal - do no change size for beta estimate */
  u=3.3653;
  scale=0.9553;
  pi=3.14159; 
  if weeks>0 then slnorm=1-probnorm((log(weeks)-u)/scale);

  /* loglogistic - note change in sign for beta estimate */
  alpha=exp(-3.3233/0.5398);
  gamma=1/0.5398; /* 1/scale */
  sllog=1/(1+alpha*weeks**gamma);
run;

symbol1 c=red v=point i=join r=1 l=1;
symbol2 c=blue v=point i=join r=1 l=3;
symbol3 c=black v=point i=join r=1 l=5;
symbol4 c=green v=point i=join r=1 l=8;
symbol5 c=black v=point i=stepjl l=2;
axis1 order=(0 to 1 by .1) minor=none label=('S(t)');
axis2 order=(0 to 200 by 50) minor=none label=('weeks');

proc gplot data=osurv2;
  plot (sexp sweibull slnorm sllog survival)*weeks/overlay vaxis=axis1 haxis=axis2 legend;
  label survival='KM' sexp='Exponential' sweibull='Weibull'
        slnorm='LogNormal' sllog='LogLogistic'
        weeks='Weeks';
  title1 'Comparison of Survival Models';
run;
quit;

/* alternate way to obtain fitted model (unadjusted) */
data brain2;
  if _N_=1 
      then do;
         weeks=.;  /* set up dummy observation with no time value, if adjusted model, set values for covariates */
         control=1;
         output;
      end;
   set sda.brain;
   control=0;
   output;
run;

proc lifereg data=brain2;
  model weeks*event(0)=/d=gamma;
  output out=brain3 quantiles=.01 to .99 by .01 std=std p=predg
         control=control;  /* obtain cdf (quantiles) for dummy observation only */
  title 'LifeReg: Overall Survival - Gamma';
run;

data brain4;
  set osurv /* KM */ brain3(in=b);
  if b /* convert quantiles and predicted values */
     then do;
             st=1-_prob_; /* gamma survival curve */
             t=predg;
             lt=log(t);
             d=-0.4929; /* shape */
             s=0.9272; /* scale */
             mu=3.1407; /* xbeta - in this case, intercept only */
             v=exp((log(t)-mu)/s);
             d2=d**2;
             vd=v**d;
             gv=abs(d)/(v*gamma(1/d2)) * (1/d2*vd)**(1/d2) * exp(-vd/d2);
             gt=v*gv/(t*s); /* gamma distribution */
             ht=gt/st; /* gamma hazard function */
		  end;
run;

symbol1 v=point c=black i=join r=1;
legend1 label=('Distribution');
proc gplot data=brain4;
  plot gt*t/ legend=legend1;
  title 'Fitted distribution - generalized gamma';
run;
quit;

symbol1 v=point c=black i=join r=1;
legend1 label=('Hazard');
axis2 order=(0 to 200 by 50) minor=none label=('weeks');
proc gplot data=brain4;
  plot ht*t/ legend=legend1 haxis=axis2;
  title 'Fitted hazard - generalized gamma';
run;
quit;

symbol1 v=point c=black i=join r=1;
legend1 label=('Survival');
proc gplot data=brain4;
  plot st*t/ legend=legend1;
  title 'Fitted survival - generalized gamma';
run;
quit;

symbol1 c=red v=point i=join r=1 l=1;
symbol2 c=black v=point i=stepjl l=2;
axis1 order=(0 to 1 by .1) minor=none label=('S(t)');
axis2 order=(0 to 200 by 50) minor=none label=('weeks');
proc gplot data=brain4;
  plot st*t survival*weeks/overlay vaxis=axis1 haxis=axis2 legend;
  label survival='KM' st='Gamma'
        weeks='Weeks';
  title1 'Comparison of Gamma and KM curve';
run;
quit;


/* add in covariates */
/* median treat=No :      50      23.570     20.570      28.000
   median treat=Yes :     50      31.500     26.290      37.000
 */

symbol1 c=red l=1;
symbol2 c=blue l=3;
proc lifetest data=sda.brain plot=(s, ls, lls) /*notable*/  outsurv=tsurv;
  time weeks*event(0);
  strata treat;
  format treat yn.;
  title 'LifeTest: Treatment group';
run;

/* median age<50 :      50      32.430    27.140      39.710
   median age>=50 :     50      21.500    19.000      27.290

 */

symbol1 c=red l=1;
symbol2 c=blue l=3;
proc lifetest data=sda.brain plot=(s, ls, lls) /*notable*/  outsurv=asurv;
  time weeks*event(0);
  strata age50;
  format age50 agef.;
  title 'LifeTest: Age group';
run;

/* Phreg 
                    Parameter      Standard                                  Hazard
Parameter    DF      Estimate         Error    Chi-Square    Pr > ChiSq       Ratio

treat         1      -0.20006       0.14011        2.0389        0.1533       0.819
age50         1       0.51040       0.14037       13.2213        0.0003       1.666

*/

proc lifetest data=sda.brain noprint  outsurv=atsurv;
  time weeks*event(0);
  strata age50 treat;
  format age50 agef. treat yn.;
  title 'LifeTest: Age & Treatment groups';
run;

proc phreg data=sda.brain;
  model weeks*event(0)=treat age50;
  title 'PHReg: Treatment & Age groups';
run;

/* note similar HR, but sign is reversed from PHReg

                           Standard   95% Confidence     Chi-
 Parameter     DF Estimate    Error       Limits       Square Pr > ChiSq

 Intercept      1   3.9439   0.1206   3.7075   4.1804 1069.08     <.0001
 treat          1   0.1825   0.1390  -0.0900   0.4549    1.72     0.1893
 age50          1  -0.5007   0.1390  -0.7732  -0.2283   12.98     0.0003
 Scale          0   1.0000   0.0000   1.0000   1.0000
 Weibull Shape  0   1.0000   0.0000   1.0000   1.0000

*/
proc lifereg data=sda.brain;
  model weeks*event(0)=treat age50/d=exponential;
  title 'LifeReg: Treatment & Age groups - Exponential';
run;

/*  treatment only
                          Standard   95% Confidence     Chi-
Parameter     DF Estimate    Error       Limits       Square Pr > ChiSq

Intercept      1   3.7220   0.0981   3.5298   3.9142 1440.73     <.0001
treat          1   0.1853   0.1390  -0.0871   0.4578    1.78     0.1825
Scale          0   1.0000   0.0000   1.0000   1.0000
Weibull Shape  0   1.0000   0.0000   1.0000   1.0000
*/

proc lifereg data=sda.brain;
  model weeks*event(0)=treat/d=exponential;
  title 'LifeReg: Treatment groups - Exponential';
run;
/* note closed form solution for hazard ratio can be calulated from this summary data (unadjusted for age)
  i.e. No: 104/4300=0.0242 (note ln(0.0242)=-3.722) and  Yes: 103/5126=0.0201
       logHR: ln(0.0201/0.0242)=-0.1856 */
proc means data=sda.brain sum maxdec=0;
  class treat;
  var event weeks;
  format treat yn.;
  title 'Data Summary';
run;

/*
                          Standard   95% Confidence     Chi-
Parameter     DF Estimate    Error       Limits       Square Pr > ChiSq

Intercept      1   4.0394   0.0985   3.8463   4.2325 1680.64     <.0001
age50          1  -0.5018   0.1390  -0.7743  -0.2293   13.03     0.0003
*/
proc lifereg data=sda.brain;
  model weeks*event(0)=age50/d=exponential;
  * if covariate data set is not provided, then predicted survival,
    if quantile not provided then predicted median survival
    output residuals; 
  output out=eout cres=cres p=predm  std=stdm;
  title 'LifeReg: Age group - Exponential';
run;

/*
                          Standard   95% Confidence     Chi-
Parameter     DF Estimate    Error       Limits       Square Pr > ChiSq

Intercept      1   4.0569   0.0933   3.8740   4.2397 1891.78     <.0001
age50          1  -0.4927   0.1303  -0.7480  -0.2374   14.31     0.0002
Scale          1   0.9356   0.0481   0.8459   1.0349
Weibull Shape  1   1.0688   0.0550   0.9663   1.1822
*/

proc lifereg data=sda.brain;
  model weeks*event(0)=age50/d=weibull;
  output out=wout cres=cres p=predm  std=stdm;
  title 'LifeReg: Age group - Weibull';
run;
 
/* importing K-M survival estimates and estimating model values at each event weeks
    followed by plots of models for each treatment group */
/* be careful of parameterization of the model; change sign */ 
data asurv2;
  set asurv;
  where _censor_ in (0,.);
  /* models from lifereg  */
  if age50=0
     then sexp=exp(-exp(-4.0394)*weeks); 
     else sexp=exp(-exp(-(4.0394-0.5018))*weeks); 
  if age50=0
     then sweibull=exp(-exp(-1.0688*4.0569)*weeks**1.0688);
     else sweibull=exp(-exp(-1.0688*(4.0569-0.4927))*weeks**1.0688);
run;

symbol1 c=black v=point i=stepjl l=1;
symbol2 c=red v=point i=join l=1;
symbol3 c=blue v=point i=join l=3;
proc gplot data=asurv2;
  where age50=0;
  plot (survival sexp sweibull)*weeks/overlay vaxis=axis1 haxis=axis2 legend;
  label survival='KM' sexp='Exponential' sweibull='Weibull' weeks='Weeks';
  title1 'Comparison of Exponential and Weibull Models-Age<50';
run;

symbol1 c=black v=point i=stepjl l=1;
symbol2 c=red v=point i=join l=1;
symbol3 c=blue v=point i=join l=3;
proc gplot data=asurv2;
  where age50=1;
  plot (survival sexp sweibull)*weeks/overlay vaxis=axis1 haxis=axis2 legend;
  label survival='KM' sexp='Exponential' sweibull='Weibull' weeks='Weeks';
  title1 'Comparison of Exponential and Weibull Models-Age>=50';
run;
quit;

/* cox-snell residuals:  age50 model*/
proc lifetest data=eout plots=(ls) notable;
  * looking for evidence that cres is exponential using the -log(S(t)) plot;
  * note that censoring value is maintained from original data set;
  time cres*event(0);
  title1 'Cox-Snell Residuals-Exponential';
run;

/* simulated exponential data */
data simexp;
  do k=1 to 500;
  * ranuni: uniform random number generator (fixed seed);
  *   obtain random number between (0,1) - i.e. the P(S>t) and solve for t;
  x1=ranuni(8283944); 
  *   second random number for independent censoring process;
  x2=ranuni(2830483);
  * choose rate of event and censoring rate;
  lambda=0.0034;
  lambdac=0.0003;
  t=-1/lambda*log(x1);
  c=-1/lambdac*log(x2);
  * t2 in Weeks, maximum number of Weeks followed up is 800,
     otherwise minimum of c and t;
  t2=min(t,c,800); 
  if (t2=c and c^=t) or t2=800
     then event=0;
     else event=1;
  output;
  end;
  format t c t2 7.1;
run;

proc freq data=simexp;
  tables event;
  title 'Simulated Exponential Data';
run;

/* simulated exponential data - are LS and LLS plots straight lines? */
proc lifetest data=simexp plot=(s,ls,lls) notable;
  time t2*event(0);
  title 'Simulated Exponential Data';
run; 

proc lifetest data=simexp plot=(hazard) notable;
  time t2*event(0);
  title 'Simulated Exponential Data - hazard';
run;

/* obtain Cox-Snell residuals from exponential model */
proc lifereg data=simexp;
  model t2*event(0)=/d=exponential;
  output out=simout cres=cres;
run; 

/* do Cox-Snell residuals have exp(1) distribution?  Is the LS plot a straight line with slope 1? */
proc lifetest data=simout plot=(ls) notable;
  time cres*event(0);
run; 


/* example models with treatment and age50 covariates */

/* create output data set with coxsnell (cres) and standardized (sres) residuals */
/* lifereg has less residuals than phreg */

/*

-2 Log Likelihood                        647.648

                          Standard   95% Confidence     Chi-
Parameter     DF Estimate    Error       Limits       Square Pr > ChiSq

Intercept      1   3.9439   0.1206   3.7075   4.1804 1069.08     <.0001
treat          1   0.1825   0.1390  -0.0900   0.4549    1.72     0.1893
age50          1  -0.5007   0.1390  -0.7732  -0.2283   12.98     0.0003
Scale          0   1.0000   0.0000   1.0000   1.0000
Weibull Shape  0   1.0000   0.0000   1.0000   1.0000
*/
proc lifereg data=sda.brain;
  model weeks*event(0)=treat age50/d=exponential;
  output out=eout cres=cres sres=sres p=predm  std=stdm;
  title 'LifeReg: Treatment & Age groups - Exponential';
run;

/*

-2 Log Likelihood                        645.784

                          Standard   95% Confidence     Chi-
Parameter     DF Estimate    Error       Limits       Square Pr > ChiSq

Intercept      1   3.9622   0.1132   3.7403   4.1842 1224.33     <.0001
treat          1   0.1825   0.1294  -0.0711   0.4361    1.99     0.1585
age50          1  -0.4904   0.1296  -0.7444  -0.2363   14.31     0.0002
Scale          1   0.9308   0.0479   0.8415   1.0296
Weibull Shape  1   1.0744   0.0553   0.9713   1.1884
*/
proc lifereg data=sda.brain;
  model weeks*event(0)=treat age50/d=weibull;
  output out=wout cres=cres sres=sres p=predm  std=stdm;
  title 'LifeReg: Treatment & Age groups - Weibull';
run;


/*

-2 Log Likelihood                        595.383

                           Standard   95% Confidence     Chi-
 Parameter     DF Estimate    Error       Limits       Square Pr > ChiSq

 Intercept      1   3.4768   0.1072   3.2667   3.6869 1051.87     <.0001
 treat          1   0.1744   0.1253  -0.0711   0.4200    1.94     0.1639
 age50          1  -0.4144   0.1254  -0.6602  -0.1686   10.92     0.0010
 Scale          1   0.9288   0.0466   0.8418   1.0247
*/
proc lifereg data=sda.brain;
  model weeks*event(0)=treat age50/d=lnormal;
  output out=lnout cres=cres sres=sres p=predm  std=stdm;
  title 'LifeReg: Treatment & Age groups - Lognormal';
run;

/*

 -2 Log Likelihood                        589.891

          Analysis of Maximum Likelihood Parameter Estimates

                          Standard   95% Confidence     Chi-
Parameter     DF Estimate    Error       Limits       Square Pr > ChiSq

Intercept      1   3.4289   0.1031   3.2268   3.6309 1105.96     <.0001
treat          1   0.2029   0.1200  -0.0323   0.4380    2.86     0.0909
age50          1  -0.4204   0.1200  -0.6555  -0.1852   12.28     0.0005
Scale          1   0.5198   0.0304   0.4635   0.5830
*/
proc lifereg data=sda.brain;
  model weeks*event(0)=treat age50/d=llogistic;
  output out=llout cres=cres sres=sres p=predm  std=stdm;
  title 'LifeReg: Treatment & Age groups - Loglogistic';
run;

/* age continuous
                          Standard   95% Confidence     Chi-
Parameter     DF Estimate    Error       Limits       Square Pr > ChiSq

Intercept      1   4.0332   0.2350   3.5726   4.4938  294.55     <.0001
treat          1   0.2220   0.1200  -0.0132   0.4573    3.42     0.0643
age            1  -0.0170   0.0046  -0.0259  -0.0080   13.71     0.0002
Scale          1   0.5187   0.0303   0.4626   0.5817

*/
proc lifereg data=sda.brain;
  model weeks*event(0)=treat age/d=llogistic;
  output out=llout2 cres=cres sres=sres p=predm  std=stdm;
  title 'LifeReg: Treatment & Age continuous - Loglogistic';
run;



/* gamma model may help to explore nested models
gamma with shape=1 is Weibull
gamma with shape=0 is LogNormal

 -2 Log Likelihood                        590.937


          Analysis of Maximum Likelihood Parameter Estimates

                          Standard   95% Confidence     Chi-
Parameter     DF Estimate    Error       Limits       Square Pr > ChiSq

Intercept      1   3.2902   0.1386   3.0184   3.5619  563.22     <.0001
treat          1   0.1579   0.1236  -0.0844   0.4002    1.63     0.2015
age50          1  -0.3511   0.1279  -0.6018  -0.1005    7.54     0.0060
Scale          1   0.9168   0.0463   0.8304   1.0122
Shape          1  -0.3699   0.1733  -0.7096  -0.0302

*/
proc lifereg data=sda.brain;
  model weeks*event(0)=treat age50/d=gamma;
  output out=gout cres=cres sres=sres p=predm  std=stdm;
  title 'LifeReg: Treatment group - Gamma';
run;


/* cox-snell residuals - see Hosmer & Lemeshow Section 8.3 */
proc lifetest data=eout plots=(ls) notable;
  * looking for evidence that cres is exponential using the -log(S(t)) plot;
  * note that censoring value is maintained from original data set;
  time cres*event(0);
  title1 'Cox-Snell Residuals - Exponential';
run;

proc lifetest data=wout plots=(ls) notable;
  time cres*event(0);
  title1 'Cox-Snell Residuals - Weibull';
run;

proc lifetest data=lnout plots=(ls) notable;
  time cres*event(0);
  title1 'Cox-Snell Residuals - LogNormal';
run;

proc lifetest data=llout plots=(ls) notable;
  time cres*event(0);
  title1 'Cox-Snell Residuals - LogLogistic';
run;

/* Normal-deviate residuals from Nardi & Schemper */
/* Martingale & deviance residuals from Klein & Moeschberger */
data nardi;
  set sda.brain;
  id=_n_;

  /* calculate estimate S(ti) or S(ci) */

  /* exponential */
  lambda=exp(-(3.9439+0.1825*treat-0.5007*age50));
  sexp=exp(-lambda*weeks);
  xbexp=3.9439+0.1825*treat-0.5007*age50;
  chexp=-log(sexp);
  martexp=event-chexp;
  devexp=sign(martexp)*(-2*(martexp+event*log(event-martexp)))**1/2;

  /* weibull */
  lambda=exp(-1.0744*(3.9622+0.1825*treat+-0.4904*age50));
  gamma=1.0744;
  sweibull=exp(-lambda*weeks**gamma);
  xbweibull=3.9622+0.1825*treat+-0.4904*age50;
  chweibull=-log(sweibull);
  martweibull=event-chweibull;
  devweibull=sign(martweibull)*(-2*(martweibull+event*log(event-martweibull)))**1/2;

  /* lognormal */
  u=3.4768+0.1744*treat-0.4144*age50;
  scale=0.9288;
  pi=3.14159; 
  if weeks>0 then slnorm=1-probnorm((log(weeks)-u)/scale);
  xblnorm=u;
  chlnorm=-log(slnorm);
  martlnorm=event-chlnorm;
  devlnorm=sign(martlnorm)*(-2*(martlnorm+event*log(event-martlnorm)))**1/2;

  /* loglogistic */
  alpha=exp(-(3.4289+0.2029*treat-0.4204*age50)/0.5198);
  gamma=1/0.5198;
  sllog=1/(1+alpha*weeks**gamma);
  xbllog=3.4289+0.2029*treat-0.4204*age50;
  chllog=-log(sllog);
  martllog=event-chllog;
  devllog=sign(martllog)*(-2*(martllog+event*log(event-martllog)))**1/2;

  /* loglogistic - age continuous */
  alpha=exp(-(4.0332+0.2220*treat-0.0170*age)/0.5187);
  gamma=1/0.5187;
  sllog2=1/(1+alpha*weeks**gamma);
  xbllog2=4.0332+0.2220*treat-0.0170*age;
  chllog2=-log(sllog2);
  martllog2=event-chllog2;
  devllog2=sign(martllog2)*(-2*(martllog2+event*log(event-martllog2)))**1/2;

  /* Nardi & Schemper normal-deviate residual */
  array res{4} ndr1-ndr4;
  if event=1  /* event, as coded in UIS data set */
     then do;
             rexp=probit(sexp);
             rweibull=probit(sweibull);
             rlnorm=probit(slnorm);
             rllog=probit(sllog);
		  end;
     else do; /* censored */
	     do i=1 to 4;
		     rn=sexp*ranuni(3289372);  /* random number from U(0,S(ci)) */
	         res{i}=probit(rn);
		 end;
         rexp=mean(of ndr1-ndr4);
	     do i=1 to 4;
		     rn=sweibull*ranuni(28177261);  /* random number from U(0,S(ci)) */
	         res{i}=probit(rn);
		 end;
         rweibull=mean(of ndr1-ndr4);
	     do i=1 to 4;
		     rn=slnorm*ranuni(482783072);  /* random number from U(0,S(ci)) */
	         res{i}=probit(rn);
		 end;
         rlnorm=mean(of ndr1-ndr4);
	     do i=1 to 4;
		     rn=sllog*ranuni(593847302);  /* random number from U(0,S(ci)) */
	         res{i}=probit(rn);
		 end;
         rllog=mean(of ndr1-ndr4);
	 end;
run;

/* Nardi & Schemper suggest censoring should not exceed 40-50 percent 
proc freq data=nardi;
  tables event;
  format event yn.;
  title 'Censoring summary';
run;

proc univariate data=nardi noprint;
  class event;
  var rexp rweibull rlnorm rllog;
  histogram rexp rweibull rlnorm rllog/normal;
  title 'Normal-deviate residuals by event/censor';
run;

proc univariate data=nardi noprint;
  var rexp rweibull rlnorm rllog;
  histogram rexp rweibull rlnorm rllog/normal;
  title 'Normal-deviate residuals';
run;
*/


/* Martingale & deviance residuals */
symbol1 c=red v=star i=none;
symbol2 c=blue v=circle i=none;
axis1 order=(0 to 200 by 50) minor=none label=('Weeks');
axis2 order=(-6 to 1 by 1) minor=none label=('MR');
proc gplot data=nardi;
  plot martweibull*weeks=event/vaxis=axis2 haxis=axis1 vref=0;
  title 'Martingale Residual Plots - Weibull Model';
run; 
quit;

proc gplot data=nardi;
  plot martllog*weeks=event/vaxis=axis2 haxis=axis1 vref=0;
  title 'Martingale Residual Plots - Log Logistic Model';
run; 
quit;

symbol1 c=red v=star i=none;
symbol2 c=blue v=circle i=none;
axis1 order=(0 to 200 by 50) minor=none label=('Weeks');
axis2 order=(-10 to 5 by 1) minor=none label=('DR');
proc gplot data=nardi;
  plot devweibull*weeks=event/vaxis=axis2 haxis=axis1 vref=0;
  title 'Deviance Residual Plots - Weibull Model';
run; 
quit;

proc gplot data=nardi;
  plot devllog*weeks=event/vaxis=axis2 haxis=axis1 vref=0;
  title 'Deviance Residual Plots - Log Logistic Model - Age Grouped';
run; 
quit;

proc gplot data=nardi;
  plot devllog2*weeks=event/vaxis=axis2 haxis=axis1 vref=0;
  title 'Deviance Residual Plots - Log Logistic Model - Age Continuous';
run; 
quit;

proc gplot data=nardi;
  plot devlnorm*weeks=event/vaxis=axis2 haxis=axis1 vref=0;
  title 'Deviance Residual Plots - Log Normal Model - Age Grouped';
run; 
quit;



/* model interpretation - code developed to create summary tables and plots */

data models;
  /* choose covariate values */
  age50=1;  /*  age>=50 */

  label 
    weeks='t (weeks)'
	age50='Age>=50'
    sexp1='Exponential S(t), Treatment=1'
    hexp1='Exponential h(t), Treatment=1'
    mexp1='Exponential Median, Treatment=1'
    sexp0='Exponential S(t), Treatment=0'
    hexp0='Exponential h(t), Treatment=0'
    mexp0='Exponential Median, Treatment=0'
    hrexp='Exponential Hazard Ratio(t)'
    beta_exp='Exponential beta(Treatment);'
    trexp='Exponential Time Ratio'
    mrexp='Exponential Median Ratio'
    sweibull1='Weibull S(t), Treatment=1'
    hweibull1='Weibull h(t), Treatment=1'
    mweibull1='Weibull Median, Treatment=1'
    sweibull0='Weibull S(t), Treatment=0'
    hweibull0='Weibull h(t), Treatment=0'
    mweibull0='Weibull Median, Treatment=0'
    hrweibull='Weibull Hazard Ratio(t)'
    beta_weibull='Weibull beta(Treatment);'
    trweibull='Weibull Time Ratio'
    mrweibull='Weibull Median Ratio'
    slnorm1='Log Normal S(t), Treatment=1'
    hlnorm1='Log Normal h(t), Treatment=1'
    mlnorm1='Log Normal Median, Treatment=1'
    slnorm0='Log Normal S(t), Treatment=0'
    hlnorm0='Log Normal h(t), Treatment=0'
    mlnorm0='Log Normal Median, Treatment=0'
    hrlnorm='Log Normal Hazard Ratio(t)'
    beta_lnorm='Log Normal beta(Treatment);'
    trlnorm='Log Normal Time Ratio'
    mrlnorm='Log Normal Median Ratio'
    sllog1='Log Logistic S(t), Treatment=1'
    hllog1='Log Logistic h(t), Treatment=1'
    mllog1='Log Logistic Median, Treatment=1'
    sllog0='Log Logistic S(t), Treatment=0'
    hllog0='Log Logistic h(t), Treatment=0'
    mllog0='Log Logistic Median, Treatment=0'
    hrllog='Log Logistic Hazard Ratio(t)'
    beta_llog='Log Logistic beta(Treatment);'
    trllog='Log Logistic Time Ratio'
    mrllog='Log Logistic Median Ratio'
	odds1='Log Logistic Odds S(t)/(1-S(t), Treatment=1'
	odds0='Log Logistic Odds S(t)/(1-S(t), Treatment=0'
	oddsratio='Log Logisitic Odds Ratio'
	alpha1='Log Logistic Alpha, Treatment=1'
	alpha0='Log Logistic Alpha, Treatment=0'
	alpharatio='Log Logisitic Alpha Ratio'
  ;

  /* generate values on the survival curves */


  do weeks=0 to 200 by 1;  /* time frame */

  /* calculate S(t) and h(t) */

  treat=1; /* treatment=yes survival curve */

  /* exponential */
  lambda=exp(-(3.9439+0.1825*treat-0.5007*age50));
  sexp1=exp(-lambda*weeks);
  hexp1=lambda;
  mexp1= -log(.5)/lambda;

  /* weibull */
  lambda=exp(-1.0744*(3.9622+0.1825*treat+-0.4904*age50));
  gamma=1.0744;
  sweibull1=exp(-lambda*weeks**gamma);
  if weeks>0 then hweibull1=gamma*lambda*(weeks**(gamma-1));
  mweibull1= (-log(.5)/lambda)**(1/gamma);  /* or (-log(.5)/lambda)**scale */

  /* lognormal */
  u=3.4768+0.1744*treat-0.4144*age50;
  scale=0.9288;
  pi=3.14159; 
  if weeks>0 then slnorm1=1-probnorm((log(weeks)-u)/scale);
  if weeks>0 then flnorm1=1/(sqrt(2*pi)*weeks*scale)*exp(-1/2*((log(weeks)-u)/scale)**2);
  hlnorm1=flnorm1/slnorm1;
  mlnorm1=exp(scale*probit(.5)+u);

  /* loglogistic */
  alpha1=exp(-(3.4289+0.2029*treat-0.4204*age50)/0.5198);
  gamma=1/0.5198;
  sllog1=1/(1+alpha1*weeks**gamma);
  if weeks>0 then fllog1=alpha1*gamma*weeks**(gamma-1)/(1+alpha1*weeks**gamma)**2;
  hllog1=fllog1/sllog1;
  mllog1=(1/alpha1)**(1/gamma);

  treat=0; /* treatment=no survival curve */

  /* exponential */
  lambda=exp(-(3.9439+0.1825*treat-0.5007*age50));
  sexp0=exp(-lambda*weeks);
  hexp0=lambda;
  mexp0= -log(.5)/lambda;

  /* weibull */
  lambda=exp(-1.0744*(3.9622+0.1825*treat+-0.4904*age50));
  gamma=1.0744;
  sweibull0=exp(-lambda*weeks**gamma);
  if weeks>0 then hweibull0=gamma*lambda*(weeks**(gamma-1));
  mweibull0= (-log(.5)/lambda)**(1/gamma);  /* or (-log(.5)/lambda)**scale */

  /* lognormal */
  u=3.4768+0.1744*treat-0.4144*age50;
  scale=0.9288;
  pi=3.14159; 
  if weeks>0 then slnorm0=1-probnorm((log(weeks)-u)/scale);
  if weeks>0 then flnorm0=1/(sqrt(2*pi)*weeks*scale)*exp(-1/2*((log(weeks)-u)/scale)**2);
  hlnorm0=flnorm0/slnorm0;
  mlnorm0=exp(scale*probit(.5)+u);

  /* loglogistic */
  alpha0=exp(-(3.4289+0.2029*treat-0.4204*age50)/0.5198);
  gamma=1/0.5198;
  sllog0=1/(1+alpha0*weeks**gamma);
  if weeks>0 then fllog0=alpha0*gamma*weeks**(gamma-1)/(1+alpha0*weeks**gamma)**2;
  hllog0=fllog0/sllog0;
  mllog0=(1/alpha0)**(1/gamma);

  /* model - treatment comparison */
  hrexp=hexp1/hexp0;
  beta_exp=-log(hrexp);  /* hrexp = exp(-beta_exp) */
  trexp=exp(beta_exp);
  mrexp=mexp1/mexp0;

  hrweibull=hweibull1/hweibull0;
  beta_weibull=-log(hrweibull)/1.0744;  /* weibull gamma; hrweibull = exp(-beta_weibull*1.0744) */
  trweibull=exp(beta_weibull);    
  mrweibull=mweibull1/mweibull0; 


  hrlnorm=hlnorm1/hlnorm0;
  beta_lnorm=0.9288*(probit(1-slnorm0) - probit(1-slnorm1));	/* log normal scale */
  trlnorm=exp(beta_lnorm);  
  mrlnorm=mlnorm1/mlnorm0;


  hrllog=hllog1/hllog0;
  if sllog1^=1 then odds1=sllog1/(1-sllog1);
  if sllog0^=1 then odds0=sllog0/(1-sllog0);
  oddsratio=odds1/odds0;
  alpharatio=alpha0/alpha1;
  beta_llog=log(oddsratio)*0.5198; /* log logistic scale; oddsratio=exp(beta_llog/0.5198) */
  trllog=exp(beta_llog); 
  mrllog=mllog1/mllog0;

  output;
  end;
  drop treat gamma lambda u pi scale flnorm1 flnorm0 fllog1 fllog0;
  format age50 yn. sexp1--alpharatio 7.4;
run;


proc print data=models label noobs;
  where weeks in (26, 52, 104);
  id weeks;
  var age50 sexp1--mrexp;
  title 'Exponential Model';
run;

proc print data=models label noobs;
  where weeks in (26, 52, 104);
  id weeks;
  var age50 sweibull1--mrweibull;
  title 'Weibull Model';
run;

proc print data=models label noobs;
  where weeks in (26, 52, 104);
  id weeks;
  var age50 slnorm1--mrlnorm;
  title 'Log normal Model';
run;

proc print data=models label noobs;
  where weeks in (26, 52, 104);
  id weeks;
  var age50 sllog1--alpharatio;
  title 'Log logistic Model';
run;


axis1 order=(1.15 to 1.25 by .05) minor=none label=('TR');
axis2 order=(0 to 200 by 50) minor=none label=('weeks');
symbol1 c=red v=point i=join r=1 l=1;
symbol2 c=blue v=point i=join r=1 l=3;
symbol3 c=black v=point i=join r=1 l=5;
symbol4 c=green v=point i=join r=1 l=8;
proc gplot data=models;
  plot (trexp trweibull trlnorm trllog)*weeks/overlay legend vaxis=axis1 haxis=axis2;
  format trexp trweibull trlnorm trllog 6.2;
  title 'Time Ratios for Treatment - All Models';
run;
quit;

data sampleplots;
  set atsurv(where=(_censor_ in (0,.) and treat=1 and age50=1))
      models(keep=weeks sexp1 sweibull1 slnorm1 sllog1);
run;

symbol1 c=red v=point i=join r=1 l=1;
symbol2 c=blue v=point i=join r=1 l=3;
symbol3 c=black v=point i=join r=1 l=5;
symbol4 c=green v=point i=join r=1 l=8;
symbol5 c=black v=point i=stepjl l=2;
axis1 order=(0 to 1 by .1) minor=none label=('S(t)');
axis2 order=(0 to 200 by 50) minor=none label=('weeks');
proc gplot data=sampleplots;
  plot (sexp1 sweibull1 slnorm1 sllog1 survival)*weeks
       /overlay legend vaxis=axis1 haxis=axis2;
  label survival='KM' sexp1='Exponential' sweibull1='Weibull'
        slnorm1='LogNormal' sllog1='LogLogistic'
        weeks='Weeks';
  title 'Survival: Treatment=Yes, Age>=50';
run;
quit;

data sampleplots;
  set atsurv(where=(_censor_ in (0,.) and treat=0 and age50=1))
      models(keep=weeks sexp0 sweibull0 slnorm0 sllog0);
run;

axis3 order=(0 to 0.05 by 0.01) minor=none label=('h(t)');
proc gplot data=models;
  plot (hexp1 hweibull1 hlnorm1 hllog1)*weeks
       /overlay legend vaxis=axis3 haxis=axis2;
  label hexp1='Exponential' hweibull1='Weibull'
        hlnorm1='LogNormal' hllog1='LogLogistic'
        weeks='Weeks';
  title 'Hazard: Treatment=Yes, Age>=50';
run;
quit;

symbol1 c=red v=point i=join r=1 l=1;
symbol2 c=blue v=point i=join r=1 l=3;
symbol3 c=black v=point i=join r=1 l=5;
symbol4 c=green v=point i=join r=1 l=8;
symbol5 c=black v=point i=stepjl l=2;
axis1 order=(0 to 1 by .1) minor=none label=('S(t)');
axis2 order=(0 to 200 by 50) minor=none label=('weeks');
proc gplot data=sampleplots;
  plot (sexp0 sweibull0 slnorm0 sllog0 survival)*weeks
       /overlay legend vaxis=axis1 haxis=axis2;
  label survival='KM' sexp0='Exponential' sweibull0='Weibull'
        slnorm0='LogNormal' sllog0='LogLogistic'
        weeks='Weeks';
  title 'Survival: Treatment=No, Age>=50';
run;
quit;

axis3 order=(0 to 0.05 by 0.01) minor=none label=('h(t)');
proc gplot data=models;
  plot (hexp0 hweibull0 hlnorm0 hllog0)*weeks
       /overlay legend vaxis=axis3 haxis=axis2;
  label hexp0='Exponential' hweibull0='Weibull'
        hlnorm0='LogNormal' hllog0='LogLogistic'
        weeks='Weeks';
  title 'Hazard: Treatment=No, Age>=50';
run;
quit;

/* alternate way to obtain fitted model (adjusted) */
proc lifereg data=sda.brain noprint;
  model weeks*event(0)=treat age50/d=gamma;
  /* output the percentile of the CDF */
  output out=gbrainadj p=pred q=0.01 to .99 by .01;
  title 'Gamma-adjusted';
run;

data atsurv2;
  /* bring in the KM estimates */
  set atsurv(where=( _censor_ in (0,.) and treat=1 and age50=1))   /* treated, 50+ */
     /* since gamma percentiles are repeated for each observation, keep the first set of percentiles with covariates of interest */
     gbrainadj(firstobs=298 obs=396 keep=_prob_ pred);
  /* calculate SDF=1-CDF */
  if _prob_^=. then sgamma=1-_PROB_;
run;

symbol1 c=red v=point i=join r=1 l=1;
symbol2 c=black v=point i=stepjl l=2;
axis1 order=(0 to 1 by .1) minor=none label=('S(t)');
axis2 order=(0 to 200 by 50) minor=none label=('weeks');
proc gplot data=atsurv2;
  plot sgamma*pred survival*weeks
       /overlay legend vaxis=axis1 haxis=axis2;
  label survival='KM' sgamma='Gamma'
        weeks='Weeks';
  title 'Survival: Treatment=Yes, Age>=50';
run;
quit;

/* changepoint model (piecewise exponential) */
data brain2(keep=id weeks event weeks2 event2 year1);
  set sda.brain;
  id=_n_;
  if weeks<=52
     then do;
	    event2=event;
		weeks2=weeks;
		year1=1;
		output;
	 end;
     else do;
	    event2=0;
		weeks2=52;
        year1=1;
		output;
	    event2=event;
		weeks2=weeks-52;
        year1=0;
		output;
	 end;
run;

proc lifereg data=brain2;
  model weeks2*event2(0)=year1/d=exponential;
  title 'Changepoint model';
run;

data brain3;
  do weeks=0 to 200 by 1;  /* time frame */
  lambda1=exp(-(4.533-.9175));
  lambda2=exp(-(4.533));
  if weeks<=52
     then sexp=exp(-lambda1*weeks);
     else sexp=exp(-lambda1*52)*exp(-lambda2*(weeks-52));
  output;
  end;
run;

data sampleplots;
  set osurv(where=(_censor_ in (0,.)))
      brain3;
run;

symbol1 c=red v=point i=join r=1 l=1;
symbol2 c=black v=point i=stepjl l=2;
axis1 order=(0 to 1 by .1) minor=none label=('S(t)');
axis2 order=(0 to 200 by 50) minor=none label=('weeks');
proc gplot data=sampleplots;
  plot (sexp survival)*weeks
       /overlay legend vaxis=axis1 haxis=axis2;
  label survival='KM' sexp='CP' 
        weeks='Weeks';
  title 'Survival: Changepoint model (Tau=52 weeks)';
run;
quit;


/* Gamel-Boag model */

/* KM estimates */
proc lifetest data=sda.brain outsurv=asurv noprint plots=none;
  time weeks*event(1);
  strata age50;
run;

/* log normal model using proc lifereg */
proc lifereg data=sda.brain;
  model weeks*event(0)=age50/d=lnormal;
  title 'LifeReg: Survival by age group (Log normal)';
run;

/* log normal model using proc nlp (SAS/OR 12.2 legacy procedure) - maximize loglikelihood */
proc nlp data=sda.brain tech=tr cov=2 stderr;
   parms int gamma sig;
   pi=3.14159; 
   u=int+gamma*age50;
   if event=1 then logl=log(1/(sqrt(2*pi)*weeks*sig)*exp(-1/2*((lweeks-u)/sig)**2));
   if event=0 then logl=log(1-probnorm((lweeks-u)/sig));
   max logl;
run;

/* Gamel-Boag cure model using proc nlp (SAS/OR), maximize modified loglikelihood, reference Frankel & Longmate */
proc nlp data=sda.brain tech=tr cov=2 stderr;
   parms intg gamma intb beta sig;
   pi=3.14159; 
   p=exp(intb+beta*age50)/(1+exp(intb+beta*age50)); /* model proportion cured by age group */
   u=intg+gamma*age50;
   if event=1 then logl=log((1-p)*(1/(sqrt(2*pi)*weeks*sig)*exp(-1/2*((lweeks-u)/sig)**2)));
   if event=0 then logl=log(p+(1-p)*(1-probnorm((lweeks-u)/sig)));
   max logl;
run;

/* create data for plotting models fitted above */
data models;
  do weeks=1 to 200 by 1;  /* time frame */
  /* lognormal, age>=50 */
  age50=1;
  u=3.5643-0.4161*age50;
  sig=0.9335;
  pi=3.14159; 
  slnorm=1-probnorm((log(weeks)-u)/sig);
  /* Gamel-Boag cure model, age>=50 */
  u=3.3833-0.2413*age50;
  p=exp(-2.3547-3.6952*age50)/(1+exp(-2.3547-3.6952*age50));
  sig=0.8419;
  pi=3.14159; 
  sgb=p+(1-p)*(1-probnorm((log(weeks)-u)/sig));
  output;
  end;
  do weeks=1 to 200 by 1;  /* time frame */
  /* lognormal, age<50 */
  age50=0;
  u=3.5643-0.4161*age50;
  sig=0.9335;
  pi=3.14159; 
  slnorm=1-probnorm((log(weeks)-u)/sig);
  /* Gamel-Boag cure model, age<50 */
  u=3.3833-0.2413*age50;
  p=exp(-2.3547-3.6952*age50)/(1+exp(-2.3547-3.6952*age50));
  sig=0.8419;
  pi=3.14159; 
  sgb=p+(1-p)*(1-probnorm((log(weeks)-u)/sig));
  output;
  end;
run;

data sampleplots;
  set asurv(where=(_censor_ in (0,.) and age50=1))
      models(keep=weeks age50 slnorm sgb where=(age50=1));
run;

symbol1 c=red v=point i=join r=1 l=1;
symbol2 c=blue v=point i=join r=1 l=3;
symbol3 c=black v=point i=stepjl l=2;
axis1 order=(0 to 1 by .1) minor=none label=('S(t)');
axis2 order=(0 to 200 by 50) minor=none label=('weeks');
proc gplot data=sampleplots;
  plot (slnorm sgb survival)*weeks
       /overlay legend vaxis=axis1 haxis=axis2;
  label survival='KM' 
        slnorm='LogNormal' 
        sgb='Gamel-Boag'
        weeks='Weeks';
  title 'Survival: Age>=50';
run;
quit;

data sampleplots;
  set asurv(where=(_censor_ in (0,.) and age50=0))
      models(keep=weeks age50 slnorm sgb where=(age50=0));
run;

symbol1 c=red v=point i=join r=1 l=1;
symbol2 c=blue v=point i=join r=1 l=3;
symbol3 c=black v=point i=stepjl l=2;
axis1 order=(0 to 1 by .1) minor=none label=('S(t)');
axis2 order=(0 to 200 by 50) minor=none label=('weeks');
proc gplot data=sampleplots;
  plot (slnorm sgb survival)*weeks
       /overlay legend vaxis=axis1 haxis=axis2;
  label survival='KM' 
        slnorm='LogNormal' 
        sgb='Gamel-Boag'
        weeks='Weeks';
  title 'Survival: Age<50';
run;
quit;

/* relate estimates of cure probabilites back to model parameters */
data survprob;
   /* age>=50 */
   p1=exp(-2.3547-3.6952)/(1+exp(-2.3547-3.6952));
   /* age<50 */
   p0=exp(-2.3547)/(1+exp(-2.3547));
   /* odds ratio(odds: cure to not cured) age>=50 to age<50 */
   or=(p1/(1-p1))/(p0/(1-p0));
   /* log odds ratio, modeled as beta above */
   lor=log(or);
run;

proc print data=survprob noobs;
  var p1 p0 or lor;
  format p1 p0 or lor 6.3;
  title 'Cure model: survival probabilities';
run;

/* overall survival probability of 0.038, i.e. not related to age */
proc nlp data=sda.brain tech=tr cov=2 stderr;
   parms intg gamma intb sig;
   pi=3.14159; 
   p=exp(intb)/(1+exp(intb)); /* model proportion cured overall */
   u=intg+gamma*age50;
   if event=1 then logl=log((1-p)*(1/(sqrt(2*pi)*weeks*sig)*exp(-1/2*((lweeks-u)/sig)**2)));
   if event=0 then logl=log(p+(1-p)*(1-probnorm((lweeks-u)/sig)));
   max logl;
run;

/* Bayesian analysis by Gibbs sampling in the location-scale models  */

/* Weibull model */
proc lifereg data=sda.brain;
  model weeks*event(0)=age50/d=weibull;
  title 'LifeReg: Age group - Weibull';
run;

/* Weibull model-bayesian analysis */
proc lifereg data=sda.brain;
  model weeks*event(0)=age50/d=weibull;
  bayes WeibullShapePrior=gamma seed=1254 outpost=postweibull;
  title 'LifeReg: Weibull,Bayesian - prior: shape parameter:gamma age:uniform';
run;

ods rtf close;
