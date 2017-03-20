* real data analyze;
LIBNAME realdata 'C:\Users\hcg0630\Desktop\MZIBR\realdata';
DATA realdata.MTPestimate;
	SET _NULL_;
RUN;
DATA realdata.MTPfit;
	SET _NULL_;
RUN;
DATA realdata.TPMestimate;
	SET _NULL_;
RUN;
DATA  realdata.TPMfit;
	SET _NULL_;
RUN;

PROC IMPORT DATAFILE = "C:\Users\hcg0630\Desktop\MZIBR\realdata\realdata.csv" OUT = realdata.mbdata REPLACE;
RUN;
data baselinedata;
set realdata.mbdata;
if time = 1;
run;

data realdata.one;
set baselinedata;
Al=0; Ba=0; Bi=0; Clo=0; col=0; cop=0; di=0; do=0; es=0; eu=0; fa=0; h=0; l=0; p=0; ro=0; ru=0; s=0; v=0;
if Alistipes>0 then Al=1;
if Bacteroides>0 then Ba=1;
if Bifidobacterium>0 then Bi=1;
if Clostridium>0 then clo=1;
if Collinsella>0 then col=1;
if Coprobacillus>0 then cop=1;
if Dialister>0 then di=1;
if Dorea>0 then do=1;
if Escherichia>0 then es=1;
if Eubacterium>0 then eu=1;
if Faecalibacterium>0 then fa=1;
if Haemophilus>0 then h=1;
if Lactobacillus>0 then l=1;
if Parabacteroides>0 then p=1;
if Roseburia>0 then ro=1;
if Ruminococcus>0 then ru=1;
if Streptococcus>0 then s=1;
if Veillonella>0 then v=1;
run;

%macro TPM(var, var2);
PROC nlmixed DATA=realdata.one METHOD=GAUSS;
*bounds vara>0;
/* log-likelihood */
teta = alpha0 + alpha2*Treatment;* + alpha1*Time  + alpha3*Baseline&var2;
expteta = exp(teta);
p = expteta / (1 + expteta);
model &var ~ binary(p);
*random a ~ normal(0, vara) subject=Subject;
ods output ParameterEstimates=est1 FitStatistics=fit1; 
run;
data two;
set realdata.one;
if &var>0;
run; 
* Part II for Beta distribution;
PROC nlmixed DATA=two METHOD=GAUSS;
* bounds varb phi>0;
bounds phi > 0;
/* log-likelihood */
	eta=gamma0 + gamma2*Treatment;* + gamma1*Time + gamma3*Baseline&var2;
	mu = exp(eta)/(1+exp(eta));
	loglik = (mu*phi-1)*log(&var2) + ((1-mu)*phi-1)*log(1-&var2) + lgamma(phi) - lgamma(mu*phi) - lgamma((1-mu)*phi);
model &var2 ~ general(loglik);
* random b ~ normal(0, varb) subject=Subject;
ods output ParameterEstimates=est2 FitStatistics=fit2; 
run;

***************************************conventional two-part model******************************************;
data initialvalue1;
set est1 est2;
run;
PROC nlmixed DATA=realdata.one METHOD=GAUSS;
/* define initial values and bounds */
parms/data=initialvalue1;
* bounds vara varb phi>0;
bounds phi > 0;
/* log-likelihood */
teta = alpha0 + alpha2*Treatment;* + alpha1*Time + alpha3*Baseline&var2;
expteta = exp(teta);
p = expteta / (1 + expteta);
/* For zero Faecalibacterium */
if &var2=0 then loglik=log(1-p);
/* For positive Faecalibacterium */
if &var2>0 then do;
	eta=gamma0 + gamma2*Treatment;* + gamma1*Time + gamma3*Baseline&var2;
	mu = exp(eta)/(1+exp(eta));
	loglik = log(p) + (mu*phi-1)*log(&var2) + ((1-mu)*phi-1)*log(1-&var2) + lgamma(phi) - lgamma(mu*phi) - lgamma((1-mu)*phi);
end;
/* fit the above modle */
model &var2 ~ general(loglik);
* random a b~ normal([0,0],[vara,0,varb]) subject=Subject;
ods output ParameterEstimates=TPM&var2._est3 FitStatistics=TPM&var2._fit3; 
run;

DATA TPM&var2._est;
	SET TPM&var2._est3;
	genera = "&var2"; * this indicates that this is results for genera var2;
RUN;
DATA TPM&var2._fit; 
	SET TPM&var2._fit3;
	genera = "&var2";
RUN;
* put the results of every simulation in the same data set: result.estimates ;
DATA realdata.TPMestimate;
	SET realdata.TPMestimate TPM&var2._est;
RUN;
DATA realdata.TPMfit;
	SET realdata.TPMfit TPM&var2._fit;
RUN;

********************************************marginalized two-part model*****************************************;
PROC NLMIXED DATA = realdata.one;
* parameters and initial values;
parms / data = TPM&var2._est3;
bounds phi > 0;
* loglik;
temp1 = alpha0+ alpha2 * treatment;* + alpha1 * time ;
p = exp(temp1) / (1 + exp(temp1));
if &var2 = 0 then 
	loglik = log(1 - p);
if &var2 > 0 then do;
	temp2 = gamma0 + gamma2 * treatment;* + gamma1 * time;
	mu = (1 + exp(-temp1)) / (1 + exp(-temp2));
	loglik = log(p) + lgamma(phi) - lgamma(mu * phi) - lgamma((1 - mu) * phi) +
			 (mu * phi - 1) * log(&var2) + ((1 - mu) * phi - 1) * log(1 - &var2);
end;
* fit model;
model &var2 ~ general(loglik);
* save results;
ODS OUTPUT ParameterEstimates = MTP&var2._est FitStatistics = MTP&var2._fit;
*ODS SELECT NONE;
RUN;
* if simuest&j is created, then combine it with result.estimates;
* IF EXIST(&var2._est) THEN DO;
DATA MTP&var2._est;
	SET MTP&var2._est;
	genera = "&var2"; * this indicates that this is results for genera var2;
RUN;
DATA MTP&var2._fit; 
	SET MTP&var2._fit;
	genera = "&var2";
RUN;
* put the results of every simulation in the same data set: result.estimates ;
DATA realdata.MTPestimate;
	SET realdata.MTPestimate MTP&var2._est;
RUN;
DATA realdata.MTPfit;
	SET realdata.MTPfit MTP&var2._fit;
RUN;
* END;
%mend;

%TPM(al, Alistipes); %TPM(ba, Bacteroides); %TPM(bi, Bifidobacterium); %TPM(clo, Clostridium); %TPM(col, Collinsella); %TPM(cop, Coprobacillus); 
%TPM(di, Dialister); %TPM(do, Dorea); %TPM(es, Escherichia); %TPM(eu, Eubacterium); %TPM(fa, Faecalibacterium); %TPM(h, Haemophilus); 
%TPM(l, Lactobacillus); %TPM(p, Parabacteroides); %TPM(ro, Roseburia); %TPM(ru, Ruminococcus); %TPM(s, Streptococcus); %TPM(v, Veillonella); 

* save the results;
PROC EXPORT DATA = realdata.MTPestimate OUTFILE = "C:\Users\hcg0630\Desktop\MZIBR\realdata\MTPrealdata_estimate.csv" REPLACE;
RUN;
PROC EXPORT DATA = realdata.MTPfit OUTFILE = "C:\Users\hcg0630\Desktop\MZIBR\realdata\MTPrealdata_fit.csv" REPLACE;
RUN;
* save the results;
PROC EXPORT DATA = realdata.TPMestimate OUTFILE = "C:\Users\hcg0630\Desktop\MZIBR\realdata\TPMrealdata_estimate.csv" REPLACE;
RUN;
PROC EXPORT DATA = realdata.TPMfit OUTFILE = "C:\Users\hcg0630\Desktop\MZIBR\realdata\TPMrealdata_fit.csv" REPLACE;
RUN;
