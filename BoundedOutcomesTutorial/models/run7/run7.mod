; PRO Tutorial Beta Regression ------------------------------------------------
$PROBLEM run7.mod
; Identifier: ID
; Dependent Variable: DV
; Time Variable: TAFD
; Subject Level: ID STDID SEX RACE MENOS ECOG DTH CHEMO SURGERY RADIO MHDEP 
;   MDANX ERS PGRS AGE LESION1 HEIGHT WEIGHT DTHDY PFSDY
; Observation Level: C TAFD TAFR SCALE FLAG DV
; Dropped Variables: none

$INPUT C ID TAFD=TIME TAFR SCALE STDID SEX RACE MENOS ECOG DTH CHEMO SURGERY 
  RADIO MHDEP MDANX ERS PGRS AGE LESION1 HEIGHT WEIGHT DTHDY PFSDY NQS FLAG DV
  
$DATA EQLQ_subscores_21OCT2025.csv
  IGNORE=C
  IGNORE=(SCALE.NE.2) ; Physical Functioning only
  IGNORE=(FLAG.NE.5) ; Absolute scores only
  
$ABBR FUNCTION BETA(VBT, 10) ; define VBT for calculating beta PDF
  
$PRED
; Fixed parameters
  DVLO = 1.0*NQS
  DVUP = 4.0*NQS
  
; Population parameters
  POPTBASE = (THETA(1)*NQS - DVLO)/(DVUP - DVLO) ; Baseline physical functioning score ([1,4] -> [0,1])
  POPLGTBASE = LOG(POPTBASE/(1 - POPTBASE)) ; Transformed baseline physical functioning score ([0,1] -> [-Inf,Inf])
  POPPMAX = THETA(2) ; Maximum placebo effect
  POPKPBO = THETA(3) ; Placebo rate constant
  POPPHI = THETA(4) ; Beta regression precision parameter
  POPZETA1 = THETA(5) ; Zero-/One-inflation intercept parameter
  POPZETA2 = THETA(6) ; Zero-/One-inflation slope parameter
  
; Individual parameters
  LGTBASE = POPLGTBASE + ETA(1)
  PMAX = POPPMAX + ETA(2)
  KPBO = POPKPBO
  PHI = POPPHI
  ZETA1 = POPZETA1
  ZETA2 = POPZETA2  
  
; Time handling
  PTIME = TIME
  IF (TIME.LE.0) PTIME = 0
  
; Define model
; Latent variable has domain [-Inf,Inf]
  PEFF = PMAX*(1 - EXP(-KPBO*PTIME)); placebo effect
  LGTMU = LGTBASE - PEFF
  MU = 1/(1 + EXP(-LGTMU))
  
; Zero- and one-inflation
  LGTP0 = -ZETA1 - ZETA2*LGTMU
  LGTP1 = -ZETA1 + ZETA2*LGTMU
  P0 = 1E-6 + 1/(1 + EXP(-LGTP0))
  P1 = 1E-6 + 1/(1 + EXP(-LGTP1))
  
; Predictions (i.e. expectation)
  IPRED = (0*P0 + MU*(1 - P0 - P1) + 1*P1)*(DVUP - DVLO) + DVLO
  IF (COMACT.EQ.1) PPRED = IPRED
  
; Transform observed data
; Score has domain [1,4], needs to be [0,1] for beta regression
  TDV = (DV - DVLO)/(DVUP - DVLO)
  
; Beta distribution
  VBT(1) = TDV ; theta - observation
  VBT(2) = MU*PHI ; alpha - shape parameter 
  VBT(3) = (1 - MU)*PHI ; beta - shape parameter  
  BTLL = -BETA(VBT) + LOG(1 - P0 - P1)
  
; Log-Likelihood
  IF (DV.GT.DVLO.AND.DV.LT.DVUP) LL = BTLL
  IF (DV.EQ.DVLO) LL = LOG(P0)   
  IF (DV.EQ.DVUP) LL = LOG(P1)
  Y = -2*LL
  
; Pearson weighted residuals
  PWRES = (TDV - MU)/(SQRT(MU*(1 - MU)/(1 + PHI)))
  
$THETA 
  (1.00001, 1.7, 3.99999) ; TVBASE	
  0.35 ; TVPMAX	
  (0.00001, 0.006) ; TVKPBO
  (0.00001, 19) ; TVPHI
  9 ; TVZETA1
  (0.00001, 4) ; TVZETA2

$OMEGA
  0.8 ; PPVBASE
  1.2 ; PPVPMAX
  
;$SIGMA
;  1 FIX ; ERRDUM
 
;$SIM(555) ONLYSIM SUBPROB=1
  
$ESTIMATION METHOD=SAEM -2LL LAPLACIAN NUMERICAL SLOW AUTO=1 CTYPE=3 SEED=555
  CALPHA=0.01 SIGL=8 NOPRIOR=1 NOABORT
$ESTIMATION METHOD=IMP EONLY=1 AUTO=1 NITER=40 CTYPE=3 
  PRINT=5 SIGL=8 NOPRIOR=1 FILE=run7.ext MSFO=run7.msf

;$ESTIMATION MAXEVAL=0 PRINT=1 METHOD=COND -2LL NUMERICAL LAPLACIAN NOABORT 
;  SIG=3 SIGL=9 FILE=run7.ext MSFO=run7.msf

$COVARIANCE PRINT=E UNCONDITIONAL

$TABLE ID TIME TAFR DV PPRED IPRED PWRES SCALE FLAG BTLL MU P0 P1 LGTMU 
  POPLGTBASE TDV
  NOPRINT NOAPPEND ONEHEADER IDFORMAT=I FILE=run7.fit

$TABLE ID ETAS(1:LAST) OBJI STDID SEX RACE MENOS ECOG DTH CHEMO SURGERY RADIO 
  MHDEP MDANX ERS PGRS AGE LESION1 HEIGHT WEIGHT DTHDY PFSDY
  NOPRINT NOAPPEND FIRSTONLY ONEHEADER IDFORMAT=I FILE=run7.eta
  
; - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 