; PRO Tutorial Beta Regression ------------------------------------------------
$PROBLEM run9.mod
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
  
$ABBR FUNCTION BETACDF(VBL, 10) ; define function BETACDF and vector VBL
$ABBR VECTOR VBU(10) ; define vector VBU
  
$PRED
; Fixed parameters
  DVLO = 1.0*NQS
  DVUP = 4.0*NQS
  
; Population parameters
  POPTBASE = (THETA(1)*NQS - DVLO)/(DVUP - DVLO) ; Baseline physical functioning score ([1,4] -> [0,3])
  POPLGTBASE = LOG(POPTBASE/(1 - POPTBASE)) ; Transformed baseline physical functioning score ([0,3] -> [0,1])
  POPPMAX = THETA(2) ; Maximum placebo effect
  POPKPBO = THETA(3) ; Placebo rate constant
  POPKDIS = THETA(4) ; Disease progression rate constant
  POPPHI = THETA(5) ; Beta regression precision parameter
  
; Individual parameters
  LGTBASE = POPLGTBASE + ETA(1)
  PMAX = POPPMAX + ETA(2)
  KPBO = POPKPBO
  KDIS = POPKDIS*EXP(ETA(3))
  PHI = POPPHI
  
; Time handling
  PTIME = TIME
  IF (TIME.LE.0) PTIME = 0
  
; Define model
; Latent variable has domain [-Inf,Inf]
  PEFF = PMAX*(1 - EXP(-KPBO*PTIME)); placebo effect
  DPRO = EXP(KDIS*TIME) - 1 ; disease progression
  LGTMU = LGTBASE - PEFF + DPRO
  MU = 1/(1 + EXP(-LGTMU))
  IF (COMACT.EQ.1) PMU = MU ; population prediction
  
; Transform observed data
; Score has domain [1,4], needs to be [0,1] for beta regression
  TDVL = (DV - DVLO)/(DVUP - DVLO + 1)
  TDVU = (DV - DVLO + 1)/(DVUP - DVLO + 1)
  
; Lower beta cumulative distribution
  VBL(1) = TDVL ; theta - observation
  VBL(2) = MU*PHI ; alpha - shape parameter 
  VBL(3) = (1 - MU)*PHI ; beta - shape parameter  
  
; Upper beta cumulative distribution
  VBU(1) = TDVU ; theta - observation
  VBU(2) = MU*PHI ; alpha - shape parameter 
  VBU(3) = (1 - MU)*PHI ; beta - shape parameter  
  
; Log-Likelihood
  DUM = 1.0E-30 ; placeholder
  BTL = BETACDF(VBU) - BETACDF(VBL)
  Y = DUM
  IF (BTL.GE.DUM) Y = BTL
  
$THETA 
  (1.00001, 1.7, 3.99999) ; TVBASE	
  0.35 ; TVPMAX	
  (0.00001, 0.006) ; TVKPBO
  (0.00001, 0.0004) ; TVKDIS
  (0.00001, 19) ; TVPHI

$OMEGA
  0.8 ; PPVBASE
  1.2 ; PPVPMAX
  3.6 ; PPVKDIS
  
;$SIGMA
;  1 FIX ; ERRDUM
 
;$SIM(555) ONLYSIM SUBPROB=1
  
$ESTIMATION METHOD=SAEM LIKELIHOOD LAPLACIAN NUMERICAL SLOW AUTO=1 CTYPE=3 SEED=555
  CALPHA=0.01 SIGL=8 NOPRIOR=1 NOABORT
$ESTIMATION METHOD=IMP EONLY=1 AUTO=1 NITER=40 CTYPE=3 
  PRINT=5 SIGL=8 NOPRIOR=1 FILE=run9.ext MSFO=run9.msf

;$ESTIMATION MAXEVAL=0 PRINT=1 METHOD=COND -2LL NUMERICAL LAPLACIAN NOABORT 
;  SIG=3 SIGL=9 FILE=run9.ext MSFO=run9.msf

$COVARIANCE PRINT=E UNCONDITIONAL

$TABLE ID TIME TAFR DV SCALE FLAG MU PMU PHI DVLO DVUP
  NOPRINT NOAPPEND ONEHEADER IDFORMAT=I FILE=run9.fit

$TABLE ID ETAS(1:LAST) OBJI STDID SEX RACE MENOS ECOG DTH CHEMO SURGERY RADIO 
  MHDEP MDANX ERS PGRS AGE LESION1 HEIGHT WEIGHT DTHDY PFSDY
  NOPRINT NOAPPEND FIRSTONLY ONEHEADER IDFORMAT=I FILE=run9.eta
  
; - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 