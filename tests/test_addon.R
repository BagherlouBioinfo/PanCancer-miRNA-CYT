# Run in RStudio: source(file.path(code_root, 'test_addon.R'))
# Synthetic tests only. No discovery data are loaded or changed.
if(!exists('code_root')) code_root <- getwd()
source(file.path(code_root,'cluster_survival.R'))
npass <- 0L
check <- function(ok,label) {if(!isTRUE(ok)) stop('[FAIL] ',label);npass <<- npass+1L;cat('[PASS] ',label,'\n',sep='')}
set.seed(8126)
n <- 1000L
z <- data.table(sample_id=sprintf('S%04d',seq_len(n)),cluster01=rbinom(n,1,.5),
                cyt_z=rnorm(n),age_z=rnorm(n),sex01=rep(0,n))
z[,stage_z:=rnorm(n)]
t <- rexp(n,exp(.8*z$cluster01+.3*z$cyt_z+.2*z$age_z)/1000)
censor <- rexp(n,1/2200)
z[,`:=`(time=pmin(t,censor),event=as.integer(t<=censor))]
a <- cs_compare(z)
check(a$status=='ESTIMATED','synthetic estimable model')
check(a$full$n==a$baseline$n && a$full$n==n,'identical baseline/full sample size')
check(!grepl('sex01',a$covariates),'constant nuisance omitted')
manual0 <- coxph(Surv(time,event)~cyt_z+age_z+stage_z,data=z,ties='efron')
manual1 <- coxph(Surv(time,event)~cyt_z+age_z+stage_z+cluster01,data=z,ties='efron')
manualp <- pchisq(2*(as.numeric(logLik(manual1))-as.numeric(logLik(manual0))),1,lower.tail=FALSE)
check(abs(a$p-manualp)<1e-12,'LRT matches independent direct model')
check(exp(coef(a$full)['cluster01'])>1,'known positive cluster effect recovered')
check(all(c('GLOBAL','cluster01')%in%cs_ph(a$full)$term),'global and cluster-specific PH diagnostics')
check(cs_compare(z[1:20])$status=='NOT_ESTIMABLE','minimum sample guard')
zz <- copy(z);zz[,event:=0L]
check(cs_compare(zz)$status=='NOT_ESTIMABLE','event-per-parameter guard')
zz <- copy(z);zz[,cluster01:=0L]
check(cs_compare(zz)$status=='NOT_ESTIMABLE','single-cluster guard')
zz <- copy(z);zz[cluster01==0,event:=0L]
check(cs_compare(zz)$status=='NOT_ESTIMABLE','zero events in one cluster rejected')
v <- cs_compare(z,tv=TRUE)
check(v$status=='ESTIMATED'&&v$df==2,'time-varying joint 2-df comparison')
check(all(c('cluster01','tt(cluster01)')%in%names(coef(v$full))),'time-varying cluster coefficient exists')
check(all(paste0('tt(',c('cyt_z','age_z','stage_z'),')')%in%names(coef(v$baseline))),
      'same time-varying nuisance adjustment in baseline')
# Preparation rejects silent join loss and duplicate records.
s <- data.table(sample_id=c('a','b','c'),OS=c(1,0,1),OS.time=c(2,3,4),
 age_at_initial_pathologic_diagnosis=c(50,60,70),gender=c('MALE','FEMALE','MALE'),
 ajcc_pathologic_tumor_stage=c('I','','II'),clinical_stage=c('','III',''))
ann <- data.table(sample_id=c('a','b','c'),cancer='AAA')
cy <- data.table(sample_id=c('a','b','c'),CYT=1:3)
mem <- data.table(sample_id=c('c','a','b'),cancer='AAA',cluster=c(2,1,2))
stage_fun <- function(x) match(x,c('I','II','III','IV'))
d <- cs_prepare(s,ann,cy,mem,stage_fun)
check(identical(d$cluster01,c(0L,1L,1L)),'ID matching independent of row order')
check(identical(d$stage_numeric,c(1L,3L,2L)),'AJCC and clinical coalescing')
check(inherits(try(cs_prepare(s,ann,cy,mem[-1],stage_fun),silent=TRUE),'try-error'),'missing sample rejected')
check(inherits(try(cs_prepare(s,ann,cy,rbind(mem,mem[1]),stage_fun),silent=TRUE),'try-error'),'duplicate sample rejected')
mem[1,cancer:='BBB']
check(inherits(try(cs_prepare(s,ann,cy,mem,stage_fun),silent=TRUE),'try-error'),'cancer disagreement rejected')
q <- p.adjust(c(.01,NA_real_,.04),'BH',n=3)
check(abs(q[1]-.03)<1e-12,'non-estimable attempted tests retained in BH denominator')
cat('\n',npass,' checks PASS\n',sep='')
