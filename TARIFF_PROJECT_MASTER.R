# TARIFF_PROJECT_MASTER.R -- consolidated local workflow, version 1.1
# SOURCE THIS FILE, THEN RUN: results_folder <- run_project()
# No computation starts just from sourcing. No downloads of the original data.
# Newly generated code; not executed in R by the assistant.
# Primary: 2017 alternative-country HHI, price response, cumulative lags 0-3.
# Other outcomes, horizons, and baseline years are exploratory sensitivities.

TP_VERSION <- "1.1"
tp_root <- function() path.expand("~/tariff_research")
tp_dependencies <- function() {
  needed <- c("haven", "data.table", "fixest", "tidyselect")
  absent <- needed[!vapply(needed, requireNamespace, logical(1), quietly=TRUE)]
  if(length(absent)) stop("Missing packages: ",paste(absent,collapse=", "),
    ". Run the package-installation block before proceeding.")
  suppressPackageStartupMessages(library(data.table))
}
tp_paths <- function(root) {
  root <- path.expand(root)
  list(root=root,
    raw=file.path(root,"extracted","rtp","data","analysis","m_flow_hs10_fm_new.dta"),
    cache=file.path(root,"master_project","cache"),
    runs=file.path(root,"master_project","runs"))
}
tp_status <- function(root=tp_root()) {
  p <- tp_paths(root)
  cat("Master version:",TP_VERSION,"\nData:",p$raw,"\nData exists:",file.exists(p$raw),"\n")
  if(dir.exists(p$runs)) print(list.dirs(p$runs,recursive=FALSE))
  invisible(p)
}
tp_sig <- function(f) {
  x <- file.info(f)
  list(version=TP_VERSION,size=x$size,mtime=as.character(x$mtime))
}
tp_write_rds <- function(x,path) {
  tmp <- paste0(path,".tmp")
  saveRDS(x,tmp,compress=FALSE)
  if(!file.rename(tmp,path)) stop("Cannot finalize cache/output: ",path)
  invisible(path)
}
tp_read <- function(P,S,force=FALSE) {
  if(!file.exists(P$raw)) stop("Data file not found: ",P$raw,
    ". Do not redownload automatically; check the project root.")
  cf <- file.path(P$cache,"selected_native_2016plus_v11.rds")
  mf <- paste0(cf,".meta")
  sig <- tp_sig(P$raw)
  if(!force&&file.exists(cf)&&file.exists(mf)&&identical(readRDS(mf),sig)) {
    S$log("Loading selected-data cache.");return(readRDS(cf))
  }
  cols <- c("id","cty_code","cty_name","hs10","hs8","hs6","year","mdate","naics_str",
    "m_val","m_duty","m_stattariff2","m_apptariff","m_china_hit","m_steel_hit",
    "m_alum_hit","m_solar_hit","m_washer_hit","lm_p","lm_pduty","lm_q1","lm_val")
  S$log("Reading selected columns from the original DTA. The first run creates a cache.")
  d <- as.data.table(haven::read_dta(P$raw,col_select=tidyselect::all_of(cols)))
  d <- d[!is.na(year)&year>=2016&!is.na(cty_code)&cty_code>=0]
  if(anyNA(d$id)||anyNA(d$mdate)||anyNA(d$hs10)) stop("Missing panel keys.")
  if(anyDuplicated(d[,.(id,mdate)])||anyDuplicated(d[,.(hs10,cty_code,mdate)])) stop("Duplicate panel keys; no aggregation was attempted.")
  if(any(d$m_val<0,na.rm=TRUE)) stop("Negative import values require review.")
  if(any(d$m_stattariff2<=-1,na.rm=TRUE)||any(d$m_apptariff<=-1,na.rm=TRUE)) stop("Invalid tariff <= -1.")
  names_cn <- unique(trimws(toupper(as.character(d[cty_code==5700,cty_name]))))
  names_hk <- unique(trimws(toupper(as.character(d[cty_code==5820,cty_name]))))
  if(!"CHINA"%in%names_cn||!"HONG KONG"%in%names_hk) stop("Country code/name verification failed.")
  d[,naics4:=substr(as.character(naics_str),1,4)]
  d[is.na(naics_str)|!nzchar(trimws(naics_str)),naics4:=NA_character_]
  setorder(d,id,mdate)
  tp_write_rds(d,cf);tp_write_rds(sig,mf)
  S$log("Native cache saved; rows=",nrow(d),"; approximate data RAM GiB=",round(as.numeric(object.size(d))/1024^3,2))
  d
}
tp_diff <- function(d,mapping) {
  setorder(d,id,mdate)
  d[,tp_consecutive:=mdate-shift(mdate)==1,by=id]
  for(nm in names(mapping)) {
    old <- unname(mapping[nm])
    d[,(nm):=get(old)-shift(get(old)),by=id]
    d[is.na(tp_consecutive)|!tp_consecutive,(nm):=NA_real_]
  }
  d[,tp_consecutive:=NULL]
  d
}
tp_formula <- function(y,x,fe,iv=NULL) {
  right <- if(length(x)) paste(x,collapse=" + ") else "1"
  as.formula(paste(y,"~",right,"|",fe,if(!is.null(iv))paste("|",iv)else""))
}
tp_subset <- function(d,vars) {
  ok <- complete.cases(d[,..vars])
  for(v in vars) if(is.numeric(d[[v]])) ok <- ok & is.finite(d[[v]])
  d[ok]
}
tp_fit <- function(name,fml,d,S,groups=NULL,diagnostics=FALSE,lead_test=FALSE) {
  S$log("MODEL: ",name)
  warns <- character()
  ans <- tryCatch(withCallingHandlers({
    m <- fixest::feols(fml,data=d,vcov=~hs8+cty_code,ssc=S$ssc,lean=FALSE,mem.clean=TRUE)
    ct <- as.data.table(fixest::coeftable(m),keep.rownames="term")
    setnames(ct,c("term","estimate","standard_error","statistic","p_value"))
    ct[,model:=name];ct[,N:=stats::nobs(m)]
    B <- stats::coef(m);V <- stats::vcov(m)
    omitted <- if(length(m$collin.var)) m$collin.var else character()
    covs <- list(default_two_way=V);eigs <- list();contrasts <- list()
    if(diagnostics) {
      forms <- list(two_way_raw=~hs8+cty_code,two_way_repaired=~hs8+cty_code,
        HS8_only=~hs8,country_only=~cty_code,HS6_only=~hs6)
      for(vn in names(forms)) {
        vv <- tryCatch(stats::vcov(m,vcov=forms[[vn]],ssc=S$ssc,
          vcov_fix=identical(vn,"two_way_repaired")),error=function(e) {
            S$log("VCOV diagnostic failed [",vn,"]: ",conditionMessage(e));NULL
          })
        if(is.null(vv)) next
        covs[[vn]] <- vv
        ev <- eigen((vv+t(vv))/2,symmetric=TRUE,only.values=TRUE)$values
        eigs[[length(eigs)+1L]] <- data.table(model=name,covariance=vn,
          min_eigenvalue=min(ev),max_eigenvalue=max(ev),negative_below_minus_1e10=sum(ev< -1e-10))
      }
    }
    for(gn in names(groups)) {
      g <- groups[[gn]]
      if(!all(names(g)%in%names(B))) {S$log("Contrast not identified: ",gn);next}
      w <- setNames(rep(0,length(B)),names(B));w[names(g)] <- g
      for(vn in names(covs)) {
        vm <- covs[[vn]][names(B),names(B),drop=FALSE]
        est <- sum(w*B);vr <- as.numeric(t(w)%*%vm%*%w)
        se <- if(is.finite(vr)&&vr>0) sqrt(vr) else NA_real_
        contrasts[[length(contrasts)+1L]] <- data.table(model=name,contrast=gn,covariance=vn,
          estimate=est,SE=se,variance=vr,normal95_lower=est-1.96*se,normal95_upper=est+1.96*se,
          normal_p=if(is.finite(se))2*pnorm(-abs(est/se))else NA_real_,N=stats::nobs(m),
          inference="Normal approximation; clustering assumptions and covariance repairs require review")
      }
    }
    own <- intersect(paste0("own",0:6),names(d))
    own_support <- rbindlist(lapply(own,function(v)data.table(model=name,term=v,input_rows=nrow(d),
      finite_nonzero=sum(is.finite(d[[v]])&abs(d[[v]])>1e-12),
      unique_finite=length(unique(d[[v]][is.finite(d[[v]])])))),fill=TRUE)
    lead <- if(lead_test) tryCatch(fixest::wald(m,keep="^lead",print=FALSE),
      error=function(e)list(error=conditionMessage(e))) else NULL
    repair <- NA_real_
    if(all(c("two_way_raw","two_way_repaired")%in%names(covs))) {
      den <- sqrt(sum(covs$two_way_raw^2))
      if(den>0) repair <- sqrt(sum((covs$two_way_raw-covs$two_way_repaired)^2))/den
    }
    list(status="completed",name=name,formula=paste(deparse(fml),collapse=" "),table=ct,
      coefficients=B,covariances=covs,N=stats::nobs(m),omitted=omitted,
      contrasts=rbindlist(contrasts,fill=TRUE),diagnostics=rbindlist(eigs,fill=TRUE),
      own_support=own_support,lead_test=lead,relative_repair_change=repair)
  },warning=function(w) {
    text <- conditionMessage(w);warns <<- c(warns,text)
    S$log("WARNING [",name,"]: ",text);invokeRestart("muffleWarning")
  }),error=function(e) {
    S$log("FAILED [",name,"]: ",conditionMessage(e))
    list(status="failed",name=name,error=conditionMessage(e))
  })
  ans$warnings <- unique(warns)
  target <- file.path(S$models,paste0(gsub("[^A-Za-z0-9]","_",name),".rds"))
  tp_write_rds(ans,target)
  S$manifest[[name]] <- data.table(model=name,status=ans$status,
    warning_count=length(ans$warnings),file=basename(target))
  fwrite(rbindlist(S$manifest,fill=TRUE),file.path(S$run,"model_manifest.csv"))
  S$log("SAVED: ",name,"; status=",ans$status)
  invisible(ans)
}
tp_replication <- function(raw,S) {
  d <- copy(raw[year>=2017,.(id,cty_code,hs10,hs8,hs6,mdate,naics4,m_stattariff2,
    m_apptariff,lm_p,lm_pduty,lm_q1,lm_val,m_val,m_duty)])
  d[,ls:=log1p(m_stattariff2)];d[,la:=log1p(m_apptariff)]
  gap <- d$lm_pduty-d$lm_p;duty <- log1p(d$m_duty/d$m_val)
  check <- function(x,label) {
    ok <- is.finite(x)&is.finite(duty)&d$m_val>0
    data.table(check=label,rows=sum(ok),median_absolute_error=if(any(ok))median(abs(x[ok]-duty[ok]))else NA_real_,
      max_absolute_error=if(any(ok))max(abs(x[ok]-duty[ok]))else NA_real_)
  }
  fwrite(rbindlist(list(check(gap,"price gap versus duties/value"),check(d$la,"applied rate versus duties/value"))),
    file.path(S$run,"accounting_levels.csv"))
  rm(gap,duty,check)
  d <- tp_diff(d,c(dp="lm_p",dpduty="lm_pduty",dq="lm_q1",dv="lm_val",ds="ls",dapp="la"))
  d[,dimplicit:=dpduty-dp]
  d[,c("m_stattariff2","m_apptariff","lm_p","lm_pduty","lm_q1","lm_val","ls","la","m_val","m_duty"):=NULL]
  fe <- "cty_code^mdate + hs10^mdate + cty_code^naics4"
  for(y in c("dp","dpduty","dq","dv")) tp_fit(paste0("REP_RF_",y),tp_formula(y,"ds",fe),d,S)
  tp_fit("REP_FIRST_STAGE_all",tp_formula("dapp","ds",fe),d,S)
  for(y in c("dp","dpduty","dq","dv")) tp_fit(paste0("REP_APPLIED_IV_",y),tp_formula(y,character(),fe,"dapp ~ ds"),d,S)
  tp_fit("REP_INVERSE_SUPPLY_IV",tp_formula("dp",character(),fe,"dq ~ ds"),d,S)
  tp_fit("REP_IMPORT_DEMAND_IV",tp_formula("dq",character(),fe,"dpduty ~ ds"),d,S)
  common <- tp_subset(d,c("dp","dpduty","dimplicit","dapp","ds","cty_code","hs10","hs8","mdate","naics4"))
  rr <- list()
  for(y in c("dp","dpduty","dimplicit","dapp")) rr[[y]] <- tp_fit(paste0("ACCOUNTING_COMMON_",y),tp_formula(y,"ds",fe),common,S)
  if(all(vapply(rr,function(x)identical(x$status,"completed"),logical(1)))) {
    ns <- vapply(rr,function(x)x$N,numeric(1))
    cf <- function(n)unname(rr[[n]]$coefficients["ds"])
    fwrite(data.table(check=c("equal fitted N","paid minus pre minus implicit","applied minus implicit"),
      value=c(as.numeric(length(unique(ns))==1),cf("dpduty")-cf("dp")-cf("dimplicit"),cf("dapp")-cf("dimplicit"))),
      file.path(S$run,"accounting_regressions.csv"))
  }
  invisible(NULL)
}
tp_exposures <- function(raw,year0) {
  if(uniqueN(raw[year==year0,mdate])!=12L)stop("Baseline year lacks 12 observed calendar months: ",year0)
  b <- raw[year==year0&is.finite(m_val),.(value=sum(m_val)),by=.(hs10,cty_code)]
  totals <- b[,.(total=sum(value),china=sum(value[cty_code==5700])),by=hs10]
  alt <- b[!cty_code%in%c(5700,5820)&value>0]
  alt[,share:=value/sum(value),by=hs10]
  hh <- alt[,.(hhi=sum(share^2),alternative_countries=.N),by=hs10]
  e <- merge(totals,hh,by="hs10",all.x=TRUE);e <- e[total>0]
  e[,china_share:=china/total]
  broad <- raw[year>=2017,.(broad=any(m_steel_hit>0|m_alum_hit>0|m_solar_hit>0|m_washer_hit>0,na.rm=TRUE)),by=hs10]
  merge(e,broad,by="hs10",all.x=TRUE)
}
tp_rivals <- function(raw,P,S,year0=2017,force=FALSE) {
  cf <- file.path(P$cache,paste0("rival_baseline",year0,"_v11.rds"));mf <- paste0(cf,".meta")
  sig <- list(source=tp_sig(P$raw),baseline=year0,maxlag=6L)
  if(!force&&file.exists(cf)&&file.exists(mf)&&identical(readRDS(mf),sig)) {
    S$log("Loading rival cache: baseline ",year0);return(readRDS(cf))
  }
  base <- tp_exposures(raw,year0)
  cn <- raw[year>=2017&cty_code==5700,.(hs10,mdate,rate=m_stattariff2)]
  setorder(cn,hs10,mdate)
  cn[,dcn:=log1p(rate)-shift(log1p(rate)),by=hs10]
  cn[,consecutive:=mdate-shift(mdate)==1,by=hs10]
  cn[is.na(consecutive)|!consecutive,dcn:=NA_real_]
  eligible <- base[!is.na(broad)&!broad&is.finite(hhi)&hs10%in%unique(cn$hs10)]
  if(!nrow(eligible))stop("No eligible baseline products.")
  mh <- mean(eligible$hhi);ms <- mean(eligible$china_share)
  hq <- as.numeric(quantile(eligible$hhi,c(.25,.75)))
  eligible[,hhi10:=(hhi-mh)/.1];eligible[,china10:=(china_share-ms)/.1]
  d <- raw[year>=2017&!cty_code%in%c(5700,5820)&hs10%in%eligible$hs10,
    .(id,cty_code,hs10,hs8,hs6,mdate,naics4,lm_p,lm_q1,lm_val,m_stattariff2)]
  d[,lo:=log1p(m_stattariff2)]
  d <- tp_diff(d,c(dp="lm_p",dq="lm_q1",dv="lm_val",down="lo"))
  d <- merge(d,eligible[,.(hs10,hhi10,china10)],by="hs10",all.x=TRUE,sort=FALSE)
  d <- merge(d,cn[,.(hs10,mdate,dcn)],by=c("hs10","mdate"),all.x=TRUE,sort=FALSE)
  missing <- sum(!is.finite(d$dcn));setorder(d,id,mdate)
  for(k in 0:6) {
    cvar <- paste0("cn",k);ovar <- paste0("own",k)
    if(k==0) {d[,(cvar):=dcn];d[,(ovar):=down]} else {
      d[,(cvar):=fifelse(mdate-shift(mdate,k)==k,shift(dcn,k),NA_real_),by=id]
      d[,(ovar):=fifelse(mdate-shift(mdate,k)==k,shift(down,k),NA_real_),by=id]
    }
    d[,(paste0("h",k)):=get(cvar)*hhi10]
    d[,(paste0("s",k)):=get(cvar)*china10]
  }
  for(k in 1:3) {
    v <- paste0("lead",k)
    d[,(v):=fifelse(shift(mdate,k,type="lead")-mdate==k,shift(dcn,k,type="lead"),NA_real_),by=id]
    d[,(paste0("leadh",k)):=get(v)*hhi10];d[,(paste0("leads",k)):=get(v)*china10]
  }
  d[,c("lm_p","lm_q1","lm_val","m_stattariff2","lo","down","dcn"):=NULL]
  result <- list(data=d,base=base,eligible=eligible,mean_hhi=mh,mean_china=ms,hq=hq,missing_cn=missing)
  tp_write_rds(result,cf);tp_write_rds(sig,mf)
  result
}
tp_rhs <- function(L=3L)c(paste0("cn",0:L),paste0("own",0:L),paste0("h",0:L),paste0("s",0:L))
tp_groups <- function(R,L=3L) {
  c <- setNames(rep(1,L+1),paste0("cn",0:L));h <- setNames(rep(1,L+1),paste0("h",0:L))
  list(cumulative_at_baseline_center=c,PRIMARY_HHI_per_0_10=h,
    at_baseline_HHI_p25=base::c(c,h*(R$hq[1]-R$mean_hhi)/.1),
    at_baseline_HHI_p75=base::c(c,h*(R$hq[2]-R$mean_hhi)/.1),p75_minus_p25=h*diff(R$hq)/.1)
}
tp_extension <- function(raw,P,S,extended=TRUE,force=FALSE,diagnostics_only=FALSE) {
  R <- tp_rivals(raw,P,S,2017,force);d <- R$data
  fwrite(R$base,file.path(S$run,"baseline2017_all_exposures.csv"))
  fwrite(R$eligible,file.path(S$run,"baseline2017_eligible_exposures.csv"))
  rhs <- tp_rhs();keys <- c("cty_code","mdate","hs6","hs8","naics4")
  a <- tp_subset(d,c("dp",rhs,keys))
  fwrite(a[,.(rows=.N,products=uniqueN(hs10),countries=uniqueN(cty_code)),by=mdate],file.path(S$run,"extension_primary_calendar.csv"))
  fwrite(data.table(rival_rows_before_restrictions=nrow(d),missing_cn_changes=R$missing_cn,
    input_primary_rows=nrow(a),input_primary_products=uniqueN(a$hs10),
    baseline_hhi_mean=R$mean_hhi,baseline_china_share_mean=R$mean_china,
    hhi_p25=R$hq[1],hhi_p75=R$hq[2]),file.path(S$run,"extension_support.csv"))
  fe <- "cty_code^mdate + hs6^mdate + cty_code^naics4";g <- tp_groups(R)
  tp_fit("EXT_PRIMARY_lag03_price",tp_formula("dp",rhs,fe),a,S,g,TRUE)
  tp_fit("EXT_SENS_pair_FE",tp_formula("dp",rhs,paste(fe,"+ id")),a,S,g,TRUE)
  if(diagnostics_only)return(invisible(NULL))
  tp_fit("EXT_DIAG_contemporaneous_same_input",tp_formula("dp",c("cn0","own0","h0","s0"),fe),a,S,diagnostics=TRUE)
  leads <- c(paste0("lead",1:3),paste0("leadh",1:3),paste0("leads",1:3))
  ap <- tp_subset(a,leads)
  if(nrow(ap)) {
    tp_fit("EXT_LEAD_restricted_without_leads",tp_formula("dp",rhs,fe),ap,S,g,TRUE)
    tp_fit("EXT_LEAD_restricted_with_leads",tp_formula("dp",c(rhs,leads),fe),ap,S,g,TRUE,TRUE)
  }
  if(extended) {
    for(y in c("dq","dv")) {
      ay <- tp_subset(d,c(y,rhs,keys))
      tp_fit(paste0("EXPLORATORY_positive_flow_",y),tp_formula(y,rhs,fe),ay,S,g,TRUE)
      rm(ay)
    }
    rhs6 <- tp_rhs(6);a6 <- tp_subset(d,c("dp",rhs6,keys))
    tp_fit("EXPLORATORY_horizon06_sample_horizon03",tp_formula("dp",rhs,fe),a6,S,g,TRUE)
    tp_fit("EXPLORATORY_horizon06",tp_formula("dp",rhs6,fe),a6,S,tp_groups(R,6),TRUE)
    rm(a6);invisible(gc())
    # Failure of this optional baseline is recorded without deleting earlier outputs.
    tryCatch({
      R16 <- tp_rivals(raw,P,S,2016,force)
      a16 <- tp_subset(R16$data,c("dp",rhs,keys))
      fwrite(R16$eligible,file.path(S$run,"baseline2016_eligible_exposures.csv"))
      tp_fit("EXPLORATORY_baseline2016",tp_formula("dp",rhs,fe),a16,S,tp_groups(R16),TRUE)
    },error=function(e) {
      S$log("Optional 2016 baseline stage failed: ",conditionMessage(e))
      S$stage_errors <- c(S$stage_errors,paste("2016 baseline:",conditionMessage(e)))
    })
  }
  invisible(NULL)
}
tp_report <- function(S) {
  ff <- list.files(S$models,pattern="[.]rds$",full.names=TRUE)
  objs <- lapply(ff,function(f)tryCatch(readRDS(f),error=function(e)list(status="failed",name=basename(f),error=conditionMessage(e))))
  done <- Filter(function(x)identical(x$status,"completed"),objs)
  collect <- function(field,name) {
    x <- lapply(done,function(z)z[[field]]);x <- Filter(function(z)!is.null(z)&&NROW(z)>0,x)
    if(length(x))fwrite(rbindlist(x,fill=TRUE),file.path(S$run,name))
  }
  collect("table","ALL_coefficients.csv");collect("contrasts","ALL_contrasts.csv")
  collect("diagnostics","ALL_covariance_diagnostics.csv");collect("own_support","ALL_own_tariff_support.csv")
  cts <- rbindlist(lapply(done,function(x)x$contrasts),fill=TRUE)
  if(nrow(cts)) {
    scen <- copy(cts)
    scen[,illustrative_tariff_from:=0];scen[,illustrative_tariff_to:=.25]
    scen[,response_ratio_change_pct:=100*expm1(estimate*log1p(.25))]
    scen[,response_ratio_lower_pct:=100*expm1(normal95_lower*log1p(.25))]
    scen[,response_ratio_upper_pct:=100*expm1(normal95_upper*log1p(.25))]
    fwrite(scen,file.path(S$run,"illustrative_25pct_tariff_scenarios.csv"))
    primary <- cts[contrast=="PRIMARY_HHI_per_0_10"&covariance=="default_two_way"]
    if(nrow(primary)&&all(is.finite(primary$normal95_lower))&&all(is.finite(primary$normal95_upper))) {
      grDevices::png(file.path(S$run,"HHI_contrasts_diagnostic.png"),width=1700,height=max(700,nrow(primary)*85),res=150)
      tryCatch({
        par(mar=c(5,22,4,2))
        plot(primary$estimate,seq_len(nrow(primary)),yaxt="n",ylab="",pch=19,
          xlim=range(c(primary$normal95_lower,primary$normal95_upper,0)),
          xlab="Cumulative interaction per 0.10 HHI; normal 95% interval",
          main="Defined specifications; covariance diagnostics still require review")
        axis(2,at=seq_len(nrow(primary)),labels=primary$model,las=2,cex.axis=.65)
        segments(primary$normal95_lower,seq_len(nrow(primary)),primary$normal95_upper,seq_len(nrow(primary)))
        abline(v=0,lty=2)
      },error=function(e)warning("Plot not created: ",conditionMessage(e)),finally=grDevices::dev.off())
    }
  }
  report <- file.path(S$run,"SEND_BACK_MASTER_REPORT.txt")
  con <- file(report,"wt");on.exit(close(con),add=TRUE)
  put <- function(x)writeLines(as.character(x),con)
  put(c("TARIFF PROJECT MASTER REPORT",paste("Version",TP_VERSION),paste("Run",S$run),
    "All primary and exploratory models retained; do not select models by statistical significance.",
    "Limits: unweighted continuing observed flows; country not firm concentration; endogenous baseline structure;",
    "incomplete observed-China tariff schedules; unit values mix prices and composition; no automatic markup/welfare claim.",
    "No missing China tariff changes were set to zero. No full zero-flow universe was manufactured.",
    "No industry PPI/lobbying/Canadian EV/full-universe PPML analysis is implemented without additional data.",
    "Normal contrasts and covariance repairs are diagnostics, not certified causal inference."))
  status <- file.path(S$run,"RUN_STATUS.txt")
  if(file.exists(status))put(readLines(status,warn=FALSE))
  if(length(S$stage_errors)){put("STAGE ERRORS");put(S$stage_errors)}
  man <- file.path(S$run,"model_manifest.csv")
  if(file.exists(man)){put("\nEXECUTION MANIFEST");put(capture.output(print(fread(man))))}
  for(f in c("accounting_levels.csv","accounting_regressions.csv","extension_support.csv")) {
    p <- file.path(S$run,f)
    if(file.exists(p)){put(paste("\n",f));put(capture.output(print(fread(p))))}
  }
  for(x in done) {
    put(paste("\nMODEL:",x$name));put(x$formula);put(capture.output(print(x$table)))
    put(paste("Omitted controls:",paste(x$omitted,collapse=", ")))
    if(length(x$warnings)){put("WARNINGS");put(x$warnings)}
    put(paste("Relative covariance repair change:",x$relative_repair_change))
    if(NROW(x$own_support))put(capture.output(print(x$own_support)))
    if(NROW(x$diagnostics))put(capture.output(print(x$diagnostics)))
    if(NROW(x$contrasts))put(capture.output(print(x$contrasts[covariance%in%c("default_two_way","two_way_raw")])) )
    if(!is.null(x$lead_test)){put("JOINT LEAD TEST");put(capture.output(print(x$lead_test)))}
  }
  for(x in Filter(function(z)identical(z$status,"failed"),objs))put(paste("FAILED MODEL:",x$name,x$error))
  put("\nSESSION");put(capture.output(sessionInfo()))
  invisible(report)
}
run_project <- function(stage=c("all","replication","extension","diagnostics","report"),
  root=tp_root(),extended=TRUE,force_rebuild=FALSE,threads=2) {
  stage <- match.arg(stage);tp_dependencies();P <- tp_paths(root)
  for(p in c(P$cache,P$runs))dir.create(p,recursive=TRUE,showWarnings=FALSE)
  if(stage=="report") {
    rr <- sort(list.dirs(P$runs,recursive=FALSE))
    if(!length(rr))stop("No saved run folders exist.")
    S <- new.env();S$run <- tail(rr,1);S$models <- file.path(S$run,"models");S$stage_errors <- character()
    efile <- file.path(S$run,"stage_errors.txt")
    if(file.exists(efile))S$stage_errors <- readLines(efile,warn=FALSE)
    tp_report(S);cat("Report refreshed:",S$run,"\n");return(invisible(S$run))
  }
  stamp <- format(Sys.time(),"%Y%m%d_%H%M%S")
  run <- file.path(P$runs,paste0(stamp,"_",stage))
  if(dir.exists(run))stop("A run folder already exists for this second. Wait before starting another run.")
  dir.create(run,recursive=TRUE,showWarnings=FALSE)
  S <- new.env();S$run <- run;S$models <- file.path(run,"models")
  S$manifest <- list();S$stage_errors <- character();dir.create(S$models)
  S$log <- function(...) {
    text <- paste0(...,collapse="");cat(text,"\n",sep="")
    cat(text,"\n",file=file.path(run,"execution.log"),append=TRUE,sep="")
  }
  S$ssc <- fixest::ssc(K.adj=TRUE,K.fixef="nonnested",G.adj=TRUE,G.df="min",t.df="min")
  fixest::setFixest_nthreads(threads)
  S$log("START master ",TP_VERSION,"; stage=",stage,"; exploratory sensitivities=",extended)
  writeLines(capture.output(sessionInfo()),file.path(run,"sessionInfo.txt"))
  saveRDS(list(version=TP_VERSION,stage=stage,extended=extended,root=root,threads=threads,
    primary="2017 baseline, country HHI, price, cumulative lags 0-3, observed broad-policy exclusions",
    input=if(file.exists(P$raw))tp_sig(P$raw)else NULL),file.path(run,"run_config.rds"))
  raw <- tryCatch(tp_read(P,S,force_rebuild),error=function(e) {
    S$stage_errors <- c(S$stage_errors,paste("Input preparation:",conditionMessage(e)));NULL
  })
  if(!is.null(raw)) {
    if(stage%in%c("all","replication"))tryCatch(tp_replication(raw,S),error=function(e) {
      S$stage_errors <- c(S$stage_errors,paste("Replication stage:",conditionMessage(e)))
      S$log("REPLICATION STAGE ERROR: ",conditionMessage(e))
    })
    if(stage%in%c("all","extension","diagnostics"))tryCatch(
      tp_extension(raw,P,S,extended,force_rebuild,stage=="diagnostics"),error=function(e) {
        S$stage_errors <- c(S$stage_errors,paste("Extension stage:",conditionMessage(e)))
        S$log("EXTENSION STAGE ERROR: ",conditionMessage(e))
      })
    rm(raw);invisible(gc())
  }
  if(length(S$stage_errors))writeLines(S$stage_errors,file.path(run,"stage_errors.txt"))
  failed <- if(length(S$manifest))sum(vapply(S$manifest,function(x)x$status!="completed",logical(1)))else 0L
  good <- !length(S$stage_errors)&&failed==0L&&length(S$manifest)>0L
  write_status <- function(ok)writeLines(if(ok)
    "Execution finished; inference/causal validity not certified." else
    "Execution incomplete; inspect log and failed models.",file.path(run,"RUN_STATUS.txt"))
  write_status(good)
  report <- tryCatch(tp_report(S),error=function(e) {
    good <<- FALSE;S$log("REPORT ERROR: ",conditionMessage(e));NA_character_
  })
  write_status(good)
  small <- list.files(run,pattern="[.](txt|csv|log|png)$",full.names=TRUE)
  bundle <- file.path(run,"SEND_BACK_RESULTS.zip");old <- getwd()
  zipped <- tryCatch({
    setwd(run);code <- utils::zip(bundle,files=basename(small),flags="-j")
    file.exists(bundle)&&file.info(bundle)$size>0
  },error=function(e)FALSE,finally=setwd(old))
  cat("\nFINISHED. Status:",if(good)"execution completed"else"incomplete; inspect output",
    "\nResults folder:\n",run,"\nReport:\n",report,"\n",sep="")
  if(zipped)cat("Results ZIP:\n",bundle,"\n",sep="")else cat("ZIP not created; share the report, manifest and execution log instead.\n")
  invisible(run)
}
cat("Master script loaded. Run: results_folder <- run_project()\n")
cat("Do not run a second batch while the first one is busy.\n")
