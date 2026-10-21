# Copy original inputs from an existing project, or explicitly download
# public inputs. Every accepted file must match the ORIGINAL allowlist.
prepare_inputs <- function(project_root='.',from=NULL,download=FALSE) {
  if(is.null(from)&&!isTRUE(download)) stop('Supply from=/old/project or download=TRUE')
  if(!requireNamespace('digest',quietly=TRUE)) stop('Install digest first')
  root <- normalizePath(project_root,winslash='/',mustWork=TRUE)
  manifest <- read.delim(file.path(root,'config/FROM_SCRATCH_ALLOWED_INPUTS.tsv'),stringsAsFactors=FALSE)
  needed <- c('data/EB++AdjustPANCAN_IlluminaHiSeq_RNASeqV2.geneExp.xena',
    'data/pancanMiRs_EBadjOnProtocolPlatformWithoutRepsWithUnCorrectMiRs_08_04_16.xena',
    'data/Survival_SupplementalTable_S1_20171025_xena_sp',
    'reference/targets/Predicted_Targets_Context_Scores.default_predictions.txt.zip',
    'reference/targets/miR_Family_Info.txt.zip','reference/mcpcounter_signatures/genes.txt',
    'reference/mcpcounter_signatures/SIGNATURE_MANIFEST.tsv')
  urls <- c(
    'https://pancanatlas.xenahubs.net/download/EB%2B%2BAdjustPANCAN_IlluminaHiSeq_RNASeqV2.geneExp.xena',
    'https://pancanatlas.xenahubs.net/download/pancanMiRs_EBadjOnProtocolPlatformWithoutRepsWithUnCorrectMiRs_08_04_16.xena',
    'https://pancanatlas.xenahubs.net/download/Survival_SupplementalTable_S1_20171025_xena_sp',
    'https://www.targetscan.org/vert_72/vert_72_data_download/Predicted_Targets_Context_Scores.default_predictions.txt.zip',
    'https://www.targetscan.org/vert_72/vert_72_data_download/miR_Family_Info.txt.zip',
    'https://raw.githubusercontent.com/ebecht/MCPcounter/master/Signatures/genes.txt',NA_character_)
  for(i in seq_along(needed)) {
    rel <- needed[i];dest <- file.path(root,rel)
    sha <- manifest$sha256[match(rel,manifest$relative_path)]
    if(length(sha)!=1L||is.na(sha)) stop('Missing original allowlist entry: ',rel)
    if(file.exists(dest)) {
      if(digest::digest(file=dest,algo='sha256')!=sha) stop('Existing input has the wrong checksum; not overwritten: ',dest)
      cat('[OK] ',rel,'\n',sep='');next
    }
    dir.create(dirname(dest),recursive=TRUE,showWarnings=FALSE)
    tmp <- paste0(dest,'.partial')
    if(!is.null(from)&&file.exists(file.path(from,rel))) {
      if(!file.copy(file.path(from,rel),tmp,overwrite=TRUE)) stop('Copy failed: ',rel)
    } else if(isTRUE(download)&&!is.na(urls[i])) {
      options(timeout=max(3600,getOption('timeout')))
      tryCatch(download.file(urls[i],tmp,mode='wb'),error=function(e)stop('Download failed for ',rel,': ',conditionMessage(e)))
    } else stop('Input unavailable: ',rel,'. Copy the original file into this location.')
    if(digest::digest(file=tmp,algo='sha256')!=sha) {
      unlink(tmp)
      stop('Downloaded/copied file differs from the original input: ',rel,'. No input hash was changed.')
    }
    if(!file.rename(tmp,dest)) stop('Cannot finalise input: ',dest)
    cat('[VERIFIED] ',rel,'\n',sep='')
  }
  invisible(TRUE)
}
if(sys.nframe()==0L) {
  a <- commandArgs(trailingOnly=TRUE)
  f <- sub('^--file=','',grep('^--file=',commandArgs(),value=TRUE)[1])
  root <- dirname(dirname(normalizePath(f,winslash='/',mustWork=TRUE)))
  src <- grep('^--from=',a,value=TRUE)
  prepare_inputs(root,from=if(length(src)) sub('^--from=','',src[1]) else NULL,download='--download'%in%a)
}
