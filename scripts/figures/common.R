`%||%` <- function(a,b) if (is.null(a)) b else a
truth <- function(x) !is.na(x) & toupper(as.character(x)) %in% c("TRUE","T","1")
need <- function(d, columns, label) {
  absent <- setdiff(columns,names(d))
  if(length(absent)) stop(label, ": missing columns: ",paste(absent,collapse=", "),call.=FALSE)
  invisible(d)
}
pick <- function(d, candidates, label) {
  hit <- intersect(candidates,names(d))
  if(!length(hit)) stop(label, ": none of [",paste(candidates,collapse=", "),
    "] found. Actual headers: ",paste(names(d),collapse=", "),call.=FALSE)
  hit[1]
}
read_run <- function(relative, columns=character()) {
  path <- file.path(CONFIG$run_root,relative)
  if(!file.exists(path)) stop("Missing input: ",path,call.=FALSE)
  d <- data.table::fread(path,na.strings=c("NA","NaN",""),check.names=FALSE,
    data.table=FALSE,showProgress=FALSE)
  need(d,columns,relative)
  if("run_id" %in% names(d) && any(!is.na(d$run_id) & d$run_id!=CONFIG$run_id))
    stop("Run identity mismatch in ",relative,call.=FALSE)
  sources[[relative]] <<- data.frame(path=relative, rows=nrow(d), columns=ncol(d))
  d
}
down <- function(name,columns=character()) read_run(paste0("tables/downstream/",name,".tsv"),columns)
audit <- function(name,columns=character()) read_run(paste0("tables/audit/",name,".tsv"),columns)
key <- function(d,cols=c("mirna_id","gene")) {
  need(d,cols,"key"); do.call(paste,c(d[cols],sep="\r"))
}
unique_key <- function(d,cols,label) {
  if(anyDuplicated(key(d,cols))) stop("Duplicate key in ",label,call.=FALSE)
}
pair_label <- function(d) paste(sub("hsa-","",d$mirna_id),d$gene,sep=" / ")
source_table <- function(d,name) {
  data.table::fwrite(d,file.path(OUT,"source_data",paste0(name,".tsv")),sep="\t",na="NA")
  d
}
theme_pub <- function() ggplot2::theme_classic(base_size=CONFIG$base_size,base_family=CONFIG$font)+
  ggplot2::theme(plot.title=ggplot2::element_blank(),plot.subtitle=ggplot2::element_blank(),
    plot.tag=ggplot2::element_text(face="bold",size=12),
    legend.position="bottom",legend.box="vertical",legend.direction="horizontal",
    legend.text=ggplot2::element_text(size=CONFIG$legend_size),
    legend.title=ggplot2::element_text(size=CONFIG$legend_size),
    legend.spacing.y=grid::unit(1,"mm"),
    axis.text=ggplot2::element_text(colour="black"),
    plot.margin=ggplot2::margin(13,10,10,14),
    plot.background=ggplot2::element_rect(fill="white",colour=NA),
    panel.background=ggplot2::element_rect(fill="white",colour=NA),
    legend.background=ggplot2::element_rect(fill="white",colour=NA))
wrap_text <- function(x,width=35) vapply(as.character(x),function(z)
  paste(strwrap(z,width=width),collapse="\n"),character(1))
save_panels <- function(name,plots,legend,layout=NULL) {
  if(!length(plots)) stop("No panels for ",name)
  size <- CONFIG$figure_sizes[[name]] %||% c(CONFIG$width_mm,CONFIG$height_mm)
  layout <- layout %||% CONFIG$layouts[[name]]
  if(is.null(layout)) layout <- cbind(seq_along(plots),seq_along(plots),1,1,1)
  if(!identical(sort(as.integer(layout[,1])),seq_along(plots))) stop("Invalid panel layout: ",name)
  nr <- max(layout[,2]+layout[,4]-1); nc <- max(layout[,3]+layout[,5]-1)
  draw <- function() {
    grid::grid.newpage();grid::pushViewport(grid::viewport(layout=grid::grid.layout(nr,nc)))
    grobs<-lapply(plots,function(p) {
      p<-p+ggplot2::labs(tag=NULL,title=NULL,subtitle=NULL)
      p$theme<-theme_pub()+p$theme
      ggplot2::ggplotGrob(p)
    })
    aligned<-CONFIG$align_panels[[name]]
    if(length(aligned)>1L) {
      if(length(unique(vapply(grobs[aligned],function(g)length(g$heights),integer(1))))!=1L)
        stop("Cannot align paired panels: incompatible grob rows in ",name)
      heights<-do.call(grid::unit.pmax,lapply(grobs[aligned],function(g)g$heights))
      for(j in aligned) grobs[[j]]$heights<-heights
    }
    for(i in seq_along(plots)) {
      a <- layout[layout[,1]==i,,drop=FALSE]
      grid::pushViewport(grid::viewport(layout.pos.row=seq(a[1,2],length.out=a[1,4]),
        layout.pos.col=seq(a[1,3],length.out=a[1,5])))
      grid::grid.draw(grobs[[i]])
      # Draw letters independently of plot guides so oversized legends cannot remove tags.
      if(length(plots)>1L) grid::grid.text(LETTERS[i],x=grid::unit(1.5,"mm"),
        y=grid::unit(1,"npc")-grid::unit(1.5,"mm"),just=c("left","top"),
        gp=grid::gpar(fontface="bold",fontsize=CONFIG$panel_letter_size,fontfamily=CONFIG$font))
      grid::popViewport()
    }
    grid::popViewport()
  }
  pdfpath <- file.path(OUT,"figures",paste0(name,".pdf"))
  if(capabilities("cairo")) grDevices::cairo_pdf(pdfpath,width=size[1]/25.4,height=size[2]/25.4,
    family=CONFIG$font,bg="white") else grDevices::pdf(pdfpath,width=size[1]/25.4,height=size[2]/25.4,bg="white")
  tryCatch(draw(),finally=grDevices::dev.off())
  grDevices::png(file.path(OUT,"figures",paste0(name,".png")),width=size[1],height=size[2],
    units="mm",res=CONFIG$dpi,bg="white")
  tryCatch(draw(),finally=grDevices::dev.off())
  saveRDS(plots,file.path(OUT,"source_data",paste0(name,"_editable_panels.rds")))
  legends[[name]] <<- legend
  invisible(name)
}
primary <- function() {
  p <- read_run("tables/primary/triplets_rho030_CPE.tsv",
    c("mirna_id","gene","qualifies","strict_recurrence","supporting_cancers"))
  p <- p[truth(p$qualifies),,drop=FALSE]
  unique_key(p,c("mirna_id","gene"),"primary set")
  if(nrow(p)!=CONFIG$expected_primary_pairs) stop("Unexpected primary count: ",nrow(p))
  p[order(-p$strict_recurrence,p$mirna_id,p$gene),,drop=FALSE]
}
task <- function(name,fun) {
  message("Building ",name)
  before<-length(table_issues)
  answer <- tryCatch({withCallingHandlers(fun(),warning=function(w) {
    build_warnings[[length(build_warnings)+1L]] <<- data.frame(module=name,message=conditionMessage(w))
    invokeRestart("muffleWarning")
  });list(status=if(length(table_issues)>before) "PARTIAL" else "OK",detail="")},error=function(e)
    list(status="FAILED",detail=conditionMessage(e)))
  build_status[[name]] <<- data.frame(module=name,status=answer$status,detail=answer$detail)
  if(answer$status=="FAILED") message("  FAILED: ",answer$detail)
}
