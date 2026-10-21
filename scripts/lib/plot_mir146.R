# Edit dimensions, palette, font sizes or significance annotations here.
plot_mir146 <- function(output_folder) {
  r<-data.table::fread(file.path(output_folder,'tables/conditional_mir146_associations.tsv'))
  targets<-c('CYT',sort(setdiff(unique(r$outcome),'CYT')))
  palette<-grDevices::colorRampPalette(c('#2166AC','#F7F7F7','#B2182B'))(101)
  draw<-function() {
    old<-par(no.readonly=TRUE);on.exit(par(old))
    layout(matrix(1:4,2,2,byrow=TRUE));par(mar=c(4.2,6,1.8,1),oma=c(2,0,0,0),cex=.85)
    i<-0
    for(setting in c('CPE','LUMP')) for(mi in c('hsa-miR-146a-5p','hsa-miR-146b-5p')) {
      i<-i+1;sel<-which(r$setting==setting & r$mirna_id==mi);z<-r[sel]
      cancers<-sort(unique(z$cancer));nc<-length(cancers);nt<-length(targets)
      plot(NA,xlim=c(.5,nc+.5),ylim=c(.5,nt+.5),xaxs='i',yaxs='i',axes=FALSE,xlab='',ylab='')
      for(j in seq_len(nc)) for(k in seq_len(nt)) {
        cell<-z[cancer==cancers[j]&outcome==targets[k]];value<-cell$partial_rho
        color<-if(length(value)&&is.finite(value)) palette[max(1,min(101,round((value+1)*50)+1))] else '#D9D9D9'
        rect(j-.5,nt-k+.5,j+.5,nt-k+1.5,col=color,border='white')
        if(nrow(cell)&&is.finite(cell$FDR_BH)&&cell$FDR_BH<.05) text(j,nt-k+1,'*',cex=1.3)
      }
      axis(1,at=seq_len(nc),labels=cancers,las=2,cex.axis=.8,tick=FALSE)
      axis(2,at=seq_len(nt),labels=rev(targets),las=2,tick=FALSE)
      box();mtext(LETTERS[i],side=3,line=.4,adj=0,font=2,cex=1.25)
    }
    mtext('Colour: conditional rank correlation (blue negative, white zero, red positive); * BH FDR < 0.05; grey not estimable',side=1,outer=TRUE,line=.3,cex=.7)
  }
  pdf(file.path(output_folder,'figures/Figure_mir146_ABCD.pdf'),width=12,height=8);draw();dev.off()
  png(file.path(output_folder,'figures/Figure_mir146_ABCD.png'),width=12,height=8,units='in',res=600);draw();dev.off()
  writeLines(c('Conditional associations of miR-146a-5p and miR-146b-5p with CYT and shared genes from the fixed primary pair set.',
    'Panels A/B use CPE adjustment for miR-146a-5p/miR-146b-5p, respectively. Panels C/D use LUMP adjustment in the same order.',
    'Each cell shows the correlation between residualised ranks of the outcome and focal miRNA, controlling for the other miRNA and purity. Colours use the same fixed -1 to 1 scale in every panel.',
    'Asterisks indicate BH FDR <0.05 for the focal coefficient HC3 Wald test; significance is not a separate correlation test. Families are separate by purity setting and outcome class, pooling both miRNAs, all eligible cancers and shared genes.',
    'Grey cells were not estimable. These exploratory conditional associations do not establish direct targeting, causality or biological independence. See tables for VIF, effect intervals, sample size and apparent delta R2.'),file.path(output_folder,'FIGURE_LEGEND.txt'))
}
