# install appropriate version of Ensembl perl API
# set PERL5LIB appropriately

#!/usr/bin/env perl

#use strict;
#use POSIX;
#use Bio::EnsEMBL::Registry;
#
#my $registry = 'Bio::EnsEMBL::Registry';
#
#$registry->load_registry_from_db(
#  -host => 'ensembldb.ensembl.org',
#  -user => 'anonymous',
#  -db_version => '66',
#  -verbose => 1,
#  -species => 'homo_sapiens'
#);
#
##my $tr_adaptor    = $registry->get_adaptor( 'Human', 'Core', 'Transcript' );
#
#my $slice_adaptor = $registry->get_adaptor( 'Human', 'Core', 'Slice' );
#
#my $slices = $slice_adaptor->fetch_all('toplevel');
#
#while ( my $sl = shift @{$slices} ) {
#  foreach my $ge ( @{ $sl->get_all_Genes } ) {
#    foreach my $tr ( @{ $ge->get_all_Transcripts } ) {
#      print $sl->chr_name() . "\t" . $tr->id . "\t" . strftime("%Y-%m-%d",localtime($tr->created_date())) . "\n";
#    }
#  }
#}


ms <- read.table("ERR030872_A.mmseq",header=1)
ensd <- read.table("ENSTdate.txt",header=FALSE,sep="\t")
ensd <- ensd[,2:3]

ms = cbind(ms, cdate=droplevels(ensd[match(ms$feature_id,ensd[,1]),2]))


x=as.Date(ms$cdate) < as.Date("2005-01-01")
y=as.Date(ms$cdate) >= as.Date("2005-01-01") & as.Date(ms$cdate) < as.Date("2008-01-01")
z=as.Date(ms$cdate) >= as.Date("2008-01-01") & as.Date(ms$cdate) < as.Date("2012-01-01")
w=as.Date(ms$cdate) >= as.Date("2012-01-01")

barplot(c(sum(ms$log_mu[x] >0)/sum(x), sum(ms$log_mu[y] >0)/sum(y),sum(ms$log_mu[z] >0)/sum(z),sum(ms$log_mu[w] >0)/sum(w)), col=1, names=c("Before 2005", "2005-2007","2008-2011","Since 2012"), ylab="Proportion of transcripts with estimated FPKM > 1", xlab="Creation date of Ensembl transcripts")

