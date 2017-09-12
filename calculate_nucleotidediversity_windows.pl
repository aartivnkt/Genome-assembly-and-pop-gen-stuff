#!/usr/bin/perl -w
# script to calculate nucleotide diversity from an alignment, in overlapping windows of 200

use strict;
die "Usage $0: <diversity-divergence file>\t<group name>\t<window size>\n" if @ARGV <3;
#my $in_file_path="/mnt/lustre/home/aartiv/ABO/Diversity_Divergence/";

my($file, $group, $window) = ($ARGV[0],$ARGV[1], $ARGV[2]);
my ($gene, $exon, $i, $j, $line, $key, $gaps, $mean_pi, $mean_divergence, $ratio, $pi, $divergence);
my (@lines, @fields);
my %gene_hash;
my ($index, $chr_num, $genemodel,$exon_num, $chr_start);

open(FH,$file) or die "cannot open $file for reading\n";
#open(FH,$file) or "die cannot open file for reading\n";

print "reading the lines of the file ...\n";
print "hashing lines for $file\n";

#my $out_file_path = "/mnt/lustre/home/aartiv/ABO/WindowAnalysis";
my $outfile = $group.$window."_Windows.out";

print "Processing $outfile\n";
open(OUT,">$outfile") or die "cannot open $outfile for writing\n";
print OUT "##Window#\tGene_name\tExon#\tWindow_Start-Window_End\tMean_Pi\tMean_Divergence\tPi/Div\n";

while($line=<FH>){
	next if $line=~/^Chromosome/ || $line=~/^#/;
	chomp $line;
	($gene,$exon)=(split "\t",$line)[0,1];
	$key = $gene."_"."$exon";
	if(!exists($gene_hash{$key})){
		#print "$key does not exist in the hash\n";
		
		foreach $key( keys %gene_hash){
			#print "checking if there are any other keys\n";
			if($key ne ""){
				#print "processing $key, num elems for this key is\t",scalar(@{$gene_hash{$key}}),"\n";
        			for($i=0;$i<=@{$gene_hash{$key}}-$window;$i++){
                			#print $gene_hash{$key}->[$i],"\n\n";
                			#print @{$gene_hash{$key}},"\n\n";exit;
                			for($j=$i;$j<$i+$window;$j++){
                        			if($j==$i){
                                			#print $lines[$j],"\n"; 
                                			@fields=split("\t", $gene_hash{$key}->[$j]);
                                			$pi = 0;
                                			$divergence = 0;
                                			$gaps = 0;
                        
                                			if ($fields[4]=~/[^ATGC]/i || $fields[5]=~/[^ATGC]/i){
                                        			$gaps++;
                                			}
                                			else{
                                        			$pi += $fields[8];
                                        			$divergence += $fields[6];
                                			}
                                			$index = $i;
                                			$genemodel = $fields[0];
                                			$exon_num = $fields[1];
                                			$chr_num = $fields[2];
                                			$chr_start = $fields[3];
						}
	
                        			elsif($j > 1 && $j < (($window+$i)-1) ){
                                			@fields=split "\t",$gene_hash{$key}->[$j];
                        			        if ($fields[4]=~/[^ATGC]/i || $fields[5]=~/[^ATGC]/i){
                                        			$gaps++;
                                			}
                                			else{
                                        			$pi += $fields[8];
                                        			$divergence += $fields[6];
                                			}
                        			}
                        			elsif($j == (($window+$i)-1) ){
                                			@fields=split "\t",$gene_hash{$key}->[$j];
                                			if($fields[4]=~/[^ATGC]/i || $fields[5]=~/[^ATGC]/i){
                                        			$gaps++;
                                			}
                                			else{
                                        			$pi += $fields[8];
                                        			$divergence += $fields[6];
                               				}
						 	if( $gaps / $window <=0.1){ # changed gap cut off to atmost 10% -- 20th feb 2012 
                                        			$mean_pi = $pi / ($window - $gaps);
                                        			$mean_divergence = $divergence / ($window - $gaps);
                                        			$ratio = ($mean_divergence > 0 ) ? ($mean_pi / $mean_divergence) : 0;
                                        			#print "Window$index\t$genemodel\t$exon_num\t$chr_num:$chr_start-\t$fields[2]:$fields[3]\t", sprintf("%.4f",$mean_pi),"\t",sprintf("%.4f",$mean_divergence),"\t",sprintf("%.4f",$ratio),"\n";
                                        			print OUT "Window$index\t$genemodel\t$exon_num\t$chr_num:$chr_start-\t$fields[2]:$fields[3]\t", sprintf("%.4f",$mean_pi),"\t",sprintf("%.4f",$mean_divergence),"\t",sprintf("%.4f",$ratio),"\n";
                                			}
                                			else{
                                        			print "Failed Window$i\t$j\t$gaps\t",($gaps/$window),"\n";
                                			}
                                			$j = $i + 1;
                                			last;
                				}

        				}
        			}
			}
			#print "empyting hash for $key\n";
			%gene_hash=();
		}
	#print "inserting $key into hash\n";
	$gene_hash{$key} = [];
	#print "pushing $line into hash\n";
	push(@{$gene_hash{$key}},$line);
	}
	else{
		push(@{$gene_hash{$key}}, $line);	
	}
}
	foreach $key( keys %gene_hash){
                        #print "checking if there are any other keys\n";
                        if($key ne ""){
                                #print "processing $key, num elems for this key is\t",scalar(@{$gene_hash{$key}}),"\n";
                                for($i=0;$i<=@{$gene_hash{$key}}-$window;$i++){
                                        #print $gene_hash{$key}->[$i],"\n\n";
                                        #print @{$gene_hash{$key}},"\n\n";exit;
                                        for($j=$i;$j<$i+$window;$j++){
                                                if($j==$i){
                                                        #print $lines[$j],"\n"; 
                                                        @fields=split("\t", $gene_hash{$key}->[$j]);
                                                        $pi = 0;
                                                        $divergence = 0;
                                                        $gaps = 0;

                                                        if ($fields[4]=~/[^ATGC]/i || $fields[5]=~/[^ATGC]/i){
                                                                $gaps++;
                                                        }
                                                        else{
                                                                $pi += $fields[8];
                                                                $divergence += $fields[6];
                                                        }
                                                        $index = $i;
                                                        $genemodel = $fields[0];
                                                        $exon_num = $fields[1];
                                                        $chr_num = $fields[2];
                                                        $chr_start = $fields[3];
                                                }

                                                elsif($j > 1 && $j < (($window+$i)-1) ){
                                                        @fields=split "\t",$gene_hash{$key}->[$j];
                                                        if ($fields[4]=~/[^ATGC]/i || $fields[5]=~/[^ATGC]/i){
                                                                $gaps++;
                                                        }
                                                        else{
                                                                $pi += $fields[8];
                                                                $divergence += $fields[6];
                                                        }
                                                }
                                                elsif($j == (($window+$i)-1) ){
                                                        @fields=split "\t",$gene_hash{$key}->[$j];
                                                        if($fields[4]=~/[^ATGC]/i || $fields[5]=~/[^ATGC]/i){
                                                                $gaps++;                                                        
							}
                                                        else{                                                                
								$pi += $fields[8];
                                                                $divergence += $fields[6];
                                                        }
                                                        if( $gaps / $window <=0.1){ # changed gap cut off as of 20th feb 2012
                                                                $mean_pi = $pi / ($window - $gaps);
                                                                $mean_divergence = $divergence / ($window - $gaps);
                                                                $ratio = ($mean_divergence > 0 ) ? ($mean_pi / $mean_divergence) : 0;
                                                        	#print "Window$index\t$genemodel\t$exon_num\t$chr_num:$chr_start-\t$fields[2]:$fields[3]\t", sprintf("%.4f",$mean_pi),"\t",sprintf("%.4f",$mean_divergence),"\t",sprintf("%.4f",$ratio),"\n";
                                                                print OUT "Window$index\t$genemodel\t$exon_num\t$chr_num:$chr_start-\t$fields[2]:$fields[3]\t", sprintf("%.4f",$mean_pi),"\t",sprintf("%.4f",$mean_divergence),"\t",sprintf
("%.4f",$ratio),"\n";
                                                        }
                                                        else{
                                                                print "Failed Window$i\t$j\t$gaps\t",($gaps/$window),"\n";
                                                        }
                                                        $j = $i + 1;
                                                        last;
                                                }

                                        }
                                }
                        }
	}

