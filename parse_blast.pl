#!/usr/bin/perl -w
#parse a blast report to generate useful summaries
#requires -m 8 generated blast output

use strict;
use Data::Dumper;

die "Usage perl $0 <tab delim blast output> <E value cut off>\n" if @ARGV <1;
my ($blast_output , $eval) = ($ARGV[0], $ARGV[1]);
open(FH, $blast_output) or die "cannot open $blast_output for reading\n";

open(VIS,">visualizeHits.out") or die "cannot open visualizeHits.out for writing\n";
open(SUM,">summary_best_Hits.out") or die "cannot open summary_best_Hits.out for writing\n"; 

print SUM "#FRAGMENT HEADER\tHSP\t#Bp_Aligned\t#Num_Hits\tStart_Coordinates\tStop_Coordinates\n";
print SUM "#SUMMARY HEADER: Fragment\tTotal_Hits\tEval_Array\tID_Array\tTotal_Aligned_bp\tGap_Count\tMax_length_frag\n";

#query	subject	%ID	aln_length	mismatches	gap_openings	query_start	query_stop	sub_start sub_stop e_val bit_score
#1-1     scaffold5114    99.67   19036   39      4       1       19036   259275  278287  0.0     3.725e+04
#1-1     scaffold5114    99.75   2804    5       1       19026   21827   278370  281173  0.0     5497

my $line;
my @fields;
my %query;
my $frag;
my $genomic_region;
my $start;
my $stop;
my $gap;
my $overlap;
my $tmp;
my $sum;
my $max;

# first we have to sort the input file

system("rm blast_report_sorted.out"); # remove any previous instances of file

print "Sorting Input ...\n";

while($line=<FH>){
        @fields=split "\t", $line;
        next if $fields[10] > $eval; 
        if(!(exists($query{$fields[0]}))){
                if(keys %query == 1){
                        open(TMP, ">tmp.txt") or die "cannot open tmp.txt for writing\n";
                        foreach $tmp (keys %query) {
				foreach (@{$query{$tmp}}){
					print TMP $_;
                        	}
			}
                        `cat tmp.txt | sort -n -k7 >> blast_report_sorted.out`;
			%query=();
                	system("rm tmp.txt");
                }
		$query{$fields[0]} = [] unless exists $query{$fields[0]};
                push(@{$query{$fields[0]}}, $line);
                
        }
        else{  
                push(@{$query{$fields[0]}},$line);
        }
}

if(keys %query == 1){
        open(TMP, ">tmp.txt") or die "cannot open tmp.txt for writing\n";
        foreach $tmp (keys %query) {
		foreach (@{$query{$tmp}}){
                	print TMP $_;
        	}
	}
	
        `cat tmp.txt | sort -n -k7 >> blast_report_sorted.out`;
	%query=();
        system("rm tmp.txt");
}

close FH;

print "Done Sorting, Now parsing blast output ...\n";

# open and parse the sorted file
open(FH,"blast_report_sorted.out") or die "cannot open blast_report_sorted.out for reading\n";

while($line=<FH>){
	chomp $line;
	@fields = split "\t", $line;
	next if $fields[10] > $eval;
	

	if(!(exists($query{$fields[0]}))){
		
		if(scalar(keys %query) == 1){
			
			# print out statistics of the key
			# clear the hash


			foreach $frag(keys %query){
				$max = 0;
				$sum = 0;
				foreach $genomic_region(keys %{$query{$frag}}){
					next if $genomic_region=~/ID/ || $genomic_region=~/count/ || $genomic_region=~/Eval/;
					print SUM "$frag\t$genomic_region\t$query{$frag}{$genomic_region}{'aligned_bp'}\t",scalar(@{$query{$frag}{$genomic_region}{'start'}}),"\t","@{$query{$frag}{$genomic_region}{'start'}}\t@{$query{$frag}{$genomic_region}{'stop'}}\n";
					$sum += $query{$frag}{$genomic_region}{'aligned_bp'};
					$tmp = pop @{$query{$frag}{$genomic_region}{'stop'}};
					$max = ($tmp > $max) ? $tmp : $max;
				}
				print SUM "SUMMARY\t",$frag,"\t",$query{$frag}{'count'},"\t@{$query{$frag}{'Eval'}}\t@{$query{$frag}{'ID'}}\t$sum\t$gap\t$max\n";
			}
		%query	= ();
		$gap	= 0;		
	
		}	

		# start of a new query

		print VIS $line,"\n";
		
		$start 		= $fields[6];
		$stop  		= $fields[7];
		$query{$fields[0]}{'count'}++;
		$query{$fields[0]}{$fields[1]}{'aligned_bp'} = $stop - $start + 1;
		
	#	$query{$fields[0]}{$fields[1]}{'mismatches'} = $fields[4];
	#	$query{$fields[0]}{$fields[1]}{'gapOpenings'} = $fields[5];
	
		$query{$fields[0]}{$fields[1]}{'start'} = [] unless exists $query{$fields[0]}{$fields[1]}{'start'};
		push(@{$query{$fields[0]}{$fields[1]}{'start'}}, $start);
		
		$query{$fields[0]}{$fields[1]}{'stop'} = [] unless exists $query{$fields[0]}{$fields[1]}{'stop'};
		push(@{$query{$fields[0]}{$fields[1]}{'stop'}}, $stop);

		$query{$fields[0]}{'Eval'}= [] unless exists $query{$fields[0]}{'Eval'};
		push(@{$query{$fields[0]}{'Eval'}}, $fields[10]);
		$query{$fields[0]}{'ID'}= [] unless exists $query{$fields[0]}{'ID'};
		push(@{$query{$fields[0]}{'ID'}}, $fields[2]);
	}
	else{
		# same scaffold as before
		if(exists($query{$fields[0]}{$fields[1]})){
			if($fields[6] < $stop){
				if($fields[7] > $stop){
					# overlapping hits
					
					print VIS "$line\n";

					$overlap = $stop - $fields[6] + 1;
					$start = $fields[6];
					$stop = $fields[7];
					$query{$fields[0]}{'count'}++;
					$query{$fields[0]}{$fields[1]}{'aligned_bp'} += ($stop - $start + 1) - $overlap;	
					#$query{$fields[0]}{$fields[1]}{'mismatches'} += $fields[4];
					#$query{$fields[0]}{$fields[1]}{'gapOpenings'} += $fields[5];
					
					push(@{$query{$fields[0]}{$fields[1]}{'start'}}, $start);
					push(@{$query{$fields[0]}{$fields[1]}{'stop'}}, $stop);
					push(@{$query{$fields[0]}{'Eval'}}, $fields[10]);
					#$query{$fields[0]}{'ID'}= [] unless exists $query{$fields[0]}{'ID'};
					push(@{$query{$fields[0]}{'ID'}}, $fields[2]);
				}
				else{
					# nested hits
					next;
				}
			}
			else{
				print VIS "$line\n";

				# gap in coverage; would need to relax Eval cut off to look for more hits
				$gap += $fields[6] - $stop;
				#print "Gap!\n$gap\t$line\n"; exit;
				$start = $fields[6];
				$stop = $fields[7];

				$query{$fields[0]}{'count'}++;
				$query{$fields[0]}{$fields[1]}{'aligned_bp'} += ($stop - $start + 1);	
				
				#$query{$fields[0]}{$fields[1]}{'mismatches'} += $fields[4];
				#$query{$fields[0]}{$fields[1]}{'gapOpenings'} += $fields[5];

				push(@{$query{$fields[0]}{$fields[1]}{'start'}}, $start);
				push(@{$query{$fields[0]}{$fields[1]}{'stop'}}, $stop);
			
				push(@{$query{$fields[0]}{'Eval'}}, $fields[10]);
				#$query{$fields[0]}{'ID'}= [] unless exists $query{$fields[0]}{'ID'};
				push(@{$query{$fields[0]}{'ID'}}, $fields[2]);
			}
		}
		else{
			if($fields[6] < $stop){
				if($fields[7] > $stop){

					print VIS "$line\n";
					
					# overlapping hits
					$query{$fields[0]}{$fields[1]}{'aligned_bp'} = 0;
					
				#	$query{$fields[0]}{$fields[1]}{'mismatches'} = $fields[4];
				#	$query{$fields[0]}{$fields[1]}{'gapOpenings'} = $fields[5];
					
					$overlap = $stop - $fields[6] + 1;
					$start = $fields[6];
					$stop = $fields[7];
					$query{$fields[0]}{'count'}++;
					$query{$fields[0]}{$fields[1]}{'aligned_bp'} += ($stop - $start + 1) - $overlap;
					push(@{$query{$fields[0]}{$fields[1]}{'start'}}, $start);
					push(@{$query{$fields[0]}{$fields[1]}{'stop'}}, $stop);
					push(@{$query{$fields[0]}{'Eval'}}, $fields[10]);
					$query{$fields[0]}{'ID'}= [] unless exists $query{$fields[0]}{'ID'};
					push(@{$query{$fields[0]}{'ID'}}, $fields[2]);
				
				}
				else{
					# nested hits
					next;
				}
			}
			else{
				print VIS "$line\n";

				# gap in coverage; would need to relax Eval cut off to look for more hits
				$gap += $fields[6] - $stop;
				$start = $fields[6];
				$stop = $fields[7];
			
				$query{$fields[0]}{'count'}++;
			
				$query{$fields[0]}{$fields[1]}{'aligned_bp'} += ($stop - $start + 1);	

				#$query{$fields[0]}{$fields[1]}{'mismatches'} += $fields[4];
				#$query{$fields[0]}{$fields[1]}{'gapOpenings'} += $fields[5];
				
				push(@{$query{$fields[0]}{$fields[1]}{'start'}}, $start);
				push(@{$query{$fields[0]}{$fields[1]}{'stop'}}, $stop);
			
				push(@{$query{$fields[0]}{'Eval'}}, $fields[10]);
				#$query{$fields[0]}{'ID'}= [] unless exists $query{$fields[0]}{'ID'};
				push(@{$query{$fields[0]}{'ID'}}, $fields[2]);
			}

		}
	}

}

#print Dumper(\%query); 

foreach $frag(keys %query){
	$max = 0;
	$sum = 0;
	foreach $genomic_region(keys %{$query{$frag}}){
		next if $genomic_region=~/ID/ || $genomic_region=~/count/ || $genomic_region=~/Eval/;
		print SUM "$frag\t$genomic_region\t$query{$frag}{$genomic_region}{'aligned_bp'}\t",scalar(@{$query{$frag}{$genomic_region}{'start'}}),"\t","@{$query{$frag}{$genomic_region}{'start'}}\t@{$query{$frag}{$genomic_region}{'stop'}}\n";
		$sum += $query{$frag}{$genomic_region}{'aligned_bp'};
		$tmp = pop @{$query{$frag}{$genomic_region}{'stop'}};
		$max = ($tmp > $max) ? $tmp : $max;
	}
	print SUM "SUMMARY\t",$frag,"\t",$query{$frag}{'count'},"\t@{$query{$frag}{'Eval'}}\t@{$query{$frag}{'ID'}}\t$sum\t$gap\t$max\n";
}
