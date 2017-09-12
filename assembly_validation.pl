#!/usr/bin/perl -w
use strict;
use Getopt::Long();
use Data::Dumper;

# Assembly validation -- assesses completeness and accurateness of a genome assembly

######## USAGE ##########

usage() if @ARGV<1;

sub usage{
        print STDERR (
	    	"\nUsage: perl $0 options [arg]\n\n"
		. "Where options and args are:\n"
		. "--asm_dir [dir name]\n" 
		. "--validation_dir [dir name]\n"
		. "--file_list [fastq file list]\n"
		. "--index [0||1]\n"
		. "--align [0||1]\n" 
		. "--sampe [0||1]\n"
		. "--global_dist_stats [0||1]\n"
		. "--local_misassembly_stats [0||1]\n"
	       	."Specify dir name from the Lemur base dir\n\n "
        );
exit();
}

my $asm_dir 			= "";
my $validation_dir		= "";
my $file_list 			= "";
my $index 			= 0;
my $align 			= 0;
my $sampe	 		= 0;
my $global_dist_stats 		= 0;
my $local_misassembly_stats 	= 0;

Getopt::Long::GetOptions(
    'asm_dir=s'   		=> \$asm_dir,
    'validation_dir=s'		=> \$validation_dir,
    'file_list=s' 		=> \$file_list,
    'index=i'     		=> \$index,
    'align=i'     		=> \$align,
    'sampe=i'     		=> \$sampe,
    'global_dist_stats=i'	=> \$global_dist_stats,
    'local_misassembly_stats'	=> \$local_misassembly_stats,
) or die "Incorrect usage!\n", usage();

######################################################

my ($command, $scaffold_seq);

if($index){

	die "Specify assembly dir name\n" if $asm_dir eq "";

	chomp $asm_dir if $asm_dir ne "";

	$scaffold_seq = `ls $asm_dir/*scafSeq`;
	chomp $scaffold_seq;
	die "cannot find scaffold sequence in $asm_dir\n" if $scaffold_seq eq "";

	#mkdir $validation_dir unless -d $validation_dir or die "cannot create $validation_dir \n";

	# index the genome using bwa index

	chdir($validation_dir); # output is getting written to assembly dir though
				# instead of validation dir.
	
	$command = "time bwa index -a bwtsw ../$scaffold_seq";
	print "Indexing genome using $command ... \n";
	#system($command);
	print "Done indexing\n";
	
	chdir("../");
}


#perl -e '@file=</bigmemscratch/aartiv/Insert_Sorted_Sequences/*gz>; push(@files,<~/Lemur/TrimmedReads_NoQuake/Quake_upto1Kb_C3/*single*>); push(@files, <~/Lemur/TrimmedReads_NoQuake/Quake_upto1Kb_C3/*176*cor.fastq.gz>); push(@files,<~/Lemur/TrimmedReads_NoQuake/Quake_upto1Kb_C3/*1kb*cor.fastq.gz>);push(@files, <~/Lemur/TrimmedReads_NoQuake/Quake_upto1Kb_C3/*510_6*gz>);foreach $f(@files){ print $f,"\n";}' | wc -l

# align reads

my @seqfiles = `ls /scratch/aartiv/Lemur/Lemur_SNP/Insert_Sorted_Sequences/*gz`; #check path for new jobs
push(@seqfiles, </scratch/aartiv/Lemur/Lemur_SNP/SOAP_NewEC_CorrectedReads_DelErr/SOAP_NewEC_CorrectedReads_DelErr_LibSort/*innies.gz>);
push(@seqfiles, </scratch/aartiv/Lemur/Lemur_SNP/SOAP_NewEC_CorrectedReads_DelErr/SOAP_NewEC_CorrectedReads_DelErr_LibSort/Lemur_176*gzsort.gz>);
#push(@seqfiles,<~/Lemur/TrimmedReads_NoQuake/Quake_upto1Kb_C3/*1kb*cor.fastq.gz>);
#push(@seqfiles, <~/Lemur/TrimmedReads_NoQuake/Quake_upto1Kb_C3/*510_6*gz>);

my ($max_jobs, $i) = (60, 0);
my ($jobs , $outfile);
my @fields;

if($align){
	while(1){
        	$jobs=`qstat | grep -c "aln"`;
        	chomp($jobs);
        	if($jobs < $max_jobs){
                	print "$jobs currently running\n";
			chomp($seqfiles[$i]);
                	@fields = split /\//,$seqfiles[$i];
			$outfile = (split /\.gz/, $fields[$#fields])[0];
			print "submitting aln job for $seqfiles[$i], outfile is $validation_dir/$outfile\n"; 
                	`qsub /mnt/lustre/home/aartiv/Lemur/Scripts/aln.sh $validation_dir/$outfile $seqfiles[$i]`;
                	$i++;
        	}
        	else {
                	sleep 3;
                	$jobs=`qstat | grep -c "aln"`;
                	print "$jobs currently running\n";
        	}
        	last if $i==@seqfiles;
	}
	print "all jobs submitted\n";

	while(1){
		$jobs = `qstat | grep -c "aln"`;
		chomp($jobs);
		if($jobs == 0){
			last;
		}
		else{
			sleep 5;
		}
	}
}

# bwa sampe

my ($line , $file , $pair1, $pair2 , $sai_file1 , $sai_file2 , $insert_size , $out_prefix );
$command = "";

if($sampe ) {
	print "aligning PE reads\n";
	die "Please specify file_list name\n" if $file_list eq "";
	
	open(FH,$file_list) or die "cannot find $file_list in this dir\n";
	while($line = <FH>){
		next if $line=~/^#/;
		chomp $line;
		($pair1, $pair2 , $insert_size , $out_prefix) = (split "\t", $line)[0,1,2,3];
		
		$sai_file1 = (split /\.gz/, $pair1)[0];
		$sai_file2 = (split /\.gz/, $pair2)[0];
		
		foreach $file (`ls $validation_dir`){
			if($file =~/$sai_file1/){
				chomp $file;
				$sai_file1 = $file;
			}
			elsif($file =~ /$sai_file2/){
				chomp $file;
				$sai_file2 = $file;
			}
		}
		print " SUBMITTING THE FOLLOWING JOB : /mnt/lustre/home/aartiv/Lemur/Scripts/align_PE.sh $insert_size $validation_dir $out_prefix $sai_file1 $sai_file2 /scratch/aartiv/Lemur/Lemur_SNP/SOAP_NewEC_CorrectedReads_DelErr/SOAP_NewEC_CorrectedReads_DelErr_LibSort/$pair1 /scratch/aartiv/Lemur/Lemur_SNP/SOAP_NewEC_CorrectedReads_DelErr/SOAP_NewEC_CorrectedReads_DelErr_LibSort/$pair2";
		`qsub /mnt/lustre/home/aartiv/Lemur/Scripts/align_PE.sh $insert_size $validation_dir $out_prefix $sai_file1 $sai_file2 /scratch/aartiv/Lemur/Lemur_SNP/SOAP_NewEC_CorrectedReads_DelErr/SOAP_NewEC_CorrectedReads_DelErr_LibSort/$pair1 /scratch/aartiv/Lemur/Lemur_SNP/SOAP_NewEC_CorrectedReads_DelErr/SOAP_NewEC_CorrectedReads_DelErr_LibSort/$pair2`;
		print "\n\n";
	}
	
	 while(1){
                $jobs = `qstat | grep -c "Align_PE"`;
                chomp($jobs);
                if($jobs == 0){
                        last;
                }
                else{
                        sleep 5;
			$jobs = `qstat | grep -c "Align_PE"`;
			print "$jobs currently running\n";
                }
        }

	print "All jobs submitted!\n";

	# merge and index the sorted bams

	print "Merging and indexing the sorted bams...\n";
        
	# QC on assembly
	$command = "sed '/^\$/d' $scaffold_seq > $asm_dir.no_blnk_lines.fasta";
	#system($command);

	$command = "samtools merge -rh $validation_dir/Lemur.rg.txt $validation_dir/Lemur.allLibs.merged.bam $validation_dir/*sort.bam";
	#system($command);
        
	$command = "samtools index $validation_dir/Lemur.allLibs.merged.bam";
	#system($command);

	$command = "samtools view $validation_dir/*merged.bam > $validation_dir/Lemur.allLibs.merged.sam";
	#system($command);
	$command = "gzip $validation_dir/*merged.sam";
	#system($command);	
	
	# run samtools mpileup to get per base covg for PE reads
	$command = "samtools mpileup -f $asm_dir.no_blnk_lines.fasta -D -q 10 -Q 20 $validation_dir/*merged.bam > $validation_dir/Coverage.mpileup";
	#system($command);

	close FH;
}

# variables for global dist; will come handy for local too.

my (%lib , %pos , %D_Of_OlpR );
my ($count , $total_size , $key1 , $key2 , $lib_name , $mean , $stddev , $ssd , $date , $median);	
my (@coverage, @sorted_pos ,@sorted_coverage, @cigar , @values);
my ($d , $A_stat);

if($global_dist_stats){
	print "Validating assembly based on mate pair distances\n";
	
	# To interpret bitwise flags from *merged.sam.gz file

	my %flag_set=(
                '1'  , 'Read has mate reads',
                '2'  , 'All mate reads mapped',
                '4'  , 'Read is unmapped',
                '8'   , 'Next mate read unmapped',
                '16'  , 'Read on reverse strand',
                '32'  , 'Next mate read on rev str',
                '64'   , 'First of all mate reads',
                '128'  , 'Last of all mate reads',
                '256'  , 'Secondary alignment',
                '512'  , 'Read fails quality checks',
                '1024' , 'Read is PCR/opt duplicate',
	);

	my %flag_unset=(
                '1'  , 'Read has NO mate reads',
                '2'  , 'Not all mate reads mapped',
                '4'  , 'Read is mapped',
                '8'   , 'Next mate read is mapped',
                '16'  , 'Read on forward strand',
                '32'  , 'Next mate read on forw str',
                '64'   , 'Not first of all mate reads',
                '128'  , 'Not last of all mate reads',
                '256'  , 'Primary alignment',
                '512'  , 'Read passes quality checks',
                '1024' , 'Read is no duplicate',
	);
	
	push @values,1;

	for (my $i=2; $i<513; $i*=2) {
        	push @values,$i;
	}

	# first get global insert size distribution for each library based on properly mapped pairs
	
	open(FH, "zcat $validation_dir/*merged.sam.gz |") or die "cannot open $validation_dir/*merged.sam.gz for reading\n";
	#open(FH, "$validation_dir/test.sam") or die "cannot open $validation_dir/test.sam for reading\n";
	
	$date = `date`;
	chomp($date);

	print "parsing $validation_dir/*merged.sam.gz file...\nStart Time is $date\n";
	
	print "Hashing insert sizes for each library, as well as overlapping reads(innies) ...\n";
	while($line=<FH>){
		chomp $line;
		@fields= split "\t",$line;
		
		# only look at scaffolds.. not singletons

		if(($fields[6] eq "=") && ($fields[2]=~/^scaffold/) && (length($fields[5]) >= 3) && (!($fields[1] & 4)) && (!($fields[1] & 8)) && ($fields[1] & 2) && ($fields[4] >= 10)){		
		
		# also trying to only consider reads which match perfectly atleast over the avg read length = 79 (avg read length of the soap error corr reads)
	
				$lib_name = (split /:/, $fields[$#fields])[2];
				$lib{$lib_name}{$fields[8]}++ if $fields[8] > 0;
				
				# additional filters such as below

				@cigar = split "", $fields[5];
				if($lib_name=~/Aarti_3kb/ || $lib_name=~/Aarti_8kb/){
					if($cigar[0].$cigar[1] >=28){
						$lib{$lib_name}{$fields[8]}++ if $fields[8] > 0;
					}
				}
				else{
					if($cigar[0].$cigar[1] >= 79){
						$lib{$lib_name}{$fields[8]}++ if $fields[8] > 0;
					}
				}
				## DEBUG DOC FOR OVERLAPPING READS FOR 510 & 1KB INSERTS, GLOBAL DIST. ########
				
				if ($fields[8] <= 172 && $fields[8] > 0 && $lib_name ne "Lemur_176.PE.sort"){
					$D_Of_OlpR{$fields[2]}{$fields[3]}++;
				}
				
				#################################################


			}
		}
	}
	close FH;

	print "Done\n";
	
	# parse the mpileup file here to get median coverage for PE reads and set cut off to identify repetitive reads	
	# below to be executed once mpileup finishes ...
	
	# print "Calculating Median coverage for properly mapped PE reads (excluding gaps) ...\n";

	open(FH,"$validation_dir/Coverage.mpileup") or die "cannot open $validation_dir/Coverage.mpileup for reading\n";

	while($line=<FH>){
        	chomp $line;
        	@fields = split "\t", $line;
        	
		if($fields[0]=~/^scaffold/ && $fields[2] =~/[ATGC]/i){
               		push(@coverage, $fields[3]);
        	}
	}
	@sorted_coverage = sort {$a <=> $b} @coverage;
	@coverage=();
	$median = $sorted_coverage[int(0.5*$#sorted_coverage)];
	
	print "Done...Median Coverage is $median\n";
	
	close FH;

	print "Reporting innies with coverage greater than 2 fold the median covg ... \n";

	open(FHOUT,">$validation_dir/D_Of_OlpR.out") or die "cannot open $validation_dir/D_Of_OlpR.out for writing\n";

	foreach $key1 (keys %D_Of_OlpR){
		foreach $key2 (keys %{$D_Of_OlpR{$key1}}){
			print FHOUT "$key1\t$key2\t$D_Of_OlpR{$key1}{$key2}\n" if $D_Of_OlpR{$key1}{$key2} > 2*$median;
		}
	}

	close FHOUT;
	%D_Of_OlpR=();

	print "Done\n";

	# sorting through the overlap hash

	print "Sorting through the innies to report suspicious assemblies ...\n";

	open(FH,"$validation_dir/D_Of_OlpR.out") or die "cannot open D_Of_OlpR.out for reading\n";
	open(FHOUT,"$validation_dir/Global_Dist.Suspicious.CollapsedRepeats.out") or die "cannot open $validation_dir/Global_Dist.Suspicious.CollapsedRepeats.out for writing\n";

	while($line=<FH>){
        	next if $line !~/^scaffold/;
        	chomp $line;
        	@fields = split "\t", $line;
        	if(!exists($pos{$fields[0]})){
                	
			# start of new scaffold or end of prev scaffold info
                	foreach $key1 (keys %pos){
                        	@sorted_pos=();
				$A_stat = 0;
                        		if($key1 ne ""){
                                		foreach(@{$pos{$key1}}){
                                        		push(@sorted_pos, (split "=>",$_)[0]);
                                		}
                                		@sorted_pos = sort {$a <=> $b} @sorted_pos;
                                		print "$key\t@sorted_pos\n";

                                		$d = ( ($sorted_pos[$#sorted_pos] - $sorted_pos[0]) / @sorted_pos);

                                		if( $d  > 0 && $d <= 5){
                                        		@sorted_pos = ();
                                        		@sorted_pos = sort {$a cmp $b} @{$pos{$key1}};
                                        		
							#############################
							
							print FHOUT "Suspicious $key1:\tA statistic\t$A_stat\t@sorted_pos\n";
                                		}
                        		}
                	}
                	%pos=();
                	$pos{$fields[0]}=[];
                	push(@{$pos{$fields[0]}}, "$fields[1]\=>$fields[2]");
        	}
        	push(@{$pos{$fields[0]}}, "$fields[1]\=>$fields[2]");
	}

	close FH;
	close FHOUT;
	%pos=();

	print "Done\n";
	################################################################
	
	$date = `date`;
	chomp($date);
	
	print "Done parsing\nEnd Time is $date\n";

	print "Calculating Mean and Stddev from global distribution \n";
	
	foreach $key1 (keys %lib){
		open(FHOUT,">$validation_dir/$key1.insert_dist.out") or die "cannot open $validation_dir/$key1.insert_dist.out for writing\n";
		
		$count 		= 0;
		$total_size 	= 0;
		#print "$key1:\n";
		foreach $key2 (keys %{$lib{$key1}}){
			$count+= $lib{$key1}{$key2};
			$total_size += $key2* $lib{$key1}{$key2};
			print FHOUT "$key2\t$lib{$key1}{$key2}\n";
		}
		close FHOUT;
		$mean = $total_size / $count;
			
		#stddev
		$ssd = 0;

		open(FH,"$validation_dir/$key1.insert_dist.out") or die "cannot open $validation_dir/$key1.insert_dist.out for reading\n";
		while($line=<FH>){
			chomp $line;
			$ssd += (((split "\t", $line)[0] - $mean)**2)* (split "\t", $line)[1];
		}
		close FH;

		$stddev = sqrt($ssd / ($count -1));

		# write global stats to a file for use with CE/Z-scores later

		open(APPEND,">>$validation_dir/Global_Distribution_stats_NoSingletons.out") or die "cannot open $validation_dir/Global_Distribution_stats_NoSingletons.out for appending\n";
		print APPEND "$key1\tMean Insert Size\t$mean\tStd dev\t$stddev\n";
		print "$key1\tMean Insert Size\t$mean\tStd dev\t$stddev\n";
		close APPEND;

		print "Done\n";
		
		print "Plotting $key1 ... \n";
		# plot dist
		open(FHOUT,">$validation_dir/plot.R") or die "cannot open $validation_dir/plot.R script for writing\n";
		print FHOUT "insert_size <- read.table(file=\"$validation_dir/$key1.insert_dist.out\", header = FALSE)\n";
		print FHOUT "names(insert_size)<-c(\"distance\", \"frequency\")\n";
		print FHOUT "pdf(\"$validation_dir\/$key1.Rplot.NoSingletons.pdf\")\n";
		print FHOUT "plot(insert_size\$distance, insert_size\$frequency, type=\"h\" , xlab=\"insert size\", ylab=\"frequency\", main = \"Distribution of insert size for $key1\", sub= \"Mean = $mean Stddev = $stddev\")\n";
		#print FHOUT "lines(lowess(insert_size\$distance, insert_size\$frequency), col=\"red\")\n";	
		print FHOUT "dev.off()\n";
		
		close FHOUT;

		$command = "R CMD BATCH $validation_dir/plot.R";
		system($command);	
	}
}
my (%lib_stats , %seen , %scaffold);
my ($paircount, $single, $num_keys , $scaff_num);

if($local_misassembly_stats){
	
	print "Testing local assemblies for quality...\n";
	$date=`date`;
	print "Start time is $date";	
	
	# parse global distributions
	open(FH,"$validation_dir/Global_Distribution_stats_NoSingletons.out") or die "cannot open $validation_dir/Global_Distribution_stats_NoSingletons.out for writing\n";
	
	while($line=<FH>){
		chomp $line;
		@fields=split "\t", $line;
		$lib_stats{$fields[0]}{"Mean"}=$fields[2];
		$lib_stats{$fields[0]}{"Stddev"}=$fields[4];
	}
	close FH;

	open(FH, "zcat $validation_dir/*merged.sam.gz |") or die "cannot open $validation_dir/*merged.sam.gz for reading\n";
	#open(FHOUT, ">>$validation_dir/GMB_Analysis_localClones.out") or die "cannot open $validation_dir/GMB_Analysis_localClones.out for writing\n";

	#print FHOUT "##Note: Bad Values\tAvg_Zscore = 999 indicates FR_count = 0; RGB = 999 indicates MP_bad = 0 or MP_good = 0\n"; 	

	while($line=<FH>){
		chomp $line;
		@fields = split "\t", $line;
		$scaff_num = (split "d", $fields[2])[1];	
		print "processing $scaff_num ...\n";
		next if $scaff_num < 67815; # for testing

		if(!(exists $seen{$fields[2]})){
			#print "$fields[2] is new\n";
			# parse out the previous scaffolds
			foreach $key1 (keys %seen){
				if($key1 ne ""){
					#print FHOUT "$key1\t";
					print "$key1\t";
					#print "Printing Data Dumper\n";
					sleep 5;
					#print Dumper(%scaffold);

					########## get clone measures #################
						
					get_clone_measures(\%scaffold, \%lib_stats);

					###############################################
				}
			}
			#print "emptying scaffold and seen hash if they have values\n";
			%scaffold=();
			%seen=();
			$seen{$fields[2]}++;

			# pop scaffold hash;
			#print "populating scaffold hash\n";
			if($fields[0]=~/#/){
				$key1 = (split /#/,$fields[0])[0];
				if($fields[4] >= 10){
					$scaffold{$key1}{$fields[0]}=[];
					push(@{$scaffold{$key1}{$fields[0]}}, $fields[1], $fields[2],$fields[3],$fields[6],$fields[8],(split /:/, $fields[$#fields])[2]);
				}
			}
			else{
				if($fields[4] >= 10){
					$scaffold{$fields[0]}{$fields[0]}=[] if(!exists($scaffold{$fields[0]}));
					push(@{$scaffold{$fields[0]}{$fields[0]}}, $fields[1], $fields[2],$fields[3],$fields[6],$fields[8],(split /:/, $fields[$#fields])[2]);
				}
			}
			
		}
		else{
			# pop scaffold hash;
			
			if($fields[0]=~/#/){
                                $key1 = (split /#/,$fields[0])[0];
				if($fields[4] >= 10){
					$scaffold{$key1}{$fields[0]}=[];
                                	push(@{$scaffold{$key1}{$fields[0]}}, $fields[1], $fields[2],$fields[3],$fields[6],$fields[8],(split /:/, $fields[$#fields])[2]);
                        	}
			}   
                        else{
				if($fields[4] >= 10){
					$scaffold{$fields[0]}{$fields[0]}=[] if(!exists($scaffold{$fields[0]}));
                                	push(@{$scaffold{$fields[0]}{$fields[0]}}, $fields[1], $fields[2],$fields[3],$fields[6],$fields[8],(split /:/, $fields[$#fields])[2]);
                        	}
			}     				
		}

	}
	#close FHOUT;
	$date = `date`;
	print "Done testing local assemblies..\nEnd time is $date";
}

sub get_clone_measures{
	
	my ($scaffold_ref , $lib_stats) = @_;
	my ($key1 , $key2 , $i , $flag , $scaffold, $pos , $mate_scaff, $isize , $lib);
	my ($Zscore , $Avg_Zscore , $RGB);

	my $Z_threshold 		= 3;
	my $FR_count			= 0;
	my $MP_good			= 0;
	my $MP_bad			= 0;
	my $MP_stretched		= 0;
	my $MP_compressed		= 0;
	my $MP_outward_facing_reads	= 0;
	my $MP_unmapped_mate		= 0;
	my $MP_diff_scaff		= 0;
	my $MP_same_strand		= 0;

	#array order is flag, scaffold, map_pos, mate_scaff, insert, lib 

	foreach $key1 ( keys %{$scaffold_ref}){
        	foreach $key2 ( keys %{$scaffold_ref->{$key1}}){
			#print $key2,"\n";
			$flag 		= $scaffold_ref->{$key1}->{$key2}->[0];
			$scaffold 	= $scaffold_ref->{$key1}->{$key2}->[1];
			$pos 		= $scaffold_ref->{$key1}->{$key2}->[2];
			$mate_scaff	= $scaffold_ref->{$key1}->{$key2}->[3];
			$isize		= $scaffold_ref->{$key1}->{$key2}->[4];
			$lib		= $scaffold_ref->{$key1}->{$key2}->[5];
                        
	################## cal # of good and bad clones for the scaffold ##############################
	
			# -> <- FR
			if( ($flag & 2) && (!($flag & 4)) && (!($flag & 8)) && ($scaffold=~/^scaffold/) ){
				if($mate_scaff eq "="){
					
					$Zscore = (abs($isize) - ${$lib_stats{$lib}}{"Mean"}) / ${$lib_stats{$lib}}{"Stddev"};
					#print "$lib\t$Zscore\n";
					
					$Avg_Zscore += $Zscore;
					$FR_count++;
					
					if( ( ($Zscore <= $Z_threshold) && ($Zscore >=0 ) ) ||  ($Zscore >= -$Z_threshold) ) {
						$MP_good++;
					}
					elsif ($Zscore > $Z_threshold){
						$MP_stretched++;
					}
					elsif ($Zscore < -$Z_threshold){
						$MP_compressed++;
					}
				}
				else{
					$MP_diff_scaff++; # inter scaffold mapping mates
				}
			}
					
			
			#<- -> ? RF ? check on this
			elsif ( (!($flag & 2)) && (!($flag & 4)) && (!($flag & 8)) && ($scaffold=~/^scaffold/) && ($mate_scaff eq "=")){
				if( ($flag & 16) && ($flag && 32)){
					$MP_same_strand++; # <-<-
				}
				elsif( (!($flag & 16)) && (!($flag & 32))){
					$MP_same_strand++; # ->->
				}
				else{
					#$sample{$lib}{"Mean"} += $isize;
					#$sample{$lib}{"Count"}++;
					$MP_outward_facing_reads++;	
				}
			}
			#case 3: unmapped mate
			elsif ( ( (!($flag & 2)) && (!($flag & 4)) && ($flag & 8) )  || ( (!($flag & 2)) && ( $flag & 4 ) && (!($flag & 8 )) ) ){
				$MP_unmapped_mate++; # ->/* or */<-
			}
		
			last; # we look at only one read of the pair
		}
	}
	
	# CE stats didnt work well, its too sensitive to little differences
	$MP_bad = $MP_stretched + $MP_compressed + $MP_diff_scaff + $MP_same_strand + $MP_outward_facing_reads + $MP_unmapped_mate;

	if ( ($MP_bad > 0) && ($FR_count > 0) && ($MP_good > 0)){
		$RGB = log($MP_good/$MP_bad);
		$Avg_Zscore /= $FR_count;
	}
	elsif(($MP_bad > 0) && ($FR_count == 0) && ($MP_good > 0)){
		$RGB = log($MP_good/$MP_bad);
		$Avg_Zscore = 999.99999;
	}
	elsif(($MP_bad == 0) && ($FR_count > 0) && ($MP_good > 0)){
		$RGB = 999.99999;
		$Avg_Zscore /= $FR_count;
	}
	elsif($MP_good == 0){
		$RGB = 999.99999;
	}
		
	
	print "MP_good\t$MP_good\tMP_stretched\t$MP_stretched\tMP_compressed\t$MP_compressed\tMP_diff_scaff\t$MP_diff_scaff\tMP_samestrand\t$MP_same_strand\tMP_outward_facing_reads\t$MP_outward_facing_reads\tMP_unmapped_mate\t$MP_unmapped_mate\tGMB\t",$MP_good-$MP_bad,"\tRGB\t",sprintf("%.5f",$RGB),"\tAvg_Zscore\t", sprintf("%.5f",$Avg_Zscore),"\n" ;
	#print FHOUT "MP_good\t$MP_good\tMP_stretched\t$MP_stretched\tMP_compressed\t$MP_compressed\tMP_diff_scaff\t$MP_diff_scaff\tMP_samestrand\t$MP_same_strand\tMP_outward_facing_reads\t$MP_outward_facing_reads\tMP_unmapped_mate\t$MP_unmapped_mate\tGMB\t",$MP_good-$MP_bad,"\tRGB\t",sprintf("%.5f",$RGB),"\tAvg_Zscore\t", sprintf("%.5f",$Avg_Zscore),"\n" ;

}
