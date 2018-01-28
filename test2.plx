#!/usr/bin/perl -w
use warnings;
use strict;
use DBI;
use Bio::EnsEMBL::Registry;;
use Data::Dumper;
use IO::Handle;
use Storable;
use File::Path;

my $start_run = time();

#############################Establish local connection###################################
my $registry = 'Bio::EnsEMBL::Registry';


$registry->load_registry_from_db(
    -host => '127.0.0.1', # establish a local connection
    -user => 'root'
    
);

#############################Organize Files, Data, and Inputs#############################
sub reverse_complement {
        my $dna = shift;

	# reverse the DNA sequence
        my $revcomp = reverse($dna);

	# complement the reversed DNA sequence
        $revcomp =~ tr/ACGTacgt/TGCAtgca/;
        return $revcomp;
}

my @chromosomes = ('X','Y','1', '2', '3', '4', '5', '6', '7', '8', '9', '10', '11', '12', 
				   '13', '14', '15', '16', '17', '18', '19', '20', '21', '22'); 
my @chromosomeProbs = (0.0501572157997191, 0.0693367226080984, 0.149852418472859, 0.228413373195488,
					   0.292380771488473, 0.354129542616568, 0.412570793406687, 0.467846277323937,
					   0.519253010591144, 0.566533137852672, 0.612149464170332, 0.655931399095017,
					   0.699542698992307, 0.74278102042759, 0.779984469841782, 0.814661707070659,
					   0.847782532775092, 0.876969925702323, 0.903198501614418, 0.928419878912112, 
					   0.947520377811253, 0.967879579243446, 0.983427030929927, 1); 

my $dir = "indexed_chrom_trinuc/";
my $dir2 = "chromosome_marker/";
					   

my %samples;
my %mutations;
my %probsCatalogue;
my @trinucleotides;



open (SAMPLES, "<", "sample_input2.txt") or die $!; # Open Sample file
open (FILE, "<", "signatures_probabilities2.txt") or die $!; # Open Signatures File
my $outputFile = "generated_mutational_catalogue_Count_10000_v7.txt"; # Create an output file
open (my $out, ">", $outputFile) or die $!;

my $firstLine = <SAMPLES>;
my @signatures = split 'Signature', $firstLine;
@signatures = @signatures[1 .. $#signatures]; # Save the signatures of interest

while (<SAMPLES>){
	next if $. < 2;
	my @cols = split(' ');
	$samples{$cols[0]} = [@cols[1 .. $#cols]]; # Save each sample and its corresponding mutation counts for each signature
}
close SAMPLES;

print {$out} "CancerSample\tGenomeAssemblyVersion\tChromosome\tPositionStart\tPositionEnd\tReference\tMutation\n";
 


while (<FILE>){
	next if $. < 2;
	my @cols = split(' ');
	my $ref = substr($cols[0], 0, 1);
	my $mut = substr($cols[0], 2, 1);
	my $trinuc = $cols[1];
	foreach my $sig (@signatures){
		my $probability = $cols[$sig+2];
		$probsCatalogue{$sig}{$trinuc}{$ref,'>',$mut} = $probability;
	}
	if (grep {$_ eq $trinuc} @trinucleotides ){next} else {push @trinucleotides, $trinuc};
}
close FILE;

#############################Construct Mutational Catalogue###############################
my $i = 0;
foreach my $tri (@trinucleotides) {
	my @RevCompBreak = retrieve ($dir2 . $tri . 'marker.data');
	my $mutations = retrieve ($dir . $tri . '.data');
# 	foreach my $chrom (@chromosomes) {
# 		my $slice_adaptor = $registry->get_adaptor( 'Human', 'Core', 'Slice' );
# 		my $slice = $slice_adaptor->fetch_by_region('chromosome', $chrom);
# 		my $sequence = $slice->seq();
# 		my $offset = 0;
# 		my $result = index($sequence, $tri, $offset);
# 		while ($result != -1) {
# 			push @{$mutations{$chrom}}, $result + 2;
# 			$offset = $result + 1;
# 			$result = index($sequence, $tri, $offset);
# 		}
# 		push @RevCompBreak, scalar( @{ $mutations{$chrom} } ); # Marks point of reverse complement entries
# 		my $triReverse = reverse_complement($tri);
# 		$offset = 0;
# 		$result = index($sequence, $triReverse, $offset);
# 		while ($result != -1) {
# 			push @{$mutations{$chrom}}, $result + 2; # Saves location of middle base
# 			$offset = $result + 1;
# 			$result = index($sequence, $triReverse, $offset);
# 		}
# 	}
	foreach my $sample (keys %{ samples}){ 
		my $sigIndex = 0;
		foreach my $sig (@signatures){
			foreach my $type (keys %{ $probsCatalogue{$sig}{$tri} } ){
				my $mutationCount = $probsCatalogue{$sig}{$tri}{substr($type,0,1),'>',substr($type,4,1)}*int($samples{$sample}[$sigIndex]);
				if ($mutationCount - int($mutationCount) >= 0.5) {$mutationCount = int($mutationCount+1)} else {$mutationCount = int($mutationCount)}
				my @randomNumberList = ();
				for (my $k=0; $k < $mutationCount; $k++){
					my $random_number_chrom = rand(1);
					my $reverseIndex;
					for (my $m=0; $m < scalar(@chromosomeProbs); $m++){
						if ($random_number_chrom < $chromosomeProbs[$m]){
							if ($m == 0) {$random_number_chrom = 'X'} elsif ($m == 1) {$random_number_chrom = 'Y'} else {$random_number_chrom = $m - 1};
							$reverseIndex = $m;
							last;
						}
					}
					my $range = scalar(@{$mutations->{$random_number_chrom}});
					my $random_number = int(rand($range));
					if ($random_number ~~ @randomNumberList){
						my $random_number = int(rand($range));
					}
					push @randomNumberList, $random_number;	
					if ($random_number < $RevCompBreak[0][$reverseIndex]) { 
						print {$out} "$sample\tGRCh37\t$random_number_chrom\t",$mutations->{$random_number_chrom}[$random_number],"\t",$mutations->{$random_number_chrom}[$random_number],"\t",substr($type,0,1),"\t",substr($type,4,1),"\t","for", "\n";
					}
					else{
						print {$out} "$sample\tGRCh37\t$random_number_chrom\t",$mutations->{$random_number_chrom}[$random_number],"\t",$mutations->{$random_number_chrom}[$random_number],"\t",reverse_complement(substr($type,0,1)),"\t",reverse_complement(substr($type,4,1)),"\t","rev", "\n";
					}
				}
 			}
 			$sigIndex += 1;
 		}
 	}
 	$out->flush();
 	$mutations = ();
 	@RevCompBreak = ();
 	$i = $i + 1;
 	print "Finished ", $i, " iterations\n";
 	
}


close $out;

system("sort -t \$'\t' -k 3,3n -k 3 $outputFile -o $outputFile");

my $end_run = time();
my $run_time = $end_run - $start_run;
print "Job took $run_time seconds\n";
# expected count 8272
#new Version:
#actual count 8430
#took 434 seconds for one sample

#old version:
#actual count: 8430
#took 3265 seconds for one sample

#502 seconds for all samples
#try rewriting the pulling script by pulling out one chromosome at a time and finding all mutations for each chrom