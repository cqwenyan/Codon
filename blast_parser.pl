#!/usr/bin/perl

=head1 Name

	blast_parser.pl -- parse the BLAST result and convert to tabular format.

=head1 Usage

  	blast_parser.pl [options] input_file
	-nohead     do not show the first instruction line.
	-tophit     integer, to set how many subjects for a query to be displayed. 
	-topmatch   integer, to set suits(results of one subject match one query) to be displayed. 
	-eval       float or exponent,to filter the alignments which worse than the E-value cutoff.
	-verbose    output verbose information to screen.
	-help       output help information to screen.

=head1 Exmple

	perl blast_parser.pl -nohead -tophit 1 -topmatch 1 -eval 1e-5  blast.out > blast.out.xls

=cut

use strict;
use Getopt::Long;
use Data::Dumper;

my ($Nohead,$Tophit,$Topmatch,$Type,$Eval);
my ($Verbose,$Help);
GetOptions(
	"nohead"=>\$Nohead,
	"tophit:i"=>\$Tophit,
	"topmatch:i"=>\$Topmatch,
	"eval:f"=>\$Eval,
	"verbose"=>\$Verbose,
	"help"=>\$Help
);
die `pod2text $0` if (@ARGV==0 || $Help);

my $blast_file = shift;


##convert blast raw result to tabular format
&parse_blast($blast_file,$Tophit,$Topmatch,$Eval,$Nohead);



##for othe applications, you can distill information from %Data;
#my %Data;
#&read_blast($blast_file,\%Data,$Tophit,$Topmatch,$Eval);
#print Dumper \%Data;



####################################################
################### Sub Routines ###################
####################################################


##parse the BLAST files, and output in tabular formats 
####################################################
sub parse_blast{
	my ($file,$tophit,$topmatch,$eval,$nohead) = @_;
	open (BLAST,"$file") || die ("Could not open the blast file.");

	print "Query_id\tQuery_length\tQuery_start\tQuery_end\tSubject_id\tSubject_length\tSubject_start\tSubject_end\t",
		"Identity\tPositive\tGap\tAlign_length\tScore\tE_value\tQuery_annotation\tSubject_annotation\n" unless(defined $nohead);

	my $type=<BLAST>;
	$type=~s/\s+.+//s;	#blast的类型
	$/="$type";	#用$type来隔开每一个query的比对结果
	
	my $database;
	while (<BLAST>) {
		next if(/No hits found/);
		my @cycle=split (/\n>/,$_); 	#利用"\n>"将Query序列和每一条subject序列比对结果信息分开

		my ($pointer,$query,$query_len,$subject,$subject_len,$query_annotation,$subject_annotation);
		if ($cycle[0]=~/Query= (\S+)\s+(.+?)\((\S+)\s+letters\)/s) {
			$query=$1;
			$query_annotation=$2;
			$query_len=$3;
			$query_annotation=~s/\n//g;
			$query_annotation=~s/\s{2,}/ /g;
			$query_annotation="--" if($query_annotation eq " ");	#信息为空时用"--"代替
		}	#提取Query id,Query length,Query annotation信息

		shift @cycle;
		for (my $i=0;$cycle[$i];$i++) {
			last if((defined $tophit) && $i>$tophit-1);
			if ($cycle[$i]=~/(\S+)\s+(.+?)\s+Length = (\d+)/s) {
				$subject=$1;
				$subject_annotation=$2;
				$subject_len=$3;
				$subject_annotation=~s/\n//g;
				$subject_annotation=~s/\s{2,}/ /g;
				$subject_annotation="--" if($subject_annotation eq " ");	#信息为空时用"--"代替
			}	#提取Subject id,Subject length,Subject annotation信息

			my @cycle_inner=split (/Score/,$cycle[$i]);	#分开同一个query和同一个subject的多个比对结果
			shift @cycle_inner;
			for (my $j=0;$cycle_inner[$j];$j++) {
				last if((defined $topmatch) && $j>$topmatch-1);
				if ($cycle_inner[$j]=~/(\S+) bits.+?Expect.+?=\s+(\S+).+?Identities\D+\d+\/(\d+)\s+\((\S+)\%\)/s) {
					$pointer->[$i][$j]{score}=$1;
					$pointer->[$i][$j]{e_value}=$2;
					$pointer->[$i][$j]{align_len}=$3;
					$pointer->[$i][$j]{identity}=$4/100;
					chop $pointer->[$i][$j]{e_value} if ($pointer->[$i][$j]{e_value}=~/,$/);
					$pointer->[$i][$j]{e_value}=~s/^e/1e/;
					last if((defined $eval) && $pointer->[$i][$j]{e_value}>$eval);
				}	#提取Score,E value,Identity,Align_len信息

				$pointer->[$i][$j]{positive}=($cycle_inner[$j]=~/Positives = \S+\s+\((\S+)\%\)/s)? $1/100 : "--";
				#BLASTN结果文件中无Positive信息,其他几种结果文件都有

				$pointer->[$i][$j]{gap}=($cycle_inner[$j]=~/Gaps = \S+\s+\((\S+)\%\)/)? $1/100 : 0;
				#Gap信息不是每一组里都有,单独考虑

				$pointer->[$i][$j]{q_start}=$1 if($cycle_inner[$j]=~/Query\D+(\d+)\D+/);
				$pointer->[$i][$j]{s_start}=$1 if($cycle_inner[$j]=~/Sbjct\D+(\d+)\D+/);
			
				$cycle_inner[$j]=~s/Database.+?$//s;	#最后一组情况中需删除文件末尾的某些多余信息
				if ($cycle_inner[$j]=~/\D+(\d+)\D+Sbjct\D+\d+\D+(\d+)\D+$/s) {
					$pointer->[$i][$j]{q_end}=$1;
					$pointer->[$i][$j]{s_end}=$2;
				}

				print "$query\t$query_len\t$pointer->[$i][$j]{q_start}\t$pointer->[$i][$j]{q_end}\t",
					"$subject\t$subject_len\t$pointer->[$i][$j]{s_start}\t$pointer->[$i][$j]{s_end}\t",
					"$pointer->[$i][$j]{identity}\t$pointer->[$i][$j]{positive}\t$pointer->[$i][$j]{gap}\t$pointer->[$i][$j]{align_len}\t",
					"$pointer->[$i][$j]{score}\t$pointer->[$i][$j]{e_value}\t$query_annotation\t$subject_annotation\n";
			}
		}
	}
	$/="\n";
	close(BLAST);
}


##parse the BLAST files, and read into data structure
##the alignments are stored in two demension array
##in the order from good to bad.
####################################################
sub read_blast{
	my ($file,$data_p,$tophit,$topmatch,$eval) = @_;
	open (BLAST,"$file") || die ("Could not open the blast file.");

	my $type=<BLAST>;
	$type=~s/\s+.+//s; 
	$/="$type";	#用$type来隔开每一个query的比对结果
	
	while (<BLAST>) {
		next if(/No hits found/);
		my @cycle=split (/\n>/,$_); 	#利用"\n>"将Query序列和每一条subject序列比对结果信息分开

		my ($pointer,$query,$query_len,$subject,$subject_len,$query_annotation,$subject_annotation);
		if ($cycle[0]=~/Query= (\S+)\s+(.+?)\((\S+)\s+letters\)/s) {
			$query=$1;
			$query_annotation=$2;
			$query_len=$3;
			$query_annotation=~s/\n//g;
			$query_annotation=~s/\s{2,}/ /g;
			$query_annotation="--" if($query_annotation eq " ");	#信息为空时用"--"代替
		}	#提取Query id,Query length,Query annotation信息

		shift @cycle;
		for (my $i=0;$cycle[$i];$i++) {
			last if((defined $tophit) && $i>$tophit-1);
			if ($cycle[$i]=~/(\S+)\s+(.+?)\s+Length = (\d+)/s) {
				$subject=$1;
				$subject_annotation=$2;
				$subject_len=$3;
				$subject_annotation=~s/\n//g;
				$subject_annotation=~s/\s{2,}/ /g;
			}	#提取Subject id,Subject length,Subject annotation信息

			my @cycle_inner=split (/Score/,$cycle[$i]);	#分开同一个query和同一个subject的多个比对结果
			shift @cycle_inner;
			for (my $j=0;$cycle_inner[$j];$j++) {
				last if((defined $topmatch) && $j>$topmatch-1);
				if ($cycle_inner[$j]=~/(\S+) bits.+?Expect.+?=\s+(\S+).+?Identities\D+\d+\/(\d+)\s+\((\S+)\%\)/s) {
					$pointer->[$i][$j]{subject_id}=$subject;
					$pointer->[$i][$j]{score}=$1;
					$pointer->[$i][$j]{e_value}=$2;
					$pointer->[$i][$j]{align_len}=$3;
					$pointer->[$i][$j]{identity}=$4/100;
					chop $pointer->[$i][$j]{e_value} if ($pointer->[$i][$j]{e_value}=~/,$/);
					$pointer->[$i][$j]{e_value}=~s/^e/1e/;
					last if((defined $eval) && $pointer->[$i][$j]{e_value}>$eval);
				}	#提取Score,E value,Identity,Align_len信息

				$pointer->[$i][$j]{positive}=($cycle_inner[$j]=~/Positives = \S+\s+\((\S+)\%\)/s)? $1/100 : "--";
				#BLASTN结果文件中无Positive信息,其他几种结果文件都有

				$pointer->[$i][$j]{gap}=($cycle_inner[$j]=~/Gaps = \S+\s+\((\S+)\%\)/)? $1/100 : 0;
				#Gap信息不是每一组里都有,单独考虑

				$pointer->[$i][$j]{q_start}=$1 if($cycle_inner[$j]=~/Query\D+(\d+)\D+/);
				$pointer->[$i][$j]{s_start}=$1 if($cycle_inner[$j]=~/Sbjct\D+(\d+)\D+/);
			
				$cycle_inner[$j]=~s/Database.+?$//s;	#最后一组情况中需删除文件末尾的某些多余信息
				if ($cycle_inner[$j]=~/\D+(\d+)\D+Sbjct\D+\d+\D+(\d+)\D+$/s) {
					$pointer->[$i][$j]{q_end}=$1;
					$pointer->[$i][$j]{s_end}=$2;
				}
				
				$subject_annotation="--" if($subject_annotation eq " ");	#信息为空时用"--"代替
				$pointer->[$i][$j]{subject_annotation} = $subject_annotation;

			}
		}

		$data_p->{$query}{query_annotation} = $query_annotation;
		$data_p->{$query}{query_length} = $query_len;
		$data_p->{$query}{alignment} = $pointer; 
	}
	$/="\n";
	close(BLAST);
}
