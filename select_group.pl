#!/usr/bin/perl -w
# 
# Copyright (c)   2016
# Writer:         xuxiong <xuxiong19880610@163.com><xiongxu\@me.com>
# Program Date:   2016.4.21
# Modifier:       xuxiong <xuxiong19880610@163.com><xuxiong\@me.com>
# Last Modified:  2016.4.22
my $ver="0.0.1";

use strict;
use Getopt::Long;
use Data::Dumper;
use FindBin qw($Bin $Script);
use FileHandle;
use File::Basename qw(basename dirname);
use List::Util qw(first max maxstr min minstr reduce shuffle sum);
use Cwd qw(abs_path getcwd realpath);

my ($infile,$outfile);
my $groupNum=12;
GetOptions(
			"help|?" =>\&USAGE,
			"i:s"=>\$infile,
			"o:s"=>\$outfile,
			"n:i"=>\$groupNum
			) or &USAGE;
&USAGE unless ($infile ) ;
###############Time_start#############################
my $Time_Start;
$Time_Start = sub_format_datetime(localtime(time()));
print STDERR "\nStart Time :[$Time_Start]\n\n";
my $BEGIN_TIME=time();
######################################################
Select_groups($infile,$outfile,$groupNum);

sub CMD{
	my ($shell,$run)=@_;
	my $CurrentTime = sub_format_datetime(localtime(time()));
	my $BEGIN_TIME=time();
	print STDERR "\n[$CurrentTime] Executing...\n$shell\n";
	my $return="";
	if ($run){
		$return=`$shell`;
		print STDERR "$return\nTime elapsed :",time()-$BEGIN_TIME," s\n";
	}
	return $return;
}

sub create_infile{
	my $infile=shift;
	if ($infile ne "-" && ! -f $infile) {
		print STDERR "Error: $infile not found \n";
		exit;
	}
	my $ii;
	if (!defined($infile) || $infile eq "-") {
		$ii = FileHandle->new();
		$ii->fdopen(fileno(STDIN), "r")or die $!;
	}else{
		$ii = FileHandle->new($infile, "r");
	}
	return $ii;
}

sub create_outfile{
	my $outfile=shift;
	my $io;
	if (!defined($outfile) || $outfile eq "-") {
		$io = FileHandle->new();
		$io->fdopen(fileno(STDOUT),"w") or die $!;
	}
	else{
		$io = FileHandle->new($outfile, "w");
	}
	return $io;
}

sub Select_groups{
	my ($infile,$outfile,$groupNum)=@_;
	my $ii=create_infile($infile);
	my @groups=();
	while (my $line=<$ii>) {
		chomp($line);
		my @unit=split(/,|\s+/,$line);
		push @groups,[@unit];
	}
	close($ii);
	@groups=shuffle(@groups);
	my $mean=sum(map {$_->[1]} @groups)/$groupNum;
	print STDERR "mean: ",$mean,"\n";
	my ($Sum,$PreIndex)=(0,0);
	my @GroupArray=();
	for(my $i=0;$i<@groups;++$i){
		$Sum+=$groups[$i][1];
		if(($i<$#groups && $Sum<=$mean && $Sum+$groups[$i+1][1]>=$mean) || $i==$#groups){
			push @GroupArray,[@groups[$PreIndex..$i]];
			$PreIndex=$i+1;
			$Sum=0;
		}
	}
	my @lastGroup=sort{$b->[1] <=> $a->[1]} @{pop(@GroupArray)};
	# map {print STDERR join("\t",@{$_}),"\n";} @lastGroup;
	@GroupArray= sort {sum(map {$_->[1]} @{$a}) <=> sum(map {$_->[1]} @{$b}) } @GroupArray;
	# for(my $i=0;$i<@GroupArray;++$i){
	# 	print STDERR $i,"\t",$GroupArray[$i][0][0],"\t",$GroupArray[$i][-1][0],"\t",sum(map {$_->[1]} @{$GroupArray[$i]}),"\n";
	# }
	# print STDERR "\n";
	for(my $i=0;$i<@lastGroup;++$i){
		push @{$GroupArray[$i]},$lastGroup[$i];
	}
	my @groupsIndex=();
	map {push @groupsIndex,-1} (0..255);
	for(my $i=0;$i<@GroupArray;++$i){
		print STDERR $i,"\t",$GroupArray[$i][0][0],"\t",$GroupArray[$i][-1][0],"\t",sum(map {$_->[1]} @{$GroupArray[$i]}),"\n";
		map {$groupsIndex[$_->[0]]=$i;} @{$GroupArray[$i]}
	}
	my $out=create_outfile($outfile);
	for (my $i=0;$i<256;$i+=16){
		print $out join(",",@groupsIndex[$i..$i+15]),",\n";
	}
	close($out);
}

###############Time_end###########
my $Time_End;
$Time_End = sub_format_datetime(localtime(time()));
print STDERR "\nEnd Time :[$Time_End]\n\n";

###############Sub_format_datetime
sub sub_format_datetime {#Time calculation subroutine
    my($sec, $min, $hour, $day, $mon, $year, $wday, $yday, $isdst) = @_;
	$wday = $yday = $isdst = 0;
    sprintf("%4d-%02d-%02d %02d:%02d:%02d", $year+1900, $mon+1, $day, $hour, $min, $sec);
}

sub USAGE {#
	my $usage=<<"USAGE";
Program: $0
Version: $ver
Contact: XiongXu <xiongxu\@me.com> <xuxiong19880610\@163.com>

Example:
  perl $0
Description:
  This program is used for define the group info from the first character of the binary safe coded string

Usage:
  -i             infile
  -o             outfile
  -n             group number

USAGE
	print STDERR $usage;
	exit;
}
