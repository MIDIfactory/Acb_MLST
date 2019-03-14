

# This pipe allows to obtain a concatenate of orthologues genes starting to be used in phylogenomic analysis

# Required inputs:
#	(1) Fasta of reference genes (for example all the orfs of a bacterium) 
#	(2) A folder conteining the ffn (multifasta of genes) files of the organisms

# For firt orthologues genes are identified using the Double Best Hist algorithm and saved into the *_genesDB folder (see into the *_tmp folder). Genes shared by all the organisms are copied into the *_genesDB_shared folder, aligned using Muscle and concatenate into the *.concatenate file.	

# Command perl Pipe_phylogenomics.pl [multifasta database] [folder]


# Currently the script is written to work with nucleotides



$db=shift;
$folder=shift;
$dbtype=shift;

if ($dbtype eq ""){$dbtype = "nucl";} #$dbtype = nucl or prot

# makeblastdb

`ls $folder/* > $folder\_file_list`;

`ls $folder | parallel \"sed -i \\\"s/>/>{}\\\&\\\&\\\&\\\&\\\&/g\\\" $folder/{}\"`;

`cat $folder/* > $folder.ffn`; 

`ls $folder/* | parallel \"makeblastdb -in {} -dbtype $dbtype\"`;
`makeblastdb -in $db -dbtype $dbtype`;

# blast

if ($dbtype eq "nucl")
{

`cat $folder\_file_list | parallel \"blastn -task blastn -query {} -db $db -evalue 0.00001 -max_target_seqs 1 -outfmt '6 qseqid sseqid' >> {}.outblast\"`; 

`cat $folder\_file_list | parallel \"blastn -task blastn -query $db -db {} -evalue 0.00001 -max_target_seqs 1 -outfmt '6 qseqid sseqid' >> {}.outblast\"`;

}

`cat $folder/*.outblast > $folder\_all_outblast`;

open(BLA,"<$folder\_all_outblast");
open(DOU,">$folder\_all_outblast.DoubleBestHit");

while(<BLA>)
{
	chomp $_;
	@s=split("\t",$_);
	$hash{$s[0]}{$s[1]}="1";
	
	if ($hash{$s[1]}{$s[0]} eq "1")
	{print DOU $s[0],"\t",$s[1],"\n";}

}

close(BLA);
close(DOU);


$DBH_tab="";
$fasta=shift;

	use Bio::Perl;
	use Bio::SeqIO;
	use Bio::Seq;

open(DBH,"<$folder\_all_outblast.DoubleBestHit");

while(<DBH>)
{
	chomp $_;
	@s=split("\t",$_);
	$hash{$s[1]}=$s[0];
}

mkdir("$folder\_genesDB");

my $seqobj  = Bio::SeqIO->new(-format => 'fasta',-file => "$folder.ffn") or die;

while($seq=$seqobj->next_seq)
{
	$name=$seq->display_id;

	if ($hash{$name} ne "")
	{
		open(GEN,">>$folder\_genesDB/$hash{$name}");
		print GEN ">",$name,"\n",$seq->seq,"\n";
		close(GEN);

		$count{$hash{$name}}=$count{$hash{$name}}+1;
	}

}

$num=`wc -l $folder\_file_list`;
@sn=split(" ",$num);

mkdir("$folder\_genesDB_shared");

for $k(keys %count)
{
	if ($count{$k} == $sn[0])
	{`cp $folder\_genesDB/$k $folder\_genesDB_shared/`;}
}

`ls $folder\_genesDB_shared/* | parallel \"muscle -in {} -out {}.aln\"`;

@aln=`ls $folder\_genesDB_shared/*.aln`;

foreach $file(@aln)
{
	my $geneobj  = Bio::SeqIO->new(-format => 'fasta',-file => $file) or die;

	while($gene=$geneobj->next_seq)
	{
		@sg=split("&&&&&",$gene->display_id);
		$conca{$sg[0]}=$conca{$sg[0]}.$gene->seq;
	}

}


open(CON,">$folder.concatenate");

foreach $org(keys %conca)
{
	print CON ">",$org,"\n",$conca{$org},"\n";
}

mkdir("$folder\_tmp");

`mv $folder\_all_outblast $folder\_tmp/`;
`mv $folder\_all_outblast.DoubleBestHit $folder\_tmp/`;
`mv $folder\_file_list $folder\_tmp/`;

`mv $folder\_genesDB $folder\_tmp/`;
`mv $folder\_genesDB_shared $folder\_tmp/`;
`mv $db.nhr $folder\_tmp/`;
`mv $db.nin $folder\_tmp/`;
`mv $db.nsq $folder\_tmp/`;








