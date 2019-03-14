

	use Bio::Perl;
	use Bio::SeqIO;
	use Bio::Seq;

$folder=shift;


@aln=`ls $folder/*`;

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

