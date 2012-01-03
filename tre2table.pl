#!/usr/bin/perl
use lib "/usr/share/perl5";
use Bio::TreeIO;
use strict;
use warnings;
use diagnostics;

if (@ARGV==0) {
print "This program takes:
1) a file containing a Newick format tree (required)
2) a file containing a list of refseq ids matching those in the tree (optional)

and outputs:
1) a treetable (for input into pangenome functions)
2) a full distance matrix (optional)
\n";
exit 0;
}

my $debug = 0;
my $printmat =0; #print out distance matrix? 0/1
my $method = 0; #0=just named subset; 1=whole subtree
$method ==1 ? print "using whole subtree\n\n" : print "using named subset\n\n" if $debug>0;

#-------------------------
# set up tree
#-------------------------

my $mytree = shift(@ARGV);

my $treeio = Bio::TreeIO->new(-format => 'newick', -file => $mytree);
my $tree = $treeio->next_tree;
$tree->move_id_to_bootstrap;
my $lca;

my @leafnodes;
my @allnodes;
my @distnodes = $tree->get_leaf_nodes;


my $blastfile = shift(@ARGV);
#if blastlist, choose only a subset of tree
if (defined($blastfile)) {
  open(FH,"<$blastfile");
  my @blastlist = <FH>;
  close(FH);

  my @usednodes;
  foreach my $id (@blastlist) {
    chomp $id;
    my $addnode = $tree->find_node(-id => $id);
    if ($addnode) {
      push(@usednodes,$addnode);
    } else {print "$id not in tree\n" if ($debug>0);}
  }

    #get last common ancestor of wanted nodes
    $lca = $tree->get_lca(-nodes => \@usednodes);
    @allnodes = $lca->get_all_Descendents;

    if ($method == 1) {
      @leafnodes = grep { $_->is_Leaf } @allnodes;
    } else {
      push(@leafnodes,@usednodes);
    }
} else {
    @allnodes = $tree->get_nodes;
    @leafnodes = $tree->get_leaf_nodes;
    $lca = $tree->get_root_node;
}


#-------------------------
# get distance tables
#-------------------------

#number nodes by how many steps to lca
#exception: number leafs first
#first count steps for every node
#built by coalescent
my %steps;
my %coalescent;
my $tmpnodeid=100000000;

#start with leaves
push (my @nodelol,[@leafnodes]);


while (scalar(@nodelol)>0) {
print "new level $tmpnodeid\n" if $debug>0;
  my $noderef = shift(@nodelol);
  my @nodelist = @$noderef;
  my @addnodes;

foreach my $node (@nodelist) {
  next unless defined($node);
  unless ($node->id) {$node->id($tmpnodeid);$tmpnodeid++;}
  next unless defined($node->ancestor);
  push(@addnodes,$node->ancestor) unless (grep { $_ eq $node->ancestor} @addnodes);

  my $step=0;
  my $tmpnode = $node;
  $coalescent{$tmpnode}++;
  until ($tmpnode eq ($lca+1)) {
    last unless defined($tmpnode->ancestor);
    $tmpnode = $tmpnode->ancestor;
    $step++;
    $coalescent{$tmpnode}++;
  }
  $steps{$node->id}{'step'} = $step;
  $steps{$node->id}{'leaf'} = $node->is_Leaf;
  $steps{$node->id}{'node'} = $node;
}
print scalar(@addnodes),"\n" if $debug>0;
push (@nodelol,[@addnodes]) if scalar(@addnodes)>0;
}
# exit 0;

#prune the subtree to only the branches that lead to our desired taxa
if ($method == 0) {
  foreach my $node (@allnodes) {
    for my $child ( $node->each_Descendent ) {
#       unless (grep { $_ eq $child} @coalescent) {
	unless ($coalescent{$child}>0) {
	$tree->splice(-remove_id=>$child->id,-preserve_lengths=>1);
	print "deleting node " . $child->id . "\n" if $debug>0;
      }
    }
  }
}

#get new list of remaining nodes
my @newnodes;
if (defined($blastfile)) {
  @newnodes = $lca->get_all_Descendents;
} else {@newnodes = $tree->get_nodes;}

#rename nodes and print out sorted list

my $newnode=0;
foreach my $nodeid (sort {$steps{$b}{'leaf'} <=> $steps{$a}{'leaf'} || $steps{$b}{'step'} <=> $steps{$a}{'step'}} keys %steps) {
my $node = $steps{$nodeid}{'node'};
  next unless defined($node);
  next unless grep {$_ eq $node} @newnodes;
  $newnode++;
  $node->id($newnode);
  print $nodeid,"\t",$steps{$nodeid}{'step'},"\t",$node->id,"\n" if $debug>0;
}



# print table for nodes in subtree
for my $node (sort {$a->id <=> $b->id} @newnodes) {
  next unless $node->ancestor;
  my @children = $node->each_Descendent;
  next if (@children==1);
  my @siblings = $node->ancestor->each_Descendent;
  print "siblings=" . @siblings . "\tchildren=" . @children . "\n" if $debug<0;
  #skip if internal node without divergence

#find ancestor in order to get distance
  my $distance;
#if node is without siblings
  if (@siblings==1) {
      print $node->id . "->" if $debug<0;
      my $ancestor = $node->ancestor;
      my @sibs = $ancestor->each_Descendent; #should be 1
    while (scalar(@sibs)<2) {
      print $ancestor->id . "->" if $debug<0;
      $ancestor = $ancestor->ancestor;
      @sibs = $ancestor->each_Descendent; #can be 1 or 2
    }
    print $ancestor->id . "(" . @children . ")\n" if $debug<0;
    #now $ancestor is the ancestor
    $distance = $tree->distance(-nodes => [$node,$ancestor]);
  } else {
    $distance = $tree->distance(-nodes => [$node,$node->ancestor]);
  }

  #get children labels
  my @kids;
  if ($node->is_Leaf) {
    @kids = (0,0);
  } else {
    foreach my $child (@children) {
      my $thischild = $child;
      #evaluate to 0 or 2 is fine
      my @newchildren = $thischild->each_Descendent;
      while (scalar(@newchildren)==1) {
	$thischild = $newchildren[0];
	@newchildren = $thischild->each_Descendent;
      }
      push(@kids,$thischild->id);
    }
    @kids = sort {$a <=> $b} @kids;
  }

  print join("\t",($node->id,@kids,sprintf("%f",$distance))) . "\n";
}

{#for lca
my @dec = $lca->each_Descendent;
my @kids;
foreach (@dec) {

  push(@kids,$_->id);
}
@kids = sort {$a <=> $b} @kids;
print join("\t",(++$newnode,@kids,0)),"\n";
}

#-------------------------
# print distance matrix
#-------------------------
if($printmat>0) {
print "printing distance matrix to $mytree.distmatrix.csv\n" if $debug>0;
open(OUT,">$mytree.distmatrix.csv") or die "$0 help!\n";
print OUT join("\t",map {$_->id} @leafnodes),"\n";
for (my $i=0;$i<@leafnodes;$i++) {
  print OUT $leafnodes[$i]->id;
    for (my $j=$i;$j<@leafnodes;$j++) {
      my $distance = $tree->distance(-nodes => [$leafnodes[$i],$leafnodes[$j]]);
      print OUT "\t$distance";
  }
  print OUT "\n";
}
close(OUT);
}
