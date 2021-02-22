import pandas as pd
from random import shuffle
import gzip
import urllib.request
import os
from Bio import SeqIO
from subprocess import run
from Bio import Align, AlignIO
from Bio.Align import substitution_matrices
from Bio.Phylo import TreeConstruction, write as ph_write
from Bio.Phylo.Consensus import bootstrap
import requests
from bs4 import BeautifulSoup
import dendropy
import re
from ete3 import Tree
import argparse


# downloading a table with Uniprot IDs for proteomes
def download_ids(outfile):
    url = 'ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/reference_proteomes/README'
    with urllib.request.urlopen(url) as response:
        text = response.read()
    text = text.decode('utf-8')
    text = 'Proteome_ID' + text.split('Proteome_ID')[1]
    text = text.split('-----------------------------')[0].strip()
    with open(outfile, 'w') as f_ids:
        f_ids.write(text)


# picking n bacterial species at random (given that the number of proteins is more than lb and less than ub)
def choose_species(prot_ids_file, species_file, n=30, lb=3500, ub=5000, write_species=None):
    proteomes = pd.read_csv(prot_ids_file, sep='\t')
    my_species = {}
    with open(species_file) as f:
        lines = f.readlines()
        shuffle(lines)
        for sp in lines:
            sp = sp.split('_')[:2]
            sp = '_'.join(sp)
            if sp[0].islower():  # not a species
                continue
            for i, row in proteomes.iterrows():
                if sp.replace('_', ' ') in row['Species Name']:
                    if row['SUPERREGNUM'] == 'bacteria' and lb < row['#(1)'] < ub:
                        url = 'https://www.uniprot.org/proteomes/' + row['Proteome_ID']
                        html_content = requests.get(url).text
                        soup = BeautifulSoup(html_content, features='html.parser').prettify()
                        if 'Chromosome' in soup.split('basket-item component-checkbox" id=')[1][:100].split('\"')[1]:
                            my_species[sp] = row['Proteome_ID'] + '_' + str(row['Tax_ID'])
                            print(sp)
                            break
            if len(my_species) == n:
                if write_species:
                    with open(write_species, 'w') as s:
                        s.write('\n'.join(list(my_species.keys())))
                return my_species


def download_proteomes(species, fams_dir):
    def download_file(url, outfile):
        with urllib.request.urlopen(url) as response:
            with gzip.GzipFile(fileobj=response) as uncompressed:
                file_content = uncompressed.read()
        with open(outfile, 'wb') as f:
            f.write(file_content)

    if not os.path.isdir(fams_dir):
        os.mkdir(fams_dir)
    for sp in species:
        url = 'ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/reference_proteomes/Bacteria/'
        download_file(url + my_species[sp] + '.fasta.gz', os.path.join(fams_dir, sp + '.fasta'))


def combine_files(fams_dir, outfile):
    with open(outfile, 'w') as f:
        for file_name in os.listdir(fams_dir):
            i = 0
            for record in SeqIO.parse(os.path.join(fams_dir, file_name), 'fasta'):
                new_name = file_name.split('.')[0]
                new_name = new_name + '_' + str(i)
                record.id = new_name
                record.description = ''
                SeqIO.write(record, f, 'fasta')
                i += 1


# MMseqs2 clustering
def mmseqs(seqs_file, name, c, cov_mode=0, cluster_mode=0, e=0.001):
    run(['mmseqs', 'easy-cluster', seqs_file, name, 'tmp', '--cov-mode', str(cov_mode), '-c', str(c), '--cluster-mode',
         str(cluster_mode), '-e', str(e)])


def get_mmseqs_clustering(name):
    clustering = []
    with open(name + '_cluster.tsv') as f:
        last_i = ''
        for line in f:
            i, j = line.split()
            i = i.strip()
            j = j.strip()
            if last_i == i:
                clustering[-1].append(j)
            else:
                clustering.append([j])
            last_i = i
    return clustering


# filtering gene clusters to get one-to-one gene-species correspondence
def write_clusters_consensus(clusters, out_dir, proteins_file, n=30, max_seqs=120):
    if not os.path.isdir(out_dir):
        os.mkdir(out_dir)
    aligner = Align.PairwiseAligner()
    aligner.substitution_matrix = substitution_matrices.load("BLOSUM62")
    aligner.mode = 'global'
    i = 0
    for cl in clusters:
        sp = ['_'.join(x.split('_')[:2]) for x in cl]
        if len(set(sp)) == n:
            if len(sp) < max_seqs:
                with open(os.path.join(out_dir, str(i) + '.fasta'), 'w') as f:
                    if len(sp) == n:
                        for record in SeqIO.parse(proteins_file, 'fasta'):
                            if record.id in cl:
                                record.id = '_'.join(record.id.split('_')[:2])
                                record.description = ''
                                SeqIO.write(record, f, 'fasta')
                    else:  # more than n proteins
                        names = {s: [] for s in sp}
                        sequences = {}
                        for record in SeqIO.parse(proteins_file, 'fasta'):
                            if record.id in cl:
                                names['_'.join(record.id.split('_')[:2])].append(record.id)
                                sequences[record.id] = record.seq
                        seqs = sequences.values()
                        final_n = []
                        rejected_seqs = set()
                        for s in names:
                            if len(names[s]) == 1:
                                final_n.append(names[s][0])
                            else:  # > 1 protein from species s
                                max_score = 0
                                for protein in names[s]:
                                    seq1 = sequences[protein]
                                    score = 0
                                    for seq2 in seqs:
                                        if seq2 not in rejected_seqs:
                                            try:
                                                score += aligner.score(seq1, seq2)
                                            except ValueError:  # sequence contains letters not in the alphabet
                                                aligner1 = Align.PairwiseAligner()
                                                aligner1.mode = 'global'
                                                score += aligner1.score(seq1, seq2)
                                    if score > max_score:
                                        max_score = score
                                        p = protein
                                final_n.append(p)
                                for protein in names[s]:
                                    if protein != p:
                                        rejected_seqs.add(sequences[protein])
                        for record in SeqIO.parse(proteins_file, 'fasta'):
                            if record.id in final_n:
                                record.id = '_'.join(record.id.split('_')[:2])
                                record.description = ''
                                SeqIO.write(record, f, 'fasta')
                print(i)
                i += 1


# b - whether to use bootstrap; boot_dir - directory for keeping bootstrap trees; b_n - number of bootstrap replicates
def build_trees(fams_dir, outfile, b=False, boot_dir='bootstraps', b_n=50):
    trees = []
    i = 0
    if b:
        if not os.path.isdir(boot_dir):
            os.mkdir(boot_dir)
    for file_name in os.listdir(fams_dir):
        print(i)
        with open('tmp.txt', 'w') as f:
            run(['mafft', os.path.join(fams_dir, file_name)], stdout=f)
        aln = AlignIO.read('tmp.txt', 'fasta')
        calculator = TreeConstruction.DistanceCalculator('blosum62')
        dm = calculator.get_distance(aln)
        constructor = TreeConstruction.DistanceTreeConstructor()
        tree = constructor.nj(dm)
        for node in tree.get_nonterminals():
            node.name = None
        trees.append(tree)
        if b:  # bootstrap
            ph_write(tree, 'tmp.nwk', 'newick')
            msas = bootstrap(aln, b_n)
            k = 0
            for msa in msas:
                dm = calculator.get_distance(msa)
                t = constructor.nj(dm)
                ph_write(t, os.path.join(boot_dir, str(k) + '.nwk'), 'newick')
                k += 1
            L1 = ['sumtrees.py', '-M', '-r', '-F', 'newick', '--suppress-annotation', '-o', 'result.nwk', '-t',
                  'tmp.nwk']
            L2 = [os.path.join(boot_dir, fn) for fn in os.listdir(boot_dir)]
            run(L1 + L2)
            with open('result.nwk') as f, open(outfile, 'a') as g:
                line = f.readline()
                line = line.split()[1]
                g.write(line)
        i += 1
    if not b:
        ph_write(trees, outfile, 'newick')


def build_maj_consensus(infile, outfile, min_freq=0.5):
    trees = dendropy.TreeList()
    trees.read_from_path(infile, 'newick')
    con_tree = trees.consensus(min_freq=min_freq)
    con_tree.write(path=outfile, schema='newick')
    with open(outfile) as f:
        line = f.readline()
        line = line.split()[1]
    with open(outfile, 'w') as f:
        f.write(line)


# without deleting paralogs, for supertree method
def build_trees_all(clusters, proteins_file, trees_outfile, max_seqs=120):
    i = 0
    for cl in clusters:
        sp = ['_'.join(x.split('_')[:2]) for x in cl]
        if len(set(sp)) > 2:
            if len(sp) < max_seqs:
                with open('tmp.fasta', 'w') as f:
                    for record in SeqIO.parse(proteins_file, 'fasta'):
                        if record.id in cl:
                            SeqIO.write(record, f, 'fasta')
                with open('tmp.txt', 'w') as f:
                    run(['mafft', '--anysymbol', 'tmp.fasta'], stdout=f)
                aln = AlignIO.read('tmp.txt', 'fasta')
                calculator = TreeConstruction.DistanceCalculator('blosum62')
                try:
                    dm = calculator.get_distance(aln)
                except ValueError:
                    calculator = TreeConstruction.DistanceCalculator('identity')
                    dm = calculator.get_distance(aln)
                constructor = TreeConstruction.DistanceTreeConstructor()
                tree = constructor.nj(dm)
                for node in tree.get_nonterminals():
                    node.name = None
                for leaf in tree.get_terminals():
                    leaf.name = '_'.join(leaf.name.split('_')[:2])
                ph_write(tree, 'tree.nwk', 'newick')
                with open('tree.nwk') as f, open(trees_outfile, 'a') as g:
                    line = f.readline()
                    g.write(re.sub(':-?[0-9]\.[0-9]*', '', line))  # usuwam długości gałęzi
                i += 1


# set del_lengths=True for using fasturec2
def super_tree(trees_file, out_super, path_to_fasturec, del_lengths=False):
    if del_lengths:
        new_trees_file = trees_file.split('.')[0] + '_no_lengths.' + trees_file.split('.')[1]
        with open(trees_file) as f, open(new_trees_file, 'a') as g:
            for line in f:
                g.write(re.sub(':-?[0-9]\.[0-9]*', '', line))
        trees_file = new_trees_file
    run(['./' + path_to_fasturec, '-Y', '-G', trees_file, '-e', 'a'])
    with open('fu.txt') as f, open(out_super, 'w') as g:
        tree = f.readline().split()[1]
        g.write(tree + ';')


# del_bad - whether to delete nodes with bad support; bad_thr - bad support threshold;
# mb - how many badly supported nodes there can be; alb - min average of nodes' supports
def good_bootstrap_trees(trees_file, outfile, del_bad=False, bad_thr=0.5, mb=8, alb=0.85):
    with open(trees_file) as f:
        trees = f.read().split(';')[:-1]
    for i in range(len(trees)):
        with open('tmp.nwk', 'w') as g:
            g.write(trees[i] + ';')
        L = []
        M = []
        t = Tree('tmp.nwk')
        for n in t.traverse():
            s = n.support
            L.append(s)
            if s < bad_thr:
                M.append(s)
        if len(M) <= mb and sum(L) / len(L) >= alb:  # dobre wartości bootstrapu
            if del_bad:
                for n in t.traverse():
                    s = n.support
                    if s < bad_thr:
                        n.delete()
            print(t.write(format=9), file=open(outfile, 'a'))


def file_exist(path):
    if not os.path.isfile(path):
        raise argparse.ArgumentTypeError(f"File {path} does not exist.")
    return path


parser = argparse.ArgumentParser()
parser.add_argument("i", type=file_exist, help="File with species names (species_name in each line)")
parser.add_argument("-n", "--n_species", type=int, default=30, help="Number of species to draw")
parser.add_argument("--min_proteins", type=int, default=3500, help="Minimum number of proteins in a proteome")
parser.add_argument("--max_proteins", type=int, default=5000, help="Maximum number of proteins in a proteome")
parser.add_argument("--my_species", type=str, default="my_species.txt", help="Name for file to save selected species")
parser.add_argument("--ids_file", type=str, default="proteome_ids.tsv",
                    help="Name for file to save the Uniprot IDs table")
parser.add_argument("--proteins_file", type=str, default="proteins.fasta",
                    help="Name for file to save all protein sequences")
parser.add_argument("--proteomes_dir", type=str, default="protdir", help="Name for directory to save the proteomes")
parser.add_argument("--clustered_dir", type=str, default="clustered",
                    help="Name for directory to save the clustered sequences")
parser.add_argument("--trees_file_c", type=str, default="trees_nj.nwk",
                    help="Name for file to save protein family trees (bijective)")
parser.add_argument("--trees_file_s", type=str, default="sups.nwk",
                    help="Name for file to save protein family trees (with paralogs)")
parser.add_argument("-c", "--c_value", type=float, default=0.5, help="c parameter for MMseqs2")
parser.add_argument("--cov_mode", type=int, default=0, help="cov_mode parameter for MMseqs2")
parser.add_argument("--cluster_mode", type=int, default=0, help="cluster_mode parameter for MMseqs2")
parser.add_argument("-e", "--e_value", type=float, default=0.001, help="e parameter for MMseqs2")
parser.add_argument("--mm_file", type=str, default="mm_05_0", help="Name for file to save MMseqs2 results")
parser.add_argument("--max_seqs", type=int, default=120, help="Maximum number of sequences in a cluster")
parser.add_argument("--min_freq", type=float, default=0.5, help="Minimal frequency of a split to add to consensus tree")
parser.add_argument("--consensus_file", type=str, default="nj_consensus_tree00.nwk",
                    help="Name for file to save consensus tree")
parser.add_argument("--supertree_file1", type=str, default="supertree_basic.nwk",
                    help="Name for file to save basic supertree")
parser.add_argument("--supertree_file2", type=str, default="supertree.nwk",
                    help="Name for file to save supertree (from trees with paralogs)")
parser.add_argument("--boot_file", type=str, default="bootstrapped.nwk",
                    help="Name for file to save trees with bootstrap supports")
parser.add_argument("--boot_dir", type=str, default="bootstraps", help="Name for directory to compute bootstrap trees")
parser.add_argument("--bn", type=int, default=50, help="Number of bootstrap replicates")
parser.add_argument("--good_bootstrap_file", type=str, default="good_bootstrap.nwk",
                    help="Name for file to keep trees with good bootstrap supports")
parser.add_argument("--b_consensus_file", type=str, default="bootstrap_consensus_tree00.nwk",
                    help="Name for file to save bootstrap consensus tree")
parser.add_argument("fasturec", type=file_exist, help="Path to fasturec program")

args = parser.parse_args()

# PART I
download_ids(args.ids_file)
my_species = choose_species(args.ids_file, args.i, args.n_species, args.min_proteins, args.max_proteins,
                            args.my_species)
download_proteomes(my_species, args.proteomes_dir)
combine_files(args.proteomes_dir, args.proteins_file)

# PART II  
mmseqs(args.proteins_file, args.mm_file, args.c_value, args.cov_mode, args.cluster_mode, args.e_value)
L = get_mmseqs_clustering(args.mm_file)

# PART III
write_clusters_consensus(L, args.clustered_dir, args.proteins_file, n=args.n_species, max_seqs=args.max_seqs)

# PART IV & V
build_trees(args.clustered_dir, args.trees_file_c, b=False)

# PART VI
build_maj_consensus(args.trees_file_c, args.consensus_file, args.min_freq)
super_tree(args.trees_file_c, args.supertree_file1, args.fasturec, del_lengths=True)

# PART Va
build_trees(args.clustered_dir, args.boot_file, True, args.boot_dir, b_n=args.bn)
good_bootstrap_trees(args.boot_file, args.good_bootstrap_file, del_bad=False)
# good_bootstrap_trees(args.boot_file, args.good_bootstrap, del_bad=True)  # with deleting badly supported nodes

build_maj_consensus(args.good_bootstrap_file, args.b_consensus_file, args.min_freq)

# PART Vb
build_trees_all(L, args.proteins_file, args.trees_file_s, max_seqs=args.max_seqs)
super_tree(args.trees_file_s, args.supertree_file2, args.fasturec)

