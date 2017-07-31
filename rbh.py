#! /usr/bin/env python

import sys
import math
import itertools
import collections
import argparse

# TODO: incorporate metric option to use either bitscore or evalue.
# TODO: bitscore can be approximated with a robust linear fit (eventually also
# for evalues)
class Hit(object):
    def __init__(self,qname,sname,pctid,length,mismatches,ngaps,
                 qstart,qend,sstart,send,evalue,bitscore,
                 qcov=0.,qlen=0,slen=0):
        self.qname = qname
        self.sname = sname
        self.pctid = float(pctid)
        self.length = int(length)
        self.mismatches = int(mismatches)
        self.ngaps = int(ngaps)
        self.qstart = int(qstart)
        self.qend = int(qend)
        self.sstart = int(sstart)
        self.send = int(send)
        self.evalue = float(evalue)
        self.bitscore = float(bitscore)
        self.qcov = float(qcov)
        self.qlen = int(qlen)
        self.slen = int(slen)
    def cmp_evalue(self,other):
        return cmp(self.evalue,other.evalue)
    def cmp_bitscore(self,other):
        return cmp(other.bitscore,self.bitscore)
    def cmp_sbjct(self,other):
        return cmp(self.sname,other.sname)
    def __cmp__(self,other):
        return self.cmp_evalue(other) or self.cmp_bitscore(other) or \
                self.cmp_sbjct(other)
    def nonidentical_tie(self,other):
        return (self.cmp_evalue(other)  == 0 and 
                self.cmp_bitscore(other) == 0 and
                self.cmp_sbjct(other) != 0)
    def __str__(self):
        return '\t'.join(
            [str(x) for x in [
                self.qname,self.sname,self.pctid,self.length,self.mismatches,
                self.ngaps,self.qstart,self.qend,self.evalue,self.bitscore,
                self.qcov]])

class HitGroup(object):
    def __init__(self):
        self._hits = []
        self._best_hits = []
        self._dirty = False
    def add(self,hit):
        self._hits.append(hit)
        self._dirty = True
    def _update(self):
        if self._dirty:
            self._hits = list(sorted(self._hits))
            self._best_hits = [self._hits[0]] + [
                h for h in self._hits if h.nonidentical_tie(self._hits[0])]
        self._dirty = False
    def _get_hits(self):
        self._update()
        return self._hits
    hits = property(_get_hits)
    def _get_best_hits(self):
        self._update()
        return self._best_hits

    best_hits = property(_get_best_hits)

class HitPair(object):
    def __init__(self,hit1,hit2=None):
        self.hit1 = hit1
        self.hit2 = hit2
    def _qname_id(self):
        return '%s_%s' % (self.hit1.qname,self.hit2.qname)
    def __hash__(self):
        return hash(self._qname_id())
    def __cmp__(self,other):
        return cmp(self._qname_id(),other._qname_id())
    def __str__(self):
        return '\t'.join(str(x) for x in [
            self.hit1.qname,self.hit2.qname,
            self.hit1.evalue,self.hit2.evalue,
            self.hit1.qlen,self.hit2.qlen,
            self.hit1.length])

class CRBFitting(dict):
    def check(self,hit):
        length = min(hit.length,len(self)-1)
        score = CRB.transform_evalue(hit.evalue)
        return self[length] <= score
    
class CRB(object):
    MIN_EVALUE=sys.float_info.min
    def __init__(self,rbh,output_training=None,output_fitting_data=None,output_fitting=None):
        CRB.MIN_EVALUE = min(hit.evalue for hit in
                             itertools.chain((pair.hit1 for pair in rbh),
                                             (pair.hit2 for pair in rbh))
                             if hit.evalue > 0.0)
        self.fitting1 = self.get_fitting(
            set(pair.hit1 for pair in rbh),
            output_training,output_fitting_data,output_fitting)
        self.fitting2 = self.get_fitting(
            set(pair.hit2 for pair in rbh))
    
    @classmethod
    def transform_evalue(cls,evalue):
        return -1*math.log(evalue or CRB.MIN_EVALUE,10)

    def get_fitting(self,training_hits,output_training=None,output_fitting_data=None,output_fitting=None):
        lengths = [h.length for h in training_hits]
        shortest_hit = min(lengths)
        longest_hit = max(lengths)
        evalue_dict = collections.defaultdict(list)
        for hit in training_hits:
            transformed_evalue = CRB.transform_evalue(hit.bitscore)
            evalue_dict[hit.length].append(transformed_evalue)
            if output_training is not None:
                output_training.write("%d\t%.3f\n" %
                                      (hit.length,transformed_evalue))

        fitting = CRBFitting()
        for l in range(shortest_hit,longest_hit+1):
            flank_len = int(min(100,max(l*.1,5)))
            evalues = 0.
            count = 0
            for flank in range(-flank_len,flank_len+1):
                if l+flank in evalue_dict:
                    evalues += sum(evalue_dict[l+flank])
                    count += len(evalue_dict[l+flank])
                    if output_fitting_data is not None:
                        for evalue in evalue_dict[l+flank]:
                            output_fitting_data.write("%d\t%.3f\n" % (l,evalue))
            if count:
                evalue_mean = evalues/count
                fitting[l] = evalue_mean
            if l-1 in fitting and (
                l not in fitting or fitting[l] < fitting[l-1]):
                    fitting[l] = fitting[l-1]
        for l in range(shortest_hit):
            fitting[l] = fitting[shortest_hit]
        if output_fitting is not None:
            for l in range(longest_hit+1):
                output_fitting.write('%d\t%.3f\n' % (l,fitting[l]))
        return fitting

    def _get_crb_pairs(self,blastA,blastB,fittingA,fittingB,exclude_pairs,pair_fac):
        log = lambda x: sys.stderr.write(x)
        pairs = []
        for qname,hit_group in blastA.items():
            for hitA in hit_group.hits:
                if fittingA.check(hitA) and hitA.sname in blastB:
                    for hitB in blastB[hitA.sname].hits:
                        if hitB.sname == qname and fittingB.check(hitB):
                            pair = pair_fac(hitA,hitB)
                            if pair not in exclude_pairs:
                                pairs.append(pair)
        return set(pairs)


    def get_crb(self,blast1,blast2,exclude_pairs=set()):
        pairs = self._get_crb_pairs(
            blast1,blast2,
            self.fitting1,self.fitting2,
            exclude_pairs,lambda a,b: HitPair(a,b)
        ).union(self._get_crb_pairs(
            blast2,blast1,
            self.fitting2,self.fitting1,
            exclude_pairs,lambda a,b: HitPair(b,a)))
        return pairs

class OrthologyGroup(object):
    def __init__(self,species_a_orthologs=[],species_b_orthologs=[]):
        self.a = set(species_a_orthologs)
        self.b = set(species_b_orthologs)
    def add_orthologs(self,species_a_orthologs=[],species_b_orthologs=[]):
        self.a = self.a.union(set(species_a_orthologs))
        self.b = self.b.union(set(species_b_orthologs))
    def update(self,other_orthology_group):
        self.add_orthologs(other_orthology_group.a,other_orthology_group.b)
    def intersect(self,species_a_orthologs=None,species_b_orthologs=None):
        if species_a_orthologs is not None:
            self.a = self.a.intersection(set(species_a_orthologs))
        if species_b_orthologs is not None:
            self.b = self.b.intersection(set(species_b_orthologs))
    def has_single_linkage(self,other):
        if not self.a.isdisjoint(other.a) or not self.b.isdisjoint(other.b):
            return True
        return False
    def all(self):
        return set(self.a).union(self.b)
    def __str__(self):
        return 'a:%s,b:%s' % (' '.join(self.a), ' '.join(self.b))
    def __len__(self):
        return len(self.a) + len(self.b)

class ReciprocalOrthologyGroup(object):
    def __init__(self):
        self.reciprocal_hits = set()
        self.reciprocal_group = OrthologyGroup()
        self.crb_group = OrthologyGroup()
    def add_reciprocal_pair(self,hit1,hit2):
        self.reciprocal_group.add_orthologs([hit1.qname],[hit2.qname])
        self.reciprocal_hits.add(HitPair(hit1,hit2))
    def add_crb_pair(self,hit1,hit2):
        if hit1.qname not in self.reciprocal_group.a:
            self.crb_group.a.add(hit1.qname)
        if hit2.qname not in self.reciprocal_group.b:
            self.crb_group.b.add(hit2.qname)
    def intersect(self,species_a_orthologs=None,species_b_orthologs=None):
        self.reciprocal_group.intersect(species_a_orthologs,species_b_orthologs)
        self.crb_group.intersect(species_a_orthologs,species_b_orthologs)
    def all(self):
        return self.reciprocal_group.all().union(self.crb_group.all())
    def __str__(self):
        return '%s\t%s\t%s\t%s' % (
            ' '.join(sorted(self.reciprocal_group.a)),
            ' '.join(sorted(self.crb_group.a)),
            ' '.join(sorted(self.reciprocal_group.b)),
            ' '.join(sorted(self.crb_group.b))
        )

class GeneSet(object):
    def __init__(self,species_name,genes):
        self.species_name = species_name
        self.total_genes = set(genes)
        self.rbh_genes = set()
        self.crb_genes = set()

class RBHSummary(object):
    def __init__(self,species_a,species_b,total_genes_a,total_genes_b,
                 ortholog_set):
        self.a = GeneSet(species_a,total_genes_a)
        self.b = GeneSet(species_b,total_genes_b)
        self.rbh_groups = 0
        self.crb_groups = 0
        for ortholog in ortholog_set:
            self.a.rbh_genes.update(ortholog.reciprocal_group.a)
            self.b.rbh_genes.update(ortholog.reciprocal_group.b)
            self.a.crb_genes.update(ortholog.crb_group.a)
            self.b.crb_genes.update(ortholog.crb_group.b)
            if ortholog.reciprocal_group.a and ortholog.reciprocal_group.b:
                self.rbh_groups += 1
            if ortholog.crb_group.a or ortholog.crb_group.b:
                self.crb_groups += 1
        self.a.crb_genes.difference_update(self.a.rbh_genes)
        self.b.crb_genes.difference_update(self.a.rbh_genes)

    def summary_header(self):
        return '\t%s Genes\t%s Genes\tOrthology Groups\n' % (self.a.species_name,self.b.species_name)
    def summary_table(self,prefix=''):
        return '%sGenes\t%d\t%d\t-\n%sRBB\t%d\t%d\t%d\n%sCRB\t%d\t%d\t%d\n' % (
            prefix,len(self.a.total_genes),len(self.b.total_genes),
            prefix,len(self.a.rbh_genes),len(self.b.rbh_genes),self.rbh_groups,
            prefix,len(self.a.crb_genes),len(self.b.crb_genes),self.crb_groups)

class RBHArgumentParser(argparse.ArgumentParser):
    def error(self,message):
        sys.stderr.write('error: %s\n' % message)
        self.print_help()
        sys.exit(2)

class RBH(object):
    def __init__(self):
        self.parse_args()
        self.blast1 = self.qualifying_hits(self.args.BLAST_FILE1)
        self.blast2 = self.qualifying_hits(self.args.BLAST_FILE2)

    '''
    def _get_summary(self,ortholog_set,total_genes1,total_genes2):
        summary = RBHSummary(self.args.SPECIES_NAME1,self.args.SPECIES_NAME2)
        summary.a.total_genes.update(total_genes1)
        summary.b.total_genes.update(total_genes2)
        for ortholog in ortholog_set:
            summary.a.rbh_genes.update(ortholog.reciprocal_group.a)
            summary.b.rbh_genes.update(ortholog.reciprocal_group.b)
            summary.a.crb_genes.update(ortholog.crb_group.a)
            summary.b.crb_genes.update(ortholog.crb_group.b)
        summary.a.crb_genes.difference_update(summary.a.rbh_genes)
        summary.b.crb_genes.difference_update(summary.b.rbh_genes)
        return summary
    '''

    def _get_rbh(self):
        orthologs = collections.defaultdict(ReciprocalOrthologyGroup)
        for qname in sorted(self.blast1.keys()):
            hits = self.blast1[qname]
            for hit in hits.best_hits:
                recip_qname = self.process_name(hit.sname)
                if recip_qname in self.blast2:
                    for recip_hit in self.blast2[recip_qname].best_hits:
                        if self.process_name(recip_hit.sname) == qname:
                            orthologs[qname].add_reciprocal_pair(hit,recip_hit)

        for q1,q2 in itertools.combinations(orthologs.keys(),r=2):
            if orthologs[q1].reciprocal_group.has_single_linkage(
                orthologs[q2].reciprocal_group):
                orthologs[q1].reciprocal_hits.update(
                    orthologs[q2].reciprocal_hits)
                orthologs[q1].reciprocal_group.update(
                    orthologs[q2].reciprocal_group)
                orthologs[q2] = orthologs[q1]

        return set(orthologs.values())

    def _get_crb(self,ortholog_set):
        reciprocal_hits = set(itertools.chain.from_iterable(
            ortholog.reciprocal_hits for ortholog in ortholog_set))
            
        crb_finder = CRB(reciprocal_hits,
                         self.args.output_evalue_table,
                         self.args.output_fit_table,
                         self.args.output_fit_values)
        crb_pairs = crb_finder.get_crb(self.blast1,self.blast2,reciprocal_hits)

        orthologs = {}
        for ortholog in ortholog_set:
            for name in ortholog.reciprocal_group.all():
                orthologs[name] = ortholog

        for pair in crb_pairs:
            if pair.hit1.qname in orthologs:
                orthologs[pair.hit1.qname].add_crb_pair(pair.hit1,pair.hit2)
            if pair.hit2.qname in orthologs:
                orthologs[pair.hit2.qname].add_crb_pair(pair.hit1,pair.hit2)
            if not any(hit.qname in orthologs for hit in [pair.hit1,pair.hit2]):
                ortholog = ReciprocalOrthologyGroup()
                ortholog.add_crb_pair(pair.hit1,pair.hit2)
                orthologs[pair.hit1.qname] = ortholog
                orthologs[pair.hit2.qname] = ortholog
        return ortholog_set

    def _filter_targets(self,ortholog_set):
        if self.args.target_list1:
            list1_orthologs = set(
                group for group in ortholog_set 
                if not group.all().isdisjoint(self.args.target_list1))
        else:
            list1_orthologs = set(ortholog_set)

        if self.args.target_list2:
            list2_orthologs = set(
                group for group in ortholog_set 
                if not group.all().isdisjoint(self.args.target_list2))
        else:
            list2_orthologs = dict(orthologs)

        if self.args.target_list1 and self.args.target_list2:
            if self.args.target_list_intersection:
                ortholog_set = list1_orthologs.intersection(list2_orthologs)
            else:
                ortholog_set = list1_orthologs.union(list2_orthologs)
        else:
            # the one without the list is just the full list -- just take
            # the intersection
            ortholog_set = list1_orthologs.intersection(list2_orthologs)


        if self.args.target_list_intersection:
            if self.args.target_list1:
                for group in ortholog_set:
                    group.intersect(species_a_orthologs=self.args.target_list1)
            if self.args.target_list2:
                for group in ortholog_set:
                    group.intersect(species_b_orthologs=self.args.target_list2)

        return ortholog_set

    def run(self):
        ortholog_set = self._get_rbh()

        if self.args.crb:
            ortholog_set = self._get_crb(ortholog_set)

        if self.args.output_summary:
            summary_prefilter = RBHSummary(
                self.args.SPECIES_NAME1,self.args.SPECIES_NAME2,
                self.blast1.keys(), self.blast2.keys(), ortholog_set)
            self.args.output_summary.write(summary_prefilter.summary_header())
            self.args.output_summary.write(summary_prefilter.summary_table())
            
        if self.args.target_list1 or self.args.target_list2:
            ortholog_set = self._filter_targets(ortholog_set)

            if self.args.output_summary:
                summary_postfilter = RBHSummary(
                    self.args.SPECIES_NAME1,self.args.SPECIES_NAME2,
                    self.args.target_list1,self.args.target_list2,ortholog_set)
                self.args.output_summary.write(
                    summary_postfilter.summary_table(prefix='Target '))

        if not self.args.silent:
            if not self.args.no_header:
                print '%s RBB Genes\t%s CRB Genes\t%s RBB Genes\t%s CRB Genes' % (
                    self.args.SPECIES_NAME1,self.args.SPECIES_NAME1,
                    self.args.SPECIES_NAME2,self.args.SPECIES_NAME2)
            print '\n'.join(str(x) for x in ortholog_set)
        return 

    def parse_args(self):
        parser = RBHArgumentParser(
            description='''Find reciprocal best BLAST hits based on BLAST 
            output of queries whole proteomes. Expected inputs are NCBI-BLAST
            tables.''')
        setfile = lambda fn: set([l.strip() for l in open(fn)])
        parser.add_argument('SPECIES_NAME1')
        parser.add_argument('SPECIES_NAME2')
        parser.add_argument('BLAST_FILE1')
        parser.add_argument('BLAST_FILE2')
        parser.add_argument('--max-evalue',metavar='E',type=float,default=1e-5,
                            help='maximum evalue per hit [default: 1e-5]')
        parser.add_argument('--min-qcov',metavar='C',type=float,default=0.,
                            help='minimum query coverage per hit [default: 0]')
        parser.add_argument('--crb',action='store_true',
                            help='Add conditional best reciprocal blast hits')
        parser.add_argument('--target-list-intersection',action='store_true',
                            help='''Report the intersection of the
                            riciprocal/conditional target hits lists rather 
                            than the union''')
        parser.add_argument('--target-list1',metavar='LIST',type=setfile,
                           help='''Output hits containing queries from
                           BLAST_FILE1 in this list''')
        parser.add_argument('--target-list2',metavar='LIST',type=setfile,
                           help='''Output hits containing queries from
                           BLAST_FILE2 in this list''')
        parser.add_argument('--output-summary',metavar='FILE',
                            type=argparse.FileType('w'),help='''Output summary
                            statistics for the number of RBB and CRB genes in a
                            tab-delimited table''')
        parser.add_argument('--silent',action='store_true',help='''Do not output
                            hits''')
        parser.add_argument('--no-header',action='store_true',help='''Do not
                            inlcude a header line in the output''')
        # hidden arguments
        '''Output all evalue and length values used for training the CRB hits'''
        parser.add_argument('--output-evalue-table',metavar='TAB',
                            type=argparse.FileType('w'),help=argparse.SUPPRESS)
        '''Output all training data including data re-binned into neighboring 
        length categories.'''
        parser.add_argument('--output-fit-table',metavar='TAB',
                            type=argparse.FileType('w'),help=argparse.SUPPRESS)
        '''Output fitting function'''
        parser.add_argument('--output-fit-values',metavar='VAL',
                            type=argparse.FileType('w'),help=argparse.SUPPRESS)

        # XXX: this does not increase RBH and is not implemented in CRB so 
        # for now we will exclude this
        #parser.add_argument('--trinity-gene-prefix',metavar='PREFIX',
                            #help='''assume genes which begin with PREFIX are
                            #trinity genes, and trim the isoform ID off when
                            #searching for reciprocal hits.''')
        self.args = parser.parse_args()

    def parse_hit_table(self, file):
        for lineno,line in enumerate(open(file)):
            if line.startswith('#'):
                continue
            fields = line.split()
            if len(fields) < 12 or len(fields) > 15:
                raise ValueError, 'Could not parse BLAST hit on line %d of %s' % (
                    lineno,file)
            try:
                yield Hit(*fields)
            except:
                print fields
                raise

    def qualifying_hit(self, hit):
        if hit.evalue <= self.args.max_evalue and \
           hit.qcov >= self.args.min_qcov:
            return True

    def process_name(self,name):
        #if self.args.trinity_gene_prefix is not None \
           #and name.startswith(self.args.trinity_gene_prefix):
            #return '_'.join(name.split('_')[:-1])
        #else:
        return name

    def qualifying_hits(self,file):
        hits = collections.defaultdict(HitGroup)
        for hit in self.parse_hit_table(file):
            if self.qualifying_hit(hit):
                hits[self.process_name(hit.qname)].add(hit)
        return hits


if __name__ == '__main__':
    rbh = RBH()
    rbh.run()
