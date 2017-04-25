#! /usr/bin/env python

import sys
import math
import itertools
import collections
import argparse

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
        self._winners = []
        self._dirty = False
    def add(self,hit):
        self._hits.append(hit)
        self._dirty = True
    def _update(self):
        if self._dirty:
            self._hits = list(sorted(self._hits))
            self._winners = [self._hits[0]] + [
                h for h in self._hits if h.nonidentical_tie(self._hits[0])]
        self._dirty = False
    def _get_hits(self):
        self._update()
        return self._hits
    hits = property(_get_hits)
    def _get_winners(self):
        self._update()
        return self._winners

    winners = property(_get_winners)

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
        return self[length] <= CRB.transform_evalue(hit.evalue)
    
class CRB(object):
    MIN_EVALUE=sys.float_info.min
    def __init__(self,rbh,output_training=None,output_fitting_data=None,output_fitting=None):
        CRB.MIN_EVALUE = min(hit.evalue for hit in
                             itertools.chain((pair.hit1 for pair in rbh),
                                             (pair.hit2 for pair in rbh))
                             if hit.evalue > 0.0)
        self.fitting1 = self.get_fitting([pair.hit1 for pair in rbh],
                                         output_training,output_fitting_data,output_fitting)
        self.fitting2 = self.get_fitting([pair.hit2 for pair in rbh])#,
                                         #output_training,output_fitting)
    
    @classmethod
    def transform_evalue(cls,evalue):
        return -1*math.log(evalue or CRB.MIN_EVALUE,10)

    def get_fitting(self,training_hits,output_training=None,output_fitting_data=None,output_fitting=None):
        lengths = [h.length for h in training_hits]
        shortest_hit = min(lengths)
        longest_hit = max(lengths)
        evalue_dict = collections.defaultdict(list)
        for hit in training_hits:
            transformed_evalue = CRB.transform_evalue(hit.evalue)
            evalue_dict[hit.length].append(transformed_evalue)
            if output_training is not None:
                output_training.write("%d\t%.3f\n" %
                                      (hit.length,transformed_evalue))

        fitting = CRBFitting()
        for l in range(shortest_hit,longest_hit+1):
            flank_len = int(max(l*.1,5))
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

    def get_crb(self,blast1,blast2,exclude_pairs=set()):
        pairs = set()
        for qname,hit_group in blast1.items():
            for hit1 in hit_group.hits:
                if self.fitting1.check(hit1) and hit1.sname in blast2:
                    for hit2 in blast2[hit1.sname].hits:
                        if hit2.sname == qname and  self.fitting2.check(hit2):
                            pair = HitPair(hit1,hit2)
                            if pair not in exclude_pairs and pair not in pairs:
                                pairs.add(pair)
        return pairs

class HitList(object):
    def __init__(self):
        self.reciprocal_pair = None
        self.crb_pairs = []
    def _get_query(self):
        if self.reciprocal_pair is not None:
            return self.reciprocal_pair.hit1.qname
        elif len(self.crb_pairs):
            return self.crb_pairs[0].hit1.qname
        else:
            return ''
    query = property(_get_query)

    def all(self):
        if self.reciprocal_pair is not None:
            return [self.reciprocal_pair] + self.crb_pairs
        else:
            return self.crb_pairs

    def __str__(self):
        result = '%s\t' % self.query
        if self.reciprocal_pair is not None:
            result += self.reciprocal_pair.hit1.sname
        if len(self.crb_pairs):
            result += '\t' + '\t'.join(pair.hit1.sname for pair in self.crb_pairs)
        return result

class RBH(object):
    def __init__(self):
        self.parse_args()
        self.blast1 = self.qualifying_hits(self.args.BLAST_FILE1)
        self.blast2 = self.qualifying_hits(self.args.BLAST_FILE2)

    def run(self):
        orthologs = collections.defaultdict(HitList)
        for qname in sorted(self.blast1.keys()):
            hits = self.blast1[qname]
            for winner in hits.winners:
                rname = self.process_name(winner.sname)
                if winner.sname in self.blast2:
                    for rwinner in self.blast2[rname].winners:
                        if self.process_name(rwinner.sname) == qname:
                            orthologs[qname].reciprocal_pair = HitPair(winner,rwinner)
        if self.args.crb:
            reciprocal_hits = set(ortholog.reciprocal_pair for ortholog in orthologs.values())
            crb_finder = CRB(reciprocal_hits,
                             self.args.output_evalue_table,
                             self.args.output_fit_table,
                             self.args.output_fit_values)
            crb_pairs = crb_finder.get_crb(self.blast1,self.blast2,reciprocal_hits)
            for pair in crb_pairs:
                orthologs[pair.hit1.qname].crb_pairs.append(pair)

        if self.args.target_list1:
            orthologs = {q:hl for q,hl in orthologs.items() 
                         if q in self.args.target_list1}
        if self.args.target_list2:
            orthologs = {q:hl for q,hl in orthologs.items()
                         if any(hit_pair.hit2.qname in self.args.target_list2
                               for hit_pair in hl.all()) }
            # XXX: filtering out the list should be more elegant than this
            for hit_list in orthologs.values():
                if hit_list.reciprocal_pair is not None and \
                   hit_list.reciprocal_pair.hit2.qname \
                   not in self.args.target_list2:
                    hit_list.reciprocal_pair = None
                hit_list.crb_pairs = [
                    hit_pair for hit_pair in hit_list.crb_pairs
                    if hit_pair.hit2.qname in self.args.target_list2]

        print '\n'.join(str(x) for x in orthologs.values())
        return 

    def parse_args(self):
        parser = argparse.ArgumentParser(
            description='''Find reciprocal best BLAST hits based on BLAST 
            output of queries whole proteomes. Expected inputs are NCBI-BLAST
            tables.''')
        listfile = lambda fn: [l.strip() for l in open(fn)]
        parser.add_argument('BLAST_FILE1')
        parser.add_argument('BLAST_FILE2')
        parser.add_argument('--max-evalue',metavar='E',type=float,default=1e-5,
                            help='maximum evalue per hit [default: 1e-5]')
        parser.add_argument('--min-qcov',metavar='C',type=float,default=0.,
                            help='minimum query coverage per hit [default: 0]')
        parser.add_argument('--crb',action='store_true',
                            help='Add conditional best reciprocal blast hits')
        parser.add_argument('--target-list1',metavar='LIST',type=listfile,
                           help='''Output hits containing queries from
                           BLAST_FILE1 in this list''')
        parser.add_argument('--target-list2',metavar='LIST',type=listfile,
                           help='''Output hits containing queries from
                           BLAST_FILE2 in this list''')
        parser.add_argument('--output-evalue-table',metavar='TAB',
                            type=argparse.FileType('w'),help='''Output all
                            evalue and length values used for training the CRB
                            hits''')
        parser.add_argument('--output-fit-table',metavar='TAB',
                            type=argparse.FileType('w'),help='''Output all
                            training data including data re-binned into
                            neighboring length categories.''')
        parser.add_argument('--output-fit-values',metavar='VAL',
                            type=argparse.FileType('w'),help='''Output fitting
                            function''')

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
