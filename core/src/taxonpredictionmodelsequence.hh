/*
taxator-tk predicts the taxon for DNA sequences based on sequence alignment.

Copyright (C) 2021 Johannes Dröge, Eik Dahms

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.

*/

#ifndef taxonpredictionmodelsequence_hh_
#define taxonpredictionmodelsequence_hh_

#include <seqan/align.h>
#include <seqan/graph_msa.h>
#include <seqan/file.h>
#include <seqan/basic.h>
#include <boost/tuple/tuple.hpp>
#include <boost/tuple/tuple_comparison.hpp>
#include <boost/format.hpp>
#include <assert.h>
#include <limits>
#include <set>
#include <ostream>
#include <cmath>
#include "taxonpredictionmodel.hh"
#include "sequencestorage.hh"
#include "profiling.hh"

// alignment object
// We represent each alignment with both, a distance and a similarity metric
// because some alignment algorithms will report a distance while others will
// report a similarity score. We normalize both to add up to the length of the
// calculated alignment (unit fractional base pair)
// TODO: use AlignmentStats object in Seqan
template <typename StringType>
struct Alignment {
  unsigned int matches; // only exact positional matches
  unsigned int mismatches;
  unsigned int gaps;
  unsigned int length;
  float distance;
  float similarity;
  seqan::Align<StringType, seqan::ArrayGaps> alignment;
};

template <typename StringType>
std::ostream& operator<<(std::ostream& os, const Alignment<StringType>& aln) {
  if(seqan::length(seqan::rows(aln.alignment)) > 0) { // ugly hack because otherwise outstream fails on empty object
    os << aln.alignment;
  }
  return os;
}

// The MyersHirschberg implementation is working but currently disabled, because
// it takes about 2 to 3 times as long as the MyersBitVector implementation
Alignment<seqan::String<seqan::Dna5>> getAlignmentDNAExact(const seqan::String<seqan::Dna5>& A, const seqan::String<seqan::Dna5>& B) {
  typedef seqan::String<seqan::Dna5> StringType;
  // typedef seqan::EditDistanceScore ScoringScheme;
  typedef seqan::MyersHirschberg AlignmentAlgorithm;
  typedef seqan::Align<StringType, seqan::ArrayGaps> TAlign;
  typedef seqan::Row<TAlign>::Type TRow;
  typedef seqan::Iterator<TRow>::Type TRowIterator;

  // instantiate static objects
  auto alignAlgo = AlignmentAlgorithm();
  // auto alignScoring = AlignmentScoring();

  // align sequence A to B
  TAlign alignAB;
  resize(rows(alignAB), 2);
  assignSource(row(alignAB, 0), A);
  assignSource(row(alignAB, 1), B);

  // alignment with traceback for alignment statistics
  int mutualscore = -seqan::globalAlignment( // make edit distance positive
    alignAB,
    // alignScoring,
    alignAlgo
  );
  assert(mutualscore >= 0 && "distance metric cannot be negative");

  // extract stats from mutual alignment
  TRow & row1 = row(alignAB, 0);
  TRow & row2 = row(alignAB, 1);

  unsigned int gap = 0;
  unsigned int match = 0;
  unsigned int mismatch = 0;

  TRowIterator itRow1 = begin(row1);
  TRowIterator itEndRow1 = end(row1);
  TRowIterator itRow2 = begin(row2);
  for (; itRow1 != itEndRow1; ++itRow1, ++itRow2) {
    if (seqan::isGap(itRow1) || seqan::isGap(itRow2)) {
      gap++;
    } else if(*itRow1 == *itRow2){
      match++;
    } else{
      mismatch++;
    }
  }

  assert(static_cast<unsigned int>(mutualscore) == gap + mismatch  && "edit distance does not equal gaps+mismatches");

  // construct return object
  Alignment<StringType> aln;
  aln.distance = static_cast<float>(mutualscore);
  aln.similarity = match;
  aln.matches = match;
  aln.mismatches = mismatch;
  aln.gaps = gap;
  aln.length = match + mismatch + gap; // TODO: remove duplicate information
  aln.alignment = alignAB;

  return aln;
}

// the MyersBitVector implementation for DNA sequences is fastest but only uses
// lower bound estimates for the number of matches (which is good enough for
// most applications)
// template<typename ContainerT, typename QStorType, typename DBStorType, typename StringType>
Alignment<seqan::String<seqan::Dna5>> getAlignmentDNA(const seqan::String<seqan::Dna5>& A, const seqan::String<seqan::Dna5>& B) {
  typedef seqan::String<seqan::Dna5> StringType;
  typedef seqan::MyersBitVector AlignmentAlgorithm;

  // instantiate static objects
  auto alignAlgo = AlignmentAlgorithm();
  // auto alignScoring = AlignmentScoring();

  // fix long and short sequence to make things easier
  const StringType* long_seq = &A;
  const StringType* short_seq = &B;
  if(seqan::length(A) < seqan::length(B)) {
    long_seq = &B;
    short_seq = &A;
  }

  // alignment without traceback using approximated alignment statistics
  int mutualscore = -seqan::globalAlignmentScore(*short_seq, *long_seq, alignAlgo); // make edit distance positive
  assert(mutualscore >= 0 && "distance metric cannot be negative");

  // approximate required statistics (lower bounds on matches)
  int gap_or_mismatch = mutualscore;
  int lendiff = seqan::length(*long_seq) - seqan::length(*short_seq);
  assert(lendiff <= gap_or_mismatch && "alignment score ignores some positions");
  int gap = lendiff; // lower bound
  int mismatch = gap_or_mismatch - lendiff; // upper bound
  int match = seqan::length(*short_seq) - mismatch; // lower bound

  // construct return object
  Alignment<StringType> aln;
  aln.distance = mutualscore;
  aln.similarity = match;
  aln.matches = match;
  aln.mismatches = mismatch;
  aln.gaps = gap;
  aln.length = match + mismatch + gap; // TODO: remove duplicate information

  return aln;
}

Alignment<seqan::String<seqan::AminoAcid>> getAlignmentProtein(const seqan::String<seqan::AminoAcid>& A, const seqan::String<seqan::AminoAcid>& B) {
  typedef seqan::String<seqan::AminoAcid> StringType;
  typedef seqan::Blosum62 AlignmentScoring;
  typedef seqan::LinearGaps AlignmentAlgorithm; // fast
  // typedef seqan::DynamicGaps AlignmentAlgorithm; // medium
  // typedef seqan::AffineGaps AlignmentAlgorithm; // slow
  // typedef seqan::EditDistanceScore ScoringScheme;
  // typedef typename seqan::MyersHirschberg AlignmentAlgorithm;
  typedef typename seqan::Align<StringType, seqan::ArrayGaps> TAlign;
  typedef typename seqan::Row<TAlign>::Type TRow;
  typedef typename seqan::Iterator<TRow>::Type TRowIterator;

  // instantiate static objects
  auto alignAlgo = AlignmentAlgorithm();
  auto alignScoring = AlignmentScoring();

  // TODO: instead use only scoring matrix to infer selfscore to save runtime
  int selfscore = seqan::globalAlignmentScore(A, A, alignScoring, alignAlgo) +
                  seqan::globalAlignmentScore(B, B, alignScoring, alignAlgo);

  // align sequence A to B
  TAlign alignAB;
  resize(rows(alignAB), 2);
  assignSource(row(alignAB, 0), A);
  assignSource(row(alignAB, 1), B);

  int mutualscore = seqan::globalAlignment(
    alignAB,
    alignScoring,
    alignAlgo
  );

  assert(selfscore >= mutualscore && "sequence self comparision must yield highest possible score");

  // extract stats from mutual alignment
  TRow & row1 = row(alignAB, 0);
  TRow & row2 = row(alignAB, 1);

  int gap = 0;
  int match = 0;
  int mismatch = 0;

  TRowIterator itRow1 = begin(row1);
  TRowIterator itEndRow1 = end(row1);
  TRowIterator itRow2 = begin(row2);
  for (; itRow1 != itEndRow1; ++itRow1, ++itRow2) {
    if (seqan::isGap(itRow1) || seqan::isGap(itRow2)) {
      gap++;
    } else if(*itRow1 == *itRow2){
      match++;
    } else{
      mismatch++;
    }
  }
  unsigned int len = (gap + match + mismatch);
  float normfactor = len/static_cast<float>(selfscore);

  // construct return object
  Alignment<StringType> aln;
  aln.distance = (selfscore - 2*mutualscore) * normfactor; // simple symmetric scoring formula
  aln.similarity = (2*mutualscore) * normfactor;
  aln.matches = match;
  aln.mismatches = mismatch;
  aln.gaps = gap;
  aln.length = len; // TODO: remove duplicate information
  aln.alignment = alignAB;

  assert(aln.distance >= 0 && "distance metric cannot be negative");
  return aln;
}

template<class T>
Alignment<T> getAlignment(const T& A, const T& B);

template<>
Alignment<seqan::String<seqan::Dna5>> getAlignment<seqan::String<seqan::Dna5>>(const seqan::String<seqan::Dna5>& A, const seqan::String<seqan::Dna5>& B) {
  // return getAlignmentDNAExact(A, B);
  return getAlignmentDNA(A, B);
}

template<>
Alignment<seqan::String<seqan::AminoAcid>> getAlignment<seqan::String<seqan::AminoAcid>>(const seqan::String<seqan::AminoAcid>& A, const seqan::String<seqan::AminoAcid>& B) {
  return getAlignmentProtein(A, B);
}

// helper class
class BandFactor {
public:
  BandFactor(TaxonomyInterface& taxinter, uint reserve = 0) :
  bandfactor_(-1),
  taxinter_(taxinter) {
    if(reserve) data_.reserve(reserve);
  }

  void addSequence(const float score, const TaxonNode* node) {
    data_.push_back(boost::make_tuple(score, node));
  }

  float getFactor() {
    if(bandfactor_ < 0) {
      sort();
      setBandFactor();
    }
    return sqrt(bandfactor_);
  }

private:
  void setBandFactor(const float min_bandfactor = 1., const float max_bandfactor = std::numeric_limits< float >::max()) { //data_ must be sorted TODO: optimize
    float bandfactor = min_bandfactor;
    float score;
    const TaxonNode* anchor;
    const TaxonNode* node;
    std::map< small_unsigned_int, float > worstscore_per_rank;
    boost::tie(score, anchor) = data_[0];
    small_unsigned_int last_rank = anchor->data->root_pathlength;
    worstscore_per_rank[ last_rank ] = score; //TODO: change name

    for (uint i = 1; i < data_.size(); ++i) {
      boost::tie(score, node) = data_[i];
      node = taxinter_.getLCA(node, anchor);
      small_unsigned_int rank = node->data->root_pathlength;
      if (rank == last_rank) {
      } else if (rank < last_rank) {
        worstscore_per_rank[rank] = score;
        last_rank = rank;
      } else { //disorder
        float refscore;
        small_unsigned_int r = rank - 1;
        do {
          auto it = worstscore_per_rank.find(r);
          if (it != worstscore_per_rank.end()) {
            refscore = it->second;
            if (refscore) bandfactor = std::max(bandfactor, score/refscore);
          }
        } while (r--);
      }
    }
    bandfactor_ = std::min(bandfactor, max_bandfactor);
  }

  void sort() { //sort by increasing order
    auto start_it = data_.begin();
    std::sort(++start_it, data_.end(), comparator_);
  }

  compareTupleFirstLT< boost::tuple< float, const TaxonNode* >, 0 > comparator_;
  typedef std::vector< boost::tuple< float, const TaxonNode* > > data_type_;
  float bandfactor_;
  data_type_ data_;
  TaxonomyInterface taxinter_;
};

// TODO: make timers thread-safe
template< typename ContainerT, typename QStorType, typename DBStorType, typename StringType>
class RPAPredictionModel : public TaxonPredictionModel< ContainerT > {
public:
  RPAPredictionModel(const Taxonomy* tax, QStorType& q_storage, const DBStorType& db_storage, float exclude_factor ,float reeval_bandwidth = .1) :
  TaxonPredictionModel< ContainerT >(tax),
  query_sequences_(q_storage),
  db_sequences_(db_storage),
  exclude_alignments_factor_(exclude_factor),
  reeval_bandwidth_factor_(1. - reeval_bandwidth),
  measure_sequence_retrieval_("sequence retrieval using index"),
  measure_pass_0_alignment_("best reference re-evaluation alignments (pass 0)"),
  measure_pass_1_alignment_("best reference anchor alignments (pass 1)"),
  measure_pass_2_alignment_("distant anchor alignments (pass 2)")
  {};

  void predict(ContainerT& recordset, PredictionRecord& prec, std::ostream& logsink) {
    this->initPredictionRecord(recordset, prec);  // set query name and length
    const std::string& qid = prec.getQueryIdentifier();
    boost::format seqname_fmt("%d:%d@%s");  // local variable because not thread-safe
    StopWatchCPUTime stopwatch_init("initializing this record");  // log overall time for this predict phase
    stopwatch_init.start();

    // push records into active_records  TODO: remove intermediate active_records?
    active_list_type_ active_records;
    {
      for(typename ContainerT::iterator rec_it = recordset.begin(); rec_it != recordset.end(); ++rec_it) {
        if(!(*rec_it)->isFiltered()) active_records.push_back(*rec_it);
      }
    }
    const uint n = active_records.size();

    // with no unmasked alignment, set to unclassified and return
    if(n==0) {  //TODO: record should not be reported at all in GFF3
      const std::string qrseqname = boost::str(seqname_fmt % -1 % -1 % qid);
      logsink << "ID" << tab << qrseqname << std::endl;
      logsink << "  NUMREF" << tab << n << std::endl << std::endl;
      logsink << "    RANGE\t" << this->taxinter_.getRoot()->data->annotation->name << tab << this->taxinter_.getRoot()->data->annotation->name << tab << this->taxinter_.getRoot()->data->annotation->name << std::endl << std::endl;
      logsink << "STATS" << tab << qrseqname << tab << n << "\t0\t0\t0\t0\t0\t0\t0\t.0" << std::endl << std::endl;

      TaxonPredictionModel< ContainerT >::setUnclassified(prec);
      return;
    }

    // with one alignment, don't align and return
    if(n==1) {
      typename ContainerT::value_type rec = active_records.front();
      large_unsigned_int qrstart = rec->getQueryStart();
      large_unsigned_int qrstop = rec->getQueryStop();
      const std::string qrseqname = boost::str(seqname_fmt % qrstart % qrstop % qid);

      logsink << "ID" << tab << qrseqname << std::endl;
      logsink << "  NUMREF" << tab << n << std::endl;
      logsink << "  RANGE\t" << rec->getReferenceNode()->data->annotation->name << tab << rec->getReferenceNode()->data->annotation->name << tab << this->taxinter_.getRoot()->data->annotation->name << std::endl << std::endl;
      logsink << "STATS" << tab << qrseqname << tab << n << "\t0\t0\t0\t0\t0\t0\t0\t.0" << std::endl << std::endl;

      prec.setQueryFeatureBegin(qrstart);
      prec.setQueryFeatureEnd(qrstop);
      prec.setInterpolationValue(1.);
      prec.setNodeRange(rec->getReferenceNode(), this->taxinter_.getRoot(), rec->getIdentities());
      prec.setBestReferenceTaxon(rec->getReferenceNode());
      return;
    }

    // n>1 -> screen alignments and determine query range
    large_unsigned_int qrstart;
    large_unsigned_int qrstop;
    {
      assert(n>1);
      typename active_list_type_::iterator rec_it = active_records.begin();
      qrstart = (*rec_it)->getQueryStart();
      qrstop = (*rec_it)->getQueryStop();
      ++rec_it;
      do {  // because n >= 2
        qrstart = std::min((*rec_it)->getQueryStart(), qrstart);
        qrstop = std::max((*rec_it)->getQueryStop(), qrstop);
      } while(++rec_it != active_records.end());
    }
    const large_unsigned_int qrlength = qrstop - qrstart + 1;

    // logging
    const std::string qrseqname = boost::str(seqname_fmt % qrstart % qrstop % qid);
    logsink << "ID" << tab << qrseqname << std::endl;
    logsink << "  NUMREF" << tab << n << std::endl;

    // sort the list by score
    sort_.filter(active_records);

    // data storage  TODO: maybe use Boost ptr containers
    const StringType qrseq = query_sequences_.getSequence(qid, qrstart, qrstop);

    std::vector< typename ContainerT::value_type > records(n);  //TODO: move below next section and do not create records if q==r_best
    {
      typename active_list_type_::iterator rec_it = active_records.begin();
      uint i = 0;
      while(rec_it != active_records.end()) {
        records[i] = *rec_it;
        ++rec_it;
        ++i;
      }
    }
    const float qmax_searchscore = records[0]->getScore();

    // n>1 and query is identical to reference, we will use local alignment scores only
    // TODO: add option to recalculate (github issue #24)
    if(records[0]->getAlignmentLength() == qrlength && records[0]->getIdentities() == qrlength) {
      float searchscore_best = records[0]->getScore();
      const TaxonNode* lnode = records[0]->getReferenceNode();
      const TaxonNode* unode = nullptr;
      uint i = 1;

      while (true) {  // 2 breaks
        if(i == n) {  // 1 break
          unode = this->taxinter_.getRoot();
          break;
        }

        float searchscore = records[i]->getScore();

        if(searchscore == searchscore_best) {
          const TaxonNode* cnode = records[i]->getReferenceNode();
          lnode = this->taxinter_.getLCA(lnode, cnode);
          logsink << "    current ref/lower node: " << "("<< searchscore <<") "<< lnode->data->annotation->name << " (+ " << cnode->data->annotation->name << " )" << std::endl;
        }
        else {  // 1 break
          float uscore = searchscore;
          unode = lnode;
          do {
            const TaxonNode* cnode = records[i]->getReferenceNode();
            unode = this->taxinter_.getLCA(unode, cnode);
            logsink << "    current upper node: " << "("<< uscore <<") "<< unode->data->annotation->name << " (+ " << cnode->data->annotation->name << " at " << static_cast<int>(this->taxinter_.getLCA(cnode, lnode)->data->root_pathlength) << " )" << std::endl;
          } while (++i < n && records[i]->getScore() == uscore);
          break;
        }
        ++i;
      }

      logsink << "  RANGE\t" << lnode->data->annotation->name << tab << lnode->data->annotation->name << tab << unode->data->annotation->name << std::endl << std::endl;
      logsink << "STATS" << tab << qrseqname << tab << n << "\t0\t0\t0\t0\t" << stopwatch_init.read() << "\t0\t0\t.0" << std::endl << std::endl;

      prec.setQueryFeatureBegin(qrstart);
      prec.setQueryFeatureEnd(qrstop);
      prec.setInterpolationValue(.0);
      prec.setNodeRange(lnode, unode, qrlength);
      prec.setBestReferenceTaxon(lnode);
      return;
    }

    std::vector< StringType > segments(n);    //TODO: don't call element constructors
    std::vector< float > querydistance(n, std::numeric_limits< float >::max()); // distance to query sequence
    std::vector< float > querysimilarity(n, .0);   // number of matches for nucleotide
    stopwatch_init.stop();

    // count number of alignment calculations in each of the three passes
    uint gcounter = 0;
    uint pass_0_counter = 0;
    uint pass_0_counter_naive = 0;
    uint pass_1_counter = 0;
    uint pass_1_counter_naive = 0;
    uint pass_2_counter = 0;
    uint pass_2_counter_naive = 0;

    StopWatchCPUTime stopwatch_seqret("retrieving sequences for this record");  // log overall time for this predict phase
    StopWatchCPUTime stopwatch_process("processing this record");  // log overall time for this predict phase
    stopwatch_process.start();

    std::set<uint> qgroup;
    large_unsigned_int anchors_support = 0;
    const TaxonNode* rtax = NULL;  // taxon of closest evolutionary neighbor(s)
    const TaxonNode* lca_allnodes = records.front()->getReferenceNode();  // used for optimization

    {   // pass 0 (re-alignment to most similar reference segments)
      logsink << std::endl << "  PASS\t0" << std::endl;
      float dbalignment_searchscore_threshold = reeval_bandwidth_factor_*qmax_searchscore;
      uint index_best = 0;

      for (uint i = 0; i < n; ++i) { // get scores for best-scoring references
        float dist;
        float sim;
        Alignment<StringType> queryalignment;

        const float qsearchscore = records[i]->getScore();
        const large_unsigned_int qsearchmatch = records[i]->getIdentities();
        const double qsearchpid = static_cast<double>(qsearchmatch)/qrlength;

        if(records[i]->getAlignmentLength() == qrlength && records[i]->getIdentities() == qrlength) {

          qgroup.insert(i);
          dist = 0;
          sim = records[i]->getIdentities();
          logsink << std::setprecision(2) << "    *ALN " << i << " <=> query" << tab  << "dist=" << dist << "; sim=" << sim << "; qsearchscore=" << qsearchscore << "; qsearchmatch=" << qsearchmatch << "; qpid=1.0" << std::endl;
          ++pass_0_counter_naive;
        } else if (records[i]->getScore() >= dbalignment_searchscore_threshold) {

          qgroup.insert(i);

          stopwatch_seqret.start();
          if(seqan::empty(segments[i])) segments[i] = getSequence(records[i]->getReferenceIdentifier(),  records[i]->getReferenceStart(), records[i]->getReferenceStop(), records[i]->getQueryStart() - qrstart, qrstop - records[i]->getQueryStop());
          stopwatch_seqret.stop();

          queryalignment = getAlignment(segments[i], qrseq);
          dist = queryalignment.distance;
          ++pass_0_counter;
          ++pass_0_counter_naive;

          sim = std::max(queryalignment.similarity, static_cast<float>(records[i]->getIdentities()));
          double qpid = static_cast<double>(sim)/qrlength;
          logsink << std::setprecision(2) << "    +ALN " << i << " <=> query" << tab  << "dist=" << dist << "; sim=" << sim << "; qsearchscore=" << qsearchscore << "; qsearchmatch=" << qsearchmatch << "; qsearchpid=" << qsearchpid << "; qpid=" << qpid << std::endl;
          logsink << queryalignment << std::endl;

        } else {  // not similar -> fill in some dummy values
          dist = std::numeric_limits< float >::max();
          sim = records[i]->getIdentities();
        }
        querydistance[i] = dist;
        querysimilarity[i] = sim;
        if (dist < querydistance[index_best]) index_best = i;
        else if (dist == querydistance[index_best]) {
          if (sim > querysimilarity[index_best]) index_best = i;
          else if (sim == querysimilarity[index_best] && qsearchscore > records[index_best]->getScore()) index_best = i;
        }
        anchors_support = std::max(anchors_support, static_cast<large_unsigned_int>(sim));  // TODO: make support a float?
        lca_allnodes = this->taxinter_.getLCA(lca_allnodes, records[i]->getReferenceNode());
      }

      // TODO: resort alignments by local score (already done), and by secondary vector (sim)
      // only keep and use the best-scoring reference sequences
      rtax = records[index_best]->getReferenceNode();
      for (std::set< uint >::iterator it = qgroup.begin(); it != qgroup.end();) {
        if (querydistance[*it] != querydistance[index_best] || querysimilarity[*it] != querysimilarity[index_best] || records[*it]->getScore() != records[index_best]->getScore()) qgroup.erase(it++);
        else {
          const TaxonNode* cnode = records[*it]->getReferenceNode();
          rtax = this->taxinter_.getLCA(rtax, cnode);
          logsink << "      current ref node: " << "("<< querydistance[*it] <<") "<< rtax->data->annotation->name << " (+ " << cnode->data->annotation->name << " )" << std::endl;
          ++it;
        }
      }
      assert(! qgroup.empty());  // TODO: only in debug mode

      logsink << "    NUMALN\t" << pass_0_counter << tab << pass_0_counter_naive - pass_0_counter << std::endl << std::endl;
    }


    float anchors_taxsig = 1.;  // a measure of tree-like scores
    float ival_global = 0.;
    const TaxonNode* lnode_global = rtax;
    const TaxonNode* unode_global = rtax;
    std::set<uint> outgroup;
    float bandfactor_max = 1.;

    {   // pass 1 (best reference alignment)
      logsink << "  PASS\t1" << std::endl;

      small_unsigned_int lca_root_dist_min = std::numeric_limits<small_unsigned_int>::max();
      do {  // for each most similar reference segment
        BandFactor bandfactor1(this->taxinter_, n);
        const uint index_anchor = *qgroup.begin();
        qgroup.erase(qgroup.begin());
        const float qdist = querydistance[index_anchor];
        const TaxonNode* rnode = records[index_anchor]->getReferenceNode();
        bandfactor1.addSequence(0, rnode);
        const TaxonNode* lnode = rtax;
        const TaxonNode* unode = NULL;
        float ldist = 0;
        float udist = std::numeric_limits<float>::max();

        std::list< boost::tuple< uint, int > > outgroup_tmp;

        // align all others <=> anchor TODO: adaptive cut-off
        logsink << std::setprecision(2) << "      query: (" << qdist << ") unknown" << std::endl;
        pass_1_counter_naive += n - 1;

        // TODO: implement heuristic cut-off
        double qpid_upper = 0.;
        double qpid_thresh_guarantee = 0.;
        double qpid_thresh_heuristic = 0.;
        int qsearchscore_thresh_heuristic = 0.;

        for(uint i = 0; lnode != this->taxinter_.getRoot() && i < n && records[i]->getScore() >= qsearchscore_thresh_heuristic; ++i) {  //TODO: break loop when qsearchscore < qsearchscore_thresh_heuristic
          const TaxonNode* cnode = records[i]->getReferenceNode();
          const large_unsigned_int qsearchmatch = records[i]->getIdentities();
          const double qsearchpid = static_cast<double>(qsearchmatch)/qrlength;  // must be used for cutoff for stability reasons
          const double qpid = static_cast<double>(querysimilarity[i])/qrlength;
          const float qsearchscore = records[i]->getScore();
          double qpid_thresh = std::max(qpid_thresh_guarantee, qpid_thresh_heuristic);

          if(qpid >= qpid_thresh) {  //TODO: implement command line option
            float dist;
            //float sim;
            Alignment<StringType> segmentalignment;

            if (i == index_anchor) dist = .0;
            else {
              // use triangle relation to avoid alignment
              if (querydistance[i] == .0) { // && querydistance[index_anchor] == 0 ) { //&& querysimilarity[i]) { // TODO: correct?
                dist = querydistance[index_anchor];
                //sim = querysimilarity[index_anchor];
              }
              else {
                stopwatch_seqret.start();
                if(seqan::empty(segments[index_anchor])) segments[index_anchor] = getSequence(records[index_anchor]->getReferenceIdentifier(),  records[index_anchor]->getReferenceStart(), records[index_anchor]->getReferenceStop(), records[index_anchor]->getQueryStart() - qrstart, qrstop - records[index_anchor]->getQueryStop());
                if(seqan::empty(segments[i])) segments[i] = getSequence(records[i]->getReferenceIdentifier(),  records[i]->getReferenceStart(), records[i]->getReferenceStop(), records[i]->getQueryStart() - qrstart, qrstop - records[i]->getQueryStop());
                stopwatch_seqret.stop();

                segmentalignment = getAlignment(segments[i],segments[index_anchor]);
                dist = segmentalignment.distance;

                ++pass_1_counter;
                float sim = segmentalignment.similarity;

                logsink << std::setprecision(2) << "    +ALN " << i << " <=> " << index_anchor << tab << "dist=" << dist << "; sim=" << sim << "; qsearchscore=" << qsearchscore << "; qsearchmatch=" << qsearchmatch << "; qsearchpid=" << qsearchpid << "; qpid=" << qpid << "; qsearchscore_cut=" << qsearchscore_thresh_heuristic << "; qpid_cutg=" << qpid_thresh_guarantee << "; qpid_cut_h=" << qpid_thresh_heuristic << std::endl;
                logsink << segmentalignment << std::endl;
              }
            }

            bandfactor1.addSequence(dist, cnode);

            // place sequence
            if (dist == .0) qgroup.erase(i);  // remove this from list of qnodes because it is sequence-identical
            else {
              if(dist <= qdist) {
                lnode = this->taxinter_.getLCA(lnode, cnode);
                if(dist > ldist) ldist = dist;
                logsink << std::setprecision(2) << "      current lower node: " << "("<< dist <<") "<<lnode->data->annotation->name << " (+ " << cnode->data->annotation->name << " at " << static_cast<int>(this->taxinter_.getLCA(cnode, rnode)->data->root_pathlength) << " )" << std::endl;
              }
              else {
                if(dist < udist) {  // true if we find a segment with a lower dist than query
                  udist = dist;
                  if(qsearchpid > qpid_upper) {
                    qpid_upper = qsearchpid;
                    qpid_thresh_guarantee = qsearchpid*2. - 1.;  // hardcoded inequation: qpid+1.-qpid_up < qpid_up
                    qpid_thresh_heuristic = qsearchpid*exclude_alignments_factor_;
                  }
                  if(!qsearchscore_thresh_heuristic) qsearchscore_thresh_heuristic = records[i]->getScore()*exclude_alignments_factor_;
                }
                outgroup_tmp.push_back(boost::make_tuple(i,dist));
              }
            }
          }
        }

        float bandfactor = bandfactor1.getFactor();  //TODO: limit and check
        bandfactor_max = std::max(bandfactor_max, bandfactor);
        float qdist_ex = qdist * bandfactor;
        float min_upper_dist = std::numeric_limits< int >::max();

        logsink << std::endl << "    EXT\tquerydist = " << qdist << "; threshold = " << qdist_ex << "; bandfactor = " << bandfactor << std::endl;
        for(std::list< boost::tuple<uint,int> >::iterator it = outgroup_tmp.begin(); it != outgroup_tmp.end();) {
          float dist = it->get<1>();

          if(dist > qdist_ex) {
            if (dist > min_upper_dist) it = outgroup_tmp.erase(it);
            else {
              if(dist < min_upper_dist) min_upper_dist = dist;
              ++it;
            }
          } else {
            if(min_upper_dist > qdist_ex) min_upper_dist = dist;
            else min_upper_dist = std::max(min_upper_dist, dist);
            ++it;
          }
        }

        // push elements from temporary to outgroup set
        if(min_upper_dist != std::numeric_limits< float >::max()) unode = lnode;
        for(std::list< boost::tuple<uint,int> >::iterator it = outgroup_tmp.begin(); it != outgroup_tmp.end(); ++it) {
          uint i;
          float dist;
          boost::tie(i, dist) = *it;
          const TaxonNode* cnode = records[i]->getReferenceNode();

          if(dist > min_upper_dist) continue;

          // add to upper node if(dist <= min_upper_dist)
          unode = this->taxinter_.getLCA(cnode, unode);
          logsink << "      current upper node: " << "("<< dist <<") "<< unode->data->annotation->name << " (+ " << cnode->data->annotation->name << " at " << static_cast<int>(this->taxinter_.getLCA(cnode, rnode)->data->root_pathlength) << " )" << std::endl;

          // curate minimal outgroup TODO: only keep dist == min_upper_dist in outgroup?
          const small_unsigned_int lca_root_dist = this->taxinter_.getLCA(cnode, rtax)->data->root_pathlength;
          if(lca_root_dist > lca_root_dist_min) continue;
          else if(lca_root_dist < lca_root_dist_min) {
            lca_root_dist_min = lca_root_dist;
            outgroup.clear();
          }
          outgroup.insert(i);
        }

        // adjust interpolation value and upper node
        float ival = 0.;                                 //  TODO: initialize
        if(!unode) {
          unode = this->taxinter_.getRoot();
          udist = -1;
          ival = 1.;
        } else if(unode != lnode && ldist < qdist) ival = (qdist - ldist)/(udist - ldist);

        logsink << std::endl << std::setprecision(2) << "    SCORE\tldist = " << ldist << "; udist = " << udist << "; querydist = " << qdist << "; querydist_ex = " << qdist_ex << "; ival = " << ival  << std::endl << std::endl;
        const float taxsig = .0;  // TODO: placer.getTaxSignal(qdist);

        ival_global = std::max(ival, ival_global);  // combine interpolation values conservatively
        anchors_taxsig = std::min(taxsig, anchors_taxsig);  // combine taxonomic signal values conservatively
        unode_global = this->taxinter_.getLCA(unode_global, unode);
        lnode_global = this->taxinter_.getLCA(lnode_global, lnode);

      } while (! qgroup.empty() && lnode_global != this->taxinter_.getRoot());

      logsink << "    NUMALN\t" << pass_1_counter << tab << pass_1_counter_naive - pass_1_counter << std::endl;
      logsink << "    NUMOUTGRP\t" << outgroup.size() << std::endl;
    }

    logsink << "    RANGE\t" << rtax->data->annotation->name << tab << lnode_global->data->annotation->name << tab << unode_global->data->annotation->name << std::endl << std::endl;

    {   // pass 2 (stable upper node estimation alignment)
      logsink << "  PASS\t2" << std::endl;
      while (! outgroup.empty()) {
        const uint index_anchor = *outgroup.begin();
        outgroup.erase(outgroup.begin());

        if( unode_global == lca_allnodes ) {
          if( querydistance[index_anchor] == std::numeric_limits<float>::max() ) pass_2_counter_naive += n;
          else pass_2_counter_naive += n - 1;
          continue;
        }

        // align all others <=> anchor TODO: heuristic
        const double qpid_anchor = static_cast<double>(querysimilarity[index_anchor])/qrlength;
        const double qpid_thresh_guarantee = qpid_anchor*2. - 1.;  // hardcoded inequation: qpid+1.-qpid_up < qpid_up;
        const double qpid_thresh_heuristic = qpid_anchor*exclude_alignments_factor_;
        const double qpid_thresh = std::max(qpid_thresh_guarantee, qpid_thresh_heuristic);
        const float qsearchscore_thresh_heuristic = records[index_anchor]->getScore()*exclude_alignments_factor_;
        ++pass_2_counter_naive; // query angainst reference alignment

        for (uint i = 0; i < n && records[i]->getScore() >= qsearchscore_thresh_heuristic; ++i) {
          const double qpid = static_cast<double>(querysimilarity[i])/qrlength;
          if(qpid >= qpid_thresh) {

            const TaxonNode* cnode = records[i]->getReferenceNode();
            const float qsearchscore = records[i]->getScore();
            const large_unsigned_int qsearchmatch = records[i]->getIdentities();
            float dist;
            float sim;
            Alignment<StringType> segmentalignment;

            if (i == index_anchor) dist = .0;
            else {
              ++pass_2_counter_naive;
              if( this->taxinter_.isParentOf(unode_global, cnode) || cnode == unode_global ) continue;
              else {
                stopwatch_seqret.start();
                if(seqan::empty(segments[index_anchor])) segments[index_anchor] = getSequence(records[index_anchor]->getReferenceIdentifier(),  records[index_anchor]->getReferenceStart(), records[index_anchor]->getReferenceStop(), records[index_anchor]->getQueryStart() - qrstart, qrstop - records[index_anchor]->getQueryStop());
                if(seqan::empty(segments[i])) segments[i] = getSequence(records[i]->getReferenceIdentifier(),  records[i]->getReferenceStart(), records[i]->getReferenceStop(), records[i]->getQueryStart() - qrstart, qrstop - records[i]->getQueryStop());
                stopwatch_seqret.stop();

                segmentalignment = getAlignment(segments[i], segments[index_anchor]);
                dist = segmentalignment.distance;
                sim = segmentalignment.similarity;

                logsink << std::setprecision(2) << "    +ALN " << i << " <=> " << index_anchor << tab << "dist=" << dist << "; sim=" << sim << "; qsearchscore=" << qsearchscore << "; qsearchmatch=" << qsearchmatch << "; qpid=" << qpid << std::endl;
                logsink << segmentalignment << std::endl;
                ++pass_2_counter;
                querydistance[i] = dist;
              }
            }

            if (dist == .0) outgroup.erase(i);
            else {
              float qdist_ex;
              if (querydistance[index_anchor] == std::numeric_limits<float>::max()) { //need to align query <=> anchor
                stopwatch_seqret.start();
                if(seqan::empty(segments[index_anchor])) segments[index_anchor] = getSequence(records[index_anchor]->getReferenceIdentifier(),  records[index_anchor]->getReferenceStart(), records[index_anchor]->getReferenceStop(), records[index_anchor]->getQueryStart() - qrstart, qrstop - records[index_anchor]->getQueryStop());
                stopwatch_seqret.stop();

                segmentalignment =  getAlignment(segments[index_anchor], qrseq);
                float dist = segmentalignment.distance;
                float sim = std::max(segmentalignment.similarity, querysimilarity[index_anchor]);

                double qpid = static_cast<double>(sim)/qrlength;
                logsink << std::setprecision(2) << "    +ALN query <=> " << index_anchor << tab << "dist=" << dist << "; sim=" << sim << "; qsearchscore=" << records[index_anchor]->getScore() << "; qsearchmatch=" << qsearchmatch << "; qpid=" << qpid << std::endl;
                logsink << segmentalignment << std::endl;
                querydistance[index_anchor] = dist;
                querysimilarity[index_anchor] = sim;
                qdist_ex = dist*bandfactor_max;
                logsink << "      query: (" << qdist_ex << ") unknown" << std::endl;
                ++pass_2_counter;
              } else qdist_ex = querydistance[index_anchor]*bandfactor_max;

              if(dist <= qdist_ex) {
                const TaxonNode* rnode = records[index_anchor]->getReferenceNode();
                unode_global = this->taxinter_.getLCA(unode_global, cnode);
                logsink << "      current upper node: " << "("<< dist <<") "<< unode_global->data->annotation->name << " (+ " << cnode->data->annotation->name << " at " << static_cast<int>(this->taxinter_.getLCA(cnode, rnode)->data->root_pathlength) << " )" << std::endl;
              }
            }
          }
        }
        logsink << std::endl;
      }
      logsink << "    NUMALN\t" << pass_2_counter << tab << pass_2_counter_naive - pass_2_counter << std::endl;
    }

    if(unode_global == lnode_global) ival_global = 1.;

    logsink << "    RANGE\t" << rtax->data->annotation->name << tab << lnode_global->data->annotation->name << tab << unode_global->data->annotation->name << std::endl << std::endl;

    prec.setSignalStrength(anchors_taxsig);
    prec.setQueryFeatureBegin(qrstart);
    prec.setQueryFeatureEnd(qrstop);
    prec.setInterpolationValue(ival_global);
    prec.setNodeRange(lnode_global, unode_global, anchors_support);
    prec.setBestReferenceTaxon(rtax);
    gcounter = pass_0_counter + pass_1_counter + pass_2_counter;
    float normalised_rt = (float)gcounter/(float)n;
    stopwatch_process.stop();
    logsink << "STATS" << tab << qrseqname << tab << n << tab << pass_0_counter << tab << pass_1_counter << tab << pass_2_counter << tab << gcounter << tab << stopwatch_init.read() << tab << stopwatch_seqret.read() << tab << stopwatch_process.read() << tab << std::setprecision(2) << std::fixed << normalised_rt << std::endl << std::endl;
  }

protected:
  typedef std::list<typename ContainerT::value_type> active_list_type_;
  QStorType& query_sequences_;
  const DBStorType& db_sequences_;
  SortFilter< active_list_type_ > sort_;
  compareTupleFirstLT< boost::tuple< int, uint >, 0 > tuple_1_cmp_le_;

private:
  const float exclude_alignments_factor_;
  const float reeval_bandwidth_factor_;
  StopWatchCPUTime measure_sequence_retrieval_;
  StopWatchCPUTime measure_pass_0_alignment_;
  StopWatchCPUTime measure_pass_1_alignment_;
  StopWatchCPUTime measure_pass_2_alignment_;

  // compile type switch for member function specialization TODO: simplify
  const StringType getSequence(const std::string& id, const large_unsigned_int start, const large_unsigned_int stop, const large_unsigned_int left_ext = 0, const large_unsigned_int right_ext = 0) {
    if(std::is_same<StringType, seqan::String<seqan::AminoAcid>>::value) getSequenceProtein(id, start, stop, left_ext, right_ext);
    // assume DNA sequence
    return getSequenceDNA(id, start, stop, left_ext, right_ext);
  }

  const seqan::String<seqan::AminoAcid> getSequenceProtein(const std::string& id, const large_unsigned_int start, const large_unsigned_int stop, const large_unsigned_int left_ext = 0, const large_unsigned_int right_ext = 0 ) {
    //static_assert(std::is_same<StringType, seqan::String<seqan::AminoAcid>>::value, "StringType mismatch");
    assert(start <= stop);
    large_unsigned_int newstart = left_ext < start ? start - left_ext : 1;
    large_unsigned_int newstop = stop + right_ext;
    return db_sequences_.getSequence(id, newstart, newstop); //TODO: can we avoid copying
  }

  const seqan::String<seqan::Dna5> getSequenceDNA(const std::string& id, const large_unsigned_int start, const large_unsigned_int stop, const large_unsigned_int left_ext = 0, const large_unsigned_int right_ext = 0 ) {
    //static_assert(std::is_same<StringType, seqan::String<seqan::Dna5>>::value, "StringType mismatch");
    if(start <= stop) {
      large_unsigned_int newstart = left_ext < start ? start - left_ext : 1;
      large_unsigned_int newstop = stop + right_ext;
      return db_sequences_.getSequence(id, newstart, newstop); //TODO: can we avoid copying
    }
    large_unsigned_int newstart = right_ext < stop ? stop - right_ext : 1;
    large_unsigned_int newstop = start + left_ext;
    return db_sequences_.getSequenceReverseComplement(id, newstart, newstop);
  }
};

#endif // taxonpredictionmodelsequence_hh_
