#pragma once
#include "common.h"
#include "DataReader.h"
#include "OrderedData.h"
#include "ChromData.h"
#include "Spline.h"
//#include <memory>
#include <fstream>
#include <stdexcept>
#include <iterator>

#include <stdlib.h>     // abs
#include <stdio.h>

#define MY_DEBUG

using point = pair<chrlen, float>;
using coviter = covmap::const_iterator;

const fraglen FragDefLEN = 200;
const coval CUTOFF_STRAND_EXT_RGN = 5;
const uint16_t ReadSplineBASE = 10;	//half-length of moving window for reads spline

// Message verbose level
static class Verb
{
public:
	enum eVerb {	// verbose level
		CRIT,	// print critical messages only
		RES,	// print vCRIT + results
		RT,		// print vRES + runtime info
		DBG		// print vPAR + additional info
	};
	static const char* ValTitles[];
	static const char* ValDescr;

	static uint8_t Size() { return 4; }

	static void Set(uint8_t level) { _level = eVerb(level); }

	static void PrintMsg(eVerb level, const char* msg) { if (Level(level))	printf("%s\n", msg); }

	// Return true if given level is allowed
	static bool Level(eVerb level) { return _level >= level; }

private:
	static eVerb _level;
} verb;

struct StrandOp {
	int Factor;							// 1 for direct reads, -1 for reversed ones
	void (*Next)(coviter&);				// decrement/increment iterator (for direct/reversed)
	coviter (*RetNext)(coviter&);				// decrement/increment iterator (for direct/reversed)
	coviter (*GetPrev)(const coviter&);	// returns prev/this iterator (for direct/reversed)
	bool (*EqLess)(chrlen pos, chrlen lim);
};

// strand-dependent operations: 0 - direct, 1 - reversed
static const StrandOp StrandOps[2]{
	{-1,	
	[](coviter& it) { it--; },
	[](coviter& it) { return --it; },
	[](const coviter& it) { return prev(it); },
	[](chrlen pos, chrlen lim) { return pos <= lim; }
	},
	{ 1, 
	[](coviter& it) { it++; }, 
	[](coviter& it) { return ++it; },
	[](const coviter& it) { return it; },
	[](chrlen pos, chrlen lim) { return pos >= lim; }
	}
};

static struct Glob {
	static readlen ReadLen;
	static fraglen FragLen;		// average fragment length
	static fraglen ROI_ext;		// Regions of interest extention	
} glob;


//===== WRITING AN OBLIQUE CURVE

// 'LinearWriter' implements the method of writing an oblique curve in FixedStep wig format
class LinearWriter : WigWriter
{
public:
	// Creates new instance for writing cover to wiggle_0 file.
	//	@param strand: strand
	//	@param fields: BED/WIG track fields
	LinearWriter(eStrand strand, const TrackFields& fields)
		: WigWriter(FT::eType::WIG_FIX, strand, fields) {}

	// Writes oblique line
	//	@param cID: chrom's ID
	//	@param start: start (left) position
	//	@param ptCnt: number of points in oblique line
	//	@param shift: shift of value at each point of the oblique line: if < 0 then direct (positive) line, otherwise reversed (negative) one
	void WriteOblique(chrid cID, chrlen start, chrlen ptCnt, float shift)
	{
		WriteFixStepLine(cID, start + 1, ptCnt, shift);	// +1 to match "0-start, half-open" coordinates used by bedGraph format
	}
};

// Ordered LinearWriter
class OLinearWriter : OrderedData<int, LinearWriter>
{
	chrid _cID = Chrom::UnID;

public:
	// Primer constructor
	//	@param cSizes: chrom sizes
	//	@param dim: number (dimension) of data
	//	@param write: if true then data sould be save to output file
	//	@param fname: name of output file
	//	@param descr: track decsription in declaration line 
	//OLinearWriter(const ChromSizes& cSizes, BYTE dim, bool write, const TrackFields& fields)
	OLinearWriter(const ChromSizes& cSizes, BYTE dim, bool write, const string& name, const char* descr)
		: OrderedData<int, LinearWriter>(cSizes, dim, write, name, descr) {}

	void SetChromID(chrid cID) { _cID = cID; }

	bool IsWriterSet() const { return _writers != nullptr; }

	// Writes oblique line (positive or negative strand)
	//	@param reverse: true if tag is reversed (neg strand)
	//	@param start: start (zero value) position
	//	@param itStop: stop (max value) iterator
	void WriteOblique(BYTE reverse, chrlen start, coviter& itStop);
};


//===== DATA

// Linear Regression result
struct LRegrResult
{
	chrlen	Pos;
	coval	PtCnt;
	float	Deriv;

	bool Valid() const { return
		Pos					// zero Pos is the result of a failed regression calculation
		&& Deriv < 1		// 1 (vertical line) is useless
		&& Deriv >= 0.01;	// too small angle is useless
	}

	void Clear() { Deriv = 0; PtCnt = Pos = 0; }
};

// temporary Reads collection
class Reads
{
	vector<Region>	_reads[2];

public:
	void Reserve(size_t capac) {
		_reads[0].reserve(capac >> 1);
		_reads[1].reserve(capac >> 1);
	}

	const vector<Region>& GetReads(bool reverse) const { return _reads[reverse]; }

	void AddRead(const Region& rgn, bool reverse) { _reads[reverse].push_back(rgn); }

	void Clear() {
		for (BYTE s : {0, 1}) { _reads[s].clear(); _reads[s].shrink_to_fit(); }
	}
};


//=== TREATED COVER & COLLECTION
struct CoverRegion;

class TreatedCover : public AccumCover
{
#ifdef MY_DEBUG
	class TxtOutFile
	{
		FILE* _file;
	public:
		TxtOutFile(const char* name) { _file = fopen(name, "w"); }
		~TxtOutFile() { fclose(_file); }

		template<typename... Args>
		void Write(const char* format, Args ... args) {	fprintf(_file, format, args...); }
	};

	shared_ptr<TxtOutFile> _fDdelim = nullptr;
#endif

public:
#ifdef MY_DEBUG
	static bool WriteDelim;

	TreatedCover() { if(WriteDelim) _fDdelim.reset(new TxtOutFile("delim.txt")); }
#endif
	coval GetMaxVal() const;

	// Returns cover mass centre for given region
	//	@param rng: potential region
	//chrlen GetRegionCentre(const CoverRegion& rng) const;

	// Calculates Linear Regression
	//	@param it: start iterator
	//	@param itStop: end iterator
	//	@param op: strand operations
	//	@result[out]: calculated results
	void LinearRegr(coviter it, const coviter& itStop, const StrandOp& op, LRegrResult& result) const;
};

// 'CombCover' keeps the chromosome covers and optionally the set of writers these covers to file.
class CombCover : public OrderedData<TreatedCover, BedGrWriter>
{
	string	_outFName;				// common name of output files

public:
	// Primer constructor
	//	@param cSizes: chrom sizes
	//	@param dim: number of data (dimension)
	//	@param write: if true then data sould be save to output file
	//	@param fname: output common file name (without extention)
	//	@param descr: 'description' field in BED/WIG track line
	CombCover(const ChromSizes& cSizes, BYTE dim, bool write, const string& fname, const char* descr)
		: _outFName(fname)
		, OrderedData<TreatedCover, BedGrWriter>(cSizes, dim, write, fname, descr)
	{}

	//coval GetMaxVal() const { return _data->StrandData(POS).GetMaxVal(); }

	// For current chrom adds extended SE tag to total coverage, and pure tag to strand coverage
	//	@param[in] read: added tag
	//	@param[in] reverse: true if tag is reversed (neg strand)
	void AddRead(const Region& read, bool reverse)
	{
		_data->StrandDataByInd(reverse).AddRegion(read);	// strand read coverage
	}

	// For current chrom adds extended SE tag to total coverage, and pure tag to strand coverage
	//	@param[in] read: added tag
	//	@param[in] reverse: true if tag is reversed (neg strand)
	//	@param[in] addTotal: true if tag is added to total cover as well
	void AddExtRead(const Region& read, bool reverse, bool addTotal);

	void Fill(const Reads& reads, bool addTotal);
};


//=== COVER REGION & COLLECTION

// 'CoverRegion' represetns potential region
struct CoverRegion
{
	coviter	itStart;
	coviter	itEnd;
	coval	value;

	CoverRegion(coviter& start, coviter& end, coval val) : itStart(start), itEnd(end), value(val) {}

	chrlen	Start()	 const { return itStart->first; }
	chrlen	End()	 const { return itEnd->first; }
	chrlen	Length() const { return End() - Start(); }
};

class CoverRegions : public vector<CoverRegion>
{
	using Iter = vector<CoverRegion>::iterator;
	using cIter = vector<CoverRegion>::const_iterator;

	friend class DataCoverRegions;

	// Fills the instance by potential regions on extended read coverage
	//	@param cover: fragment coverage
	//	@param capacity: capacity to reserve the instance
	//	@param cutoff: fragment coverage cut off value
	void SetPotentialRegions(const TreatedCover& cover, size_t capacity, coval cutoff);

public:
	// methods used in EliminateNonOverlaps()
	static chrlen	Start	(cIter it) { return it->Start(); }
	static chrlen	End		(cIter it) { return it->End(); }
	static chrlen	Length	(cIter it) { return it->Length(); }
	static bool		IsWeak	(cIter it) { return it->value <= CUTOFF_STRAND_EXT_RGN; }
	static bool	IsNotEmpty	(cIter it) { return it->value; }
	static void	MarkAsEmpty	(Iter it)  { it->value = 0; }
};

class DataCoverRegions : public DataSet<CoverRegions>
{
public:
	// Fills the instance by potential regions of extended SE read coverage
	//	@param cover: fragment coverage
	//	@param cLen: length of chromosome
	//	@param cutoff: fragment coverage cut off value
	//	@param noMultiOverl: if true then eliminate multi overlaps regions as well
	void SetPotentialRegions(const DataSet<TreatedCover>& cover, chrlen cLen, coval cutoff, bool noMultiOverl);

	// Fills the instance by potential regions of fragment coverage
	//	@param cover: fragment coverage
	//	@param cLen: length of chromosome
	//	@param cutoff: fragment coverage cut off value
	void SetPotentialRegions(const DataSet<TreatedCover>& cover, chrlen cLen, coval cutoff);

	// Returns mean fragment length based on cover regions mass centre comparison
	//	@param cover: fragment coverage
	//fraglen GetFragMean(const DataSet<TreatedCover>& cover) const;

	void Clear();
};


//=== BINDING SITES DATA 

struct BS_PosVal
{
	BYTE	Reverse;
	chrlen	GrpNumb;
	coval	PointCount;
	float	Score;
	BYTE	Unreliable = 0;
	bool	BS_NegWidth = false;

	//BS_PosVal(BYTE reverse, chrlen grpNumb, coval pointCount, float deriv, float splineVal) :
	//	Reverse(reverse), GrpNumb(grpNumb), PointCount(pointCount), Score(deriv* splineVal) {}

	BS_PosVal(BYTE reverse, chrlen grpNumb, const LRegrResult& result, float splineVal, BYTE unreliable) :
		Reverse(reverse), GrpNumb(grpNumb), PointCount(result.PtCnt), Score(result.Deriv * splineVal),
		Unreliable(unreliable) {}
};

class BS_Map : public map<chrlen, BS_PosVal>
{
	//using iter = map<chrlen, BS_PosVal>::iterator;

public:
	using iter = map<chrlen, BS_PosVal>::iterator;
	using citer = map<chrlen, BS_PosVal>::const_iterator;

	struct ValPos {
		iter	Iter;
		float	Val;

		ValPos(iter it, float val) : Iter(it), Val(val) {}
	};

	//void AddPos(BYTE reverse, chrlen grpNumb, chrlen pos, coval pointCount, float deriv, float splineVal) {
	//	emplace(pos, BS_PosVal(reverse, grpNumb, pointCount, deriv, splineVal));	// !!! use hint?
	//}
	void AddPos(BYTE reverse, chrlen grpNumb, const LRegrResult& result, float splineVal, BYTE unreliable) {
		emplace(
			result.Pos + StrandOps[reverse].Factor * Glob::ReadLen,
			BS_PosVal(reverse, grpNumb, result, splineVal, unreliable)
		);	// !!! use hint?
	}

	// Brings the instance to canonical order of placing BS borders
	void Refine();

	// Applies lambda to each group of binding sites, passing borders collection
	template<typename F>
	void Do(F&& lambda)
	{
		vector<ValPos> VP[2];	// 0 - direct, 1 - reversed

		VP[0].reserve(4), VP[1].reserve(4);
		for (auto it = begin(); it != end(); it++)
			if (it->second.Score) {
				if (it->second.Reverse && VP[0].size() && VP[1].size()) {
					lambda(VP);
					VP[0].clear(), VP[1].clear();
				}
				VP[it->second.Reverse].emplace_back(it, it->second.Score);
			}

		// last element
		if (VP[0].size() && VP[1].size())
			lambda(VP);
	}

	// Applies lambda to each group denotes the binding site, passing start-end iterators
	template<typename F>
	void DoBasic(F&& lambda) const
	{
		for (auto it0 = begin(), it = next(it0); it != end(); it0++, it++)
			if (it->second.Score && !it->second.Reverse && it0->second.Reverse)
				lambda(it0, it);
	}

	void NormalizeScore();

	void NormalizeBSwidth();

	void PrintStat() const;

#ifdef MY_DEBUG
	// Prints "unsorted" scores
	void CheckScoreHierarchy();

	// Prints positions
	//	@param real: if true then prints position with unzero score 
	//	@param stopPos: max printed position or all by default 
	void Print(chrid cID, bool real, chrlen stopPos = 0) const;
#endif
};

// 'BedWriter' implements methods for writing in ordinary bed format
class BedWriter : public RegionWriter
{
	const BYTE MIN_BS_WIDTH = 5;
	static bool rankScore;	// true if scores in the main bed should be normalized relative to 1000

	void LineAddIntScore	(float score, bool delim) { LineAddInt(int(round(score * 1000)), delim); }
	void LineAddFloatScore	(float score, bool delim) { LineAddFloat(score, delim); }

	typedef void (BedWriter::* tAddScore)(float, bool);
	static tAddScore fLineAddScore;

	void LineAddScore(float score, bool delim) { assert(fLineAddScore); (this->*fLineAddScore)(score, delim); }

	// Adds BS start/end positions and number, correcting 'negative' BS width
	void LineAddRegion(
		//const BS_Map::iter& itStart, const BS_Map::iter& itEnd, chrlen bsNumb, bool addDelim, fraglen expLength = 0);
		BS_Map::citer itStart, BS_Map::citer itEnd, chrlen bsNumb, bool addDelim, fraglen expLength = 0);

	
public:
	static bool RankScore() { return rankScore; }

	static void SetRankScore(bool rank) {
		fLineAddScore = (rankScore = rank) ? &BedWriter::LineAddIntScore : &BedWriter::LineAddFloatScore;
	}

	// Creates new instance for writing.
	//	@param strand: strand
	//	@param fields: BED/WIG track fields
	BedWriter(eStrand strand, const TrackFields& fields)
		: RegionWriter(FT::eType::BED, strand, fields) {}

	// Writes regions as lines containing 4 or 9-field feature; adds a gray color to feature with zero score
	//	@param cID: chrom ID
	//	@param rgns: range values
	void WriteChromData(chrid cID, const CoverRegions& rgns);

	// Writes BSs as lines containing 7 or 11-fields feature; adds a gray color to feature with zero score
	//	@param cID: chrom ID
	//	@param bss: binding sites (const in fact)
	void WriteChromData(chrid cID, BS_Map& bss);

	// Writes extended BSs; adds a gray color to feature with zero score
	//	@param cID: chrom ID
	//	@param bss: binding sites (const in fact)
	void WriteChromExtData(chrid cID, BS_Map& bss);

	// Writes Regions of interest (extended binding sites) as a 4-fields features
	//	@param cID: chrom ID
	//	@param bss: binding sites
	void WriteChromROI(chrid cID, const BS_Map& bss);
};

class BedWriters
{
	BedWriter _main;
	BedWriter _ext;
	BedWriter _roi;

public:
	BedWriters(eStrand strand, const TrackFields& fields, const string* commLine = nullptr)
		: _main(strand, TrackFields(fields, false, BedWriter::RankScore(), "Black"))
		, _ext (strand, TrackFields(fields, "_ext", "extended called binding sites", true, false, "Black"))
		, _roi (strand, TrackFields(fields, ".ROI", "IGV Region Navigator list", false, false))
	{
		_main.SetFloatFractDigits(3);
	}

	void WriteChromData(chrid cID, BS_Map& bss)
	{
		_main.WriteChromData(cID, bss);
		_ext.WriteChromExtData(cID, bss);
		_roi.WriteChromROI(cID, bss);
	}
};

using OBS_Map = OrderedData<BS_Map, BedWriters>;
//class OBS_Map : public OrderedData<BS_Map, BedWriter>
//{
//public:
//	// Primer constructor
//	//	@param cSizes: chrom sizes
//	//	@param dim: number (dimension) of data
//	//	@param write: if true then data sould be save to output file
//	//	@param fname: name of output file
//	//	@param descr: track decsription in declaration line
//	OBS_Map(const ChromSizes& cSizes, BYTE dim, bool write, const string& fname, const char* descr)
//		: OrderedData<BS_Map, BedWriter>(cSizes, dim, write, TrackFields(fname, descr)) {}
//};


//=== SPLINE & COLLECTION

// Ñollection of sequential float values
class Values : public vector<float>
{
	float _maxVal;

public:
	chrlen GrpNumb = 0;

	Values() noexcept : _maxVal(0) { Reserve(); }
	Values(Values&& rvals) noexcept : _maxVal(rvals._maxVal), GrpNumb(0), vector<float>(move(rvals)) { rvals._maxVal = 0; }
	Values(const Values& rvals) = default;

	// Returns region length
	fraglen Length() const { return fraglen(size()); }

	// Returns value for given region relative position
	float Val(chrlen relPos) const { return (*this)[relPos]; }

	// Returns maximum region value
	float MaxVal() const { return _maxVal; }

	void MarkAsEmpty() { _maxVal = 0; }

	void Reserve() { reserve(40); }

	void Clear() { _maxVal = 0; clear(); }

	// Adds calue to the instance
	void AddVal(float val);

	// Adds positions with local maximum value to the vector pos
	void GetMaxValPos(chrlen startPos, vector<chrlen>& pos) const;

#ifdef MY_DEBUG
	void Print() const;
#endif
};

using tValuesMap = map<chrlen, Values>;

// Maped collections of float values
class ValuesMap : public tValuesMap
{
	using Iter = map<chrlen, Values>::iterator;
	using cIter = map<chrlen, Values>::const_iterator;

	float _maxVal = 0;

#ifdef MY_DEBUG
	friend class DataValuesMap;		// access to this->Print

	void Print(chrid cID, BYTE reverse, chrlen stopNumb) const;
#endif

	// Filters and fill spline curve by read cover within potential region
	//	@param cover: raw read cover
	//	@param rgn: potential region
	//	@param redifRgns: if true then redifine regions position
	//	@param splineBase: half-length of spliner moving window
	void BuildRegionSpline(const TreatedCover& cover, const CoverRegion& rgn, bool redifRgns, fraglen splineBase);

	// Adds values with given position
	//	@param pos: starting position of values
	//	@param vals: values (will be moved)
	void AddRegion(chrlen pos, Values& vals);

public:
	// methods used in EliminateNonOverlaps()
	static chrlen	Start (cIter it) { return it->first; }
	static chrlen	End	  (cIter it) { return it->first + it->second.Length(); }
	static chrlen	Length(cIter it) { return it->second.Length(); }
	static bool		IsWeak(cIter)	 { return false; }
	static bool		IsNotEmpty(cIter it) { return it->second.MaxVal(); }
	static void		MarkAsEmpty(Iter it) { it->second.MarkAsEmpty(); }

	chrlen	Start() const { return Start(begin()); }
	chrlen	End()	const { return End(begin()); }
	chrlen	Length()const { return Length(begin()); }

	float MaxVal() const { return _maxVal; }

	// Filters and fill spline curve by read cover within potential regions
	//	@param cover: raw read/frag cover
	//	@param rgns: potential regions
	//	@param redifRgns: if true then redifine regions position
	//	@param splineBase: half-length of spliner moving window
	void BuildSpline(const TreatedCover& cover, const CoverRegions& rgns, bool redifRgns, fraglen splineBase);

	// Resets non overlapping spline value
	void EliminateNonOverlaps();

	// Sets a single consecutive number (group number) to each overlapped not empty spline
	void Numerate();

	// Prints potential regions before and after selection
	void PrintStat(chrlen clen) const;
};

// Datset of maped collections of float values
class DataValuesMap : public DataSet<ValuesMap>
{
public:
	// Filters and fill spline curve by [extended] SE read cover within potential regions
	//	@param cover: raw read/frag cover
	//	@param rgns: potential regions
	//	@param redifineRgns: if true then redifine regions position
	//	@param splineBase: half-length of spliner moving window
	void BuildSplineSE(const DataSet<TreatedCover>& cover, const DataCoverRegions& rgns, bool redifineRgns, fraglen splineBase);// = SSpliner::ReadSplineBase);

	// Filters and fill spline curve by PE read cover within potential regions
	//	@param cover: raw read/frag cover
	//	@param rgns: potential regions
	//	@param redifineRgns: if true then redifine regions position
	//	@param splineBase: half-length of spliner moving window
	void BuildSplinePE(const DataSet<TreatedCover>& cover, const DataCoverRegions& rgns, bool redifineRgns, fraglen splineBase);// = SSpliner::ReadSplineBase);

	//void BuildSpline(const DataSet<TreatedCover>& cover, const DataFeatureVals& rgns, fraglen splineBase = SSpliner::ReadSplineBase);
	
	// Returns fragment's mean calculated by minimizing the difference between peak positions
	fraglen GetFragMean() const;

	void Clear();

#ifdef MY_DEBUG
	void Print(chrid cID, chrlen stopNumb = 0) const;
#endif
};


//=== DERIVATIVES & COLLECTION

// Ñollection of float values with its associated position and minimum/maximum value
class BoundValues : public Values
{
	chrlen	_pos;		// region start position
	float	_val[2];	// region min, max values
public:

	BoundValues(chrlen startPos, float valMin, float valMax, Values& vals) 
		: _pos(startPos), _val{valMin,valMax}, Values(move(vals))
	{}

	chrlen	Start()	const { return _pos; }
	chrlen	End()	const { return _pos + Length(); }
	// Copies min, max values to the external pair
	void	GetValues(float(&val)[2]) const { memcpy(val, _val, 2 * sizeof(float)); }
};

// BoundValues collection representing 
// rising (for the right side of the peak) or falling (for the left side of the peak) derivatives
class BoundsValues : public vector<BoundValues>
{
	using citer = vector<BoundValues>::const_iterator;

	// Start/end positions and min/max values with storage of previous ones
	class TracedPosVal
	{
		chrlen	_pos[2], _pos0[2]{ 0,0 };	// current, previous start/end positions
		float	_val[2], _val0[2]{ 0,0 };	// current, previous min/max values 

	public:
		// Returns current position
		//	@param reverse: 0 - start, 1 - end
		chrlen	Pos(BYTE reverse) const { return _pos[reverse]; }

		// Returns previous position
		//	@param reverse: 0 - start, 1 - end
		chrlen	PrevPos(BYTE reverse) const { return _pos0[reverse]; }

		// Returns current value
		//	@param lim: 0 - min, 1 - max
		float	Val(BYTE lim) const { return _val[lim]; }

		// Set start/end positions and values by cover iterator, and merge successive raises
		//	@param reverse: 0 for firect, 1 for reverse
		//	@param it: direct/reverse cover iterator
		template<typename T>
		void Set(BYTE reverse, T it)
		{
			// set positions & values
			_pos[!reverse] = it->Start();	// start: right for direct, left for reverse
			_pos[reverse] = it->End();		// end: left for direct, right for reverse
			it->GetValues(_val);
			// merge successive raises
			if (_val0[reverse] && _val[!reverse] / _val0[reverse] >= 0.9)
				_pos[0] = _pos0[0];
		}

		// Saves previous positions/values (by swapping)
		void Trace()
		{
			std::swap(_pos, _pos0);
			std::swap(_val, _val0);
		}
	};

	float _maxVal = 0;
	chrlen _grpNumb;

	// Calculates cover linear regression and adds BS position between start and stop
	//	@param reverse: true if tag is reversed (neg strand)
	//	@param posVal: traced positions and values
	//	@param cover: coverage on which the regression is calculated
	//	@param bs: BS collection to which the new position is added
	//	@param lwriter: writer to save linear regression representation
	void AddBSpos(
		BYTE reverse,
		const TracedPosVal& posVal,
		chrlen& prevStart,
		chrlen& prevBSpos,
		const TreatedCover& cover,
		BS_Map& bs,
		OLinearWriter& lwriter
	) const;

public:
	float	MaxVal()	const { return _maxVal; }
	//chrlen	GroupNumb()	const { return _grpNumb; }

	//using cIter = vector<ValuesMap>::const_iterator;
	//using rIter = vector<ValuesMap>::reverse_iterator;
	//typedef vector<ValuesMap>::iterator iter_type;

	void SetDirectBSpos	(const TreatedCover& cover, BS_Map& bs, OLinearWriter& lwriter) const;

	void SetReverseBSpos(const TreatedCover& cover, BS_Map& bs, OLinearWriter& lwriter) const;

	// Adds derivative values for given spline relative position
	//	@param spline: given spline
	//	@param relPos: spline relative position
	//	@param deriv: derivative values (will be moved)
	void AddValues(const tValuesMap::value_type& spline, chrlen relPos, Values& deriv);
};

class BoundsValuesMap : public map<chrlen, BoundsValues>
{
public:
	using Iter = map<chrlen, BoundsValues>::iterator;
	using cIter = map<chrlen, BoundsValues>::const_iterator;

	// Adds values with given position
	//	@param pos: starting position of values
	//	@param vals: values (will be moved)
	void AddRegions(chrlen pos, BoundsValues& vals) { emplace(pos, move(vals)); }

	// Fill instance with derivatives of splines
	//	@param splines: splines from which derivatives are built
	//void BuildDerivs(BYTE reverse, const ValuesMap& splines);
	void BuildDerivs(int factor, const ValuesMap& splines);

	void SetBSpos(BYTE reverse, const TreatedCover& cover, BS_Map& bs, OLinearWriter& lwriter) const;

#ifdef MY_DEBUG
	void Print(eStrand strand, chrlen stopPos = 0) const;
#endif
};

class DataBoundsValuesMap : public DataSet<BoundsValuesMap>
{
public:
	void BuildDerivs(const DataValuesMap& splines)
	{
		StrandData(POS).BuildDerivs(StrandOps[0].Factor, splines.StrandData(POS));
		StrandData(NEG).BuildDerivs(StrandOps[1].Factor, splines.StrandData(NEG));
	}

	void SetBSpos(const DataSet<TreatedCover>& rCover, BS_Map& bss, OLinearWriter& lwriter)
	{
		StrandData(POS).SetBSpos(0, rCover.StrandData(POS), bss, lwriter);
		StrandData(NEG).SetBSpos(1, rCover.StrandData(NEG), bss, lwriter);
	}

#ifdef MY_DEBUG
	void Print(chrlen stopPos = 0) const
	{
		StrandData(POS).Print(POS, stopPos);
		StrandData(NEG).Print(NEG, stopPos);
	}
#endif

};


//===== WIG WRITERS

// 'FixWigWriter' implements methods for writing in wiggle_0 (fixed step) format
class FixWigWriter : WigWriter
{
public:
	// Creates new instance for writing cover to wiggle_0 file.
	//	@param strand: strand
	//	@param fields: BED/WIG track fields
	FixWigWriter(eStrand strand, const TrackFields& fields)
		: WigWriter(FT::eType::WIG_FIX, strand, fields) {}

	void WriteChromData(chrid cID, const ValuesMap& vals);
};

class FixWigWriterSet : WigWriter
{
public:
	// Creates new instance for writing cover to wiggle_0 file.
	//	@param strand: strand
	//	@param fields: BED/WIG track fields
	FixWigWriterSet(eStrand strand, const TrackFields& fields)
		: WigWriter(FT::eType::WIG_FIX, strand, fields) {}

	void WriteChromData(chrid cID, const BoundsValuesMap& set);
};

//=====  ORDERED DATA

// Ordered FeatureVals
class OCoverRegions : public OrderedData<CoverRegions, BedWriter>
{
public:
	// Primer constructor
	//	@param cSizes: chrom sizes
	//	@param dim: number (dimension) of data
	//	@param write: if true then data sould be save to output file
	//	@param fname: name of output file
	//	@param descr: track decsription in declaration line
	OCoverRegions(const ChromSizes& cSizes, BYTE dim, bool write, const string& fname, const char* descr)
		: OrderedData<CoverRegions, BedWriter>(cSizes, dim, write, fname, descr) {}

#ifdef MY_DEBUG
	void Print() const {
		for (BYTE i = 0; i <= 1; i++) {
			printf("%s regions:\n", sStrandTITLES[i]);
			const auto& data = _data->StrandDataByInd(i);
			for (auto it = data.cbegin(); it != data.cend(); it++)
				printf("%9d%9d%5d\n", it->Start(), it->End(), it->value);
		}
	}
#endif
};

// Ordered ValuesMap (splines)
using OValuesMap = OrderedData<ValuesMap, FixWigWriter>;
//class OValuesMap : public OrderedData<ValuesMap, FixWigWriter>
//{
//public:
//	// Primer constructor
//	//	@param cSizes: chrom sizes
//	//	@param dim: number (dimension) of data
//	//	@param write: if true then data sould be save to output file
//	//	@param fname: name of output file
//	//	@param descr: track decsription in declaration line
//	OValuesMap(const ChromSizes& cSizes, BYTE dim, bool write, const string& fname, const char* descr)
//		: OrderedData<ValuesMap, FixWigWriter>(cSizes, dim, write, TrackFields(fname, descr)) {}
//};

// Ordered BoundsValuesMap (derivatives)
using OBoundsValuesMap = OrderedData<BoundsValuesMap, FixWigWriterSet>;
//class OBoundsValuesMap : public OrderedData<BoundsValuesMap, FixWigWriterSet>
//{
//public:
//	// Primer constructor
//	//	@param cSizes: chrom sizes
//	//	@param dim: number (dimension) of data
//	//	@param write: if true then data sould be save to output file
//	//	@param fname: name of output file
//	//	@param descr: track decsription in declaration line
//	OBoundsValuesMap(const ChromSizes& cSizes, BYTE dim, bool write, const string& fname, const char* descr)
//		: OrderedData<BoundsValuesMap, FixWigWriterSet>(cSizes, dim, write, TrackFields(fname, descr)) {}
//	//OBoundsValuesMap(const ChromSizes& cSizes, BYTE dim, bool write, const TrackFields& fields)
//	//	: OrderedData<BoundsValuesMap, FixWigWriterSet>(cSizes, dim, write, fields) {}
//};
