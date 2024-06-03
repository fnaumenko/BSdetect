/**********************************************************
Treatment.h
Provides support for binding sites discovery
Fedor Naumenko (fedor.naumenko@gmail.com)
Last modified: 06/02/2024
***********************************************************/
#pragma once
#include "common.h"
#include "DataReader.h"
#include "OrderedData.h"
#include "ChromData.h"
#include "Spline.h"
#include <fstream>
#include <stdexcept>
#include <iterator>

#include <stdlib.h>     // abs
#include <stdio.h>

#define MY_DEBUG

using point = pair<chrlen, float>;
using coviter = covmap::const_iterator;

const fraglen	FragDefLEN = 200;
const coval		CUTOFF_STRAND_EXT_RGN = 5;
const uint16_t	ReadSplineBASE = 10;	//half-length of moving window for reads spline

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

	static void PrintMsg(eVerb level, const char* msg = NULL);

	//static void PrintMsgVar(eVerb level, const char* format, ...);

	template<typename... Args>
	static void PrintMsgVar(eVerb level, const char* format, Args ... args) { 
		if (Level(level)) printf(format, args...);
	}

	// Return true if given level is allowed
	static bool Level(eVerb level) { return _level >= level; }

	// Return true if given level is set
	static bool StrictLevel(eVerb level) { return _level == level; }

private:
	static eVerb _level;
} verb;

struct StrandOp {
	int Factor;							// 1 for forward reads, -1 for reversed ones
	void (*Next)(coviter&);				// decrement/increment iterator (for forward/reversed)
	coviter (*RetNext)(coviter&);		// decrement/increment iterator (for forward/reversed)
	coviter (*GetPrev)(const coviter&);	// returns prev/this iterator (for forward/reversed)
	bool (*EqLess)(chrlen, chrlen);
	bool (*Less)(chrlen, chrlen);
};

// strand-dependent operations: 0 - forward, 1 - reversed
static const StrandOp StrandOps[2]{
	{-1	
	,[](coviter& it) { it--; }
	,[](coviter& it) { return --it; }
	,[](const coviter& it) { return prev(it); }
	,[](chrlen pos, chrlen lim) { return pos <= lim; }
	,[](chrlen pos, chrlen lim) { return pos < lim; }
	},
	{ 1
	,[](coviter& it) { it++; }
	,[](coviter& it) { return ++it; }
	,[](const coviter& it) { return it; }
	,[](chrlen pos, chrlen lim) { return pos >= lim; }
	,[](chrlen pos, chrlen lim) { return pos > lim; }
	}
};

static struct Glob {
	static bool		IsPE;		// is sequense paired-end
	static bool		FragLenUndef;
	static readlen	ReadLen;	// length of read
	static fraglen	FragLen;	// average fragment length
	static fraglen	ROI_ext;	// Regions of interest extention

	static void SetPE(bool isPE) { if ((IsPE = isPE)) FragLenUndef = false; }

	static void SetFragLen(int fragLen)
	{
		if (fragLen != vUNDEF) {
			FragLenUndef = false;
			FragLen = fragLen;
		}
	}
} glob;


// Sequential float values
class Values : public vector<float>
{
	float _maxVal;

public:
	chrlen RgnNumb = 0;

	Values() noexcept : _maxVal(0) { Reserve(); }
	Values(Values&& rvals) noexcept : _maxVal(rvals._maxVal), RgnNumb(0), vector<float>(move(rvals)) { rvals._maxVal = 0; }
	Values(const Values& rvals) = default;

	// Returns values length
	fraglen Length() const { return fraglen(size()); }

	// Returns value for given relative position within sequential values
	float Value(chrlen relPos) const { return (*this)[relPos]; }

	// Returns maximum value within sequential values
	float MaxVal() const { return _maxVal; }

	void MarkAsEmpty() { _maxVal = 0; }

	void Reserve() { reserve(40); }

	void Clear() { _maxVal = 0; clear(); }

	void AddValue(float val);

	// Adds sequential values to the instance
	void AddValues(const Values& vals);

	// Adds positions with local maximum value
	//	@param startPos[in]: start position to search
	//	@oaram pos[out]: vector of maximum value positions to which the value is added
	void GetMaxValPos(chrlen startPos, vector<chrlen>& pos) const;

#ifdef MY_DEBUG
	void Print(bool prValues = false) const;
#endif
};

#ifdef MY_DEBUG
//===== SPECIAL WRITERS

// 'SpecialWriter' implements the method of debug writing for the specified chromosome in FixedStep wig format
class SpecialWriter : WigWriter
{
public:
	// Creates new instance for writing cover to wiggle_0 file.
	//	@param strand: strand
	//	@param fields: BED/WIG track fields
	SpecialWriter(eStrand strand, const TrackFields& fields)
		: WigWriter(FT::eType::WIG_FIX, strand, fields)
	{
		SetFloatFractDigits(3);
	}

	// Writes inclined line
	//	@param cID: chrom's ID
	//	@param start: start (left) position
	//	@param ptCnt: number of points in incline
	//	@param shift: shift of value at each point of the incline: if < 0 then forward (positive) line, otherwise reversed (negative) one
	void WriteIncline(chrid cID, chrlen start, chrlen ptCnt, float shift)
	{
		WriteFixStepLine(cID, start + 1, ptCnt, shift);	// +1 to match "0-start, half-open" coordinates used by bedGraph format
	}

	// Writes curve
	//	@param cID: chrom's ID
	//	@param start: start (left) position
	//	@param vals: values to write
	void WriteChromData(chrid cID, chrlen start, const Values& vals)
	{
		WriteFixStepRange(cID, start, vals);
	}
};

// Ordered SpecialWriter
class OSpecialWriter : OrderedData<int, SpecialWriter>
{
	chrid _cID = Chrom::UnID;

public:
	// Primer constructor
	//	@param cSizes: chrom sizes
	//	@param dim: number (dimension) of data
	//	@param write: if true then data sould be save to output file
	//	@param fname: name of output file
	//	@param descr: track decsription in declaration line 
	OSpecialWriter(const ChromSizes& cSizes, BYTE dim, bool write, const string& name, const char* descr)
		: OrderedData<int, SpecialWriter>(cSizes, dim, write, name, descr) {}

	void SetChromID(chrid cID) { _cID = cID; }

	bool IsWriterSet() const { return _writers != nullptr; }

	// Writes oblique line (positive or negative strand)
	//	@param reverse: true if tag is reversed (neg strand)
	//	@param start: start (zero value) position
	//	@param itStop: stop (max value) iterator
	void WriteIncline(BYTE reverse, chrlen start, coviter& itStop);

	// Writes curve for TOTAL
	//	@param start: start (left) position
	//	@param vals: values to write
	void WriteChromData(chrlen start, const Values& vals)
	{
		(_writers->_files)[TOTAL]->WriteChromData(_cID, start, vals);
	}
};

//===== TXT FILE AS DUMP

class TxtOutFile
{
	FILE* _file;
public:
	TxtOutFile(const char* name) { _file = fopen(name, "w"); }
	~TxtOutFile() { fclose(_file); }

	template<typename... Args>
	void Write(const char* format, Args ... args) { fprintf(_file, format, args...); }
};

#endif	// MY_DEBUG

//===== DATA

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

// Inclined line represents average derivative of the area of read's coverage growth/decline,
// as a resul of Linear Regression
struct Incline
{
	chrlen	Pos;		// position of the incline on the x-axis (chromosome's position)
	chrlen	TopPos;		// position of the top point of the incline
	float	Deriv;		// derivative (tangent of the angle of incline)

	bool Valid() const {
		return
			Pos					// zero Pos is the result of a failed regression calculation
			&& Deriv < 1		// 1 (vertical line) is useless
			&& Deriv >= 0.01;	// too small angle is useless
	}

	void Clear() { Deriv = 0; Pos = 0; }

	bool operator != (const Incline& incl) const { return Pos != incl.Pos || TopPos != incl.TopPos; }

#ifdef MY_DEBUG
	void Print() const { printf("%d %d, Deriv: %-2.2f\n", Pos, TopPos, Deriv); }

private:
	static shared_ptr<TxtOutFile> OutFile;
	static chrid cID;

public:
	static void SetOutFile(chrid cid, const char* fname) { cID = cid;  OutFile.reset(new TxtOutFile(fname)); }

	void Write();
#endif
};

class TreatedCover : public AccumCover
{
#ifdef MY_DEBUG
	static OSpecialWriter* SplineWriter;	// spline writer to save local splines
#endif

public:
#ifdef MY_DEBUG
	static void SetSpecialWriter(OSpecialWriter& splineWriter) { SplineWriter = &splineWriter; }
	static bool WriteDelim;
#endif
	coval GetMaxVal() const;

	// Returns cover mass centre for given region
	//	@param rng: potential region
	//chrlen GetRegionCentre(const CoverRegion& rng) const;

	// Calculates Linear Regression
	//	@param it: start iterator
	//	@param itStop: end iterator
	//	@param op: strand operations
	//	@incline[out]: resulting inclined line
	void LinearRegr(coviter it, const coviter& itStop, const StrandOp& op, Incline& incline) const;

	// Sets spline of the instance between start-end positions
	//	@param spliner[in]: spliner that does the work
	//	@param startPos[in,out]: start position; returns modified (real) position
	//	@param endPos[in]: end position
	//	@param vals[out]: step (one bp) splineed values
	void SetLocalSpline(
		SSpliner<coval>& spliner,
		chrlen& startPos,
		chrlen endPos,
		Values& vals
	) const;
};

// 'CombCover' keeps the chromosome covers and optionally the set of writers these covers to file.
class CombCover : public OrderedData<TreatedCover, BedGrWriter>
{
	string	_outFName;				// common name of output files
	UniBedReader* _file;

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
	void AddExtRead(const Region& read, bool reverse);

	void Fill(const Reads& reads);
};

using tChromsFreq = map<chrid, BYTE>;

// Wrapper for initializing a coverage from a file
class CombCoverReader
{
	CombCover&	 _cover;
	tChromsFreq& _chrFreq;
	eStrand		 _strand;
	UniBedReader _file;

public:
	// Constructor-initializer
	//	@param fName[in]: file name
	//	@param cSizes[in]: chrom sizes
	//	@param cover[out]: coverage to initialize
	//	@param chrFreq[out]: chrom reading frequency
	//	@param strand[in]: strand
	CombCoverReader(
		const char* fName,
		ChromSizes& cSizes,
		CombCover& cover,
		tChromsFreq& chrFreq,
		eStrand strand
	)
		: _cover(cover), _chrFreq(chrFreq), _strand(strand)
		, _file(fName, FT::eType::BGRAPH, &cSizes, 4, 0, eOInfo::LAC, Verb::Level(Verb::DBG), false, true)
	{
		Verb::PrintMsg(Verb::DBG);
		_file.Pass(*this);
	}

	// treats current item
	//	@returns: true if item is accepted
	bool operator()() {
		_cover.AddNextRegion(_strand, _file.ItemRegion(), coval(_file.ItemValue()));
		return true;
	}

	// Closes current chrom, open next one
	//	@param cID: current chrom ID
	//	@param cLen: chrom length
	//	@param cnt: current chrom items count
	//	@param nextcID: next chrom ID
	void operator()(chrid cID, chrlen cLen, size_t cnt, chrid nextcID) {
		_cover.SetChrom(nextcID);
		_chrFreq[nextcID]++;
	}

	// Closes last chrom
	//	@param cID: last chrom ID
	//	@param cLen: chrom length
	//	@param cnt: last chrom items count
	//	@param tCnt: total items count
	void operator()(chrid cID, chrlen cLen, size_t cnt, size_t)	{}
};

//=== COVER REGION & COLLECTION

// 'CoverRegion' represents potential region
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
	void SetPotentialRegions(const TreatedCover& cover, chrlen capacity, coval cutoff);

#ifdef MY_DEBUG
	// Prints frequency distribution of potentail regions value ('score')
	void PrintScoreDistrib() const;
#endif

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
	//	@returns: true in case of empty instance (and print message), false otherwise
	bool SetPotentialRegions(const DataSet<TreatedCover>& cover, chrlen cLen, coval cutoff, bool noMultiOverl = false);

	// Returns mean fragment length based on cover regions mass centre comparison
	//	@param cover: fragment coverage
	//fraglen GetFragMean(const DataSet<TreatedCover>& cover) const;

	void Clear()
	{
		StrandData(POS).clear();
		StrandData(NEG).clear();
	}
};


//=== SPLINE & COLLECTION

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
	//	@param redifineRgns: if true then redifine regions position
	//	@param splineBase: half-length of spliner moving window
	void BuildRegionSpline(const TreatedCover& cover, const CoverRegion& rgn, bool redifineRgns, fraglen splineBase);

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

	// Sets a single consecutive potential region number to each overlapped not empty spline
	void Numerate();

	// Prints potential regions before and after selection
	void PrintStat(chrlen clen) const;
};

// Datset of maped collections of float values
class DataValuesMap : public DataSet<ValuesMap>
{
public:
	// Filters and fill spline curve by PE read cover within potential regions
	//	@param cover: raw read/frag cover
	//	@param rgns: potential regions
	//	@param redifineRgns: if true then redifine regions position
	//	@param splineBase: half-length of spliner moving window
	void BuildSpline(const DataSet<TreatedCover>& cover, const DataCoverRegions& rgns, bool redifineRgns, fraglen splineBase = ReadSplineBASE);

	// Calculates the deviation from the default average fragment length
	float GetPeakPosDiff() const;

	// Resets non overlapping spline value
	void EliminateNonOverlaps() { Data()->EliminateNonOverlaps(); }

	// Sets a single consecutive potential region number to each overlapped not empty spline
	void Numerate() { Data()->Numerate(); }

	// Prints potential regions before and after selection
	void PrintStat(chrlen clen) const { Data()->PrintStat(clen); }

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
#ifdef MY_DEBUG
	static OSpecialWriter* LineWriter;	// line writer to save inclined
#endif
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
		//	@param it: forward/reverse cover iterator
		template<typename T>
		void Set(BYTE reverse, T it)
		{
			// set positions & values
			_pos[!reverse] = it->Start();	// start: right for forward, left for reverse
			_pos[reverse] = it->End();		// end: left for forward, right for reverse
			it->GetValues(_val);
			// merge successive raises
			if (_val0[reverse] && _val[!reverse] / _val0[reverse] >= 0.9)
				_pos[0] = _pos0[0];
		}

		// Saves previous positions/values (by swapping)
		void Retain()
		{
			//memcpy(_pos0, _pos, sizeof(_pos));
			//memcpy(_val0, _val, sizeof(_val));
			std::swap(_pos, _pos0);
			std::swap(_val, _val0);
		}
	};

	float _maxVal = 0;
	chrlen _rgnNumb;	// potential region number

	void PushIncline(
		chrlen rgnNumb,
		BYTE reverse,
		const TracedPosVal& posVal,
		const TreatedCover& cover,
		vector<Incline>& inclines
	) const;

public:
#ifdef MY_DEBUG
	static void SetSpecialWriter(OSpecialWriter& lineWriter) { LineWriter = &lineWriter; }
#endif
	float	MaxVal()	const { return _maxVal; }
	// Returns potential region number
	chrlen	RgnNumb()	const { return _rgnNumb; }

	//using cIter = vector<ValuesMap>::const_iterator;
	//using rIter = vector<ValuesMap>::reverse_iterator;
	//typedef vector<ValuesMap>::iterator iter_type;

	// Collects forward inclined lines
	//	@param rCover[in]: read coverage
	//	@param inclines[out]: filled collection of forward inclined lines
	void CollectForwardInclines	(const TreatedCover& cover, vector<Incline>& inclines) const;

	// Collects reversed inclined lines
	//	@param rCover[in]: read coverage
	//	@param inclines[out]: filled collection of reversed inclined lines
	void CollectReverseInclines(const TreatedCover& cover, vector<Incline>& inclines) const;

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
	//	@param factor: 1 for forward reads, -1 for reversed ones
	//	@param splines: splines from which derivatives are built
	void BuildDerivs(int factor, const ValuesMap& splines);

#ifdef MY_DEBUG
	void Print(eStrand strand, chrlen stopPos = 0) const;
#endif
};

class DataBoundsValuesMap : public DataSet<BoundsValuesMap>
{
public:
	// Fills the instance with spline derivatives
	//	@param splines: Read coverage splines
	void Set(const DataValuesMap& splines)
	{
		StrandData(POS).BuildDerivs(StrandOps[0].Factor, splines.StrandData(POS));
		StrandData(NEG).BuildDerivs(StrandOps[1].Factor, splines.StrandData(NEG));
	}

#ifdef MY_DEBUG
	void Print(chrlen stopPos = 0) const
	{
		StrandData(POS).Print(POS, stopPos);
		StrandData(NEG).Print(NEG, stopPos);
	}
#endif
};

//=== BINDING SITES DATA 

struct BS_PosVal
{
	BYTE	Reverse;
	chrlen	RgnNumb;	// potential region number
	float	Score = 1;

	// Constructor
	//	@param reverse: 0 for firect, 1 for reverse
	//	@param rgnNumb: potential region number
	//	@param incln: inclined line
	BS_PosVal(BYTE reverse, chrlen rgnNumb) : Reverse(reverse), RgnNumb(rgnNumb) {}
};

class BS_map : public map<chrlen, BS_PosVal>
{
public:
	using iter = map<chrlen, BS_PosVal>::iterator;
	using citer = map<chrlen, BS_PosVal>::const_iterator;

private:
	iter _lastIt;	// last inserted iterator (for the left bounds only)

	// Inserts BS position (bound)
	//	@param reverse: 0 for firect (right bounds), 1 for reverse (left bounds)
	//	@param rgnNumb: potential region number
	//	@param incl: inclined line
	void AddPos(BYTE reverse, chrlen rgnNumb, const Incline& incl);

	// Inserts BS positions (left/right bounds)
	//	@param reverse: 0 for firect (right bounds), 1 for reverse (left bounds)
	//	@param rgnNumb: potential region number
	//	@param inclines: forward/reversed (right/left) inclined lines
	void AddBounds(BYTE reverse, chrlen rgnNumb, vector<Incline>& inclines);

	// Fills the instance with recognized left/right BS positions (bounds)
	//	@param reverse[in]: 0 for firect (right bounds), 1 for reverse (left bounds)
	//	@param derivs[in]: derivatives
	//	@param rCover[in]: read coverage
	void SetBounds(BYTE reverse, const BoundsValuesMap& derivs, const TreatedCover& rCover);

public:
	// positioned value
	struct PosValue {
		iter	Iter;
		float	Val;

		PosValue(iter it, float val) : Iter(it), Val(val) {}
	};

	// Fills the instance with recognized binding sites
	//	@param derivs: derivatives
	//	@param rCover: read coverage
	void Set(const DataBoundsValuesMap& derivs, const DataSet<TreatedCover>& rCover)
	{
		SetBounds(0, derivs.StrandData(POS), rCover.StrandData(POS));
		_lastIt = begin();
		SetBounds(1, derivs.StrandData(NEG), rCover.StrandData(NEG));
	}

	// Brings the instance to canonical order of placing BS borders
	void Refine();

	// Sets score for each binding sites and normalizes it
	//	@param fragCovers: fragment total coverage
	void SetScore(const DataSet<TreatedCover>& fragCovers);

	void NormalizeScore();

	// Expands too narrow BSs to the minimum acceptable width
	void NormalizeBSwidth();

	void PrintStat() const;

	// Applies lambda to each potential region of binding sites, passing borders collection
	template<typename F>
	void DoExtend(F&& lambda)
	{
		vector<PosValue> VP[2];	// 0 - forward, 1 - reversed

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

	// Applies lambda to each potential region denotes the binding site, passing start-end iterators
	template<typename F>
	void DoBasic(F&& lambda) const
	{
		for (auto it0 = begin(), it = next(it0); it != end(); it0++, it++)
			if (it->second.Score && !it->second.Reverse && it0->second.Reverse)
				lambda(it0, it);
	}
	template<typename F>
	void DoBasic(F&& lambda)
	{
		for (auto it0 = begin(), it = next(it0); it != end(); it0++, it++)
			if (it->second.Score && !it->second.Reverse && it0->second.Reverse)
				lambda(it0, it);
	}

#ifdef MY_DEBUG
	// Prints "unsorted" scores
	void CheckScoreHierarchy();

	// Prints BS width distribution
	void PrintWidthDistrib() const;

	// Prints positions
	//	@param outFName: name of file to print, or NULL to print to console
	//	@param selected: if true then prints position with unzero score 
	//	@param stopPos: max printed position or all by default 
	void Print(chrid cID, const char* outFName, bool selected, chrlen stopPos = 0) const;
#endif
};

// 'BedWriter' implements methods for writing in ordinary bed format
class BedWriter : public RegionWriter
{
	static bool rankScore;	// true if scores in the main bed should be normalized relative to 1000

	void LineAddIntScore	(float score, bool delim) { LineAddInt(int(round(score * 1000)), delim); }
	void LineAddFloatScore	(float score, bool delim) { LineAddFloat(score, delim); }

	typedef void (BedWriter::* tAddScore)(float, bool);
	static tAddScore fLineAddScore;

	void LineAddScore(float score, bool delim) { assert(fLineAddScore); (this->*fLineAddScore)(score, delim); }

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
	void WriteChromData(chrid cID, BS_map& bss);

	// Writes extended BSs; adds a gray color to feature with zero score
	//	@param cID: chrom ID
	//	@param bss: binding sites (const in fact)
	void WriteChromExtData(chrid cID, BS_map& bss);

	// Writes Regions of interest (extended binding sites) as a 4-fields features
	//	@param cID: chrom ID
	//	@param bss: binding sites
	void WriteChromROI(chrid cID, const BS_map& bss);
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

	void WriteChromData(chrid cID, BS_map& bss)
	{
		_main.WriteChromData(cID, bss);
		_ext.WriteChromExtData(cID, bss);
		_roi.WriteChromROI(cID, bss);
	}
};

using OBS_Map = OrderedData<BS_map, BedWriters>;
//class OBS_Map : public OrderedData<BS_map, BedWriter>
//{
//public:
//	// Primer constructor
//	//	@param cSizes: chrom sizes
//	//	@param dim: number (dimension) of data
//	//	@param write: if true then data sould be save to output file
//	//	@param fname: name of output file
//	//	@param descr: track decsription in declaration line
//	OBS_Map(const ChromSizes& cSizes, BYTE dim, bool write, const string& fname, const char* descr)
//		: OrderedData<BS_map, BedWriter>(cSizes, dim, write, TrackFields(fname, descr)) {}
//};


//===== WIG WRITERS

// 'FixWigWriter' implements methods for writing in wiggle_0 (fixed step) format
class FixWigWriter : WigWriter
{
public:
	// Creates new instance for writing cover to wiggle_0 file.
	//	@param strand: strand
	//	@param fields: BED/WIG track fields
	FixWigWriter(eStrand strand, const TrackFields& fields)
		: WigWriter(FT::eType::WIG_FIX, strand, fields, 3) {}

	void WriteChromData(chrid cID, const ValuesMap& vals);
};

class FixWigWriterSet : WigWriter
{
public:
	// Creates new instance for writing cover to wiggle_0 file.
	//	@param strand: strand
	//	@param fields: BED/WIG track fields
	FixWigWriterSet(eStrand strand, const TrackFields& fields)
		: WigWriter(FT::eType::WIG_FIX, strand, fields, 3) {}

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
