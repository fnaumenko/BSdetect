#include "Treatment.h"
#include <algorithm>

//const float PI = 3.14159265F;

//const eCurveType CurveTYPE = eCurveType::SMOOTH;
const eCurveType CurveTYPE = eCurveType::ROUGH;

bool Glob::IsPE = false;
bool Glob::FragLenUndef = true;
readlen Glob::ReadLen = 0;
fraglen Glob::FragLen = FragDefLEN;
fraglen Glob::ROI_ext = 500;
 
//===== IGVlocus

#ifdef MY_DEBUG
class IGVlocus
{
	const string _chrom;
	mutable char _buf[2 * 10 + 5 + 2];	// 2 * max position length + Chrom::MaxAbbrNameLength + 2 separators

	// Prints IGV locus to inner buffer
	//	@return: the total number of characters written
	chrlen NPrint(chrlen start, chrlen end) const
	{
		return sprintf(_buf, "%s:%d-%d", _chrom.c_str(), start - Glob::ROI_ext, end + Glob::ROI_ext);
	}

public:
	IGVlocus(chrid cID) : _chrom(Chrom::AbbrName(cID)) {}

	// Returns inner buffer
	//const char* Buff() const { return _buf; }

	// Prints IGV locus to inner buffer
	//	@return: inner buffer
	const char* Print(chrlen start, chrlen end) const { NPrint(start, end); return _buf; }

	const char* Print(chrlen pos) const { return Print(pos, pos); }
};
#endif

//===== Verb

const char* Verb::ValTitles[] = { "SL","RES","RT","DBG" };
const char* Verb::ValDescr = "set verbose level:\n?  -\tsilent mode (show critical messages only)\n? -\tshow result summary\n?  -\tshow run-time information\n? -\tshow debug messages";
Verb::eVerb Verb::_level;

void Verb::PrintMsg(eVerb level, const char* msg)
{
	if (Level(level))
		if (msg)	printf("%s\n", msg);
		else		printf("\n");
}

//void Verb::PrintMsgVar(eVerb level, const char* format, ...)
//{
//	if (Level(level)) {
//		va_list argptr;
//		va_start(argptr, format);
//		vfprintf(stdout, format, argptr);
//		va_end(argptr);
//	}
//}

#ifdef MY_DEBUG
//===== OSpecialWriter

void OSpecialWriter::WriteIncline(BYTE reverse, chrlen start, coviter& itStop)
{
	assert(_cID != Chrom::UnID);

	const chrlen stop = itStop->first;
	const int ptCnt = int(stop) - int(start);
	int8_t k = 1;
	if (!reverse) {
		k = -1;
		start = stop;
		--itStop;
	}
	(_writers->_files)[reverse]->WriteIncline(_cID, start, k * ptCnt, float(itStop->second) / ptCnt);
}

//===== Incline

shared_ptr<TxtOutFile> Incline::OutFile = nullptr;
chrid Incline::cID;

void Incline::Write()
{
	if (OutFile) {
		IGVlocus locus(cID);
		OutFile->Write("%1.3f\t%d\t%s\n", Deriv, Pos, locus.Print(Pos));
	}
}
#endif

//===== TreatedCover

#ifdef MY_DEBUG
OSpecialWriter* TreatedCover::SplineWriter = nullptr;
bool TreatedCover::WriteDelim = false;
#endif

coval TreatedCover::GetMaxVal() const
{
	coval val = 2;

	for (const auto& item : *this)
		if (val < item.second)
			val = item.second;
	return val;
}

//chrlen TreatedCover::GetRegionCentre(const CoverRegion& rng) const
//{
//	//const chrlen centre = itStart->first + (itEnd->first - itStart->first) / 2;
//	const chrlen centre = rng.itStart->first + rng.Length() / 2;
//	int32_t area = 0, estimHalfArea = 0;
//	coviter itC;	// centre iterator
//
//	// ** define area and estimated halfArea
//	for (coviter it0 = rng.itStart, it = next(it0); it0 != rng.itEnd; it0++, it++) {
//		area += it0->second * (it->first - it0->first);
//		if (!estimHalfArea && it0->first >= centre) {
//			estimHalfArea = area;	// assigned just once
//			itC = it0;
//		}
//	}
//	// ** find mass centre
//	area /= 2;			// from now on half the area
//	if (estimHalfArea < area) {
//		for (coviter it = next(itC); estimHalfArea < area; itC++, it++)
//			estimHalfArea += itC->second * (it->first - itC->first);
//
//		if (estimHalfArea == area || itC->first - prev(itC)->first == 1)
//			return itC->first;
//		// precise fit
//		int32_t d = 0;
//		for (coval val = prev(itC)->second; estimHalfArea > area; d++)
//			estimHalfArea -= val;	// step for one pos
//		return itC->first - d;
//	}
//	else {
//		for (coviter it = prev(itC); estimHalfArea > area; itC--, it--)
//			estimHalfArea -= it->second * (itC->first - it->first);
//
//		if (estimHalfArea == area || itC->first - prev(itC)->first == 1)
//			return itC->first;
//		// precise fit
//		int32_t d = 0;
//		for (coval val = itC->second; estimHalfArea < area; d++)
//			estimHalfArea += val;	// step for one pos
//		return itC->first + d;
//	}
//}

void TreatedCover::LinearRegr(coviter it, const coviter& itStop, const StrandOp& op, Incline& incline) const
{
	const int shift = op.Factor * itStop->first;
	float sumX = 0, sumX2 = 0, sumY = 0, sumXY = 0;
	coval PtCnt = 0;	// number of points along which the incline was drawn

	incline.Clear();
	for (; it != itStop; op.Next(it), PtCnt++) {
		const int x = shift - op.Factor * it->first;
		const coval y = op.GetPrev(it)->second;
		sumX += x;
		sumX2 += x * x;
		sumY += y;
		sumXY += x * y;
	}
	if (PtCnt < 2)	return;

	incline.Deriv = 
		-(PtCnt * sumXY - sumX * sumY) /	// numerator
		(PtCnt * sumX2 - sumX * sumX);		// denominator
	//angl = coeff * 180 / PI;	if (angl < 0) angl = -angl;

	float x = (sumY + incline.Deriv * sumX) / (PtCnt * incline.Deriv);
	incline.Pos = itStop->first - op.Factor * chrlen(round(x));
	incline.TopPos = itStop->first;

#ifdef MY_DEBUG
	incline.Write();
#endif
}

void TreatedCover::SetLocalSpline(SSpliner<coval>& spliner, chrlen& startPos, chrlen endPos, Values& vals) const
{
	startPos -= spliner.SilentLength();
	endPos += spliner.SilentLength();

	coviter itCov = prev(upper_bound(startPos));	// no check for end() since cover is well-defined in this region
	coviter itCovEnd = itCov;

	// set itCovEnd
	for (itCovEnd++; itCovEnd->first <= endPos; itCovEnd++);

	// build spline
	chrlen pos = itCov->first + 1;
	for (auto it = next(itCov); itCov != itCovEnd; itCov++, it++)	// loop through the cover
		for (; pos <= it->first; pos++) {							// loop through positions between iterators
			float val = spliner.Push(itCov->second);
			if (val) {
				if (!vals.MaxVal())
					startPos = spliner.CorrectX(pos);	// just once
				vals.AddValue(val);
			}
		}
#ifdef MY_DEBUG
	if (SplineWriter && SplineWriter->IsWriterSet())
		SplineWriter->WriteChromData(startPos, vals);
#endif
	spliner.Clear();
}


//===== CombCover

void CombCover::AddExtRead(const Region& read, bool reverse)
{
	Region frag(read, Glob::FragLen, reverse);

	_data->TotalData().AddRegion(frag);			// total frag coverage
	AddRead(frag, reverse);
}

void CombCover::Fill(const Reads& reads)
{
	for (BYTE s : {0, 1})
		for (auto& rd : reads.GetReads(s))
			AddExtRead(rd, s);
}


// marks as empty non-overlapping and weak-overlapping regions
template<typename T>
void EliminateNonOverlapsRegions(T rgns[2], fraglen minOverlapLen)
{
	const typename T::iterator itEnd[2]{ rgns[0].end(), rgns[1].end() };
	typename T::iterator it[2]{ rgns[0].begin(), rgns[1].begin() };
	BYTE s, suspend = 0;	// if 1st or 2nd bit is set then pos or neg region is suspended and analyzed in the next pass

	auto start	= [&it](BYTE s) -> chrlen { return T::Start(it[s]); };
	auto end	= [&it](BYTE s) -> chrlen { return T::End(it[s]); };
	auto isWeak	= [&it](BYTE s) -> bool   { return T::IsWeak(it[s]); };
	auto markAsEmpty = [&it](BYTE s) { T::MarkAsEmpty(it[s]); };

	// unconditional close of the region; always applied to the left region
	auto closeRgn = [&](BYTE s) {
		markAsEmpty(s);
		suspend &= ~(1 << s);	// reset suspend s
		it[s]++;
	};
	auto closeLastRgns = [&it, &itEnd, &closeRgn](BYTE s) {
		while (it[s] != itEnd[s])
			closeRgn(s);
	};

	while (it[0] != itEnd[0] && it[1] != itEnd[1]) {
		s = start(0) > start(1);	// s denotes the index of the left started region: 0 - pos, 1 - neg
		if (end(s) > start(!s) + minOverlapLen) {	// strong intersection
			/***************************************
				-----    ?????	more regions?		left overlaps
			+++++++++++++		suspended...		left overlaps
				---------		suspended...		right overlaps
			+++++++	  ?????		more regions?		right overlaps
			***************************************/
			bool valid = true;
			if (isWeak(s))	closeRgn(s),  valid = false;
			if (isWeak(!s))	closeRgn(!s), valid = false;
			if (valid) {
				if (end(s) < end(!s))	s = !s;
				it[!s]++;
				suspend |= 1 << s;	// set suspend s
			}
		}
		else {							// weak intersection or no one
			// conditional close of the region
			if (suspend)	suspend = 0;
			else			markAsEmpty(s);
			// close remaining complementary regions
			if (++it[s] == rgns[s].end()) {
				closeLastRgns(!s);
				break;			// no need to check in while()
			}
		}
	}

	// close 'out of scope' regions
	if ((s = it[1] != itEnd[1]) || it[0] != itEnd[0]) {
		it[s]++;
		closeLastRgns(s);
	}
}

// marks as empty multi-overlapping (more than 2) regions
template<typename T>
void EliminateMultiOverlapsRegions(T rgns[2])
{
	const typename T::iterator End[]{ rgns[0].end(), rgns[1].end() };
	typename T::iterator it0[2]{ rgns[0].begin(), rgns[1].begin() };
	typename T::iterator it[2] { next(it0[0]), next(it0[1]) };

	auto start	= [&it] (BYTE s) -> chrlen	{ return T::Start(it[s]); };
	auto end	= [&it0](BYTE s) -> chrlen	{ return T::End(it0[s]); };
	auto isWeak = [&it0](BYTE s) -> bool	{ return T::IsWeak(it0[s]); };
	auto markAsEmpty = [&it, &rgns](BYTE s) { T::MarkAsEmpty(it[s]); };
	auto markBothAsEmpty = [&it0, &rgns]()	{
		T::MarkAsEmpty(it0[0]);
		T::MarkAsEmpty(it0[1]);
	};
	auto jumpOverWeaks = [&](BYTE s) {
		if (isWeak(s)) {
			it0[s]++;
			if (isWeak(s))	it0[s]++;
			it[s] = it0[s] != End[s] ? next(it0[s]) : it0[s];
			return true;
		}
		return false;
	};

	while (it[0] != End[0] && it[1] != End[1]) {
		if (jumpOverWeaks(0))	continue;
		if (jumpOverWeaks(1))	continue;

		BYTE s = end(0) > end(1);	// s denotes the index of the left ended region: 0 - pos, 1 - neg
		if (start(s) <= end(!s)) {	// in this case it[s] (next) region is not weak
			markBothAsEmpty();
			markAsEmpty(s);
			it0[s]++;
		}
		it0[s]++;	it[s] = next(it0[s]);
		it0[!s]++;	it[!s]++;
	}
}

// prints regions before and after selection
template<typename T>
void PrintRegionStats(const T* rgns, chrlen chrLen, bool strands = true)
{
	const char* format[] = {
		"%s: %4d (%2.2f%%) %4d (%2.2f%%)\n",	// features
		"%s: %4d (%1.3f%%) %4d (%1.3f%%)\n"	// splines
	};
	const char* titles[] = {
		"POTENTIAL REGIONS:",
		"SPLINED REGIONS:"
	};
	const bool isFeatures = is_same<T, ValuesMap>::value;

	const char* title = "strand      received      selected";
	printf("\n%s\n%s\n", titles[isFeatures], title);
	PrintSolidLine(USHORT(strlen(title) + 2));
	for (BYTE s = 0; s < 1 + strands; s++) {
		chrlen rawLen = 0, refineLen = 0;
		const auto& rgn = rgns[s];
		const auto itEnd = rgn.end();
		size_t realCnt = 0;

		for (auto it = rgn.begin(); it != itEnd; it++) {
			const chrlen len = T::Length(it);
			rawLen += len;
			if (T::IsNotEmpty(it)) refineLen += len, realCnt++;
		}
		printf(format[isFeatures], sStrandTITLES[s + strands],
			rgn.size(), Percent(rawLen, chrLen),
			realCnt, Percent(refineLen, chrLen));
	}
}

//===== CoverRegions

void CoverRegions::SetPotentialRegions(const TreatedCover& cover, chrlen capacity, coval cutoff)
{
	const fraglen minLen = 3 * (Glob::FragLen / 2);		//1.5 fragment lengths
	chrlen	start = 0, end = 0;
	coviter itStart, itEnd;

	this->reserve(capacity);
	for (auto it0 = cover.cbegin(), it = it0; it != cover.end(); it0 = it++)
		if (it->second >= cutoff || it0->second >= cutoff) {	// look for summit
			if (!start)
				start = (itStart = it)->first;	// set after prev region processed only
			end = (itEnd = it)->first;
		}
		else {
			if (start && end - start > minLen) {
				// find raw summit
				coval val = itStart->second;
				for (auto it = next(itStart); it != itEnd; it++)
					if (val < it->second)
						val = it->second;
				this->emplace_back(itStart, itEnd, val);
			}
			start = 0;
		}
}

#ifdef MY_DEBUG
void CoverRegions::PrintScoreDistrib() const
{
	map<coval, chrlen> freq;

	for (const auto& rgn : *this)
		freq[rgn.value]++;
	printf("\nREGIONS LENGTH FREQUENCY\n");
	printf("length freq\n");
	for (const auto& item : freq)
		printf("%5d %d\n", item.first, item.second);
}
#endif

//===== DataCoverRegions

bool DataCoverRegions::SetPotentialRegions(const DataSet<TreatedCover>& cover, chrlen cLen, coval cutoff, bool noMultiOverl)
{
	chrlen capacity = cLen / (Glob::FragLen * 100);
	if(Glob::IsPE)
		TotalData().SetPotentialRegions(cover.TotalData(), capacity, cutoff);
	else {
		StrandData(POS).SetPotentialRegions(cover.StrandData(POS), capacity, cutoff);
		StrandData(NEG).SetPotentialRegions(cover.StrandData(NEG), capacity, cutoff);

		EliminateNonOverlapsRegions<CoverRegions>(Data(), Glob::FragLen);
		if (noMultiOverl)
			EliminateMultiOverlapsRegions<CoverRegions>(Data());
		if (Verb::Level(Verb::DBG))
			PrintRegionStats<CoverRegions>(Data(), cLen);
	}
	if (Empty()) { Verb::PrintMsg(Verb::CRIT, "No enriched regions found"); return true; }
	return false;
}

//fraglen DataCoverRegions::GetFragMean(const DataSet<TreatedCover>& cover) const
//{
//	auto& pData = StrandData(POS);
//	auto& nData = StrandData(NEG);
//	auto& pCover = cover.StrandData(POS);
//	auto& nCover = cover.StrandData(NEG);
//	uint32_t negCnt = 0;
//	vector<int16_t> diffs;
//	IGVlocus locus(0);
//
//	// get the difference of the splines maximums 
//	diffs.reserve(pData.size());
//	for (auto itP = pData.begin(), itN = nData.begin(); itP != pData.end() && itN != nData.end(); ) {
//		if (!itP->value) { itP++; continue; }
//		if (!itN->value) { itN++; continue; }
//
//		chrlen pCentre = pCover.GetRegionCentre(*itP);
//		chrlen nCentre = nCover.GetRegionCentre(*itN);
//		int16_t diff = pCentre - nCentre;
//
//		negCnt += diff < 0;
//		diffs.push_back(diff);
//		printf("%d\t%d\t%d\t%s\n", diff, pCentre, nCentre, locus.Print(pCentre));
//
//		itP++, itN++;
//	}
//
//	// average the difference, get frag length
//	bool mostPositive = negCnt < diffs.size() / 2;
//	auto CompareVal = mostPositive ? &PositiveVal : &NegativeVal;
//	int sum = 0;
//
//	for (auto diff : diffs)
//		if (CompareVal(diff))
//			sum += diff;
//
//	return FragDefLEN - round(float(sum) / (mostPositive ? (diffs.size() - negCnt) : negCnt));
//}


//===== Values

Values::Values(Values& vals, const Region& rgn) 
	: vector<float>(vals.begin() + rgn.End, vals.end())
	, _maxVal(vals._maxVal)		// shared _maxVal doesn't matter
{
	vals.resize(rgn.Start);
}

void Values::AddValue(float val)
{
	if (_maxVal < val)	_maxVal = val;
	push_back(val);
}

void Values::AddValues(const Values& vals)
{
	if (_maxVal < vals._maxVal)	_maxVal = vals._maxVal;
	insert(end(), vals.begin(), vals.end());
}

void Values::GetMaxValPos(chrlen startPos, vector<chrlen>& pos) const
{
	float val0 = front();
	bool increase = false;
	USHORT equalCnt = 0;	// equal value counter

	for (auto it = next(begin()); it != end(); val0 = *it, it++, startPos++)
		if (*it == val0)
			equalCnt++;
		else {
			if (*it > val0)			// increase value
				increase = true;
			else if (increase) {	// decrease value
				pos.push_back(startPos - equalCnt / 2);	// correct startPos for the 'flat' summit
				increase = false;
			}
			equalCnt = 0;
		}
}

#ifdef MY_DEBUG
void Values::Print(bool prValues) const
{
	printf("MaxVal: %2.2f\n", _maxVal);
	if (prValues) {
		printf("Values: ");
		if (size())
			for (auto v : *this)	printf("%2.2f ", v);
		else printf("none");
		printf("\n");
	}
}
#endif


//===== ValuesMap

#ifdef MY_DEBUG
void ValuesMap::Print(chrid cID, BYTE reverse, chrlen stopNumb) const
{
	IGVlocus locus(cID);

	printf("SPLINES %s\n", sStrandTITLES[reverse + 1]);
	printf(" N start\tend\tval\tIGV view\n");
	for (const auto& x : *this) {
		if (stopNumb && x.second.RgnNumb > stopNumb)	break;
		if (x.second.MaxVal()) {
			chrlen end = x.first + x.second.Length();
			printf("%2d %d\t%d\t%2.2f\t%s\n",
				x.second.RgnNumb, x.first, end, x.second.MaxVal(), locus.Print(x.first, end));
		}
	}
}
#endif

void ValuesMap::BuildRegionSpline(const TreatedCover& cover, const CoverRegion& rgn, bool redifineRgns, fraglen splineBase)
{
	assert(Glob::ReadLen);
	coviter it0;	// at the beginning the start it, then used as a variable
	coviter itEnd;	// the end it

	// set up start/end iterators
	if (redifineRgns) {
		it0 = cover.upper_bound(rgn.itStart->first);
		itEnd = cover.lower_bound(rgn.itEnd->first);		// !!! optimize

		if (--it0 == cover.end())	it0 = cover.begin();
		if (itEnd == cover.end())				itEnd--;
		else if (next(itEnd) != cover.end())	itEnd++;
	}
	else {
		it0 = rgn.itStart;
		itEnd = next(rgn.itEnd) != cover.end() ? next(rgn.itEnd) : rgn.itEnd;
	}

	// *** spline via covmap local copy, filtering unsignificant splines
	chrlen pos = it0->first + 1, newPos = 0;
	bool isZeroBefore = true;
	SSpliner<coval> spliner(CurveTYPE, splineBase);
	Values vals;

	// adds spline, eliminating unsignificant one
	auto addDecentRgn = [this,&vals,&newPos]() {
		if (vals.MaxVal())
			if (vals.MaxVal() > 2.)
				this->AddRegion(newPos, vals);
			else
				vals.Clear();
	};

	for (auto it = next(it0); it != itEnd; it0++, it++) {	// loop through the cover

		if (isZeroBefore)
			// skip 2 duplicated reads or standalone read or contiguous reads
			if (!it->second) {
				if (it0->second <= 2) {
					if (++it == itEnd)	break;
					it0++;
				}
			}
			// skip 2 overlapping reads
			else if (it->second == 2) {
				auto it1 = next(it); if (it1 == itEnd)	break;
				if (it1->second == 1) {
					if (++it1 == itEnd)		break;
					if (!it1->second) {
						it = move(it1);		// !!! swap?
						if(++it == itEnd)	break;
						advance(it0, 3);
					}
				}
			}

		// treat other reads
		for (; pos <= it->first; pos++) {			// loop through positions between iterators
			float val = spliner.Push(it0->second);
			if (val) {
				if (!vals.MaxVal())
					newPos = spliner.CorrectX(pos);	// start new spline
				vals.AddValue(val);
			}
			else 
				addDecentRgn();						// end new spline
		}

		isZeroBefore = !it0->second;
	}
	addDecentRgn();
}

void ValuesMap::AddRegion(chrlen pos, Values& vals)
{
	_maxVal = vals.MaxVal();	// duplicated because _maxVal is also used as 'empty' sign
	emplace(pos, move(vals));
	vals.Reserve();
}

void ValuesMap::BuildSpline(const TreatedCover& cover, const CoverRegions& rgns, bool redifRgns, fraglen splineBase)
{
	for (const auto& rgn : rgns)
		if (rgn.value)
			BuildRegionSpline(cover, rgn, redifRgns, splineBase);
}

void ValuesMap::EliminateNonOverlaps()
{
	EliminateNonOverlapsRegions<ValuesMap>(this, SSpliner<coval>::SilentLength(CurveTYPE, ReadSplineBASE));
}

void ValuesMap::Numerate()
{
	const fraglen overLen = SSpliner<coval>::SilentLength(CurveTYPE, ReadSplineBASE) + 10;	// minimum overlap length; 10 - min BS length !!! - constanta?
	ValuesMap::Iter it[2]	{ this[0].begin(),	this[1].begin() };
	ValuesMap::Iter itEnd[2]{ this[0].end(),	this[1].end()	};
	chrlen numb = 1;

	auto NextStart = [this, &it, &itEnd](BYTE s) { 
		auto it1 = next(it[s]);
		return it1 != itEnd[s] ? ValuesMap::Start(it1) : CHRLEN_MAX;
	};

	while (it[0] != itEnd[0] && it[1] != itEnd[1]) {
		if (!it[0]->second.MaxVal()) { it[0]++; continue; }
		if (!it[1]->second.MaxVal()) { it[1]++; continue; }

		it[0]->second.RgnNumb = it[1]->second.RgnNumb = numb;

		for(bool overlap = true; overlap; )
			if (NextStart(0) + overLen < End(it[1]))
				overlap = (++it[0])->second.RgnNumb = numb;
			else if (NextStart(1) + overLen < End(it[0]))
				overlap = (++it[1])->second.RgnNumb = numb;
			else 
				overlap = false;

		numb++;
		it[0]++, it[1]++;
	}
}

void ValuesMap::PrintStat(chrlen clen) const
{
	PrintRegionStats<ValuesMap>(this, clen);
}

//===== DataValuesMap

void DataValuesMap::BuildSpline(
	const DataSet<TreatedCover>& cover, const DataCoverRegions& rgns, bool redifineRgns, fraglen splineBase)
{
	BYTE strand = !Glob::IsPE;	// TOTAL for PE or POS for SE
	StrandData(POS).BuildSpline(cover.StrandData(POS), rgns.StrandData(eStrand(  strand)), redifineRgns, splineBase);
	StrandData(NEG).BuildSpline(cover.StrandData(NEG), rgns.StrandData(eStrand(2*strand)), redifineRgns, splineBase);
}

float DataValuesMap::GetPeakPosDiff() const
{
	auto& pData = StrandData(POS);
	auto& nData = StrandData(NEG);
	assert(pData.size() == nData.size());
	USHORT missed = 0;
	vector<SHORT> diffs;
	vector<chrlen> pPos, nPos;	// max positions in a positive, negative splines

	// get the difference of the splines maximums
	diffs.reserve(pData.size());	// suppose one summit for the region's spline
	pPos.reserve(2);
	nPos.reserve(2);
	for (auto itP = pData.begin(), itN = nData.begin(); itP != pData.end() && itN != nData.end(); itP++, itN++) {

		itP->second.GetMaxValPos(itP->first, pPos);
		itN->second.GetMaxValPos(itN->first, nPos);
		if (pPos.size() == nPos.size())
			for (auto itP = pPos.begin(), itN = nPos.begin(); itP != pPos.end(); itP++, itN++)
				diffs.push_back(SHORT(*itP - *itN));
		else
			missed++;		// just ignore splines with different max positions
		pPos.clear();
		nPos.clear();
	}
	if (missed)
		Verb::PrintMsgVar(Verb::DBG, "%4.1f%% regions were rejected while determining the fragment length\n", Percent(missed, pData.size()));

	int sum0 = 0;
	for (auto diff : diffs)
		sum0 += diff;

	return float(sum0) / diffs.size();
}

void DataValuesMap::Clear()
{
	StrandData(POS).clear();
	StrandData(NEG).clear();
}

#ifdef MY_DEBUG
void DataValuesMap::Print(chrid cID, chrlen stopNumb) const
{
	printf("\n");
	StrandData(POS).Print(cID, 0, stopNumb);
	StrandData(NEG).Print(cID, 1, stopNumb);
}
#endif

//===== BoundsValues
#ifdef MY_DEBUG
OSpecialWriter* BoundsValues::LineWriter = nullptr;
#endif

const BYTE L = 1;	// left bound, synonym for 'reverse' (is formed by reverse reads)
const BYTE R = 0;	// right bound, synonym for 'forward' (is formed by forward reads)

void BoundsValues::PushIncline(
	chrlen rgnNumb,
	BYTE reverse,
	const TracedPosVal& posVal,
	const TreatedCover& cover,
	vector<Incline>& inclines

) const
{
	const StrandOp& opDir = StrandOps[reverse];		// forward operations
	const StrandOp& opInv = StrandOps[!reverse];	// inversed operations
	auto itStart = opDir.GetPrev(cover.upper_bound(posVal.Pos(R)));	// min val: right for forward, left for reverse
	auto itStop = itStart;									// max val: left for forward, right for reverse

	// *** set itStop
	if (reverse)	for (itStop++; itStop->first < posVal.Pos(L); itStop++);
	else			for (itStop--; itStop->first > posVal.Pos(L); itStop--);

	// checking if there is an ordinary tag pool before the current start, masked by a spline
	if (!opDir.GetPrev(itStart)->second)
		opDir.Next(itStart);		// skip the gap after ordinary tag pool

	// *** it stop correction: calc the best BS position by successive linear regression approximations
	Incline incline;
	cover.LinearRegr(itStart, itStop, opDir, incline);
	if (!incline.Valid())	return;

	coviter it0 = itStop;
	Incline incline1{};
	auto compare = [&](coviter& it) {
		cover.LinearRegr(itStart, it, opDir, incline1);
		if (!incline1.Valid() || opDir.EqLess(incline.Pos, incline1.Pos))	return true;
		std::swap(incline, incline1);
		itStop = it;
		return false;
	};

	 //iterate itStop opposed to start
	for (auto it = it0;
		opDir.RetNext(it)->second > itStop->second || opDir.RetNext(it)->second > itStop->second;)
	{
		if (compare(it))	break;
	}

	// iterate itStop towards start
	bool nextStep = false;
	for (auto it = it0;
		opInv.RetNext(it)->second > itStop->second || (nextStep = opInv.RetNext(it)->second >= itStop->second);)
	{
		if (reverse && nextStep) { nextStep = false; it++; }
		if (compare(it))	break;
	}

	if (inclines.size()) {
		if (inclines.back() != incline)		// skip duplicate incline
			inclines.push_back(incline);
	}
	else
		inclines.push_back(incline);

#ifdef MY_DEBUG
	if (LineWriter && LineWriter->IsWriterSet()) {
		opInv.Next(itStop);
		LineWriter->WriteIncline(reverse, incline.Pos, itStop);
	}
#endif
}

void BoundsValues::CollectForwardInclines(const TreatedCover& rCover, vector<Incline>& inclines) const
{
	TracedPosVal posVal;
	for (auto it = rbegin(); it != rend(); it++) {		// loop through one derivative region
		posVal.Set(R, it);
		PushIncline(_rgnNumb, R, posVal, rCover, inclines);
		posVal.Retain();
	}
}

void BoundsValues::CollectReverseInclines(const TreatedCover& rCover, vector<Incline>& inclines) const
{
	TracedPosVal posVal;
	for (auto it = begin(); it != end(); it++) {		// loop through one derivative region
		posVal.Set(L, it);
		PushIncline(_rgnNumb, L, posVal, rCover, inclines);
		posVal.Retain();
	}
}

void BoundsValues::AddSignifValues(const tValuesMap::value_type& spline, chrlen relPos, Values& deriv)
{
#ifdef MY_DEBUG
	if (_maxVal < deriv.MaxVal())	_maxVal = deriv.MaxVal();
#endif
	_rgnNumb = spline.second.RgnNumb;
	emplace_back(
		spline.first + relPos,
		spline.second.Value(relPos),
		spline.second.Value(relPos + deriv.Length()),
		deriv
	);
}

void BoundsValues::SepSignifValues(
	const vector<Region>& rgns,
	BYTE rgnInd,
	const tValuesMap::value_type& spline,
	chrlen relPos,
	Values& deriv
)
{
	if (rgnInd == rgns.size())
		AddSignifValues(spline, relPos, deriv);
	else {
		auto rgn = rgns[rgnInd];
		USHORT rgnShift = rgnInd ? rgns[rgnInd - 1].End : 0;
		rgn -= rgnShift;
		Values sepDeriv(deriv, rgn);

		AddSignifValues(spline, relPos, deriv);
		SepSignifValues(rgns, ++rgnInd, spline, relPos + rgn.End, sepDeriv);
	}
}

void BoundsValues::AddValues(const tValuesMap::value_type& spline, chrlen relPos, Values& deriv)
{
	const float minDeriv = 0.052f;	// tangent 3 degrees
	chrlen	pos = 0;
	chrlen	start = 0;
	bool	nextPit = false;
	vector<Region> negligRgns;	// negligible derivator regions (inner only); typical capacity is 1

	// search for negligible regions
	for (auto& val : deriv) {
		if (val < minDeriv) {
			if (nextPit && !start)	// ignore starting negligible derivator region
				start = pos;
		}
		else
			if (nextPit) {
				if (start) {
					negligRgns.emplace_back(start, pos);
					start = 0;
				}
			}
			else
				nextPit = true;
		pos++;
	}

	// add derivators
	if (negligRgns.size())
		SepSignifValues(negligRgns, 0, spline, relPos, deriv);
	else
		AddSignifValues(spline, relPos, deriv);
}


//===== BoundsValuesMap

void BoundsValuesMap::BuildDerivs(int factor, const ValuesMap& splines)
{
	for (const auto& spline : splines) {
		if (!spline.second.MaxVal())	continue;

		Values deriv;			deriv.Reserve();
		BoundsValues derivSet;	derivSet.reserve(4);
		fraglen pos = 0;	// relative position
		const auto& itEnd = spline.second.end();

		// loop through spline
		for (auto it0 = spline.second.begin(), it = next(it0); it != itEnd; it0++, it++) {
			float tang = (*it - *it0) * factor;		// flip sign of pos strand delta
			if (tang > 0)							// cut off tangent of a falling spline
				deriv.AddValue(tang);
			else {
				pos++;
				if (deriv.MaxVal()) {
					fraglen len = deriv.Length();	// deriv will be cleared

					derivSet.AddValues(spline, pos, deriv);
					pos += len;
					deriv.Reserve();
				}
			}
		}

		// last region
		if (deriv.MaxVal())
			if (derivSet.size() && derivSet.back().End() == spline.first + pos)
				derivSet.back().AddValues(deriv);	// add adjacent region
			else
				derivSet.AddValues(spline, pos, deriv);
		AddRegions(spline.first, derivSet);
	}
}

#ifdef MY_DEBUG
void BoundsValuesMap::Print(eStrand strand, chrlen stopPos) const
{
	//printf("\nDERIVS %s: %2.2f\n", sStrandTITLES[strand - 1], MaxVal);
	printf("\nDERIVS %s\n", sStrandTITLES[strand]);
	for (const auto& rvss : *this) {
		if (stopPos && rvss.first > stopPos)
			break;
		printf("%d: %2.2f\n", rvss.first, rvss.second.MaxVal());
		for (const auto& rvs : rvss.second)
			printf("  %2.2f:\t%d %d\t%d\n", rvs.MaxVal(), rvs.RgnNumb, rvs.Start(), rvs.Length());
	}
}
#endif


//===== BS_map

void BS_map::AddPos(BYTE reverse, chrlen rgnNumb, const Incline& incl)
{
	// all positions are added in ascending order
	chrlen pos = incl.Pos + StrandOps[reverse].Factor * Glob::ReadLen;
	if (reverse) {		// insert left position; all right positions are already inserted
		iter lastIt = _lastIt;

		// find the iterator that will follow the inserted element
		for (; lastIt != end() && lastIt->first < pos; lastIt++);

		// the left inserted position can duplicate an already inserted right one; reduce it by 1
		if (lastIt != end() && lastIt->first == pos) {
			printf(">> pos %d, numb %d, dupl last\n", pos, rgnNumb);
			pos--;
		}
		else if (lastIt != begin() && prev(lastIt)->first == pos) {
			printf(">> pos %d, numb %d, dupl prev\n", pos, rgnNumb);
			pos--;
			lastIt--;
		}

		_lastIt = emplace_hint(lastIt, pos, BS_PosVal(1, rgnNumb));
		_lastIt->second.RefPos = _lastIt->first;
	}
	else {
		auto it = emplace_hint(end(), pos, BS_PosVal(0, rgnNumb));
		it->second.RefPos = it->first;
	}
}

void BS_map::AddBounds(BYTE reverse, chrlen rgnNumb, vector<Incline>& inclines)
{
	if (!inclines.size())	return;

	// sort inclines by ascending (for forward) / descending (for reverde) reference positions
	sort(inclines.begin(), inclines.end(),
		[&reverse](const Incline& i1, const Incline& i2) {
			return StrandOps[reverse].Less(i1.Pos, i2.Pos);
		}
	);

	// find permissible derivative limit
	float limDeriv = 0;
	for (const auto& i : inclines)
		limDeriv += i.Deriv;
	limDeriv /= inclines.size();
	limDeriv *= 0.15F;

	// *** Intersection Filter & Insignificant Incline Filter:
	// insert only the best bounds (formed by the steepest inclines), and with a significant derivative
	auto it = inclines.cbegin();
	chrlen addedPos = it->TopPos;

	AddPos(reverse, rgnNumb, *it);		// definitely add the first, steepest (tightest) incline
	for (it++; it != inclines.cend(); it++)
		if (StrandOps[!reverse].EqLess(it->TopPos, addedPos)
		&& it->Deriv > limDeriv)	
		{
			AddPos(reverse, rgnNumb, *it);
			addedPos = it->TopPos;
		}
}

void BS_map::SetBounds(BYTE reverse, const BoundsValuesMap& derivs, const TreatedCover& rCover)
{
	using tSetBSpos = void(BoundsValues::*)(const TreatedCover&, vector<Incline>& inclines) const;
	tSetBSpos fcollectInclines = reverse ?
		&BoundsValues::CollectReverseInclines :
		&BoundsValues::CollectForwardInclines;
	vector<Incline> inclines;
	inclines.reserve(4);

	for (const auto& d : derivs) {		// loop through the derivative regions
		inclines.clear();
		(d.second.*fcollectInclines)(rCover, inclines);
		AddBounds(reverse, d.second.RgnNumb(), inclines);
	}
}

const short MIN_BS_WIDTH = 5;

// Fits the BS width to a minimum by adjusting the BS reference positions
//	@param start: BS iterator pointed to base left bound 
//	@param end: BS iterator pointed to base right bound
//	@param len: basic BS width
void FitToMinWidth(BS_map::iter& start, BS_map::iter& end, short len)
{
	const auto diff = MIN_BS_WIDTH - len;
	auto Expand = [](const BS_map::iter& it, chrlen newPos) {
		if (it->first != newPos)	it->second.RefPos = newPos;
	};

	Expand(start, start->first - diff / 2);
	Expand(end, end->first + diff / 2 + diff % 2);
};

void BS_map::ExtendNarrowWidths()
{
	DoExtend([&](vector<PosValue>* VP) {
		// basic start, end
		auto& start = VP[L].back().Iter;
		auto& end	= VP[R].front().Iter;
		if (start->second.RgnNumb != end->second.RgnNumb)	return;

		const auto len = short(end->first - start->first);
		if (len < MIN_BS_WIDTH) {
			FitToMinWidth(start, end, len);

			// reset boundaries that are covered by extented base boundaries
			auto resetCoveredItems = [](BS_PosVal& pval, BS_map::iter& itBase) {
				if (pval.RefPos > itBase->second.RefPos)	return;
				pval.Score = 0;
			};

			// left BS extentions
			auto& vp = VP[L];
			const auto extLen = BYTE(vp.size() - 1);
			for (BYTE i = 0; i < extLen; i++)			// left to right
				resetCoveredItems(vp[i].Iter->second, start);

			// right BS extentions
			vp = VP[R];
			for (BYTE i = BYTE(vp.size() - 1); i; i--)	// right to left
				resetCoveredItems(vp[i].Iter->second, end);
		}
		});
}

void BS_map::PrintWidthDistrib() const
{
	map<fraglen, vector<chrlen>> freq;
	chrlen	totalLen = 0, bsNumb = 0;

	// collect numbers
	DoBasic([&](citer& start, citer& end) {
		auto len = fraglen(end->second.RefPos - start->second.RefPos);
		totalLen += len;
		freq[len].push_back(++bsNumb);
		}
	);
	// print numbers
	const char* length = "width";
	printf("\nBS WIDTH FREQUENCY:\n");
	printf("%s cnt  numbers\n", length);
	for (const auto& item : freq) {
		printf("%4d %4u  ", item.first, UINT(item.second.size()));
		auto it = item.second.begin();
		printf("%d", *it);
		for (it++; it != item.second.end(); it++)
			printf(",%d", *it);
		printf("\n");
	}
	printf("average %s: %.2f\n", length, float(totalLen) / bsNumb);
}

void BS_map::Refine()
{
	/*
	BS Left entry (bound): is formed by reverse reads; BS Right entry (bound): is formed by forward reads

	Options for placing bounds in a region:
	canonical:					[[L] [R]]
	extra right/left bounds:	[R] [[L] [R]] [L]
	'negative' BS width:		[R] [L]

	Method brings the instance to canonical form, resetting the score of all extra elements to zero,
	and marking BSs with 'negative' width.
	*/
	auto lastExtRight_it = end();	// iterator pointing to the last Right entry that starts the region 
	uint16_t extRightCnt = 0;		// count of Right entries (that starts the region)
	uint16_t extLeftCnt = 0;		// count of Left entries (that ends the region)
	bool newBS = true;				// if true then new BS in the region is registered
	bool someBS = false;			// if true then at least one BS in the region is registered
	chrlen rgnNumb = 1;

	// resets left and right extra entries
	auto ResetAllExtEntries = [&](iter it) {
		// resets left or right extra ('false') entries
		//	@param it: iterator pointing to the first extra entry
		//	@param entryCnt: count of extra entries which should be reset; becomes zero
		//	@returns: iterator pointing to the first non-extra entry
		auto ResetExtEntries = [&](iter& it, uint16_t& entryCnt) {
			for (; entryCnt && it != end(); entryCnt--, --it)
				it->second.Score = 0;
			return it;
		};

		// 1) reset extra right entries
		if (lastExtRight_it != end()) {
			if (someBS)
				ResetExtEntries(lastExtRight_it, extRightCnt);
			else {						// ** 'negative' BS width
				lastExtRight_it->second.Reverse = true;		// change 'right' bound to 'left'
				if (extRightCnt > 0)
					ResetExtEntries(--lastExtRight_it, --extRightCnt);	// reset other extra entries
			}
			lastExtRight_it = end();
		}
		// 2) reset extra left entries
		if (someBS)
			ResetExtEntries(--it, extLeftCnt);
		else if(extLeftCnt > 0) {		// ** 'negative' BS width
			auto itR = ResetExtEntries(--it, --extLeftCnt);	// reset other extra entries
			if (itR != end()) {
				itR->second.Reverse = false;				// change 'left' bound to 'right

				// decrease the width of the updated BS if needed
				auto itL = prev(itR);
				auto len = short(itR->first - itL->first);
				if (len > MIN_BS_WIDTH)
					FitToMinWidth(itL, itR, len);
			}
		}
	};

	// *** refine
	for (auto it = begin(); it != end(); it++)
	{
		const bool newRgn = rgnNumb != it->second.RgnNumb;
		if (newRgn) {
			// 'close' previous region
			ResetAllExtEntries(it);
			// reset current region
			newBS = rgnNumb = it->second.RgnNumb;
			someBS = extLeftCnt = extRightCnt = 0;
		}

		if (it->second.Reverse) {	// Left bound
			extLeftCnt++;
			newBS = false;
		}
		else {						// Right bound
			if (newBS && !extLeftCnt)
				if (!extRightCnt || !newRgn)
					lastExtRight_it = it;
			if (extLeftCnt)
				someBS = true;
			else
				extRightCnt++;
			extLeftCnt = 0;
		}
	}
	// 'close' last region
	ResetAllExtEntries(end());

	// *** extend narrows
	ExtendNarrowWidths();
}

void BS_map::SetScore(const DataSet<TreatedCover>& fragCovers)
{
	float maxScore = 0;
	auto& cover = fragCovers.TotalData();
	SSpliner<coval> spliner(eCurveType::ROUGH, 5);
	Values	vals;			// spline values
	vals.Reserve();

	DoExtend([&](vector<PosValue>* VP) {
		// extended start, end
		auto& start = VP[L].front().Iter;
		auto& end	= VP[R].back().Iter;
		if (start->second.RgnNumb != end->second.RgnNumb)	return;

		chrlen	startPos = start->second.RefPos;
		cover.SetLocalSpline(spliner, startPos, end->second.RefPos, vals);

		// *** set score

		// set score for the BS extension (excluding base boundaries)
		auto setExtScore = [&vals](BYTE reverse, const vector<PosValue>& vp, USHORT offset, BYTE ind) {
			/*
			fill the score from left to right (for reverse)
			or from rigth to left (for direct) extended boundaries
			*/
			float score = 0;
			auto len = vp[ind + reverse].Iter->second.RefPos - vp[ind - !reverse].Iter->second.RefPos;

			assert(len + offset < vals.size());
			for (chrlen i = 0; i < len; i++, offset++)
				score += vals[offset];
			vp[ind].Iter->second.Score = score / len;
		};

		USHORT offset;
		{	// left BS extensions
			const auto& vp = VP[L];
			const auto extLen = BYTE(vp.size() - 1);
			offset = USHORT(vp.front().Iter->second.RefPos - startPos);
			for (BYTE i = 0; i < extLen; i++)			// left to right
				setExtScore(L, vp, offset, i);
		}
		{	// right BS extensions
			const auto& vp = VP[R];
			const auto extLen = BYTE(vp.size() - 1);
			if (extLen) {
				offset = USHORT(vp[extLen - 1].Iter->second.RefPos - startPos);	// from the penultimate boundary
				for (BYTE i = extLen; i; i--)	// right to left
					setExtScore(R, vp, offset, i);
			}
		}
		// base BS
		const auto& itStartBS = VP[L].back().Iter;
		const auto& itEndBS = VP[R].front().Iter;
		offset = itStartBS->second.RefPos - startPos;
		const auto len = USHORT(itEndBS->second.RefPos - itStartBS->second.RefPos);
		float score = 0;

		for (USHORT i = 0; i < len; i++, offset++)
			score += vals[offset];
		score /= len;
		if (maxScore < score)	maxScore = score;
		itStartBS->second.Score = itEndBS->second.Score = score;

		vals.Clear();
		});

	// *** normalize score
	for (auto& x : *this)
		if (x.second.Score)
			x.second.Score /= maxScore;
}

void BS_map::PrintStat() const
{
	if (Verb::StrictLevel(Verb::CRIT))	return;

	const bool stat = Verb::Level(Verb::DBG);	// collect and print statistics
	chrlen	bsNumb = 0;
	//chrlen	minNegNumb, maxNegNumb, minPosNumb, maxPosNumb = 0;
	//float	minNegRatio, minPosRatio, maxNegRatio = 0, maxPosRatio = 0;
	chrlen	minScoreNumb;
	float	minScore = 1.f;
	//minNegRatio = minPosRatio = 1000;

	// define min/max length, minScore, min***Ratio/max***Ratio
	DoBasic([&](citer& start, citer& end) {
		++bsNumb;
		if (stat) {
			const fraglen len = end->second.RefPos - start->second.RefPos;
			float score = start->second.Score;
			if (minScore > score)		minScore = score, minScoreNumb = bsNumb;

			//if (score < 1) {
			//	if (minNegRatio > score)		minNegRatio = score, minNegNumb = bsNumb;
			//	else if (maxNegRatio < score)	maxNegRatio = score, maxNegNumb = bsNumb;
			//}
			//else
			//	if (minPosRatio > score)		minPosRatio = score, minPosNumb = bsNumb;
			//	else if (maxPosRatio < score)	maxPosRatio = score, maxPosNumb = bsNumb;
		}
		});
	if (stat) {
		//auto ptTableTitle = [](const char* title) {
		//	auto len = USHORT(strlen(title) + 1);
		//	PrintSolidLine(len);
		//	printf("%s\n", title);
		//	PrintSolidLine(len);
		//};

		PrintWidthDistrib();
		// score
		printf("\nmin score: %2.2f (%d)\n", minScore, minScoreNumb);
		//ptTableTitle("RATIO:\tmin  (cnt)   max  (cnt)");
		//printf("reverse\t%2.2f (%3d)   %2.2f (%3d)\n", 1 / maxNegRatio, minNegNumb, 1 / minNegRatio, maxNegNumb);
		//printf("forward\t%2.2f (%3d)   %2.2f (%3d)\n", minPosRatio, minPosNumb, maxPosRatio, maxPosNumb);

		printf("\n");
	}

	printf("BS count: %d\n\n", bsNumb);
}

#ifdef MY_DEBUG
void BS_map::CheckScoreHierarchy()
{
	if (!Verb::Level(Verb::DBG))	return;

	chrlen bsNumb = 0;
	bool issues = false;

	printf("\nCHECK SCORES HIERARCHY\n");
	DoExtend([&](const vector<PosValue>* VP) {
		bool prNeg = false, prPos = false;
		float val = 0;
		auto prVals = [](bool sep, char sign, const vector<PosValue>& VP) {
			if (sep) 	printf(" ");
			printf("%c", sign);
			for (const PosValue& vp : VP)
				printf("%2.3f ", vp.Val);
		};

		bsNumb++;
		// set prNeg; check all
		for (const PosValue& vp : VP[L])
			if (prNeg = (val > vp.Val))		break;
			else val = vp.Val;
		// set prPos; check all except first element
		const BYTE vpLen = BYTE(VP[R].size() - 1);
		const auto& vp = VP[R];
		val = 1000;
		for (BYTE i = 1; i < vpLen; i++)
			if (prPos = (val < vp[i].Val))	break;
			else val = vp[i].Val;
		// print of excess over basic score
		if (prNeg || prPos) {
			issues = printf("%4d  %d\t", bsNumb, VP[L].back().Iter->second.RefPos);	// start position
			if (prNeg)	prVals(false, '-', VP[L]);
			if (prPos)	prVals(prNeg, '+', VP[R]);
			printf("\n");
		}
		});

	if (!issues)	printf("OK\n");
}

void BS_map::Print(chrid cID, const string& outFName, bool selected, chrlen stopPos) const
{
	string format = "%8d % 4d  %c %8d %5.2f   %s\n";
	const char bound[]{ 'R','L' };
	IGVlocus locus(cID);

	TxtOutFile file(outFName.c_str());
	file.Write(" pos     rgn bnd  ref pos score   IGV view\n");
	for (const auto& x : *this) {
		if (stopPos && x.first > stopPos)	break;
		if (selected && !x.second.Score)	continue;
		format[5] = x.second.RgnNumb % 2 ? '-' : SPACE;	// odd numbers are aligned to the left, even numbers to the right
		format[20] = x.second.Score ? '2' : '0';		// zero score without fraction
		file.Write(format.c_str(),
			x.first,
			x.second.RgnNumb,
			bound[x.second.Reverse],
			x.first != x.second.RefPos ? x.second.RefPos : 0,
			x.second.Score,
			locus.Print(x.first)
		);
	}
}
#endif // MY_DEBUG

//===== BedWriter

bool BedWriter::rankScore = false;
BedWriter::tAddScore BedWriter::fLineAddScore = nullptr;

void BedWriter::WriteChromData(chrid cID, const CoverRegions& rgns)
{
	const reclen colorLen = reclen(strlen(sGRAY));
	const reclen offset = AddChromToLine(cID);

	for (const auto& rgn : rgns) {
		LineAddUInts(rgn.Start(), rgn.End(), rgn.value, false);
		if (!rgn.value) {		// discarded items
			LineAddChars("\t.\t.\t", 5, false);
			LineAddInts(rgn.Start(), rgn.End(), true);
			LineAddChars(sGRAY, colorLen, false);
		}
		LineToIOBuff(offset);
	}
}

void BedWriter::WriteChromData(chrid cID, BS_map& bss)
{
	const reclen offset = AddChromToLine(cID);
	bool lastSep[]{ false, false };
	chrlen bsNumb = 0;

	bss.DoExtend([&](const vector<BS_map::PosValue>* VP) {
		// *** save basic info
		auto& start = VP[L].back().Iter;
		auto& end	= VP[R].front().Iter;

		LineAddUInts(start->second.RefPos, end->second.RefPos, ++bsNumb, true);	// 3 basic fields
		LineAddScore(start->second.Score, true);					// BS score

		// *** extended boudaries info
		LineAddChar(DOT, true);
		lastSep[1] = VP[R].size() - 1;
		LineAddFloat(end->second.Score, VP[1].size() - 1 || lastSep[1]);	// ratio

		// *** extra deviations info
		for (BYTE s : {L, R}) {
			const BYTE vpLen = BYTE(VP[s].size() - 1);
			if (!vpLen) continue;

			auto& vp = VP[s];
			function<void(BYTE, char)> saveExtraPos = [this, &vp](BYTE i, char specChar) {
				LineAddArgs("%c%d", specChar, vp[i].Iter->first);
			};
			function<void(BYTE, char)> saveExtraVal = [this, &vp](BYTE i, char specChar) {
				if (specChar)	LineAddChar(specChar);
				LineAddScore(vp[i].Iter->second.Score, false);
			};
			auto saveExtraFields = [&](function<void(BYTE, char)>& fn, char specChar, bool delim) {
				BYTE i = !s;
				fn(i, specChar);
				for (i++; i < vpLen; i++)	fn(i, COMMA);
				if (delim)	LineAddChar(TAB);
			};

			saveExtraFields(saveExtraPos, Read::Strands[s], true);
			saveExtraFields(saveExtraVal, 0, lastSep[s]);
		}

		LineToIOBuff(offset);
		});
}

void BedWriter::WriteChromExtData(chrid cID, BS_map& bss)
{
	chrlen bsNumb = 0;
	const reclen offset = AddChromToLine(cID);

	bss.DoExtend([&](const vector<BS_map::PosValue>* VP) {
		const BYTE COLORS_CNT = 4;
		static const string colors[]{
			// color		 ind	feature_score/BS_score
			"155,233,168",	// 1	>=0.2	light light green
			"64,196,99",	// 2	>=0.4	light green
			"48,161,78",	// 3	>=0.6	green
			"33,110,57",	// 4	>=0.8	dark green

			//"0,190,255",	// 1	>=0.2	light blue
			//"0,160,230",	// 2	>=0.4	blue
			//"0,130,205",	// 3	>=0.6	dark blue
			//"0,100,180",	// 4	>=0.8	dark dark blue

			//"140,30,30",	// 0	>=0		dark red
			//"140,85,30",	// 1	>=0.2	dark orange
			//"180,140,30",	// 2	>=0.4	dark yellow
			//"130,140,30",	// 3	>=0.6	dark yellow-green
			//"30,100,30",	// 4	>=0.8	dark green
		};
		const auto& start = VP[L].back().Iter;	// basic feature start
		const float score = VP[L].back().Val;	// basic feature score

		auto addExtraLines = [=, &bsNumb](const vector<BS_map::PosValue>& vp) {	// &bsNumb is essential, otherwise bsNumb goes out of sync
			if (vp.size() == 1)	return;
			for (auto it0 = vp.begin(), it = next(it0); it != vp.end(); it0++, it++) {
				LineAddUInts(it0->Iter->second.RefPos, it->Iter->second.RefPos, bsNumb, true);
				LineAddFloat(it0->Val, true);	// BS score
				LineAddChar(DOT, true);
				LineAddInts(it0->Iter->second.RefPos, it->Iter->second.RefPos, true);
				// colors
				auto ind = BYTE(10 * it0->Val / score) / 2;
				if (ind > COLORS_CNT - 1)	ind = COLORS_CNT - 1;
				LineAddStr(colors[ind], false);
				LineToIOBuff(offset);
			}
		};

		const auto& end = VP[R].front().Iter;
		const char* delims = ".\t.\t.\t.\t";

		++bsNumb;

		addExtraLines(VP[L]);
		// *** add basic feature
		LineAddUInts(start->second.RefPos, end->second.RefPos, bsNumb, true);
		LineAddFloat(start->second.Score, true);		// BS score
		LineAddChars(delims, reclen(strlen(delims)), false);
		LineAddFloat(end->second.Score, false);			// reverse/forward ratio
		LineToIOBuff(offset);

		addExtraLines(VP[R]);
		});
}

void BedWriter::WriteChromROI(chrid cID, const BS_map& bss)
{
	const reclen offset = AddChromToLine(cID);
	chrlen bsNumb = 0;

	bss.DoBasic([&](BS_map::citer& start, BS_map::citer& end) {
		LineAddUInts(
			start->second.RefPos - Glob::ROI_ext,
			end->second.RefPos + Glob::ROI_ext,
			++bsNumb,
			false
		);
		LineToIOBuff(offset);
		});
}


//===== FixWigWriter

void FixWigWriter::WriteChromData(chrid cID, const ValuesMap& vals)
{
	for (const auto& val : vals)
		if (val.second.MaxVal())
			WriteFixStepRange(cID, val.first, val.second);
}


//===== FixWigWriterSet

void FixWigWriterSet::WriteChromData(chrid cID, const BoundsValuesMap& set)
{
	for (const auto& rvss : set)
		for (const auto& rvs : rvss.second)
			if (rvs.MaxVal())
				WriteFixStepRange(cID, rvs.Start(), rvs);
}
