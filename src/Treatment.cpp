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

void TreatedCover::SetLocalSpline(SSpliner<coval>& spliner, chrlen startPos, chrlen endPos, Values& vals) const
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
	//const auto minLen = fraglen(1.6f * Glob::FragLen);
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
void CoverRegions::PrintScoreDistrib(const string& fname) const
{
	map<coval, chrlen> freq;

	for (const auto& rgn : *this)
		freq[rgn.value]++;
	TxtOutFile file(fname.c_str());

	//printf("\nREGIONS LENGTH FREQUENCY\n");
	file.Write("score\tfreq\n");
	for (const auto& item : freq)
		file.Write("%d\t%d\n", item.first, item.second);
}
#endif

//===== DataCoverRegions

bool DataCoverRegions::SetPotentialRegions(const DataSet<TreatedCover>& cover, chrlen cLen, coval cutoff, bool noMultiOverl)
{
	chrlen capacity = cLen / (Glob::FragLen * 100);
	if (Glob::IsPE) {
		TotalData().SetPotentialRegions(cover.TotalData(), capacity, cutoff * 2);
		//TotalData().PrintScoreDistrib("RGNS.dist");
	}
	else {
		StrandData(FWD).SetPotentialRegions(cover.StrandData(FWD), capacity, cutoff);
		StrandData(RVS).SetPotentialRegions(cover.StrandData(RVS), capacity, cutoff);

		//StrandData(FWD).PrintScoreDistrib("RGNS.pos.dist");
		//StrandData(RVS).PrintScoreDistrib("RGNS.neg.dist");

		EliminateNonOverlapsRegions<CoverRegions>(Data(), Glob::FragLen);
		if (noMultiOverl)
			EliminateMultiOverlapsRegions<CoverRegions>(Data());
		if (Verb::Level(Verb::DBG))
			PrintRegionStats<CoverRegions>(Data(), cLen);
	}
	if (!Empty())	return false;
	Verb::PrintMsg(Verb::CRIT, "No enriched regions found");
	return true;
}

#ifdef MY_DEBUG
void DataCoverRegions::PrintScoreDistrib(const string& fname) const
{
	const string ext = ".dist";

	if (Glob::IsPE)
		TotalData().PrintScoreDistrib(fname + ext);
	else {
		StrandData(FWD).PrintScoreDistrib(fname + sStrandEXT[FWD] + ext);
		StrandData(RVS).PrintScoreDistrib(fname + sStrandEXT[RVS] + ext);
	}
}
#endif

//fraglen DataCoverRegions::GetFragMean(const DataSet<TreatedCover>& cover) const
//{
//	auto& pData = StrandData(FWD);
//	auto& nData = StrandData(RVS);
//	auto& pCover = cover.StrandData(FWD);
//	auto& nCover = cover.StrandData(RVS);
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

float Values::AvrScoreInRange(int32_t& offset, int32_t len, int8_t factor) const
{
#ifdef MY_DEBUG
	if (len <= 0) {
		printf(">> %d offset: %d len: %d\n", factor, offset, len);
		return -1;
	}
#endif
	float score = 0;
	for (int32_t i = 0; i < len; i++, offset += factor)
		score += (*this)[offset];
	return score / len;
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
	printf(" N  start\tend\tval\tIGV view\n");
	for (const auto& x : *this) {
		if (stopNumb && x.second.GrpNumb > stopNumb)	break;
		if (x.second.MaxVal()) {
			chrlen end = x.first + x.second.Length();
			printf("%3d %d\t%d\t%2.2f\t%s\n",
				x.second.GrpNumb, x.first, end, x.second.MaxVal(), locus.Print(x.first, end));
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

void ValuesMap::NumberGroups()
{
	// minimum overlap length: 10 - empirical addition
	const fraglen minOverlap = SSpliner<coval>::SilentLength(CurveTYPE, ReadSplineBASE) + 10;
	ValuesMap::Iter it[2]	{ this[0].begin(),	this[1].begin() };
	ValuesMap::Iter itEnd[2]{ this[0].end(),	this[1].end()	};
	chrlen numb = 1;

	auto isOverlap = [&](BYTE s) {
		auto it1 = next(it[s]);
		auto nextStart = it1 != itEnd[s] ? ValuesMap::Start(it1) : CHRLEN_MAX;
		if (nextStart >= End(it[!s])) 
			return false;		// non-overlapping
		(++it[s])->second.GrpNumb = numb;
		if (nextStart + minOverlap >= End(it[!s]))
			it[s]->second.MarkAsEmpty();	// overlapping is insufficient
		return true;			// overlapping
	};

	while (it[0] != itEnd[0] && it[1] != itEnd[1]) {
		if (!it[0]->second.MaxVal()) { it[0]++; continue; }
		if (!it[1]->second.MaxVal()) { it[1]++; continue; }

		it[0]->second.GrpNumb = it[1]->second.GrpNumb = numb;
		for (bool overlap = true; overlap; )
			if (!(overlap = isOverlap(0)))
				overlap = isOverlap(1);
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
	BYTE strand = !Glob::IsPE;	// TOTAL for PE or FWD for SE
	StrandData(FWD).BuildSpline(cover.StrandData(FWD), rgns.StrandData(eStrand(  strand)), redifineRgns, splineBase);
	StrandData(RVS).BuildSpline(cover.StrandData(RVS), rgns.StrandData(eStrand(2*strand)), redifineRgns, splineBase);
}

float DataValuesMap::GetPeakPosDiff() const
{
	auto& pData = StrandData(FWD);
	auto& nData = StrandData(RVS);
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
		Verb::PrintMsgVar(Verb::DBG, "%4.1f%% rejected regions;\t", Percent(missed, pData.size()));

	int sum0 = 0;
	for (auto diff : diffs)
		sum0 += diff;

	return float(sum0) / diffs.size();
}

void DataValuesMap::Clear()
{
	StrandData(FWD).clear();
	StrandData(RVS).clear();
}

#ifdef MY_DEBUG
void DataValuesMap::Print(chrid cID, chrlen stopNumb) const
{
	printf("\n");
	StrandData(FWD).Print(cID, 0, stopNumb);
	StrandData(RVS).Print(cID, 1, stopNumb);
}
#endif

//===== BoundsValues
#ifdef MY_DEBUG
OSpecialWriter* BoundsValues::LineWriter = nullptr;
#endif

void BoundsValues::PushIncline(
	chrlen grpNumb,
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
		PushIncline(_grpNumb, R, posVal, rCover, inclines);
		posVal.Retain();
	}
}

void BoundsValues::CollectReverseInclines(const TreatedCover& rCover, vector<Incline>& inclines) const
{
	TracedPosVal posVal;
	for (auto it = begin(); it != end(); it++) {		// loop through one derivative region
		posVal.Set(L, it);
		PushIncline(_grpNumb, L, posVal, rCover, inclines);
		posVal.Retain();
	}
}

void BoundsValues::AddSignifValues(const tValuesMap::value_type& spline, chrlen relPos, Values& deriv)
{
#ifdef MY_DEBUG
	if (_maxVal < deriv.MaxVal())	_maxVal = deriv.MaxVal();
#endif
	_grpNumb = spline.second.GrpNumb;
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
			printf("  %2.2f:\t%d %d\t%d\n", rvs.MaxVal(), rvs.GrpNumb, rvs.Start(), rvs.Length());
	}
}
#endif


//===== BS_map

#define	POS(it)	(it)->second.RefPos
#define LEN(itStart,itEnd)	short(POS(itEnd) - POS(itStart))
#define	SCORE(it)	(it)->second.Score

void BS_map::AddPos(BYTE reverse, chrlen grpNumb, const Incline& incl)
{
	// all positions are added in ascending order
	chrlen pos = incl.Pos + StrandOps[reverse].Factor * Glob::ReadLen;
	if (reverse) {		// insert left position; all right positions are already inserted
		iter lastIt = _lastIt;

		// find the iterator that will follow the inserted element
		for (; lastIt != end() && lastIt->first < pos; lastIt++);

		// the left inserted position can duplicate an already inserted right one; reduce it by 1
		if (lastIt != end() && lastIt->first == pos) {
			//printf(">> pos %d, numb %d, dupl last\n", pos, grpNumb);
			pos--;
		}
		else if (lastIt != begin() && prev(lastIt)->first == pos) {
			//printf(">> pos %d, numb %d, dupl prev\n", pos, grpNumb);
			pos--;
			lastIt--;
		}

		_lastIt = emplace_hint(lastIt, pos, BS_PosVal(1, grpNumb));
		POS(_lastIt) = _lastIt->first;
	}
	else {
		auto it = emplace_hint(end(), pos, BS_PosVal(0, grpNumb));
		POS(it) = it->first;
	}
}

void BS_map::AddBounds(BYTE reverse, chrlen grpNumb, vector<Incline>& inclines)
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

	AddPos(reverse, grpNumb, *it);		// definitely add the first, steepest (tightest) incline
	for (it++; it != inclines.cend(); it++)
		if (StrandOps[!reverse].EqLess(it->TopPos, addedPos)
		&& it->Deriv > limDeriv)	
		{
			AddPos(reverse, grpNumb, *it);
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
		AddBounds(reverse, d.second.GrpNumb(), inclines);
	}
}

const short MIN_BS_WIDTH = 5;

// Fits the BS width to a minimum by adjusting the BS reference positions
//	@param start: BS iterator pointed to base left bound
//	@param end: BS iterator pointed to base right bound
//	@param isLess: if true then the width should be less than the minimum
//	@returns: true if the condition is met and the width is adjusted
bool FitToMinWidth(BS_map::iter& start, BS_map::iter& end, bool isLess)
{
	const auto diff = MIN_BS_WIDTH - LEN(start,end);
	if ((diff <= 0) == isLess)	return false;

	auto Expand = [](const BS_map::iter& it, chrlen newPos) {
		if (it->first != newPos)	POS(it) = newPos;
	};

	Expand(start, start->first - diff / 2);
	Expand(end, end->first + diff / 2 + diff % 2);
	return true;
};

void BS_map::ExtendNarrowBS(iter& itL, iter& itR)
{
	if (FitToMinWidth(itL, itR, true)) {
		//check left adjacent bounds
		for (auto it = prev(itL); it != end(); it--)
			if (IsValid(it))
				if (POS(it) < POS(itL))	break;
				else					SetInvalid(it);
		// check right adjacent bounds
		for (auto it = next(itR); it != end(); it++)
			if (IsValid(it))
				if (POS(it) > POS(itR))	break;
				else					SetInvalid(it);
	}
}

void BS_map::ExtendSingleNarrowBS(iter& start, iter& end)
{
	iter itEnd = next(end) == this->end() ? this->end() : next(end);

	for (auto it = start; it != itEnd; it++)
		if (IsValid(it) && !it->second.Reverse) {
			auto itL = prev(it);
			ExtendNarrowBS(itL, it);
			break;
		}
}

void BS_map::ExtendNarrowBSsInGroup(iter& start, iter& stop, bool narrowBS, bool closeProx)
{
	bool lastLeft = true;
	iter itEnd = next(stop) == end() ? end() : next(stop);
	vector<pair<iter, iter>> bss;	// BS start-end collection

	bss.reserve(4);
	// collect BS
	for (auto it = start; it != itEnd; it++)
		if (IsValid(it))
			if (it->second.Reverse)
				lastLeft = true;
			else if (lastLeft) {
				bss.emplace_back(prev(it), it);
				lastLeft = false;
			}

	// 'merge' close proximities
	if (closeProx) {
		iter	lastR = end();		// last BS right bound
		vector<pair<iter, iter>> newBss;

		newBss.reserve(bss.size() - 1);
		for (const auto& bs : bss) {
			if (lastR != end())
				if (LEN(lastR, bs.first) <= MIN_BS_WIDTH) {
					for (auto it = lastR; it != bs.second; it++)
						SetInvalid(it);
					newBss.emplace_back(prev(lastR), bs.second);	// save previous & current
				}
				else
					newBss.emplace_back(prev(lastR), lastR);		// save previous
			lastR = bs.second;
		}
		if (POS(newBss.back().first) == POS(bss.back().first))		// the same lasts
			newBss.push_back(newBss.back());						// save last

		bss.swap(newBss);
	}
	// extend narrow BS if there are any left
	if (narrowBS)
		for (auto& bs : bss)
			ExtendNarrowBS(bs.first, bs.second);
}

void BS_map::ExtendNarrowBSs()
{
	/*
	* here the groups are 'canonical', e.g. 
	* L L L R R L L R L R R
	*     |_|     |_| |_|
	*/
	bool	lastLeft = true;	// true if there was BS left bound
	bool	narrowBS = false;	// true if there's at least one narrow BS
	bool	closeProx = false;	// true if there's at least one close proximity, i.e. too close BSs
	iter	itStart = end();	// group start iterator
	iter	itEnd = end();		// group end iterator
	chrlen	grpNumb = 1;
	chrlen	lastRpos = 0;		// last BS right bound position
	BYTE	bsCnt = 0;			// count of BSs

	// draft common bypass
	for (auto it = begin(); it != end(); it++) {
		if (!IsValid(it))	continue;

		if (grpNumb != it->second.GrpNumb) {
			if (narrowBS || closeProx) {
				// rigorous group bypass
				if (bsCnt == 1)
					ExtendSingleNarrowBS(itStart, itEnd);
				else
					ExtendNarrowBSsInGroup(itStart, itEnd, narrowBS, closeProx);
				narrowBS = closeProx = false;
			}
			lastRpos = bsCnt = 0;
			itStart = end();
			grpNumb = it->second.GrpNumb;
		}

		if (it->second.Reverse)
			lastLeft = true;
		else
			if (lastLeft) {		// BS right boud
				if (LEN(prev(it), it) < MIN_BS_WIDTH)
					narrowBS = true;
				if (POS(prev(it)) - lastRpos < MIN_BS_WIDTH)
					closeProx = true;
				bsCnt++;
				lastLeft = false;
				lastRpos = POS(it);
			}

		if (itStart == end())
			itStart = it;
		itEnd = it;
	}
}

void BS_map::Refine()
{
	/*
	* BS Left entry (bound): is formed by reverse reads; BS Right entry (bound): is formed by forward reads
	*
	* Options for placing bounds in a region:
	* canonical:					[[L] [R]]
	* adjacent right/left bounds:	[R] [[L] [R]] [L]
	* 'negative' BS width:			[R] [L]
	*
	* Method brings the instance to canonical form, resetting the score of all adjacent elements to zero,
	* and marking BSs with 'negative' width.
	*/
	auto lastExtRight_it = end();	// iterator pointing to the last Right entry that starts the region 
	uint16_t extRightCnt = 0;		// count of adjacent Right entries only (that starts the region)
	uint16_t extLeftCnt = 0;		// count of adjacent Left entries only (that ends the region)
	bool newBS = true;				// if true then new BS in the region is registered
	bool someBS = false;			// if true then at least one BS in the region is registered
	chrlen grpNumb = 1;

	// resets left and right adjacent entries
	auto ResetAllExtEntries = [&](iter it) {
		// resets left or right adjacent ('false') entries
		//	@param it: iterator pointing to the first adjacent entry
		//	@param entryCnt: count of adjacent entries which should be reset; becomes zero
		//	@returns: iterator pointing to the first non-adjacent entry
		auto ResetExtEntries = [&](iter& it, uint16_t& entryCnt) {
			for (; entryCnt && it != end(); entryCnt--, --it)
				SetInvalid(it);
			return it;
		};
		bool someRights = lastExtRight_it != end();	// there're some right bounds

		// 1) reset adjacent right bounds
		if (someRights) {
			if (someBS || newBS)	// newBS is true when there're only right bounds
				ResetExtEntries(lastExtRight_it, extRightCnt);
			else {						// ** 'negative' BS width
				lastExtRight_it->second.Reverse = true;		// change 'right' bound to 'left'
				if (extRightCnt > 0)
					ResetExtEntries(--lastExtRight_it, --extRightCnt);	// reset other adjacent bounds
			}
			lastExtRight_it = end();
		}
		// 2) reset adjacent left bounds
		if (someBS || !someRights)	// someRights is false when there're only left bounds
			ResetExtEntries(--it, extLeftCnt);
		else if(extLeftCnt > 0) {		// ** 'negative' BS width
			auto itR = ResetExtEntries(--it, --extLeftCnt);	// reset other adjacent entries
			if (itR != end()) {
				itR->second.Reverse = false;				// change 'left' bound to 'right
				// decrease the width of the updated BS if needed
				auto itL = prev(itR);
				FitToMinWidth(itL, itR, false);
			}
		}
	};

	// *** refine
	for (auto it = begin(); it != end(); it++)
	{
		const bool newRgn = grpNumb != it->second.GrpNumb;
		if (newRgn) {
			// 'close' previous region
			ResetAllExtEntries(it);
			// reset current region
			newBS = grpNumb = it->second.GrpNumb;
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
				extRightCnt += !someBS;
			extLeftCnt = 0;
		}
	}
	// 'close' last region
	ResetAllExtEntries(end());

	// *** extend narrows
	ExtendNarrowBSs();
}


// Sets scores for each BS within the group according to fragment coverage spline
//	@param VP: iterators for the left/right BS bounds
//	@param spline: local fragment coverage spline
//	@maxScore[in,out]: cumulative maximum score 
void SetBSscores(const vector<BS_map::iter>* VP, const Values& spline, float& maxScore)
{
	/*
	fill the score from left to right (for reverse)
	or from rigth to left (for direct) extended bounds
	*/
	chrlen	startPos = POS(VP[L].front());
	int32_t offset;
	BYTE	extLen;

	{	// left BS extensions
		const auto& vp = VP[L];
		extLen = BYTE(vp.size() - 1);
		offset = int32_t(POS(vp.front()) - startPos);
		for (BYTE i = 0; i < extLen; i++) {		// left to right, increasing offset
			if (BS_map::IsValid(vp[i])) {
				auto len = int(LEN(vp[i], vp[i + 1]));
				if (len <= 0)
					printf("+> %d  len: %d  numb: %d  score: %.3f\n", vp[i]->first, len, vp[i]->second.GrpNumb, SCORE(vp[i]));
				SCORE(vp[i]) = spline.AvrScoreInRange(offset, LEN(vp[i], vp[i + 1]));
			}
		}
	}
	{	// right BS extensions
		const auto& vp = VP[R];
		extLen = BYTE(vp.size() - 1);
		offset = int32_t(POS(vp.back()) - startPos);
		for (BYTE i = extLen; i; i--) {			// right to left, decreasing offset
			if (BS_map::IsValid(vp[i])) {
				auto len = int(LEN(vp[i - 1], vp[i]));
				if (len <= 0)
					printf("-> %d  len: %d  numb: %d  score: %.3f\n", vp[i]->first, len, vp[i]->second.GrpNumb, SCORE(vp[i]));
				SCORE(vp[i]) = spline.AvrScoreInRange(offset, LEN(vp[i - 1], vp[i]), -1);
			}
		}
	}
	// BS
	auto& start = VP[L].back()->second;
	auto& end = VP[R].front()->second;
	offset = int32_t(start.RefPos - startPos);
	auto len = int(end.RefPos - start.RefPos);
	if (len <= 0)
		printf("!> %d  len: %d  numb: %d  score: %.3f %.3f\n", VP[L].back()->first, len, start.GrpNumb, start.Score, end.Score);
	float score = spline.AvrScoreInRange(offset, end.RefPos - start.RefPos);

	if (maxScore < score)	maxScore = score;
	start.Score = end.Score = score;
}

void BS_map::SetGroupScores(iter& itStart, iter& itEnd, const Values& spline, float& maxScore)
{
	vector<iter> VP[2];	// 0 - forward, 1 - reversed

	VP[R].reserve(4), VP[L].reserve(4);
	for (auto& it = itStart; it != itEnd; it++)
		if (IsValid(it)) {
			if (it->second.Reverse && VP[R].size() && VP[L].size()) {
				SetBSscores(VP, spline, maxScore);
				VP[R].clear(), VP[L].clear();
			}
			VP[it->second.Reverse].push_back(it);
		}
	// last BS
	if (VP[R].size() && VP[L].size())
		SetBSscores(VP, spline, maxScore);
}

void BS_map::SetScore(const DataSet<TreatedCover>& fragCovers)
{
	SSpliner<coval> spliner(eCurveType::ROUGH, 5);
	auto& cover = fragCovers.TotalData();
	Values	spline;			// fragment coverage spline
	float	maxScore = 0;
	iter	itStart = end();
	iter	itEnd = end();
	chrlen	grpNumb = 1;

	spline.Reserve(Glob::FragLen * 3);
	// *** set scores
	for (auto it = begin(); it != end(); it++) {
		if (!IsValid(it))	continue;

		if (grpNumb != it->second.GrpNumb) {
			// build single fragment coverage spline for the whole group
			cover.SetLocalSpline(spliner, POS(itStart), POS(itEnd), spline);
			SetGroupScores(itStart, ++itEnd, spline, maxScore);

			spline.clear();
			itStart = end();
			grpNumb = it->second.GrpNumb;
		}
		if (itStart == end())
			itStart = it;
		itEnd = it;
	}
	// last group
	cover.SetLocalSpline(spliner, POS(itStart), POS(itEnd), spline);
	SetGroupScores(itStart, ++itEnd, spline, maxScore);

	// *** normalize scores
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
			//const fraglen len = LEN(end, start);
			float score = SCORE(start);
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
void BS_map::PrintWidthDistrib(const string& fName) const
{
	map<fraglen, vector<chrlen>> freq;
	chrlen	totalLen = 0, bsNumb = 0;
	TxtOutFile file(fName.c_str());

	// collect numbers
	DoBasic([&](citer& start, citer& end) {
		auto len = fraglen(LEN(start, end));
		totalLen += len;
		freq[len].push_back(++bsNumb);
		}
	);
	// print numbers
	const char* length = "width";
	file.Write("BS WIDTH FREQUENCY:\n");
	file.Write("%s cnt  numbers\n", length);
	for (const auto& item : freq) {
		file.Write("%4d %4u  ", item.first, UINT(item.second.size()));
		auto it = item.second.begin();
		file.Write("%d", *it);
		for (it++; it != item.second.end(); it++)
			file.Write(",%d", *it);
		file.Write("\n");
	}
	auto avrWidth = float(totalLen) / bsNumb;
	file.Write("average %s: %.2f\n", length, avrWidth);
	printf("\naverage %s: %.2f\n", length, avrWidth);
}

void BS_map::Print(chrid cID, const string& fName, bool selected, chrlen stopPos) const
{
	string format = "%9d % 5d  %c %8d %5.2f   %s\n";
	const char bound[]{ 'R','L' };
	IGVlocus locus(cID);

	TxtOutFile file(fName.c_str());
	file.Write("  pos    numb bnd  ref pos score   IGV view\n");
	for (const auto& x : *this) {
		if (stopPos && x.first > stopPos)	break;
		if (selected && !x.second.Score)	continue;
		format[5] = x.second.GrpNumb % 2 ? '-' : SPACE;	// odd numbers are aligned to the left, even numbers to the right
		format[20] = x.second.Score ? '2' : '0';		// zero score without fraction
		file.Write(format.c_str(),
			x.first,
			x.second.GrpNumb,
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

		LineAddUInts(POS(start), POS(end), ++bsNumb, true);	// 3 basic fields
		LineAddScore(SCORE(start), true);

		// *** extended boudaries info
		LineAddChar(DOT, true);
		lastSep[1] = VP[R].size() - 1;
		LineAddFloat(SCORE(end), VP[1].size() - 1 || lastSep[1]);	// ratio

		// *** adjacent deviations info
		for (BYTE s : {L, R}) {
			const BYTE vpLen = BYTE(VP[s].size() - 1);
			if (!vpLen) continue;

			auto& vp = VP[s];
			function<void(BYTE, char)> saveExtraPos = [this, &vp](BYTE i, char specChar) {
				LineAddArgs("%c%d", specChar, vp[i].Iter->first);
			};
			function<void(BYTE, char)> saveExtraVal = [this, &vp](BYTE i, char specChar) {
				if (specChar)	LineAddChar(specChar);
				LineAddScore(SCORE(vp[i].Iter), false);
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
				LineAddUInts(POS(it0->Iter), POS(it->Iter), bsNumb, true);
				LineAddFloat(it0->Val, true);	// BS score
				LineAddChar(DOT, true);
				LineAddInts(POS(it0->Iter), POS(it->Iter), true);
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
		LineAddUInts(POS(start), POS(end), bsNumb, true);
		LineAddFloat(SCORE(start), true);		// BS score
		LineAddChars(delims, reclen(strlen(delims)), false);
		LineAddFloat(SCORE(end), false);		// reverse/forward ratio
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
			POS(start) - Glob::ROI_ext,
			POS(end) + Glob::ROI_ext,
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
