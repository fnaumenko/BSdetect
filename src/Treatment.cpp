#include "Treatment.h"
#include "Distrib.h"
#include <algorithm>

#define PRINT

//const float PI = 3.14159265F;

//const eCurveType CurveTYPE = eCurveType::SMOOTH;
const eCurveType CurveTYPE = eCurveType::ROUGH;

bool Glob::IsPE = false;
bool Glob::FragLenUndef = true;
readlen Glob::ReadLen = 0;
fraglen Glob::FragLen = FragDefLEN;
#ifdef MY_DEBUG
chrid	Glob::CurrChrom = Chrom::UnID;
#endif
BYTE	Glob::BinWidth = 1;

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

//===== OSpecialWriter
#ifdef MY_DEBUG

void OSpecialWriter::WriteIncline(BYTE reverse, chrlen start, coviter& itStop)
{
	assert(Glob::CurrChrom != Chrom::UnID);

	const chrlen stop = itStop->first;
	const int ptCnt = int(stop) - int(start);
	int8_t k = 1;
	if (!reverse) {
		k = -1;
		start = stop;
		--itStop;
	}
	(_writers->_files)[reverse]->WriteIncline(Glob::CurrChrom, start, k * ptCnt, float(itStop->second) / ptCnt);
}
#endif

//===== Incline
#ifdef MY_DEBUG
OSpecialWriter* Incline::LineWriter = nullptr;
shared_ptr<FormWriter> Incline::OutFile = nullptr;
chrid Incline::cID;

void Incline::WriteLine(BYTE reverse, coviter& itStop) const
{
	if (LineWriter && LineWriter->IsWriterSet()) {
		StrandOps[!reverse].Next(itStop);
		LineWriter->WriteIncline(reverse, Pos, itStop);
	}
}

void Incline::WriteLocus() const
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

void TreatedCover::LinearRegr(coviter it, const coviter& itStop, const StrandOp& op, Incline& incline) const
{
	const int shift = op.Factor * itStop->first;
	float sumX2 = 0, sumXY = 0;
	coval sumX = 0, sumY = 0;
	coval PtCnt = 0;	// number of points along which the incline was build

	incline.Clear();
	for (; it != itStop; op.Next(it), PtCnt++) {
		const chrlen x = shift - op.Factor * it->first;	// X-distance between current it and itStop
		const coval y = op.GetPrev(it)->second;
		sumX += x;
		sumX2 += x * x;
		sumY += y;
		sumXY += x * y;
	}
	if (PtCnt < 2)	return;

	incline.Deriv =
		-(sumXY * PtCnt - sumX * sumY) /	// numerator
		(sumX2 * PtCnt - sumX * sumX);			// denominator
	//angl = coeff * 180 / PI;	if (angl < 0) angl = -angl;

	float x = (incline.Deriv * sumX + sumY) / (incline.Deriv * PtCnt);
	incline.Pos = itStop->first - op.Factor * chrlen(round(x));
	incline.TopPos = itStop->first;

#ifdef MY_DEBUG
	incline.WriteLocus();
#endif
}

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

void TreatedCover::PushIncline(BYTE reverse, const TracedPosVal& posVal, Inclines& inclines) const
{
	const StrandOp& opDir = StrandOps[reverse];		// forward operations
	const StrandOp& opInv = StrandOps[!reverse];	// inversed operations
	auto itStart = opDir.GetPrev(upper_bound(posVal.Pos(0)));	// min val: right for forward, left for reverse
	auto itStop = itStart;										// max val: left for forward, right for reverse

	// *** set itStop
	if (reverse)	for (itStop++; itStop->first < posVal.Pos(1); itStop++);
	else			for (itStop--; itStop->first > posVal.Pos(1); itStop--);

	// checking if there is an ordinary tag pool before the current start, masked by a spline
	if (!opDir.GetPrev(itStart)->second)
		opDir.Next(itStart);		// skip the gap after ordinary tag pool

	// *** it stop correction: calc the best BS position by successive linear regression approximations
	Incline incline;
	LinearRegr(itStart, itStop, opDir, incline);
	if (!incline.Valid())	return;

	coviter it0 = itStop;
	Incline incline1{};
	auto compare = [&](coviter& it) {
		LinearRegr(itStart, it, opDir, incline1);
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

#ifdef MY_DEBUG
	incline.WriteLine(reverse, itStop);
#endif
	inclines.Add(incline);
	inclines.SpreadTopCover(reverse, posVal.Val(reverse));
}


//===== CombCover

void CombCover::SetUnsortedInput()
{
	_data->TotalData().SetUnsortedInput();	// if no total data is defined, call method twice for the forward data
	_data->StrandData(FWD).SetUnsortedInput();
	_data->StrandData(RVS).SetUnsortedInput();
}

void CombCover::AddExtRead(const Region& read, bool reverse)
{
	Region frag(read, Glob::FragLen, reverse);

	_data->TotalData().AddRegionByCond(frag);	// total frag coverage
	AddRead(frag, reverse);						// strand frag coverage
}

void CombCover::FillExtRead(const Reads& reads)
{
	for (BYTE s : {0, 1}) {
		if(s)
			_data->TotalData().SetUnsortedInput();	// set unsorted because it's being refilled
		for (auto& rd : reads.GetReads(s))
			AddExtRead(rd, s);
	}
}

//===== template


template<typename T>
/// <summary>
/// Recursively marks as invalid forward/reversed regions that intersect not strictly pairwise,
/// i.e. which are a linked chain
/// </summary>
/// <typeparam name="T">class over which polymorphic methods for working with a collection of regions are defined</typeparam>
/// <param name="it">forward/reversed region iterators</param>
/// <param name="itEnd">forward/reversed region and iterators</param>
/// <returns>true if if forward or reversed end iterator is reached</returns>
bool DiscardOverlapChain_(typename T::iterator it[2], typename T::const_iterator itEnd[2])
{
	BYTE s = T::End(it[R]) > T::End(it[L]);	// index of the left ended region: R - forward, L - reverse

	if (T::Start(it[!s]) > T::End(it[s])) {
		T::Discard(it[s]); it[s]++; return false;	// single region
	}
	if (++it[s] == itEnd[s])		return true;
	if (!T::Accepted(it[s]))		return false;
	if (T::End(it[!s]) < T::Start(it[s])) { it[!s]++; return false; }
	T::Discard(prev(it[s]));
	T::Discard(it[s]);
	T::Discard(it[!s]);
	return DiscardOverlapChain_<T>(it, itEnd);
}

template<typename T>
/// <summary>
/// Marks as invalid forward/reversed regions that intersect not strictly pairwise,
/// i.e. which are a linked chain
/// </summary>
/// <typeparam name="T">class over which polymorphic methods for working with a collection of regions are defined</typeparam>
/// <param name="rgns">forward/reversed regions</param>
void DiscardMultiOverlapRegions(T rgns[2])
{
	typename T::const_iterator itEnd[2]{ rgns[R].end(), rgns[L].end() };
	typename T::iterator it[2]{ rgns[R].begin(), rgns[L].begin() };

	DiscardOverlapChain_<T>(it, itEnd);
	while (it[R] != itEnd[R] && it[L] != itEnd[L])
		if (!T::Accepted(it[R]))		it[R]++;
		else if (!T::Accepted(it[L]))	it[L]++;
		else if (DiscardOverlapChain_<T>(it, itEnd))	break;
}

#ifdef MY_DEBUG
template<typename T>
void CheckSingleOverlapping(const T rgns[2], fraglen minOverlapLen)
{
	chrlen numb = 0;
	bool done = true;
	typename T::const_iterator itEnd[2]{ rgns[R].end(), rgns[L].end() };
	typename T::const_iterator it[2]{ rgns[R].begin(), rgns[L].begin() };

	printf("Overlapping check:");
	//printf("\n  N strand  start  -   end val accept");
	while (it[R] != itEnd[R] && it[L] != itEnd[L]) {
		BYTE s = T::End(it[R]) > T::End(it[L]);	// index of the left ended region: 0 - direct, 1 - reverse

		if (T::Start(it[!s]) + minOverlapLen > T::End(it[s])) { it[s]++; continue; }
		if (!T::Accepted(it[R]) ^ !T::Accepted(it[L])) {
			printf("\n%3d FVD: %d-%d  %d, RVS: %d-%d  %d", ++numb,
				T::Start(it[R]), T::End(it[R]), T::Accepted(it[R]),
				T::Start(it[L]), T::End(it[L]), T::Accepted(it[L])
			);
			done = false;
		}
		it[R]++, it[L]++;
	}
	if (done)	printf(" done");
	printf("\n");
}
#endif

//===== CoverRegions

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
	printf("%s\n%s\n", titles[isFeatures], title);
	PrintSolidLine(USHORT(strlen(title) + 2));
	for (BYTE s = 0; s < 1 + strands; s++) {
		chrlen rawLen = 0, refineLen = 0;
		const auto& rgn = rgns[s];
		const auto itEnd = rgn.end();
		size_t realCnt = 0;

		for (auto it = rgn.begin(); it != itEnd; it++) {
			const chrlen len = T::Length(it);
			rawLen += len;
			if (T::Accepted(it)) refineLen += len, realCnt++;
		}
		printf(format[isFeatures], sStrandTITLES[s + strands],
			rgn.size(), Percent(rawLen, chrLen),
			realCnt, Percent(refineLen, chrLen));
	}
	printf("\n");
}

//===== CoverRegions

void CoverRegions::SetPotentialRegions(const TreatedCover& cover, chrlen capacity, coval cutoff)
{
	const auto minLen = Glob::FragLen;	// empirical minLen obtained in tests
	chrlen	start = 0, end = 0;
	coviter itStart, itEnd;

	coval maxVal = 0;
	chrlen maxPos = 0;

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
				if (maxVal < val) {	maxVal = val; maxPos = itStart->first; }
			}
			start = 0;
		}
	IGVlocus locus(Glob::CurrChrom);
	printf("MAX: VAL: %d POS: %d\t%s\n", maxVal, maxPos, locus.Print(maxPos));
}

bool CoverRegions::SetTopPeakRegions(const DataSet<TreatedCover>& fragCover, coval cutoff)
{
	const USHORT TOP_RGNS_CNT = 1000;
	const USHORT THRECHOLD_COEFF = 100;
	const BYTE VAL_BAR_WIDTH = 1;
	auto& cover = fragCover.TotalData();
	const auto minLen = Glob::FragLen;	// empirical minLen obtained in tests
	chrlen	start = 0, end = 0;
	coviter itStart, itEnd;

	coval maxVal = 0;

	multimap<coval, coviter> peakFreq;	// peak frequency: number of peaks with the same max coverage

	// *** fill peak frequency
	for (auto it0 = cover.cbegin(), it = it0; it != cover.end(); it0 = it++)
		if (it->second >= cutoff || it0->second >= cutoff) {	// look for summit
			if (!start)
				start = (itStart = it)->first;	// set after prev region processed only
			end = (itEnd = it)->first;
		}
		else {
			if (start && end - start > minLen) {
				// find raw summit
				coviter itPeak = itStart;
				for (auto it = next(itStart); it != itEnd; it++)
					if (itPeak->second < it->second)
						itPeak = it;
				//this->emplace_back(itStart, itEnd, val);
				peakFreq.insert(pair<coval, coviter>(itPeak->second / VAL_BAR_WIDTH, itPeak));
				if (maxVal < itPeak->second) { maxVal = itPeak->second; }
			}
			start = 0;
		}

	if(Glob::BinWidth)
	{
		printf("ROW PEAK FREQUENCY\n");
		coval val = CHRLEN_MAX;
		for (const auto& el : peakFreq)
			if (val != el.first)	val = el.first, printf("%d\t%zu\n",val, peakFreq.count(val));
	}

	// *** set max peak frequency
	UINT maxFreq = 0;
	coval val = CHRLEN_MAX;
	for (const auto& el : peakFreq)
		if (val != el.first) {
			auto freq = UINT(peakFreq.count(val = el.first));
			if (maxFreq < freq)	
				maxFreq = freq;
			else break;
		}


	// *** fill regions of interest
	auto maxfreq = maxFreq / THRECHOLD_COEFF;
	this->reserve(TOP_RGNS_CNT + 1);
	//const auto ext = 2 * Glob::FragLen / 3;
	const auto ext = Glob::FragLen;
	coval lowLimVal = 0;
	UINT	lowLimFreq = 0;

	for (auto it = prev(peakFreq.cend()); it != peakFreq.cbegin(); it--) {
		if ((peakFreq.count(it->first) > maxfreq || size() > TOP_RGNS_CNT)	// strong cutting off - for reach coverage
			&& (size() > TOP_RGNS_CNT / 2))									// weak cutting off - for peru coverage
		{
			if (Verb::Level(Verb::DBG)) {
				it++;
				lowLimVal = it->second->second;
				lowLimFreq = UINT(peakFreq.count(it->first));
			}
			break;
		}
		// expand summit to the region of interest
		auto& itPeak = it->second;
		auto itStart = prev(itPeak);
		for (auto start = itPeak->first - ext; itStart->first > start; itStart--);	// itStart
		auto itEnd = next(itPeak);
		for (auto end = itPeak->first + ext; itEnd->first < end; itEnd++);			// itEnd

		emplace_back(itStart, itEnd, it->first);
	}

	sort(this->begin(), this->end(),
		[](const CoverRegion& r1, const CoverRegion& r2) { return r1.Start() < r2.Start(); }
	);
	Verb::PrintMsgVar(Verb::DBG,
		//"Top Peak Regions: count: %zu  cutoffVal: %d, maxVal: %d, cutoffFreq %u, maxFreq: %u\n", 
		//size(), lowLimVal, maxVal, lowLimFreq, maxFreq);
		"Top Peak Regions: count %zu,  ValRatio %.3f, valMax %d, FreqRatio %.3f, FreqMax %u\n",
		size(), float(lowLimVal)/maxVal, maxVal, float(lowLimFreq)/ maxFreq, maxFreq);
	return empty();
}

#ifdef MY_DEBUG

const string distExt = ".dist";

void CoverRegions::PrintScoreDistrib(const string& fname, bool all) const
{
	map<coval, chrlen> freq;

	for (const auto& rgn : *this)
		if(all || !rgn.Accepted())
			freq[rgn.value]++;
	
	FormWriter file((fname + distExt).c_str());
	file.Write("score\tfreq\n");
	for (const auto& item : freq)
		file.Write("%d\t%d\n", item.first, item.second);
}
#endif

//===== DataCoverRegions

bool DataCoverRegions::SetPotentialRegions(const DataSet<TreatedCover>& cover, chrlen cLen, coval cutoff, bool noMultiOverl)
{
	chrlen capacity = cLen / (Glob::FragLen * 100);
	if (Glob::IsPE)
		TotalData().SetPotentialRegions(cover.TotalData(), capacity, cutoff * 2);
	else {
		StrandData(FWD).SetPotentialRegions(cover.StrandData(FWD), capacity, cutoff);
		StrandData(RVS).SetPotentialRegions(cover.StrandData(RVS), capacity, cutoff);

		auto minOverlap = fraglen(0.7f * Glob::FragLen); // empirical minOverlap obtained in tests
		// In the working version, only strand regions are defined.
		// However, for debugging purposes, total regions can also be defined.
		// So we should call StrandData() instead of Data()
		auto data = StrandData();
		DiscardNonOverlapRegions<CoverRegions>(data, minOverlap);
		if (noMultiOverl) {
			DiscardMultiOverlapRegions<CoverRegions>(data);
#ifdef MY_DEBUG
			CheckSingleOverlapping<CoverRegions>(data, minOverlap);
#endif
		}
		if (Verb::Level(Verb::DBG))
			PrintRegionStats<CoverRegions>(data, cLen);
	}
	if (!Empty())	return false;
	Verb::PrintMsg(Verb::CRIT, "No enriched regions found");
	return true;
}

#ifdef MY_DEBUG
void DataCoverRegions::PrintScoreDistrib(const string& fname, bool all) const
{

	if (Glob::IsPE)
		TotalData().PrintScoreDistrib(fname, all);
	else {
		StrandData(FWD).PrintScoreDistrib(fname + sStrandEXT[FWD], all);
		StrandData(RVS).PrintScoreDistrib(fname + sStrandEXT[RVS], all);
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
	/*
	If the starting position is already a maximum, 
	or ending position is a maximum (the "maximum-cut" case), it is ignored.
	*/
	float val0 = front();
	bool increase = false;
	USHORT equalCnt = 0;	// count of equal values

	for (auto it = next(begin()); it != end(); val0 = *it, it++, startPos++)
		if (*it == val0)
			equalCnt++;
		else {
			if (*it > val0)			// increase value
				increase = true;
			else if (increase) {	// decrease value
				if(equalCnt < 20)
					/*
					The length of flat summit is more than 10 means either a read anomaly
					or an insignificant read clustering. Ignored.
					*/
					pos.push_back(startPos - equalCnt / 2);	// the centre of the 'flat' summit
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
	chrlen i = 0;
	for (const auto& x : *this) {
		//if (stopNumb && x.second.GrpNumb > stopNumb)	break;
		if (x.second.MaxVal()) {
			chrlen end = x.first + x.second.Length();
			printf("%3d %d\t%d\t%2.2f\t%s\n",
				x.second.GrpNumb, x.first, end, x.second.MaxVal(), locus.Print(x.first, end));
		}
	}
}
#endif

void ValuesMap::BuildRegionSpline(bool reverse, const TreatedCover& rCover, const CoverRegion& rgn, fraglen splineBase)
{
	/*
	Both forward & reverse splines are built from left to right
	*/
	assert(Glob::ReadLen);
	coviter it0;	// at the beginning the start it, then used as a variable
	coviter itEnd;	// the end it
	SSpliner<coval> spliner(CurveTYPE, splineBase);
	chrlen pos = rgn.itStart->first - spliner.SilentLength() / 2;	// ??? half of SilentLength is a good git to the region start

	/*
	incrementing it0 (for forward cover) | itEnd (for reversed cover) is required 
	to smoothly start|complete the spline in case where there are no more reads in the potential region
	after the last significant coverage, 'emptiness'
	*/

	// *** set it0
	it0 = rCover.upper_bound(pos);		// find cover iterator for region start position
	if (it0 != rCover.begin())	it0--;	// decrement it0 to smooth start the spline
	//if (reverse) {
	//	if (!it0->second)	it0--;
	//	if (it0 != rCover.begin())	it0--;
	//}

	// *** set itEnd
	pos = rgn.itEnd->first + spliner.SilentLength();
	for (itEnd = it0, advance(itEnd, 10); itEnd != rCover.end() && itEnd->first < pos; itEnd++);	// 10 iterators is not enough for significant coverage in any case
	if (itEnd == rCover.end())
		itEnd--;
	else {
		// increment itEnd to smooth complete the spline
		//if (itEnd->second)	itEnd++;	// itEnd->second != 0 means that is not the last iterator
		if (next(itEnd) != rCover.end())	itEnd++;
		if (next(itEnd) != rCover.end())	itEnd++;
	}

	// *** spline via covmap local copy, filtering unsignificant splines
	chrlen newPos = 0;
	chrlen lastZeroPos = 0;
	Values vals;

	// adds spline, eliminating unsignificant one
	auto addDecentRgn = [this,&vals,&newPos]() {
		if (vals.MaxVal() > 2.)
			this->AddRegion(newPos, vals);
		else
			vals.Clear();
	};

	chrlen currPos = it0->first;// +1;
	for (auto it = next(it0); /*it->first <= pos &&*/ it != itEnd; it0++, it++) {	// loop through the cover

		// skip reads that do not form a continuous coverage with the main heap
		if (lastZeroPos) {
			if (it0->first - lastZeroPos >= spliner.SilentLength()) {
				// skip 2 duplicated reads or standalone read or contiguous single reads
				if (!it->second) {
					if (it0->second <= 2) {
						if (++it == itEnd)	
							break;
						it0++;
					}
				}
				// skip 2 overlapping reads
				else if (it->second == 2) {
					auto it1 = next(it);
					if (it1 == itEnd)	
						break;
					if (it1->second == 1) {
						if (++it1 == itEnd)		
							break;
						if (!it1->second) {
							it = it1;
							if (++it == itEnd)	
								break;
							advance(it0, 3);
						}
					}
				}
			}
			lastZeroPos = 0;
		}

		// treat other reads
		for (; currPos <= it->first; currPos++) {		// loop through positions between iterators
			float val = spliner.Push(it0->second);
			if (val) {
				if (!vals.Length())
					newPos = spliner.CorrectX(currPos);	// start new spline
				vals.AddValue(val);
			}
			else 
				if (vals.Length())
					addDecentRgn();	// end new spline
		}

		if (spliner.CorrectX(currPos) > rgn.itEnd->first)// + spliner.SilentLength())
			break;
		if (!it0->second)
			lastZeroPos = it0->first;
	}
	if (vals.Length())
		addDecentRgn();
}

void ValuesMap::AddRegion(chrlen pos, Values& vals)
{
	_maxVal = vals.MaxVal();	// duplicated because _maxVal is also used as 'empty' sign
	emplace(pos, move(vals));
	vals.Reserve();
}

void ValuesMap::BuildSpline(bool reverse, const TreatedCover& rCover, const CoverRegions& rgns, fraglen splineBase)
{
	for (const auto& rgn : rgns)
		if (rgn.Accepted())
			BuildRegionSpline(reverse, rCover, rgn, splineBase);
}

void ValuesMap::DiscardNonOverlaps()
{
	DiscardNonOverlapRegions<ValuesMap>(this, SSpliner<coval>::SilentLength(CurveTYPE, ReadSplineBASE));
}

void ValuesMap::NumberGroups()
{
	// minimum overlap length: 10 - empirical addition
	const fraglen minOverlap = SSpliner<coval>::SilentLength(CurveTYPE, ReadSplineBASE) + 10;
	ValuesMap::Iter it[2]	{ this[R].begin(),	this[L].begin() };
	ValuesMap::Iter itEnd[2]{ this[R].end(),	this[L].end()	};
	chrlen numb = 1;

	auto isOverlap = [&](BYTE s) {
		auto it1 = next(it[s]);
		auto nextStart = it1 != itEnd[s] ? ValuesMap::Start(it1) : CHRLEN_MAX;
		if (nextStart >= End(it[!s])) 
			return false;		// non-overlapping
		(++it[s])->second.GrpNumb = numb;
		if (nextStart + minOverlap >= End(it[!s]))
			it[s]->second.Discard();	// overlapping is insufficient
		return true;			// overlapping
	};

	while (it[R] != itEnd[R] && it[L] != itEnd[L]) {
		if (!it[R]->second.MaxVal()) { it[R]++; continue; }
		if (!it[L]->second.MaxVal()) { it[L]++; continue; }

		it[R]->second.GrpNumb = it[L]->second.GrpNumb = numb;
		for (bool overlap = true; overlap; )
			if (!(overlap = isOverlap(0)))
				overlap = isOverlap(1);
		numb++;
		it[R]++, it[L]++;
	}
}

void ValuesMap::PrintStat(chrlen clen) const
{
	PrintRegionStats<ValuesMap>(this, clen);
}

//===== DataValuesMap

void DataValuesMap::BuildSpline(
	const DataSet<TreatedCover>& rCover, const DataCoverRegions& rgns, fraglen splineBase)
{
	const BYTE strand = !Glob::IsPE;	// TOTAL for PE or FWD for SE
	StrandData(FWD).BuildSpline(false, rCover.StrandData(FWD), rgns.StrandData(eStrand(strand)), splineBase);
	StrandData(RVS).BuildSpline(true, rCover.StrandData(RVS), rgns.StrandData(eStrand(2*strand)), splineBase);
}

float DataValuesMap::GetPeakPosDiff() 
{
	/*
	It is assumed that the regions listed have one peak of each strand.
	However, a check is made to see if this is the case.
	If not, the region is rejected for calculation.
	*/
	DiscardMultiOverlapRegions<ValuesMap>(StrandData());
#ifdef MY_DEBUG
	//CheckSingleOverlapping<ValuesMap>(StrandData(), 0);
	IGVlocus locus(Glob::CurrChrom);
#endif

	auto& pData = StrandData(FWD);	// "positive"
	auto& nData = StrandData(RVS);	// "negative"
	USHORT rejected = 0;
	vector<SHORT> diffs;
	vector<chrlen> pPos, nPos;	// max positions in a positive, negative splines

	// get the difference of the splines maximum
	diffs.reserve(pData.size());	// suppose one summit for the region's spline
	pPos.reserve(2);
	nPos.reserve(2);
	for (auto itP = pData.begin(), itN = nData.begin(); itP != pData.end() && itN != nData.end(); itP++, itN++) {
		while (!pData.Accepted(itP)) if (++itP == pData.end()) goto out; 
		while (!nData.Accepted(itN)) if (++itN == nData.end()) goto out;
		itP->second.GetMaxValPos(itP->first, pPos);
		itN->second.GetMaxValPos(itN->first, nPos);
		// peaks within each region
		if (pPos.size() == nPos.size())
			for (auto itP = pPos.begin(), itN = nPos.begin(); itP != pPos.end(); itP++, itN++) {
				diffs.push_back(SHORT(*itN - *itP));
				//if (/*abs*/(SHORT(*itN - *itP)) <= 0)
				//	printf("%d\t%d\t%4d\t%s\n", *itN, *itP, SHORT(*itN - *itP), locus.Print(*itN));
			}
		else {
			rejected++;		// reject region with different peak count
#ifdef MY_DEBUG
			//printf("FWD %zu %s\t", pPos.size(), pPos.size() ? locus.Print(pPos[0]) : "\t\t\t");
			//printf("RVS %zu %s\n", nPos.size(), nPos.size() ? locus.Print(nPos[0]) : "");
#endif
		}
		pPos.clear();
		nPos.clear();
	}
out:if(Verb::Level(Verb::DBG))
		if (rejected)
			printf("%3d (%.0f%%) rejected regions;\t", rejected, Percent(rejected, pData.size()));
		else
			printf("\t\t\t\t");

	// *** find most frequent value
	float mode = 0;
	for(BYTE binW : {5, 9, 15, 21})
	//for (BYTE binW : {3,5,9,15,21,31})
	{
		static const int8_t factors[]{ -1,1 };
		// *** fill differences frequent value
		Distrib freq;
		for (auto diff : diffs) {
			// signbit()?
			short bin = (diff / binW) * binW + factors[diff > 0] * binW / 2;	// position in the middle of the bin
			freq.IncrFreq(bin);
		}
#ifdef PRINT
		printf("\n>>> BIN WIDTH %d\n", int(binW));
		printf("DIFFS FREQUENCY DISTRIBUTION  size: %zu\n", freq.Size());
		//freq.Print(cout);
#endif
		freq.CalcADParams(Distrib::LNORM, Distrib::INTERPOL);
		freq.PrintADParams(cout, false, true);
	}
	return mode;
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

//===== Inclines

void Inclines::Add(const Incline& incline)
{
	if (!size() || back() != incline)
		push_back(incline);
}

void Inclines::SpreadTopCover(BYTE reverse, float topCover)
{
	auto sz = USHORT(size());
	if (sz) {
		auto& lastItem = back();
		lastItem.TopCover = topCover;
		if (sz > 1) {
			/*
			we have to use indices instead of iterators
			because of possible capacity overrun between this method calls,
			and as a consequence, displacement of the internal array and changing begin()|end() iterators
			*/
			const auto lastInd = USHORT(sz - 1);
			/*
			if lastItem.Pos and prevItem.TopPos are equal,
			the inclines are already considered to belong to different bunches
			*/
			if (StrandOps[reverse].EqLess(lastItem.Pos, at(lastInd - 1).TopPos)) {
				// spread TopCover
				for (auto i = _firstInd; i < lastInd; i++) {
					auto& incline = at(i);
					incline.TopCover = _topCover;	// for one of the inclines the value will be rewritten to the same; never mind
					incline.Bunched = lastInd - _firstInd > 1;
				}
				_topCover = topCover;
				_firstInd = lastInd;
				return;
			}
		}
		if (_topCover < topCover)
			_topCover = topCover;
	}
}


//===== BoundsValues

void BoundsValues::CollectForwardInclines(const TreatedCover& rCover, Inclines& inclines) const
{
	TracedPosVal posVal;

	for (auto it = rbegin(); it != rend(); it++) {	// loop through one derivative region
		posVal.Set(R, it);
		cout << this->TopCover() << TAB << posVal.Val(R) << LF;
		rCover.PushIncline(R, posVal, inclines);
		posVal.Retain();
	}
}

void BoundsValues::CollectReverseInclines(const TreatedCover& rCover, Inclines& inclines) const
{
	TracedPosVal posVal;

	for (auto it = begin(); it != end(); it++) {		// loop through one derivative region
		posVal.Set(L, it);
		rCover.PushIncline(L, posVal, inclines);
		posVal.Retain();
	}
}

void BoundsValues::AddSignifValues(const tValuesMap::value_type& spline, chrlen relPos, Values& deriv)
{
	if (_maxVal < deriv.MaxVal())	_maxVal = deriv.MaxVal();
	auto& vals = spline.second;
	_grpNumb = vals.GrpNumb;
	_topCover = vals.MaxVal();
	emplace_back(
		spline.first + relPos,
		vals.Value(relPos),
		vals.Value(relPos + deriv.Length()),
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

	float topVal = 0;

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

	// set relative top values
	for (auto& deriv : *this) {
		//deriv.SetRelValue(_maxVal);
		//cout << _grpNumb << TAB << deriv._pos << TAB<< _maxVal << TAB << deriv.MaxVal() << TAB << deriv.RelValue() << LF;
	}
}


//===== BoundsValuesMap

void BoundsValuesMap::BuildDerivs(int factor, const ValuesMap& splines)
{
	for (const auto& spline : splines) {
		if (!spline.second.MaxVal())	continue;

		Values deriv;			deriv.Reserve();
		BoundsValues derivSet;	derivSet.reserve(4);
		fraglen pos = 0;	// start relative position
		const auto& itEnd = spline.second.end();

		float topVal = 0;

		// loop through spline
		for (auto it0 = spline.second.begin(), it = next(it0); it != itEnd; it0++, it++) {
			float tang = (*it - *it0) * factor;		// flip sign of pos strand delta
			if (tang > 0) {							// cut off tangent of a falling spline
				deriv.AddValue(tang);
				if (topVal < *it)
					topVal = *it;
			}
			else {
				pos++;
				if (deriv.MaxVal()) {
					fraglen len = deriv.Length();	// deriv will be cleared

					derivSet.AddValues(spline, pos, deriv);
					pos += len;
					topVal = 0;
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
			printf("  %2.2f  %2.2f:\t%d %d\t%d\n", rvs.MaxVal(), rvs.RelValue(), rvs.GrpNumb, rvs.Start(), rvs.Length());
	}
}
#endif


//===== BS_map

#define	POS(it)	(it)->second.RefPos
#define LEN(itStart,itEnd)	short(POS(itEnd) - POS(itStart))
#define	SCORE(it)	(it)->second.Score
#define	GrpNUMB(it)	(it)->second.GrpNumb

#define	TopCOVER(it)	(it)->second.TopCover
#define REAL(it)		(it)->second.Real

void BS_map::AddPos(BYTE reverse, chrlen grpNumb, const Incline& incline)
{
	// all positions are added in ascending order
	chrlen pos = incline.Pos + StrandOps[reverse].Factor * Glob::ReadLen;
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

		_lastIt = emplace_hint(lastIt, pos, move(BS_bound(1, grpNumb, incline)));
		POS(_lastIt) = _lastIt->first;
	}
	else {
		auto it = emplace_hint(end(), pos, move(BS_bound(0, grpNumb, incline)));
		POS(it) = it->first;
	}
}

void BS_map::AddBounds(BYTE reverse, chrlen grpNumb, /*float topCover,*/ Inclines& inclines)
{
	// sort inclines by ascending (for forward) / descending (for reverde) reference positions
	sort(inclines.begin(), inclines.end(),
		[&reverse](const Incline& i1, const Incline& i2) {
			return StrandOps[reverse].Less(i1.Pos, i2.Pos);
		}
	);

	// find permissible derivative limit
	float minDeriv = 0;
	for (const auto& i : inclines)
		minDeriv += i.Deriv;
	minDeriv /= inclines.size();	// arithmetic mean
	minDeriv *= 0.15F;				// 0.15: empirical rate

	// *** Intersection Filter & Insignificant Incline Filter:
	// insert only the best bounds (formed by the steepest inclines), and with a significant derivative
	auto it = inclines.cbegin();
	chrlen addedPos = it->TopPos;

	AddPos(reverse, grpNumb, *it);		// definitely add the first, steepest (tightest) incline
	for (it++; it != inclines.cend(); it++)
		if (StrandOps[!reverse].EqLess(it->TopPos, addedPos)
		&& it->Deriv > minDeriv)	
		{
			AddPos(reverse, grpNumb, *it);
			addedPos = it->TopPos;
		}
}

void BS_map::SetBounds(BYTE reverse, const BoundsValuesMap& derivs, const TreatedCover& rCover)
{
	using tSetBSpos = void(BoundsValues::*)(const TreatedCover&, Inclines& inclines) const;
	tSetBSpos fcollectInclines = reverse ?
		&BoundsValues::CollectReverseInclines :
		&BoundsValues::CollectForwardInclines;
	Inclines inclines;
	inclines.reserve(6);

	for (const auto& d : derivs) {		// loop through the derivative regions
		inclines.Clear();
		(d.second.*fcollectInclines)(rCover, inclines);
		if (inclines.size()) {
			inclines.SetTopCover(d.second.TopCover());
			AddBounds(reverse, d.second.GrpNumb(),/* d.second.MaxVal(),*/ inclines);
		}
	}
}

const short MIN_BS_WIDTH = 5;

// Fits the BS length to a minimum by adjusting the BS reference positions
//	@param start: BS iterator pointed to base left bound
//	@param end: BS iterator pointed to base right bound
//	@param isLess: if true then the width should be less than the minimum
//	@returns: true if the condition is met and the width is adjusted
bool FitToMinLength(BS_map::iter& start, BS_map::iter& end, bool isLess)
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
	if (FitToMinLength(itL, itR, true)) {
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

void BS_map::ExtendSingleNarrowBS(iter& start, const iter& end)
{
	const iter itEnd = next(end) == this->end() ? this->end() : next(end);

	for (auto& it = start; it != itEnd; it++)
		if (IsValid(it) && !REVERSE(it)) {
			auto itL = prev(it);
			ExtendNarrowBS(itL, it);
			break;
		}
}

void BS_map::ExtendNarrowBSsInGroup(iter& start, const iter& stop, bool narrowBS, bool closeProx)
{
	bool lastLeft = true;
	iter itEnd = next(stop) == end() ? end() : next(stop);
	vector<pair<iter, iter>> bss;	// BS start-end collection

	bss.reserve(4);
	// collect BS
	for (auto& it = start; it != itEnd; it++)
		if (IsValid(it))
			if (REVERSE(it))
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

void BS_map::ExtendNarrowBSs0()
{
	/*
	here the groups are 'canonical', e.g. 
	L L L R R L L R L R R
	    |_|     |_| |_|
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

		if (grpNumb != GrpNUMB(it)) {
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
			grpNumb = GrpNUMB(it);
		}

		if (REVERSE(it))
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

void BS_map::ExtendNarrowBSs()
{
	/*
	here the groups are 'canonical', e.g.
	L L L R R L L R L R R
		|_|     |_| |_|
	*/

	for (auto it = begin(); it != end(); it++) {
		if (REAL(it)) {	// always reversed (left) bound
			//auto pos = POS(it);
			auto nextIt = next(it);
			if (LEN(it, nextIt) < MIN_BS_WIDTH)
				FitToMinLength(it, nextIt, true);
			it++;
		}
	}
}

void BS_map::Set(const DataBoundsValuesMap& derivs, const DataSet<TreatedCover>& rCover)
{
	SetBounds(R, derivs.StrandData(FWD), rCover.StrandData(FWD));
	_lastIt = begin();
	SetBounds(L, derivs.StrandData(RVS), rCover.StrandData(RVS));
}

void BS_map::Refine0()
{
	/*
	BS Left entry (bound): is formed by reverse reads; BS Right entry (bound): is formed by forward reads

	Options for placing bounds in a region:
	canonical:					[[L] [R]]
	adjacent right/left bounds:	[R] [[L] [R]] [L]
	'negative' BS width:			[R] [L]

	Method brings the instance to canonical form, resetting the score of all adjacent elements to zero,
	and marking BSs with 'negative' width.
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
				REVERSE(lastExtRight_it) = true;	// change 'right' bound to 'left'
				if (extRightCnt > 0)
					ResetExtEntries(--lastExtRight_it, --extRightCnt);	// reset other adjacent bounds
			}
			lastExtRight_it = end();
		}
		// 2) reset adjacent left bounds
		if (someBS || !someRights)	// someRights is false when there're only left bounds
			ResetExtEntries(--it, extLeftCnt);
		else if (extLeftCnt > 0) {		// ** 'negative' BS width
			auto itR = ResetExtEntries(--it, --extLeftCnt);	// reset other adjacent entries
			if (itR != end()) {
				REVERSE(itR) = false;				// change 'left' bound to 'right
				// decrease the length of the updated BS if needed
				auto itL = prev(itR);
				FitToMinLength(itL, itR, false);
			}
		}
	};

	// *** refine
	for (auto it = begin(); it != end(); it++) {
		const bool newRgn = grpNumb != GrpNUMB(it);

		if (newRgn) {
			// 'close' previous region
			ResetAllExtEntries(it);
			// reset current region
			newBS = grpNumb = GrpNUMB(it);	// true
			someBS = extLeftCnt = extRightCnt = 0;
		}

		if (REVERSE(it)) {	// Left bound
			extLeftCnt++;
			newBS = false;
		}
		else {				// Right bound
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
	ExtendNarrowBSs0();
}

void BS_map::Ñandidates::SetReal()
{
	auto sz = size();
	if (!sz) return;

	auto setReal0 = [this](iter& itL, iter& itR, iterator& it) {
		REAL(itL) = REAL(itR) = true;
		if (REVERSE(itR)) {
			bool turnOver = true;
			if (it->adjLeft) {
				auto len = LEN(itL, itR);
				auto len0 = LEN(prev(itL), itL);
				USHORT ind = it->index;
				bool checkScore = ind ?
					it->avrCover < 2 * at(--ind).avrCover:
					true;

				if (len > 2 * len0 && checkScore)			// empirical ratio 2
					turnOver = false;
			}
			if (turnOver && it->adjRight) {
				auto len = LEN(itL, itR);
				auto len1 = LEN(itR, next(itR));
				USHORT ind = it->index;
				bool checkScore = ind < size() - 1 ?
					it->avrCover < 2 * at(++ind).avrCover :
					true;

				if (len > 2 * len1 && checkScore)		// empirical ratio 2
					turnOver = false;
			}
			if(turnOver)
				REVERSE(itR) = false,
				REVERSE(itL) = true;
		}
	};

	auto setRealWithCheckBounds = [this](iter& itL, iter& itR, iterator& it) {
		REAL(itL) = REAL(itR) = true;
		{
			int ind = it->index;
			for (auto itLL = prev(itL); ind >= 0; ind--, itLL--)
				if (!REVERSE(itLL)) {
					itLL->second.Related = false;
					//break;
				}
		}
		{
			int ind = it->index;
			for (auto itRR = next(itR); ind < size(); ind++, itRR++)
				if (REVERSE(itRR)) {
					itRR->second.Related = false;
					//break;
				}
		}
	};

	auto setReal = [=](iter& itL, iter& itR, iterator& it) {

		if (REVERSE(itR)) {
			bool turnOver = true;
			if (it->adjLeft) {
				USHORT ind = it->index;
				if (ind == 0 || it->relCover < 2 * at(--ind).relCover)
					turnOver = false;
			}
			if (turnOver && it->adjRight) {
				USHORT ind = it->index;
				if (ind == size() - 1 || it->relCover < 2 * at(++ind).relCover)
					turnOver = false;
			}
			if (turnOver) {
				REVERSE(itR) = false;
				REVERSE(itL) = true;
				setRealWithCheckBounds(itL, itR, it);
				//REAL(itL) = REAL(itR) = true;
			}
			else {
				//REAL(itL) = true;
				//REAL(itLL) = true;
				auto itLL = prev(itL);
				setRealWithCheckBounds(itLL, itL, it);
			}
		}
		else
			//REAL(itL) = REAL(itR) = true;
			setRealWithCheckBounds(itL, itR, it);
	};

	if (sz > 1) {
		// mark adjasted candidates
		for (auto it0 = begin(), it = next(it0); it != end(); it0++, it++) {
			auto pos = POS(it0->lastIt);	// for debug
			if (it0->lastIt == prev(it->lastIt)) {	// adjacent candidates?
				it0->adjRight = it->adjLeft = true;
			}
		}
		for (auto it = begin(); it != end(); it++) {
			if (it->adjLeft || it->adjRight) {
				auto& itR = it->lastIt;
				auto itL = prev(itR);
				auto len = LEN(itL, itR);
				it->relCover = it->avrCover / len;
			}
		}

		Ñandidates cands = *this;
		// sort in descending avrCover
		sort(cands.begin(), cands.end(),
			//[](const Ñandidate& c1, const Ñandidate& c2) { return c1.avrCover > c2.avrCover; }
			[](const Ñandidate& c1, const Ñandidate& c2) { return c1.relCover > c2.relCover; }
		);

		const fraglen minDistance = Glob::FragLen / 2;
		auto it = cands.begin();
		auto& itR = it->lastIt;
		auto itL = prev(itR);
		Region rgn0{ POS(itL), POS(itR) };
 		setReal(itL, itR, it);

		for (it++; it != cands.end(); it++) {
			auto& itR = it->lastIt;
			auto itL = prev(itR);
			Region rgn{ POS(itL), POS(itR) };

			if (rgn.ToTheLeft(rgn0, minDistance) || rgn.ToTheRight(rgn0, minDistance))
				setReal(itL, itR, it);

			rgn = rgn0;
		}
	}
	else {
		auto it = begin();
		auto& itR = it->lastIt;
		auto itL = prev(itR);
		REAL(itL) = REAL(itR) = true;
		if (REVERSE(itR)) {
			REVERSE(itR) = false;
			REVERSE(itL) = true;
		}
		//setRealWithCheckBounds(itL, itR, ++it);
	}
}

void BS_map::RefineRgn(Ñandidates& cands, citer endIt)
{
	cands.SetReal();
}

void BS_map::Refine()
{
	/*
	BS Left entry (bound): is formed by reverse reads; BS Right entry (bound): is formed by forward reads

	Options for placing bounds in a region:
	canonical:					[[L] [R]]
	adjacent right/left bounds:	[R] [[L] [R]] [L]
	'negative' BS width:			[R] [L]

	Method brings the instance to canonical form, resetting the score of all adjacent elements to zero,
	and marking BSs with 'negative' width.
	*/

	chrlen grpNumb = 1;
	bool firstInRgn = true;
	BYTE lastReverse = 0;
	USHORT ind = 0;
	Ñandidates cands;	cands.reserve(4);

	// *** refine
	for (auto it = begin(); it != end(); it++) {
		if (grpNumb != GrpNUMB(it)) {
			// treat previous region
			auto pos = POS(it);
			RefineRgn(cands, it);
			// reset current region
			firstInRgn = grpNumb = GrpNUMB(it);	// true
			_lastIt = it;
			cands.clear();
			ind = 0;
		}

		if (firstInRgn)	_lastIt = it;
		else if (lastReverse != REVERSE(it)) {
			auto pos = POS(it);
			auto it0 = prev(it);
			float topCover0 = TopCOVER(prev(it));
			float topCover1 = TopCOVER(it);
			//float avrCover = topCover1 > topCover0 ? topCover0 / topCover1 : topCover1 / topCover0;
			float avrCover = (topCover0 + topCover1) / 2;
			cands.emplace_back(ind++, it, avrCover);
		}
		lastReverse = REVERSE(it);
		firstInRgn = false;
	}
	// treat last region
	RefineRgn(cands, end());

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
	* fill the score from left to right (for reverse)
	* or from rigth to left (for direct) extended bounds
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
				//auto len = int(LEN(vp[i], vp[i + 1]));
				//if (len <= 0)
				//	printf("+> %d  len: %d  numb: %d  score: %.3f\n", vp[i]->first, len, GrpNUMB(vp[i]), SCORE(vp[i]));
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
				//auto len = int(LEN(vp[i - 1], vp[i]));
				//if (len <= 0)
				//	printf("-> %d  len: %d  numb: %d  score: %.3f\n", vp[i]->first, len, GrpNUMB(vp[i]), SCORE(vp[i]));
				SCORE(vp[i]) = spline.AvrScoreInRange(offset, LEN(vp[i - 1], vp[i]), -1);
			}
		}
	}
	// BS
	auto& start = VP[L].back()->second;
	auto& end = VP[R].front()->second;
	offset = int32_t(start.RefPos - startPos);
	//auto len = int(end.RefPos - start.RefPos);
	//if (len <= 0)
	//	printf("!> %d  len: %d  numb: %d  score: %.3f %.3f\n", VP[L].back()->first, len, start.GrpNumb, start.Score, end.Score);
	float score = spline.AvrScoreInRange(offset, end.RefPos - start.RefPos);

	if (maxScore < score)	maxScore = score;
	start.Score = end.Score = score;
}

void BS_map::SetGroupScores(iter& start, const iter& end, const Values& spline, float& maxScore)
{
	vector<iter> VP[2];	// 0 - forward, 1 - reversed

	VP[R].reserve(4), VP[L].reserve(4);
	for (auto& it = start; it != end; it++)
		if (IsValid(it)) {
			if (REVERSE(it) && VP[R].size() && VP[L].size()) {
				SetBSscores(VP, spline, maxScore);
				VP[R].clear(), VP[L].clear();
			}
			VP[REVERSE(it)].push_back(it);
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

		if (grpNumb != GrpNUMB(it)) {
			// build single fragment coverage spline for the whole group
			cover.SetLocalSpline(spliner, POS(itStart), POS(itEnd), spline);
			SetGroupScores(itStart, ++itEnd, spline, maxScore);

			spline.clear();
			itStart = end();
			grpNumb = GrpNUMB(it);
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
	chrlen	minScoreNumb;
	float	minScore = 1.f;

	DoBasic([&](citer& start, citer& end) {
		++bsNumb;
		if (stat) {
			float score = SCORE(start);
			if (minScore > score)		minScore = score, minScoreNumb = bsNumb;
		}
	});
	if (stat)
		printf("min score: %2.2f (%d)\n\n", minScore, minScoreNumb);

	printf("BS count: %d\n\n", bsNumb);
}

#ifdef MY_DEBUG
void BS_map::PrintDistrib(const string& fName, const char* title, function<USHORT(const citer&, const citer&, float&)> func) const
{
	map<USHORT, vector<chrlen>> freq;
	float	total = 0;
	chrlen	cnt = 0;

	// collect distribution
	DoExtend([&](const vector<citer>* VP) {
		auto val = func(VP[L].back(), VP[R].front(), total);
		freq[val].push_back(++cnt);
		}
	);

	// print distribution
	FormWriter file((fName + distExt).c_str());
	string stitle(title);

	transform(stitle.begin(), stitle.end(), stitle.begin(),	::toupper);
	file.Write("BS %s FREQUENCY:\n", stitle.c_str());
	file.Write("%s\tcnt\tnumbers\n", title);
	for (const auto& item : freq) {
		file.Write("%4d\t%u\t", item.first, UINT(item.second.size()));
		auto it = item.second.begin();
		file.Write("%d", *it);
		for (it++; it != item.second.end(); it++)
			file.Write(",%d", *it);
		file.Write("\n");
	}
	auto avr = total / cnt;
	file.Write("average %s: %.2f\n", title, avr);
	printf("average %s: %.2f\n", title, avr);
}

void BS_map::PrintWidthDistrib(const string& fName) const
{
	PrintDistrib(fName, "width", 
		[](const citer& start, const citer& end, float& total) {
			auto len = LEN(start, end);
			total += len;
			return len;
		}
	);
}

void BS_map::PrintScoreDistrib(const string& fName) const
{
	PrintDistrib(fName, "score",
		[](const citer& start, const citer&, float& total) {
			auto fscore = SCORE(start);
			total += fscore;
			return USHORT(round(fscore * 1000)) / 10;
		}
	);
}

void BS_map::Print(chrid cID, const string& fName, bool selected, chrlen stopPos) const
{
	string format = "%9d % 5d  %c %8d %5.2f %6.1f%5d%4d\t%s\n";
	const char bound[]{ 'R','L' };
	IGVlocus locus(cID);

	FormWriter file(fName.c_str());
	file.Write("  pos     numb bnd  ref_pos score topCvr real bun\tlocus\n");
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
			x.second.TopCover,
			x.second.Real,
			x.second.Related,
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
		if (!rgn.Accepted()) {		// discarded items
			LineAddChars("\t.\t.\t", 5, false);
			LineAddInts(rgn.Start(), rgn.End(), true);
			LineAddChars(sGRAY, colorLen, false);
		}
		LineToIOBuff(offset);
	}
}

void BedWriter::WriteChromData(chrid cID, const BS_map& bss)
{
	const reclen offset = AddChromToLine(cID);
	bool lastSep[]{ false, false };
	chrlen bsNumb = 0;

	bss.DoExtend([&](const vector<BS_map::citer>* VP) {
		// *** save basic info
		auto& start = VP[L].back();
		auto& end	= VP[R].front();

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
				LineAddArgs("%c%d", specChar, vp[i]->first);
			};
			function<void(BYTE, char)> saveExtraVal = [this, &vp](BYTE i, char specChar) {
				if (specChar)	LineAddChar(specChar);
				LineAddScore(SCORE(vp[i]), false);
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

void BedWriter::WriteChromExtData(chrid cID, const BS_map& bss)
{
	chrlen bsNumb = 0;
	const reclen offset = AddChromToLine(cID);

	bss.DoExtend([&](const vector<BS_map::citer>* VP) {
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
		const auto& start = VP[L].back();	// basic feature start
		const float score = SCORE(start);	// basic feature score

		auto addExtraLines = [=, &bsNumb](const vector<BS_map::citer>& vp) {	// &bsNumb is essential, otherwise bsNumb goes out of sync
			if (vp.size() == 1)	return;
			for (auto it0 = vp.begin(), it = next(it0); it != vp.end(); it0++, it++) {
				LineAddUInts(POS(*it0), POS(*it), bsNumb, true);
				LineAddFloat(SCORE(*it0), true);
				LineAddChar(DOT, true);
				LineAddInts(POS(*it0), POS(*it), true);
				// colors
				auto ind = BYTE(10 * SCORE(*it0) / score) / 2;
				if (ind > COLORS_CNT - 1)	ind = COLORS_CNT - 1;
				LineAddStr(colors[ind], false);
				LineToIOBuff(offset);
			}
		};

		const auto& end = VP[R].front();
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
		}
	);
}

void BedWriter::WriteChromROI(chrid cID, const BS_map& bss)
{
	const fraglen ROI_ext = 500;
	const reclen offset = AddChromToLine(cID);
	chrlen bsNumb = 0;

	bss.DoBasic([&](BS_map::citer& start, BS_map::citer& end) {
		LineAddUInts(POS(start) - ROI_ext, POS(end) + ROI_ext, ++bsNumb, false);
		LineToIOBuff(offset);
		}
	);
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
