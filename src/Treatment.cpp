#include "Treatment.h"
#include <algorithm>

//const float PI = 3.14159265F;

const eCurveType CurveTYPE = SPIKED;
const char* Verb::ValTitles[] = { "SL","RES","RT","DBG" };
const char* Verb::ValDescr = "set verbose level:\n?  -\tsilent mode (show critical messages only)\n? -\tshow result summary\n?  -\tshow run-time information\n? -\tshow debug messages";
Verb::eVerb Verb::_level;

readlen Glob::ReadLen = 0;
fraglen Glob::FragLen = FragDefLEN;
fraglen Glob::ROI_ext = 500;

//===== IGVlocus

#ifdef MY_DEBUG
class IGVlocus
{
	const string _chrom;
	mutable char _buf[2 * 10 + 5 + 2];	// 2 * max position length + Chrom::MaxAbbrNameLength + 2 separators

public:
	IGVlocus(chrid cID) : _chrom(Chrom::AbbrName(cID)) {}

	// Returns inner buffer
	const char* Buff() const { return _buf; }

	// Prints IGV locus to inner buffer
	//	@return: the total number of characters written
	chrlen NPrint(chrlen start, chrlen end) const
	{
		return sprintf(_buf, "%s:%d-%d", _chrom.c_str(), start - Glob::ROI_ext, end + Glob::ROI_ext);
	}

	// Prints IGV locus to inner buffer
	//	@return: inner buffer
	const char* Print(chrlen start, chrlen end) const { NPrint(start, end); return _buf; }

	const char* Print(chrlen pos) const { return Print(pos, pos); }
};
#endif


//===== OLinearWriter

void OLinearWriter::WriteOblique(BYTE reverse, chrlen start, coviter& itStop)
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
	(_writers->_files)[reverse]->WriteOblique(_cID, start, k * ptCnt, float(itStop->second) / ptCnt);
}


//===== TreatedCover

#ifdef MY_DEBUG
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

void TreatedCover::LinearRegr(coviter it, const coviter& itStop, const StrandOp& op, LRegrResult& result) const
{
	const int shift = op.Factor * itStop->first;
	float sumX = 0, sumX2 = 0, sumY = 0, sumXY = 0;

	result.Clear();
	for (; it != itStop; op.Next(it), result.PtCnt++) {
		const int x = shift - op.Factor * it->first;
		const coval y = op.GetPrev(it)->second;
		sumX += x;
		sumX2 += x * x;
		sumY += y;
		sumXY += x * y;
	}
	if (result.PtCnt < 2)	return;

	result.Deriv = 
		-(result.PtCnt * sumXY - sumX * sumY) /		// numerator
		(result.PtCnt * sumX2 - sumX * sumX);		// denominator
	//angl = coeff * 180 / PI;	if (angl < 0) angl = -angl;

	float x = (sumY + result.Deriv * sumX) / (result.PtCnt * result.Deriv);
	result.Pos = itStop->first - op.Factor * chrlen(round(x));

#ifdef MY_DEBUG
	if (_fDdelim) {
		IGVlocus locus(0);
		_fDdelim->Write("%1.3f\t%d\t%s\n", result.Deriv, result.Pos, locus.Print(result.Pos));
	}
#endif
}


//===== CombCover

void CombCover::AddExtRead(const Region& read, bool reverse, bool addTotal)
{
	Region frag(read, Glob::FragLen, reverse);

	if(addTotal)
		_data->DataByInd().AddRegion(frag);			// total frag coverage
	AddRead(frag, reverse);
}

void CombCover::Fill(const Reads& reads, bool addTotal)
{
	for (BYTE s : {0, 1})
		for (auto& rd : reads.GetReads(s))
			AddExtRead(rd, s, addTotal);
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
	const typename T::iterator itEnd0 = rgns[0].end();
	const typename T::iterator itEnd1 = rgns[1].end();
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
			it[s] = next(it0[s]);
			return true;
		}
		return false;
	};

	while (it[0] != itEnd0 && it[1] != itEnd1) {
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
void PrintRegionStats(const T* rgns, chrlen chrLen, bool isSE = true)
{
	const char* format[] = {
		"%s:  %d (%2.2f%%)   %d (%2.2f%%)\n",	// features
		"%s:  %d (%1.3f%%)   %d (%1.3f%%)\n"	// splines
	};
	const char* title[] = {
		"POTENTIAL REGIONS:",
		"SPLINED REGIONS:"
	};
	const bool isFeatures = is_same<T, ValuesMap>::value;

	printf("\n%s\nstrand     received        selected\n", title[isFeatures]);
	printf("------------------------------------\n");
	for (BYTE s = 0; s < 1 + isSE; s++) {
		chrlen rawLen = 0, refineLen = 0;
		const auto& rgn = rgns[s];
		const auto itEnd = rgn.end();
		size_t realCnt = 0;

		for (auto it = rgn.begin(); it != itEnd; it++) {
			const chrlen len = T::Length(it);
			rawLen += len;
			if (T::IsNotEmpty(it)) refineLen += len, realCnt++;
		}
		std::printf(format[isFeatures], sStrandTITLES[s + isSE],
			rgn.size(), Percent(rawLen, chrLen),
			realCnt, Percent(refineLen, chrLen));
	}
}

//===== CoverRegions

void CoverRegions::SetPotentialRegions(const TreatedCover& cover, size_t capacity, coval cutoff)
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

//===== DataCoverRegions

bool PositiveVal(int16_t val) { return val > 0; }
bool NegativeVal(int16_t val) { return val < 0; }

void DataCoverRegions::SetPotentialRegions(const DataSet<TreatedCover>& cover, chrlen cLen, coval cutoff, bool noMultiOverl)
{
	size_t capacity = cLen / (Glob::FragLen * 100);
	StrandData(POS).SetPotentialRegions(cover.StrandData(POS), capacity, cutoff);
	StrandData(NEG).SetPotentialRegions(cover.StrandData(NEG), capacity, cutoff);

	EliminateNonOverlapsRegions<CoverRegions>(Data(), Glob::FragLen);
	if (noMultiOverl)
		EliminateMultiOverlapsRegions<CoverRegions>(Data());
	if (Verb::Level(Verb::DBG))
		PrintRegionStats<CoverRegions>(Data(), cLen);
}

void DataCoverRegions::SetPotentialRegions(const DataSet<TreatedCover>& cover, chrlen cLen, coval cutoff)
{
	size_t capacity = cLen / (Glob::FragLen * 100);
	DataByInd().SetPotentialRegions(cover.DataByInd(), capacity, cutoff);

	if (Verb::Level(Verb::DBG))
		PrintRegionStats<CoverRegions>(Data(), cLen, false);
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

void DataCoverRegions::Clear()
{
	StrandData(POS).clear();
	StrandData(NEG).clear();
}


//===== BS_Map
void BS_Map::Refine()
{
	bool openedGroup = true;
	BYTE posCnt = 0, negCnt = 0;
	chrlen grpNumb = 1;

	// 'close' unrelated positions
	for (auto it = begin(); it != end(); it++)
	{
		auto MarkSubseqAsEmpty = [&it](BYTE& entryCnt) {	// resets entryCnt
			for (auto it1 = it; entryCnt; entryCnt--)
				(--it1)->second.Score = 0;
		};
		const bool newGroup = grpNumb != it->second.GrpNumb;

		if (newGroup) {
			if (!openedGroup)
				openedGroup = !(posCnt = 0);	// set to true
			// mark last reverse positions in the previous group as empty
			MarkSubseqAsEmpty(negCnt);
			grpNumb = it->second.GrpNumb;
		}
		if (it->second.Reverse) {
			// mark first direct positions in the current group as empty
			if (openedGroup) {
				if (!newGroup)
					MarkSubseqAsEmpty(posCnt);
				openedGroup = false;
			}
			posCnt = 0;
			negCnt++;
		}
		else {
			negCnt = 0;
			posCnt++;
		}
	}
}

void BS_Map::Normalize()
{
	float maxScore = 0;

	// *** set direct&reversed position scores as sum of each other
	Do([&maxScore](vector<ValPos>* VP) {
		const float val[2]{ VP[0].front().Val ,VP[1].back().Val };	// first direct, last reversed
		const float ratio = val[0] / val[1];						// direct/reversed score ratio

		auto setScore = [&maxScore](const vector<ValPos>& vp, BYTE i, float val) {
			float score = (vp[i].Iter->second.Score += val);
			if (maxScore < score)	maxScore = score;
		};

		setScore(VP[1], BYTE(VP[1].size() - 1), val[0]);	// set basic reverse score to sum of direct & reversed
		VP[0].front().Iter->second.Score = ratio;	// set basic direct score to direct/reversed score ratio

		for (BYTE s : {0, 1}) {
			const BYTE vpLen = BYTE(VP[s].size() - 1);
			const auto& vp = VP[s];
			for (BYTE i = !s + 1; i < vpLen; i++)
				setScore(vp, i, val[!s]);
		}
		});

	// *** normalize scores relative to maximum
	bool reversed = false;
	for (auto& x : *this) {
		if (!x.second.Score)	continue;
		if(x.second.Reverse)
			reversed = true;
		else if (reversed) {
			reversed = false;
			continue;		// skip basic direct score as it keeps direct/reversed ratio
		}
		x.second.Score /= maxScore;
	}
}

void BS_Map::PrintStat() const
{
	if (!Verb::Level(Verb::DBG))	return;

	chrlen	bsNumb = 0;
	chrlen	minNumb, maxNumb, minLen = CHRLEN_MAX, maxLen = 0;
	chrlen	minNegNumb, maxNegNumb, minPosNumb, maxPosNumb, minScoreNumb;
	float	minNegRatio, minPosRatio, minScore, maxNegRatio = 0, maxPosRatio = 0;
	minNegRatio = minPosRatio = minScore = 1000;

	DoBasic([&](citer start, citer end) {
		const chrlen len = end->first - start->first;
		++bsNumb;
		if (minLen > len)		minLen = len, minNumb = bsNumb;
		else if (maxLen < len)	maxLen = len, maxNumb = bsNumb;

		float val = start->second.Score;
		if (minScore > val)		minScore = val, minScoreNumb = bsNumb;

		val = end->second.Score;
		if (val < 1) {
			if (minNegRatio > val)		minNegRatio = val, minNegNumb = bsNumb;
			else if (maxNegRatio < val)	maxNegRatio = val, maxNegNumb = bsNumb;
		}
		else
			if (minPosRatio > val)		minPosRatio = val, minPosNumb = bsNumb;
			else if (maxPosRatio < val)	maxPosRatio = val, maxPosNumb = bsNumb;
		});

	printf("\nBS count: %d\n", bsNumb);
	if (!bsNumb) return;
	printf("shortest: %d (%d),  longest: %d (%d)\n", minNumb, minLen, maxNumb, maxLen);
	printf("min score: %2.2f (%d)\n", minScore, minScoreNumb);
	printf("----------------------------\n");
	printf("RATIO:\tmin  (N)   max  (N)\n");
	printf("----------------------------\n");
	printf("reverse\t%2.2f (%d)  %2.2f (%d)\n", 1 / maxNegRatio, minNegNumb, 1 / minNegRatio, maxNegNumb);
	printf("direct\t%2.2f (%d)  %2.2f (%d)\n", minPosRatio, minPosNumb, maxPosRatio, maxPosNumb);
}

#ifdef MY_DEBUG
void BS_Map::CheckScoreHierarchy()
{
	chrlen bsNumb = 0;
	bool issues = false;

	printf("\nCHECK SCORES HIERARCHY\n");
	Do([&](const vector<ValPos>* VP) {
		bool prNeg = false, prPos = false;
		float val = 0;
		auto prVals = [](bool sep, char sign, const vector<ValPos>& VP) {
			if(sep) 	printf(" ");
			printf("%c", sign);
			for (const ValPos& vp : VP)
				printf("%2.3f ", vp.Val);
		};
		
		bsNumb++;
		// set prNeg; check all
		for (const ValPos& vp : VP[1])
			if (prNeg = (val > vp.Val))		break;
			else val = vp.Val;
		// set prPos; check all except first element
		const BYTE vpLen = BYTE(VP[0].size() - 1);
		const auto& vp = VP[0];
		val = 1000;
		for (BYTE i = 1; i < vpLen; i++)
			if (prPos = (val < vp[i].Val))	break;
			else val = vp[i].Val;
		// print of excess over basic score
		if (prNeg || prPos) {
			issues = printf("%4d  %d\t", bsNumb, VP[1].back().Iter->first);	// start position
			if (prNeg)	prVals(false, '-', VP[1]);
			if (prPos)	prVals(prNeg, '+', VP[0]);
			printf("\n");
		}
		});

	if (!issues)	printf("OK\n");
}

void BS_Map::Print(chrid cID, bool real, chrlen stopPos) const
{
	const char* format[]{
		"%d %4d  %d  %1.3f %4d   %s\n",	// odd group numb to left
		"%d %-4d  %d  %1.3f %4d   %s\n"	// even group numb to right
		//"%d %4d  %d  %1.3f  %s\n",	// odd group numb to left
		//"%d %-4d  %d  %1.3f  %s\n"	// even group numb to right
	};
	IGVlocus locus(cID);

	printf("\npos\tgrp   R  score  pCnt  IGV view\n");
	//printf("\npos\tgrp   R  score  IGV view\n");
	for (const auto& x : *this) {
		if (stopPos && x.first > stopPos)	break;
		if (!real || x.second.Score)
			printf(format[x.second.GrpNumb % 2],
				x.first, x.second.GrpNumb, int(x.second.Reverse),
				x.second.Score, x.second.PointCount, locus.Print(x.first));
	}
}
#endif


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

void BedWriter::WriteChromData(chrid cID, BS_Map& bss)
{
	const reclen offset = AddChromToLine(cID);
	bool lastSep[]{ false, false };
	chrlen bsNumb = 0;

	bss.Do([&](const vector<BS_Map::ValPos>* VP) {
		// *** save basic info
		const auto& itStart = VP[1].back().Iter;
		const auto& itEnd = VP[0].front().Iter;

		LineAddUInts(itStart->first, itEnd->first, ++bsNumb, true);	// 3 basic fields
		LineAddScore(itStart->second.Score, true);					// BS score
		
		// *** additional regions info
		LineAddChar(DOT, true);
		lastSep[1] = VP[0].size() - 1;
		LineAddFloat(itEnd->second.Score, VP[1].size() - 1 || lastSep[1]);	// ratio

		// *** extra deviations info
		for (BYTE s : {1, 0}) {
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

void BedWriter::WriteChromExtData(chrid cID, BS_Map& bss)
{
	chrlen bsNumb = 0;
	const reclen offset = AddChromToLine(cID);

	bss.Do([&](const vector<BS_Map::ValPos>* VP) {
		static const string colors[]{
			// color		 ind	feature_score/BS_score
			"200,200,200",	// 0	>=0
			"160,160,160",	// 1	>=0.2
			"120,120,120",	// 2	>=0.4
			"80,80,80",		// 3	>=0.6
			"40,40,40",		// 4	>=0.8
			"Black",		// 5	>=1
		};
		const auto& start = VP[1].back().Iter;
		const float score = VP[1].back().Val;

		auto addExtraLines = [=, &bsNumb](const vector<BS_Map::ValPos>& vp) {	// &bsNumb is essential, otherwise bsNumb goes out of sync
			if (vp.size() == 1)	return;
			for (auto it0 = vp.begin(), it = next(it0); it != vp.end(); it0++, it++) {
				LineAddUInts(it0->Iter->first, it->Iter->first, bsNumb, true);
				LineAddFloat(start->second.Score, true);	// BS score
				LineAddChar(DOT, true);
				LineAddInts(it0->Iter->first, it->Iter->first, true);
				int ind = int(10 * it0->Val / score) / 2;
				if (ind > 5)	ind = 5;
				LineAddStr(colors[ind], false);
				LineToIOBuff(offset);
			}
		};

		const auto& end = VP[0].front().Iter;
		const char* delims = ".\t.\t.\t.\t";
		const char* color = NULL;
		const char* darkRed = "150,15,15";

		++bsNumb;

		addExtraLines(VP[1]);
		// *** add basic feature
		LineAddUInts(start->first, end->first, bsNumb, true);
		LineAddFloat(start->second.Score, true);		// BS score
		if (end->second.Score > 2.0)		color = darkRed;
		else if (end->second.Score < 0.5)	color = "15,15,150";
		if (color) {
			LineAddChar(DOT, true);
			LineAddInts(start->first, end->first, true);
			LineAddChars(color, reclen(strlen(darkRed)), true);
		}
		else
			LineAddChars(delims, reclen(strlen(delims)), false);
		LineAddFloat(end->second.Score, false);			// reverse/direct ratio
		LineToIOBuff(offset);

		addExtraLines(VP[0]);
		});
}

void BedWriter::WriteChromROI(chrid cID, const BS_Map& bss)
{
	const reclen offset = AddChromToLine(cID);
	chrlen bsNumb = 0;

	bss.DoBasic([&](BS_Map::citer start, BS_Map::citer end) {
		LineAddUInts(
			start->first - Glob::ROI_ext,
			end->first + Glob::ROI_ext,
			++bsNumb,
			false
		);
		LineToIOBuff(offset);
		});
}


//===== BedWriters

//===== Values

void Values::AddVal(float val)
{
	if (_maxVal < val)	_maxVal = val;
	push_back(val);
}

void Values::GetMaxValPos(chrlen startPos, vector<chrlen>& pos) const
{
	float val0 = front();
	chrlen currPos = startPos;
	bool increase = false;

	for (auto it = next(begin()); it != end(); val0 = *it, it++, currPos++) {
		if (*it > val0) {		// increase
			increase = true;
		}
		else {
			if (increase) {
				pos.push_back(currPos);
				increase = false;
			}
		}
	}
}

#ifdef MY_DEBUG
void Values::Print() const
{
	printf("MaxVal: %2.2f\n", _maxVal);
	printf("Values: ");
	if (size())
		for (auto v : *this)	printf("%2.2f ", v);
	else printf("none");
	printf("\n");
}
#endif


//===== ValuesMap

#ifdef MY_DEBUG
void ValuesMap::Print(chrid cID, BYTE reverse, chrlen stopNumb) const
{
	IGVlocus locus(cID);

	printf("SPLINES %s\n", sStrandTITLES[reverse]);
	printf(" N start\tend\tval\tIGV view\n");
	for (const auto& x : *this) {
		if (stopNumb && x.second.GrpNumb > stopNumb)	break;
		if (x.second.MaxVal()) {
			chrlen end = x.first + x.second.Length();
			printf("%2d %d\t%d\t%2.2f\t%s\n",
				x.second.GrpNumb, x.first, end, x.second.MaxVal(), locus.Print(x.first, end));
		}
	}
}
#endif

void ValuesMap::BuildRegionSpline(const TreatedCover& cover, const CoverRegion& rgn, bool redifRgns, fraglen splineBase)
{
	assert(Glob::ReadLen);
	coviter it0, itEnd;

	// set up start/end iterators
	if (redifRgns) {
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

	//for (auto it = next(it0); it0 != itEnd; it0++, it++) {	// loop through the cover
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
				vals.AddVal(val);
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

		it[0]->second.GrpNumb = it[1]->second.GrpNumb = numb;

		for(bool overlap = true; overlap; )
			if (NextStart(0) + overLen < End(it[1]))
				overlap = (++it[0])->second.GrpNumb = numb;
			else if (NextStart(1) + overLen < End(it[0]))
				overlap = (++it[1])->second.GrpNumb = numb;
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

void DataValuesMap::BuildSplineSE(
	const DataSet<TreatedCover>& cover, const DataCoverRegions& rgns, bool redifineRgns, fraglen splineBase)
{
	StrandData(POS).BuildSpline(cover.StrandData(POS), rgns.StrandData(POS), redifineRgns, splineBase);
	StrandData(NEG).BuildSpline(cover.StrandData(NEG), rgns.StrandData(NEG), redifineRgns, splineBase);
}

void DataValuesMap::BuildSplinePE(
	const DataSet<TreatedCover>& cover, const DataCoverRegions& rgns, bool redifineRgns, fraglen splineBase)
{
	StrandData(POS).BuildSpline(cover.StrandData(POS), rgns.DataByInd(), redifineRgns, splineBase);
	StrandData(NEG).BuildSpline(cover.StrandData(NEG), rgns.DataByInd(), redifineRgns, splineBase);
}

fraglen DataValuesMap::GetFragMean() const
{
	auto& pData = StrandData(POS);
	auto& nData = StrandData(NEG);
	assert(pData.size() == nData.size());
	uint32_t missed = 0, negCnt = 0;
	vector<int16_t> diffs;
	vector<chrlen> pPos, nPos;	// number of max positions in one spline
	//IGVlocus locus(0);

	// get the difference of the splines maximums 
	pPos.reserve(2);
	nPos.reserve(2);
	for (auto itP = pData.begin(), itN = nData.begin(); itP != pData.end() && itN != nData.end(); itP++, itN++) {

		itP->second.GetMaxValPos(itP->first, pPos);
		itN->second.GetMaxValPos(itN->first, nPos);
		if (pPos.size() == nPos.size())
			for (auto itP = pPos.begin(), itN = nPos.begin(); itP != pPos.end(); itP++, itN++) {
				int16_t diff = *itP - *itN;
				negCnt += diff < 0;
				diffs.push_back(diff);
				//printf("%d\t%d\t%s\n", diff, pPos.size(), locus.Print(*itP));
			}
		else
			missed++;		// just ignore splines with different max positions
		pPos.clear();
		nPos.clear();
	}
	if (Verb::Level(Verb::DBG) && missed)
		printf("\n%2.1f%% regions were rejected while determining the fragment length\n", Percent(missed, pData.size()));

	// average the difference, get frag length
	bool mostPositive = negCnt < diffs.size() / 2;
	auto CompareVal = mostPositive ? &PositiveVal : &NegativeVal;
	int sum = 0;

	for (auto diff : diffs)
		if (CompareVal(diff))
			sum += diff;

	return FragDefLEN - fraglen(round(float(sum) / (mostPositive ? (diffs.size() - negCnt) : negCnt)));
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

//===== RegionsValues

void RegionsValues::AddBSpos(BYTE reverse, const PosVal& posVal, const TreatedCover& cover, BS_Map& bs, OLinearWriter& lwriter) const
{
	//chrlen pos[2];	// 0 - start, 1 - end
	//pos[!reverse] = it->Start();	// right for direct, left for reverse
	//pos[reverse]  = it->End();		// left for direct, right for reverse

	const StrandOp& op = StrandOps[reverse];
	const StrandOp& opInv = StrandOps[!reverse];	// inversed operations
	auto itStart = op.GetPrev(cover.upper_bound(posVal.Pos(0)));	// min val: right for direct, left for reverse
	auto itStop = itStart;									// max val: left for direct, right for reverse

	// *** set itStop
	if (reverse)	for (itStop++; itStop->first < posVal.Pos(1); itStop++);
	else			for (itStop--; itStop->first > posVal.Pos(1); itStop--);

	// checking if there is an ordinary tag pool before the current start, masked by a spline
	if (!op.GetPrev(itStart)->second)
		op.Next(itStart);		// skip the gap after ordinary tag pool

	// *** it stop correction: calc the best BS position by successive linear regression approximations
	LRegrResult result;
	cover.LinearRegr(itStart, itStop, op, result);
	if (!result.Valid())	return;

	coviter it0 = itStop;
	LRegrResult result1{};
	auto compare = [&](coviter& it) {
		cover.LinearRegr(itStart, it, op, result1);
		if (!result1.Valid() || op.EqLess(result.Pos, result1.Pos))	return true;
		std::swap(result, result1);
		itStop = it;
		return false;
	};

	// iterate itStop opposed to start
	for (auto it = it0; op.RetNext(it)->second > itStop->second || op.RetNext(it)->second > itStop->second;)
		if (compare(it))	break;
	// iterate itStop towards start
	bool nextStep = false;
	//auto it1 = opInv.RetNext(it0);
	//auto it2 = it1;
	//if(it1 != cover.end())
	//	it2 = opInv.RetNext(it1);
	for (auto it = it0;
		opInv.RetNext(it)->second > itStop->second 
		|| (nextStep = opInv.RetNext(it)->second >= itStop->second);) 
	{
		if (reverse && nextStep) { nextStep = false; it++; }
		if (compare(it))	break;
	}

	// *** add BS
	if (lwriter.IsWriterSet()) {
		opInv.Next(itStop);
		lwriter.WriteOblique(reverse, result.Pos, itStop);
	}
	bs.AddPos(reverse, _grpNumb, result, posVal.Val(reverse));
}

void RegionsValues::AddRegion(const tValuesMap::value_type& spline, chrlen relPos, Values& deriv)
{
	if (_maxVal < deriv.MaxVal())	_maxVal = deriv.MaxVal();
	_grpNumb = spline.second.GrpNumb;
	emplace_back(
		spline.first + relPos, 
		spline.second.Val(relPos),
		spline.second.Val(relPos + deriv.Length()),
		deriv
	);
}

void RegionsValues::SetDirectBSpos(const TreatedCover& cover, BS_Map& bs, OLinearWriter& lwriter) const
{
	PosVal posVal;

	for (auto it = rbegin(); it != rend(); it++) {		// loop through derivatives
		posVal.Set(0, it);
		AddBSpos(0, posVal, cover, bs, lwriter);
		posVal.Swap();
	}
}

void RegionsValues::SetReverseBSpos(const TreatedCover& cover, BS_Map& bs, OLinearWriter& lwriter) const
{
	PosVal posVal;

	for (auto it = begin(); it != end(); it++) {		// loop through derivatives
		posVal.Set(1, it);
		AddBSpos(1, posVal, cover, bs, lwriter);
		posVal.Swap();
	}
}


//===== RegionsValuesMap

void RegionsValuesMap::BuildDerivs(int factor, const ValuesMap& splines)
{
	for (const auto& spline : splines) {
		if (!spline.second.MaxVal())	continue;

		Values deriv;
		RegionsValues derivSet;	derivSet.reserve(4);
		fraglen pos = 0;

		auto addRngsVals = [&]() {
			if (!deriv.MaxVal())	return;
			fraglen len = deriv.Length();	// because deriv will be cleared

			derivSet.AddRegion(spline, pos, deriv);
			pos += len;
			deriv.Reserve();
		};
		const auto& itEnd = spline.second.end();

		for (auto it0 = spline.second.cbegin(), it = next(it0); 
			it != itEnd; it0++, it++) 
		{
			float tang = (*it - *it0) * factor;				// flip sign of pos strand delta

			if (tang > 0)		// cut off negative tangent
				deriv.AddVal(tang);
			else {
				pos++;
				addRngsVals();
			}
		}
		// last region
		addRngsVals();
		this->AddRegions(spline.first, derivSet);
	}
}

void RegionsValuesMap::SetBSpos(BYTE reverse, const TreatedCover& cover, BS_Map& bs, OLinearWriter& lwriter) const
{
	using tSetBSpos = void(RegionsValues::*)(const TreatedCover&, BS_Map&, OLinearWriter&) const;
	tSetBSpos SetBSpos = reverse ?
		&RegionsValues::SetReverseBSpos	:
		&RegionsValues::SetDirectBSpos	;

	for (auto it = begin(); it != end(); it++) {		// loop through the derivative groups
		(it->second.*SetBSpos)(cover, bs, lwriter);
	}
}

#ifdef MY_DEBUG
void RegionsValuesMap::Print(eStrand strand, chrlen stopPos) const
{
	//printf("\nDERIVS %s: %2.2f\n", sStrandTITLES[strand - 1], MaxVal);
	printf("\nDERIVS %s\n", sStrandTITLES[strand - 1]);
	for (const auto& rvss : *this) {
		if (stopPos && rvss.first > stopPos)	break;
		printf("%d: %2.2f\n", rvss.first, rvss.second.MaxVal());
		for (const auto& rvs : rvss.second)
			printf("  %2.2f:\t%d %d\t%d\n", rvs.MaxVal(), rvs.GrpNumb, rvs.Start(), rvs.Length());
	}
}
#endif


//===== FixWigWriter

void FixWigWriter::WriteChromData(chrid cID, const ValuesMap& vals)
{
	for (const auto& val : vals)
		if (val.second.MaxVal())
			WriteFixStepRange(cID, val.first, val.second);
}


//===== FixWigWriterSet

void FixWigWriterSet::WriteChromData(chrid cID, const RegionsValuesMap& set)
{
	for (const auto& rvss : set)
		for (const auto& rvs : rvss.second)
			if (rvs.MaxVal())
				WriteFixStepRange(cID, rvs.Start(), rvs);
}
