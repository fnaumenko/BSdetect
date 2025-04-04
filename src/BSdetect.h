/**********************************************************
callDist.h (c) 2021 Fedor Naumenko (fedor.naumenko@gmail.com)
All rights reserved.
-------------------------
Last modified: 04/04/2025
-------------------------
Provides main functionality
***********************************************************/

#pragma once
#include "Treatment.h"

//#define TIMING

enum optValue {		// options id
	oBIN,
	oGEN,
	oCHROM,
	oDUP_LVL,
	oFRAG_LEN,
	oSAVE_COVER,
	oSAVE_INTER,
	oSAVE_SPLINE,
	oALARM,
	oREAD_LEN,
	oSINGLE_CHROM,
	oRANK_SCORE,
	oOUTFILE,
	oTIME,
	oVERB,
	oVERSION,
	oHELP,
	oHHELP,
};

const char* FragCoverDescr	= "fragment coverage";
const char* ReadCoverDescr	= "read coverage";
const char* FragSplineDescr = "fragment coverage spline";
const char* ReadSplineDescr = "read coverage spline";
const char* RegionsDescr	= "potential regions";
const char* DerivDescr		= "derivative of read coverage spline";
const char* RegressionDescr = "linear regression";
const char* BS_Descr		= "called binding sites";

const string FragSplineExt	= ".FR_SPLINE";
const string ReadSplineExt	= ".SPLINE";
const string RegionsExt		= ".RGNS";
const string DerivExt		= ".DERIV";
const string RegressionExt	= ".LINE";
const string BS_Ext			= ".BSs";

// Class for constructing a spline of an input Bedgraph file
class Spliner
{
	ChromSizes&		_cSizes;
	CombCover		_fragCovers;			// extended reads cover to find frag Mean
	OCoverRegions	_regions;
	OValuesMap		_splines;

	void BuildSpline(chrid cID)
	{
		DataSet<TreatedCover>& fragCovers = _fragCovers.ChromData(cID);
		DataCoverRegions& regions = static_cast<DataCoverRegions&>(_regions.ChromData(cID));

		if (regions.SetPotentialRegions(fragCovers, _cSizes[cID], 3))
			return;

		(static_cast<DataValuesMap&>(_splines.ChromData(cID))).
			BuildSpline(fragCovers, regions, 30);
		_splines.WriteChrom(cID);
	}

public:
	Spliner(
		const char* inFName,
		const string& outFName,
		ChromSizes& cSizes,
		bool saveInter
	)
		: _cSizes(cSizes)
		, _fragCovers(cSizes, 1, false, strEmpty, NULL)
		, _regions(cSizes, 1, saveInter, outFName + RegionsExt, RegionsDescr)
		, _splines(cSizes, 1, true, outFName + ReadSplineExt, FragSplineDescr)
	{
		// *** preparing coverage data
		{
			tChromsOccurrs	chrReadOccurs;
			CombCoverReader cvr(inFName, cSizes, _fragCovers, chrReadOccurs, TOTAL);
			_cSizes.SetTreated(chrReadOccurs, 1);
		}
		// *** treatment
		for (const auto& c : _cSizes)
			if (c.second.Treated)
				BuildSpline(c.first);
	}
};

// BS detector
class Detector
{
	const string FNameFragExt = "_frag";
	const string FNameReadExt = "_read";

	bool			 _saveCover;
	bool			 _unsortNotSet = true;	// false if unsorting is detected while input reading
	ChromSizes&		 _cSizes;
	CombCover		 _fragCovers;			// extended reads cover to find frag Mean
	CombCover		 _readCovers;
	OCoverRegions	 _regions;
	OValuesMap		 _splines;
	OBoundsValuesMap _derivs;
#ifdef MY_DEBUG
	OSpecialWriter	 _lineWriter;
	OSpecialWriter	 _splineWriter;
	string			 _outFName;
#endif
	OBS_Map			 _bss;

	RBedReader* _file;		// needs only for input reading
	FragIdent	_fIdent;	// needs only for input reading
	Reads		_reads;		// may be filled for the first chromosome only, if fragment len is not defined
	Timer		_timer;

	// Calculates the deviation from the default average fragment length
	//	@param cID: current chromosome's ID
	//	@returns: deviation from the default average fragment length
	float GetPeakPosDiff(chrid cID);

	void CallBS(chrid cID);

public:
	// Basic constructor
	//	@param file: input BAM/BED file
	//	@param outFName: common output file name
	//	@param cSizes: chrom sizes
	//	@param saveCover: if true then save covered fragments and reads to files
	//	@param saveInter: if true then save intermediate data to files
	Detector(RBedReader& file, const string& outFName, ChromSizes& cSizes, bool saveCover, bool saveInter)
		: _cSizes(cSizes)
		, _saveCover(saveCover)
		, _fragCovers(cSizes, 3-2*Glob::IsPE, saveCover, outFName + FNameFragExt, FragCoverDescr)
		, _readCovers(cSizes, 2, saveCover, outFName + FNameReadExt, ReadCoverDescr)

		//, _regions(cSizes, 2-Glob::IsPE, saveInter, outFName + RegionsExt, RegionsDescr)
		, _regions(cSizes, 3, saveInter, outFName + RegionsExt, RegionsDescr)
		//, _splines(cSizes, 2, saveInter, outFName + ReadSplineExt, ReadSplineDescr)
		, _splines(cSizes, 3, saveInter, outFName + ReadSplineExt, ReadSplineDescr)

		, _derivs(cSizes, 2, saveInter, outFName + DerivExt, DerivDescr)
#ifdef MY_DEBUG
		, _lineWriter(cSizes, 2, saveInter, outFName + RegressionExt, RegressionDescr, DARK)
		, _splineWriter(cSizes, 1, saveInter, outFName + FragSplineExt, FragSplineDescr)
		, _outFName(outFName)
#endif
		, _bss(cSizes, 1, false, outFName + BS_Ext, BS_Descr)
		, _fIdent(true)
	{
		if (Verb::Level(Verb::RT))
			printf("%s-end sequence\n", Glob::IsPE ? "paired" : "single");
		_file = &file;

#ifndef TIMING
		if (Glob::FragLenUndef)		// no need for _reads if fragment is defined by user
#endif
		{
			auto capacity = file.EstItemCount();
			if (!Options::GetBVal(oSINGLE_CHROM))
				capacity /= 10;		// about the size of first chrom in multi-chrom case
			_reads.Reserve(capacity);
		}
		file.Pass(*this);
		_file = nullptr;
	}


	// Pre-covered data constructor
	//	@param inFName: fragment coverage file name
	//	@param outFName: common output file name
	//	@param cSizes: chrom sizes
	//	@param saveInter: if true then save intermediate data to files
	Detector(
		const char* inFName,
		const string& outFName,
		ChromSizes& cSizes, 
		bool saveInter
	)
		: _cSizes(cSizes)
		, _saveCover(false)
		, _fragCovers(cSizes,3-2*Glob::IsPE, false, outFName + FNameFragExt,	FragCoverDescr)
		, _readCovers(cSizes,	2,			 false,	outFName + FNameReadExt,	ReadCoverDescr)
		, _regions	 (cSizes,2-Glob::IsPE,saveInter,outFName + RegionsExt,		RegionsDescr)
		, _splines	 (cSizes,	2,	saveInter,		outFName + ReadSplineExt,	ReadSplineDescr)
		, _derivs	 (cSizes,	2,	saveInter,		outFName + DerivExt,		DerivDescr)
#ifdef MY_DEBUG
		, _lineWriter(cSizes,	2,	saveInter,	outFName + RegressionExt, RegressionDescr, DARK)
		, _splineWriter(cSizes,	1,	saveInter,	outFName + FragSplineExt, FragSplineDescr)
		, _outFName(outFName)
#endif
		, _bss		 (cSizes,	1,	true, outFName + BS_Ext, BS_Descr)
		, _fIdent(true)
	{
		// *** preparing coverage data
		{
			const char* pattName = strrchr(inFName, USCORE);	// the presence of USCORE has already been checked
			// check inFName for the presence of '_frag'
			if (string(pattName, strrchr(inFName, DOT) - pattName) != FNameFragExt)
				Err("only fragment coverage file is permissible", inFName).Throw();

			string baseName(inFName,  pattName - inFName);
			tChromsOccurrs chrReadOccurs;	// chromosome reading occurrences
			BYTE occursCnt = 3;				// count of reading operations
			// Fills read/frag coverage by Bedgraph file
			//	@param cover: filled read/frag coverage
			//	@param baseName: common parth of Bedgraph file's name
			auto fillStrandsCover = [this, &chrReadOccurs](CombCover& cover, const string& baseName)
			{
				auto fillStrandCover = [this, &chrReadOccurs](CombCover& cover, const string& baseName, eStrand strand)
				{
					CombCoverReader ccr(
						FS::CheckedFileName((baseName + sStrandEXT[strand] + FT::Ext(FT::BGRAPH)).c_str()),
						_cSizes, cover, chrReadOccurs, strand);
				};

				fillStrandCover(cover, baseName, FWD);
				fillStrandCover(cover, baseName, RVS);
			};

			_timer.Start();
			CombCoverReader cvr(inFName, cSizes, _fragCovers, chrReadOccurs, TOTAL);	// fill total _fragCovers
			if (!Glob::IsPE) {
				fillStrandsCover(_fragCovers, baseName + FNameFragExt);
				occursCnt += 2;
			}
			fillStrandsCover(_readCovers, baseName + FNameReadExt);
			
			// Generally speaking, the input data are independent of each other,
			// so they may contain mismatched chromosomes.
			// We set as treated only the chromosomes common to all input data.
			_cSizes.SetTreated(chrReadOccurs, occursCnt);

			_timer.Stop("Reading coverage: "); cout << LF;
		}
		// *** treatment
		for(const auto& c : _cSizes)
			if(c.second.Treated)
				CallBS(c.first);
	}

	// treats current item
	//	@param unsorted: true if unsorted input is detected
	//	@returns: true if item is accepted
	bool operator()(bool unsorted) {
		auto& rgn = _file->ItemRegion();
		bool reverse = !_file->ItemStrand();

		if (/*unsorted && */_unsortNotSet) {
			_fragCovers.SetUnsortedInput();
			_readCovers.SetUnsortedInput();
			_unsortNotSet = false;
		}
#ifdef TIMING
		_reads.AddRead(rgn, reverse);
		return true;
#endif
		if (Glob::IsPE) {
			Region frag;
			const Read read(*_file);

			if (_fIdent(read, _file->ReadLength(), frag))
				_fragCovers.AddFrag(frag);
		}
		else {
			_fragCovers.AddExtRead(rgn, reverse);
			if (Glob::FragLenUndef)
				_reads.AddRead(rgn, reverse);
		}
		_readCovers.AddRead(rgn, reverse);
		return true;
	}

	// Closes current chrom, open next one
	//	@param cID: current chrom ID
	//	@param cLen: chrom length
	//	@param cnt: current chrom items count
	//	@param nextcID: next chrom ID
	void operator()(chrid cID, chrlen cLen, size_t cnt, chrid nextcID) {
		Verb::PrintMsgVar(Verb::RT, "%s", Chrom::ShortName(nextcID).c_str(), cnt);
		if (cnt)		// not the first readed chrom
			CallBS(cID);
		_fragCovers.SetChrom(nextcID);
		_readCovers.SetChrom(nextcID);
	}

	// Closes last chrom
	//	@param cID: last chrom ID
	//	@param cLen: chrom length
	//	@param cnt: last chrom items count
	//	@param tCnt: total items count
	void operator()(chrid cID, chrlen cLen, size_t cnt, size_t)
	{ 
		Verb::PrintMsgVar(Verb::RT, ": %zu reads\n", cnt);
		_timer.Stop("Reading alignment: ");	cout << LF;
		if (cnt)
			CallBS(cID);
	}
};
