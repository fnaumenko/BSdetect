/**********************************************************
callDist.h (c) 2021 Fedor Naumenko (fedor.naumenko@gmail.com)
All rights reserved.
-------------------------
Last modified: 03/22/2025
-------------------------
Provides main functionality
***********************************************************/

#pragma once
#include "Treatment.h"

enum optValue {		// options id
	oBIN,
	oGEN,
	oCHROM,
	oDUP_LVL,
	oFRAG_LEN,
	oSAVE_COVER,
	oSAVE_INTER,
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

//#define TIMING

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

	// Fills read/frag coverage by Bedgraph file
	//	@param cover: filled read/frag coverage
	//	@param baseName: common parth of Bedgraph file's name
	//	@param chrFreq: chrom frequency counter
	//	@param strand: strand
	void FillStrandCover(CombCover& cover, const string& baseName, tChromsFreq& chrFreq, eStrand strand)
	{
		CombCoverReader ccr(
			FS::CheckedFileName((baseName + sStrandEXT[strand] + FT::Ext(FT::BGRAPH)).c_str()),
			_cSizes, cover, chrFreq, strand);
	}

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
		, _fragCovers(cSizes, 3-2*Glob::IsPE, saveCover, outFName + FNameFragExt, "fragment coverage")
		, _readCovers(cSizes, 2, saveCover, outFName + FNameReadExt, "read coverage")

		//, _regions(cSizes, 2-Glob::IsPE, saveInter, outFName + ".RGNS", "potential regions")
		, _regions(cSizes, 3, saveInter, outFName + ".RGNS", "potential regions")
		//, _splines(cSizes, 2, saveInter, outFName + ".SPLINE", "read coverage spline")
		, _splines(cSizes, 3, saveInter, outFName + ".SPLINE", "read coverage spline")

		, _derivs(cSizes, 2, saveInter, outFName + ".DERIV", "derivative of read coverage spline")
#ifdef MY_DEBUG
		, _lineWriter(cSizes, 2, saveInter, outFName + ".LINE", "linear regression", DARK)
		, _splineWriter(cSizes, 1, saveInter, outFName + ".FR_SPLINE", "fragment coverage spline")
		, _outFName(outFName)
#endif
		//, _bss(cSizes, 1, true, outFName + ".BSs", "called binding sites")
		, _bss(cSizes, 1, false, outFName + ".BSs", "called binding sites")
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
		, _fragCovers(cSizes,3-2*Glob::IsPE, false, outFName + FNameFragExt, "fragment coverage")
		, _readCovers(cSizes,	2,			 false,	outFName + FNameReadExt, "read coverage")
		, _regions	 (cSizes,2-Glob::IsPE,saveInter,outFName + ".RGNS"	, "potential regions")
		, _splines	 (cSizes,	2,	saveInter,	outFName + ".SPLINE", "read coverage spline")
		, _derivs	 (cSizes,	2,	saveInter,	outFName + ".DERIV"	, "derivative of read coverage spline")
#ifdef MY_DEBUG
		, _lineWriter(cSizes,	2,	saveInter,	outFName + ".LINE"	, "linear regression", DARK)
		, _splineWriter(cSizes,	1,	saveInter,	outFName + ".FR_SPLINE", "fragment coverage spline")
		, _outFName(outFName)
#endif
		, _bss		 (cSizes,	1,	true,		outFName + ".BSs"	, "called binding sites")
		, _fIdent(true)
	{
		// *** preparing coverage data
		{
			const char* pattName = strrchr(inFName, USCORE);	// the presence of USCORE has already been checked
			// check inFName for the presence of '_frag'
			if (string(pattName, strrchr(inFName, DOT) - pattName) != FNameFragExt)
				Err("only fragment coverage file is permissible", inFName).Throw();

			string baseName(inFName,  pattName - inFName);
			tChromsFreq	chrFreq;
			BYTE dataCnt = 3;

			_timer.Start();
			CombCoverReader cvr(inFName, cSizes, _fragCovers, chrFreq, TOTAL);	// fill total _fragCovers
			if (!Glob::IsPE) {
				const string extBaseName = baseName + FNameFragExt;
				FillStrandCover(_fragCovers, extBaseName, chrFreq, FWD);
				FillStrandCover(_fragCovers, extBaseName, chrFreq, RVS);
				dataCnt += 2;
			}
			baseName += FNameReadExt;
			FillStrandCover(_readCovers, baseName, chrFreq, FWD);
			FillStrandCover(_readCovers, baseName, chrFreq, RVS);
			
			// ** set treated chroms
			_cSizes.TreatedAll(false);
			for (const auto& c : chrFreq)	// there're chroms represented in input Bedgraph only 
				_cSizes.TreatedChrom(c.first, c.second == dataCnt);	// all active readers

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
