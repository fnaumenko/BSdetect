/**********************************************************
callDist.h (c) 2021 Fedor Naumenko (fedor.naumenko@gmail.com)
All rights reserved.
-------------------------
Last modified: 06/18/2024
-------------------------
Provides main functionality
***********************************************************/

#pragma once
#include "Treatment.h"

enum optValue {		// options id
	oGEN,
	oCHROM,
	oDUP_LVL,
	oFRAG_LEN,
	oSAVE_COVER,
	oSAVE_INTER,
	oALARM,
	oRANK_SCORE,
	oOUTFILE,
	oTIME,
	oVERB,
	oVERSION,
	oHELP,

	oREAD_LEN
};

// BS detector
class Detector
{
	const string FNameFragExt = "_frag";
	const string FNameReadExt = "_read";

	ChromSizes&		 _cSizes;
	bool			 _saveCover;
	CombCover		 _frag—overs;		// extended reads cover to find frag Mean
	CombCover		 _read—overs;
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
	FragIdent	_fIdent;		// needs only for input reading
	Reads		_reads;
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
		, _frag—overs(cSizes, 3-2*Glob::IsPE, saveCover, outFName + FNameFragExt, "fragment coverage")
		, _read—overs(cSizes, 2, saveCover, outFName + FNameReadExt, "read coverage")
		, _regions(cSizes, 2-Glob::IsPE, saveInter, outFName + ".RGNS", "potential regions")
		, _splines(cSizes, 2, saveInter, outFName + ".SPLINE", "read coverage spline")
		, _derivs(cSizes, 2, saveInter, outFName + ".DERIV", "derivative of read coverage spline")
#ifdef MY_DEBUG
		, _lineWriter(cSizes, 2, saveInter, outFName + ".LINE", "linear regression", DARK)
		, _splineWriter(cSizes, 1, saveInter, outFName + ".FR_SPLINE", "fragment coverage spline")
		, _outFName(outFName)
#endif
		, _bss		(cSizes,	1,	true,		outFName + ".BSs"	, "called binding sites")
		, _fIdent(true)
	{
		if (Verb::Level(Verb::RT))
			printf("%s-end sequence\n", Glob::IsPE ? "paired" : "single");
		_file = &file;
		_reads.Reserve(file.EstItemCount() / 10);	// about the size of first chrom in common case
		file.Pass(*this);
		_file = nullptr;
		_timer.Start();
	}

	// Pre-covered data constructor. For PE only!
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
		, _frag—overs(cSizes,3-2*Glob::IsPE, false, outFName + FNameFragExt, "fragment coverage")
		, _read—overs(cSizes,	2,			 false,	outFName + FNameReadExt, "read coverage")
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
		// preparing coverage data
		{
			const char* pattName = strrchr(inFName, '_');
			{
				// check inFName
				string fragExt(pattName, strrchr(inFName, DOT) - pattName);
				if (fragExt != FNameFragExt)
					Err("only fragment coverage file is permissible", inFName).Throw();
			}
			string baseName(inFName,  pattName - inFName);

			tChromsFreq	chrFreq;
			_timer.Start();
			// initialize covered data
			CombCoverReader a(inFName, cSizes, _frag—overs, chrFreq, TOTAL);
			if (!Glob::IsPE) {
				string fCoverPos = baseName + FNameFragExt + sStrandEXT[POS] + FT::Ext(FT::BGRAPH);
				string fCoverNeg = baseName + FNameFragExt + sStrandEXT[NEG] + FT::Ext(FT::BGRAPH);

				CombCoverReader b(FS::CheckedFileName(fCoverPos.c_str()), cSizes, _frag—overs, chrFreq, POS);
				CombCoverReader c(FS::CheckedFileName(fCoverNeg.c_str()), cSizes, _frag—overs, chrFreq, NEG);
			}
			baseName += FNameReadExt;
			CombCoverReader b(
				FS::CheckedFileName((baseName + sStrandEXT[POS] + FT::Ext(FT::BGRAPH)).c_str()),
				cSizes, _read—overs, chrFreq, POS);
			CombCoverReader c(
				FS::CheckedFileName((baseName + sStrandEXT[NEG] + FT::Ext(FT::BGRAPH)).c_str()),
				cSizes, _read—overs, chrFreq, NEG);

			cSizes.TreateAll(false);
			BYTE dataCnt = Glob::IsPE ? 3 : 5;
			for (const auto& c : chrFreq)
				_cSizes.TreateChrom(c.first, c.second == dataCnt);	// all readers worked

			_timer.Stop("Reading coverage: "); cout << LF;
		}
		// treatment
		for(const auto& c : _cSizes)
			if(c.second.Treated)
				CallBS(c.first);
	}

	// treats current item
	//	@returns: true if item is accepted
	bool operator()() {
		auto& rgn = _file->ItemRegion();
		bool reverse = !_file->ItemStrand();

		if (Glob::IsPE) {
			Region frag;
			const Read read(*_file);

			if (_fIdent(read, frag))
				_frag—overs.AddFrag(frag);
		}
		else {
			_frag—overs.AddExtRead(rgn, reverse);
			if (Glob::FragLenUndef)
				_reads.AddRead(rgn, reverse);
		}
		_read—overs.AddRead(rgn, reverse);
		return true;
	}

	// Closes current chrom, open next one
	//	@param cID: current chrom ID
	//	@param cLen: chrom length
	//	@param cnt: current chrom items count
	//	@param nextcID: next chrom ID
	void operator()(chrid cID, chrlen cLen, size_t cnt, chrid nextcID) {
		//Verb::PrintMsg(Verb::RT, Chrom::ShortName(nextcID).c_str());
		Verb::PrintMsgVar(Verb::RT, "%s", Chrom::ShortName(nextcID).c_str(), cnt);
		if (cnt)		// not the first readed chrom
			CallBS(cID);
		_frag—overs.SetChrom(nextcID);
		_read—overs.SetChrom(nextcID);
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
