/**********************************************************
BSdetect.h
main class Detector
2021 Fedor Naumenko (fedor.naumenko@gmail.com)
Last modified: 11/25/23
***********************************************************/
#pragma once

#include "Treatment.h"
#include "DataReader.h"

enum optValue {		// options id
	oGEN,
	oCHROM,
	oDUP_LVL,
	oCOVER,
	oINTERM,
	oALARM,
	oRANK_SCORE,
	oOUTFILE,
	oTIME,
	oVERB,
	oVERSION,
	oHELP,
};

bool IsFragMeanUnset = true;

// BS detector
class Detector
{
	const ChromSizes& _cSizes;
	bool			  _saveCover;
	CombCover		  _frag—overs;		// extended reads cover to find frag Mean
	CombCover		  _read—overs;
	OCoverRegions	  _rgns;
	OValuesMap		  _splines;
	ORegionsValuesMap _derivs;
	OLinearWriter	  _lineWriter;
	OBS_Map			  _bss;

	RBedReader* _file;			// needs only in constructor
	Reads	_reads;
	Timer	_timer;

	void CallBS(chrid cID);

public:
	Detector(const char* inFName, const string& outFName, ChromSizes& cSizes, bool saveCover, bool saveInter, eOInfo info)
		: _cSizes(cSizes)
		, _saveCover (saveCover)
		, _frag—overs(cSizes, 2 + saveCover, saveCover, outFName + "_frag", "fragment coverage")
		, _read—overs(cSizes, 2, saveCover, outFName + "_read", "read coverage")
		, _rgns		 (cSizes, 2, saveInter, outFName + ".RGNS", "potential regions")
		, _splines	 (cSizes, 2, true, outFName + ".SPLINE", "read coverage spline")
		, _derivs	 (cSizes, 2, saveInter, outFName + ".DERIV", "derivative of read coverage spline")
		, _lineWriter(cSizes, 2, saveInter, outFName + ".LINE", "linear regression")
		, _bss		 (cSizes, 1, false,	outFName + ".BSs", "called binding sites")
	{
		if (Verb::Level(Verb::RT))
			dout << inFName << '\n';
		RBedReader file(inFName, &cSizes, Options::GetIVal(oDUP_LVL), info, false, false);
		_file = &file;
		_reads.Reserve(file.EstItemCount() / 10);	// about the size of first chrom in common case
		file.Pass(*this);
		_file = nullptr;
		_timer.Start();
	}

	// treats current item
	//	@return: true if item is accepted
	bool operator()() {
		auto& rgn = _file->ItemRegion();
		bool reverse = !_file->ItemStrand();

		_frag—overs.AddExtRead(rgn, reverse, _saveCover && !IsFragMeanUnset);
		_read—overs.AddRead(rgn, reverse);
		if (IsFragMeanUnset)
			_reads.AddRead(rgn, reverse);
		return true;
	}

	// Closes current chrom, open next one
	//	@cID: current chrom ID
	//	@cLen: chrom length
	//	@cnt: current chrom items count
	//	@nextcID: next chrom ID
	void operator()(chrid cID, chrlen cLen, size_t cnt, chrid nextcID) {
		if (Verb::Level(Verb::RT))
			dout << Chrom::ShortName(nextcID) << LF;
		if (cnt) {		// not the first readed chrom
			CallBS(cID);
		}
		_frag—overs.SetChrom(nextcID);
		_read—overs.SetChrom(nextcID);
	}

	// Closes last chrom
	//	@cID: last chrom ID
	//	@cLen: chrom length
	//	@cnt: last chrom items count
	//	@tCnt: total items count
	void operator()(chrid cID, chrlen cLen, size_t cnt, ULONG tCnt)
	{ 
		if (cnt) {
			_timer.Stop("Reading alignment: ");
			CallBS(cID);
		}
	}
};
