/**********************************************************
callDist.h (c) 2021 Fedor Naumenko (fedor.naumenko@gmail.com)
All rights reserved.
-------------------------
Last modified: 05/092024
-------------------------
Provides main functionality
***********************************************************/

#pragma once
#include "Treatment.h"

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

bool IsFragMeanUnset = true;	// common to the entire genome

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

	RBedReader* _file;		// needs only in constructor
	FragIdent _fIdent;		// needs only in constructor
	Reads	_reads;
	Timer	_timer;

	void CallBS(chrid cID);

public:
	static bool IsPEReads;			// common to the entire genome

	Detector(RBedReader& file, const string& outFName, ChromSizes& cSizes, bool saveCover, bool saveInter)
		: _cSizes(cSizes)
		, _saveCover(saveCover)
		, _frag—overs(cSizes, IsPEReads ? 1 : 2+saveCover, saveCover, outFName + "_frag", "fragment coverage")
		, _read—overs(cSizes,	2,	saveCover,	outFName + "_read"	, "read coverage")
		, _rgns(cSizes,2-IsPEReads,	saveInter,	outFName + ".RGNS"	, "potential regions")
		, _splines	(cSizes,	2,	saveInter,	outFName + ".SPLINE", "read coverage spline")
		, _derivs	(cSizes,	2,	saveInter,	outFName + ".DERIV"	, "derivative of read coverage spline")
		, _lineWriter(cSizes,	2,	saveInter,	outFName + ".LINE"	, "linear regression")
		, _bss		(cSizes,	1,	true,		outFName + ".BSs"	, "called binding sites")
		, _fIdent(true)
	{
		if (Verb::Level(Verb::DBG))
			printf("%s-end sequencing\n", IsPEReads ? "paired" : "single");
		_file = &file;
		_reads.Reserve(file.EstItemCount() / 10);	// about the size of first chrom in common case
		file.Pass(*this);
		_file = nullptr;
		_timer.Start();
	}

	// treats current item
	//	@returns: true if item is accepted
	bool operator()() {
		auto& rgn = _file->ItemRegion();
		bool reverse = !_file->ItemStrand();

		if (IsPEReads) {
			Region frag;
			const Read read(*_file);

			if (_fIdent(read, frag))
				_frag—overs.AddFrag(frag);
		}
		else {
			_frag—overs.AddExtRead(rgn, reverse, _saveCover && !IsFragMeanUnset);
			if (IsFragMeanUnset)
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
		if (Verb::Level(Verb::RT))
			dout << Chrom::ShortName(nextcID) << LF;
		if (cnt) {		// not the first readed chrom
			CallBS(cID);
		}
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
		if (cnt) {
			_timer.Stop("Reading alignment: ");	cout << LF;
			CallBS(cID);
		}
	}
};
