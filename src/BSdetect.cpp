/************************************************************************************
BSdetect is designed to deconvolve real Binding Sites in NGS alignment

Copyright (C) 2021 Fedor Naumenko (fedor.naumenko@gmail.com)
-------------------------
Last modified: 06/10/2024
-------------------------

This program is free software. It is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY;
without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
See GNU General Public License for more details.
************************************************************************************/

#include "Options.h"
#include "BSdetect.h"
#include "Treatment.h"

using namespace std;

const string Product::Title = "BSdetect";
const string Product::Version = "1.0";
const string Product::Descr = "binding sites detector";

const char* ProgParam = "<in-file>";	// program parameter tip
const string OutFileExt = FT::Ext(FT::eType::DIST);

//const char* inputs[] = { "FRAG","READ" };	// input option; corresponds to Inp
const char* smodes[] = { "SE","PE" };				// corresponds to DataWriter::eMode

// *** Options definition

enum eOptGroup { gOTHER };						// the only member - no gropus in help
const char* Options::OptGroups[] = { NULL };	// no gropus in help
const BYTE Options::GroupCount = ArrCnt(Options::OptGroups);

const BYTE Options::Option::IndentInTabs = 3;

// { char, str, Signs (8: hidden), type, group, 
//	defVal (if NO_DEF then no default value printed),
//	minVal (if NO_VAL then value is prohibited), maxVal, strVal, descr, addDescr }
Options::Option Options::List[] = {
	{ 'g', sGen,	tOpt::NONE,	tNAME,	gOTHER, vUNDEF, 0, 0, NULL, "chromosome sizes file", NULL },
	{ 'c',Chrom::Abbr,tOpt::NONE,tNAME,	gOTHER,	NO_DEF, 0, 0, NULL, "treat specified chromosome only", NULL },
	{ 'd',"dup-lvl",tOpt::OBLIG, tINT,	gOTHER,	1, 0, 2, NULL,
	"duplicate reads rejection level:\n-1 - keep all duplicates, 1 - keep one among duplicates, 2 - keep two among duplicates", NULL },
	{ 'f',"fr-len",	tOpt::NONE,	tINT,	gOTHER, 0, 50, 1000, NULL, "mean fragment length for SE sequence [AUTO]", NULL },
	{ 's',"save-cover",tOpt::NONE,tENUM,gOTHER,	FALSE,	NO_VAL,	0, NULL, "save coverage", NULL },
	{ 'i',"save-inter",tOpt::HIDDEN,tENUM,	gOTHER,	FALSE,	NO_VAL,	0, NULL, "save intermediate data", NULL },
	{ 'w', "warn",	tOpt::HIDDEN,tENUM,	gOTHER, FALSE,	NO_VAL, 0, NULL,
	"print each read ambiguity, if they exist" },
	{ 'r',"rank-score",	tOpt::NONE,tENUM,gOTHER, TRUE, 0, 2, (char*)Booleans,
	"turn on/off rendering the main result score in greyscale", NULL },
	{ 'O', sOutput,	tOpt::NONE,	tNAME,	gOTHER,	NO_DEF,	0,	0, NULL, "output files common name", NULL },
	{ 't',	sTime,	tOpt::NONE,	tENUM,	gOTHER,	FALSE,	NO_VAL, 0, NULL, sPrTime, NULL },
	{ 'V',"verbose",tOpt::NONE,	tENUM,	gOTHER, Verb::RT, Verb::CRIT, float(Verb::Size()), (char*)Verb::ValTitles, Verb::ValDescr, NULL },
	{ 'v',	sVers,	tOpt::NONE,	tVERS,	gOTHER,	NO_DEF, NO_VAL, 0, NULL, sPrVersion, NULL },
	{ 'h',	sHelp,	tOpt::NONE,	tHELP,	gOTHER,	NO_DEF, NO_VAL, 0, NULL, sPrUsage, NULL },

	{ 'R',"rd-len",	tOpt::NONE,	tINT,	gOTHER, 50, 20, 1000, NULL,
	"fixed length of output read, or minimum length of variable reads", NULL },
};
const BYTE Options::OptCount = ArrCnt(Options::List);

const Options::Usage Options::Usages[] = {	// content of 'Usage' variants in help
	{ vUNDEF, ProgParam, true, "bla-bla" }
};
const BYTE Options::UsageCount = ArrCnt(Options::Usages);

/*****************************************/
int main(int argc, char* argv[])
{
	int fileInd = Options::Parse(argc, argv, ProgParam);
	if (fileInd < 0)	return 1;		// wrong option or tip output
	int ret = 0;						// main() return code

	Chrom::SetUserChrom(Options::GetSVal(oCHROM));
	Mutex::Init(false);
	Timer::Enabled = Options::GetBVal(oTIME);
	Timer timer;
	try {
		const char* iName = FS::CheckedFileName(argv[fileInd]);	// input name
		const char* oName = Options::GetSVal(oOUTFILE);			// output name
		const char* gName = Options::GetSVal(oGEN);				// chrom sizes

#ifdef MY_DEBUG
		TreatedCover::WriteDelim = true;
#endif
		BedWriter::SetRankScore(Options::GetBVal(oRANK_SCORE));
		Verb::Set(Options::GetUIVal(oVERB));
		Glob::SetFragLen(Options::GetIVal(oFRAG_LEN));

		auto ftype = FT::GetType(iName);
		if (!gName && ftype != FT::BAM)
			Err(Options::OptionToStr(oGEN) + " is required while input file is not BAM", iName).Throw();

		ChromSizes cSizes(gName, true);

		// pre-covered data mode
		if (ftype == FT::BGRAPH)
		{
			Glob::ReadLen = Options::GetUIVal(oREAD_LEN);
			{	// check iName for pattern match
				const char* msg = "invalid fragment coverage file name";
				const char* pattName = strchr(iName, '_');
				if (pattName) {
					pattName++;
					if(*pattName != 'P' && *pattName != 'S')
						Err(msg, iName).Throw();
					Glob::SetPE(*pattName == 'P');
				}
				else
					Err(msg, iName).Throw();
			}

			Detector bsd(
				iName,
				FS::ComposeFileName(Options::GetSVal(oOUTFILE), iName),
				cSizes,
				Options::GetBVal(oSAVE_INTER)
			);
		}
		// main mode
		else
		{
			// pre-read first item to check for PE sequence
			RBedReader file(iName, &cSizes, Options::GetIVal(oDUP_LVL), eOInfo::LAC, Verb::Level(Verb::DBG), false, true, true);
			Verb::PrintMsg(Verb::DBG);
			file.GetNextItem();		// no need to check for empty sequence
			Glob::SetPE(file.IsPaired());

			// detect BS
			Detector bsd(
				file,
				FS::ComposeFileName(Options::GetSVal(oOUTFILE), iName),
				cSizes,
				Options::GetBVal(oSAVE_COVER),
				Options::GetBVal(oSAVE_INTER)
			);
		}
	}
	catch (const Err& e) { ret = 1; cerr << e.what() << endl; }
	catch (const exception& e) { ret = 1; cerr << e.what() << endl; }
	timer.Stop();
	return ret;
}

//===== Detector

float Detector::GetPeakPosDiff(chrid cID)
{
	DataSet<TreatedCover>& fragCovers = _frag�overs.ChromData(cID);
	DataSet<TreatedCover>& readCovers = _read�overs.ChromData(cID);
	DataCoverRegions& regions = static_cast<DataCoverRegions&>(_regions.ChromData(cID));
	DataValuesMap& splines	  = static_cast<DataValuesMap&>(_splines.ChromData(cID));

	Verb::PrintMsg(Verb::DBG, "Determine mean fragment length");
	_timer.Start();
	coval maxVal = readCovers.StrandData(POS).GetMaxVal();
	Verb::PrintMsgVar(Verb::DBG, "Max cover: %d;  cutoff: %d\n", maxVal, maxVal / 4);
	if (regions.SetPotentialRegions(fragCovers, _cSizes[cID], maxVal / 4, true))	return false;
	// calculate mean difference as the average of three attempts
	float peakDiff = 0;
	BYTE cnt = 0;

	Verb::PrintMsg(Verb::DBG);
	for (BYTE splineBase = 15; splineBase <= 85; splineBase += 35) {
		splines.BuildSpline(fragCovers, regions, false, splineBase);
		auto diff = splines.GetPeakPosDiff();
		peakDiff += diff;
		Verb::PrintMsgVar(Verb::DBG, "spline base: %d  diff: %.2f\n", splineBase, diff);
		splines.Clear();
		cnt++;
	}
	_timer.Stop(0, false, true);	cout << LF;
	return peakDiff / cnt;
}

void Detector::CallBS(chrid cID)
{
	DataSet<TreatedCover>& fragCovers = _frag�overs.ChromData(cID);
	DataSet<TreatedCover>& readCovers = _read�overs.ChromData(cID);
	DataCoverRegions& regions		= static_cast<DataCoverRegions&>(_regions.ChromData(cID));
	DataValuesMap& splines		= static_cast<DataValuesMap&>(_splines.ChromData(cID));
	DataBoundsValuesMap& derivs = static_cast<DataBoundsValuesMap&>(_derivs.ChromData(cID));
	BS_map& bss = *_bss.ChromData(cID).Data();
	const chrlen cLen = _cSizes[cID];
	bool	resetCover = true;

	if (!Glob::ReadLen)	Glob::ReadLen = _file->ReadLength();

#ifdef MY_DEBUG
	_lineWriter.SetChromID(cID);	BoundsValues::SetSpecialWriter(_lineWriter);
	_splineWriter.SetChromID(cID);	TreatedCover::SetSpecialWriter(_splineWriter);
	//Incline::SetOutFile(cID, "incline.txt");
#endif

	if (Glob::FragLenUndef) {		// can be true for SE sequence only
		auto peakDiff = short(round(GetPeakPosDiff(cID)));
		Verb::PrintMsgVar(Verb::RT, "Mean fragment length: %d\n", FragDefLEN - peakDiff);
		if (peakDiff > 10) {			// significant difference
			Verb::PrintMsgVar(Verb::RT, "Rebuild coverages\n");
			Glob::FragLen -= peakDiff;
			_frag�overs.Clear();
			_frag�overs.Fill(_reads);
			regions.Clear();
		}
		else
			resetCover = false;
		_reads.Clear();
		Glob::FragLenUndef = false;
	}
	Verb::PrintMsg(Verb::RT, "Locate binding sites");
	_timer.Start();
	if (resetCover && regions.SetPotentialRegions(fragCovers, cLen, 3))	return;

	splines.BuildSpline(readCovers, regions, true);	_regions.WriteChrom(cID);
	splines.EliminateNonOverlaps();
	if (Verb::Level(Verb::DBG))		splines.PrintStat(cLen);
	splines.Numerate();
	derivs.Set(splines);			_splines.WriteChrom(cID);
	
	bss.Set(derivs, readCovers);	_read�overs.WriteChrom(cID); _derivs.WriteChrom(cID);
	//bss.Print(cID, _outFName + ".BSS_dump0.txt", false);
	bss.Refine();
	//bss.Print(cID, _outFName + ".BSS_dump1.txt", false);
	bss.SetScore(fragCovers);		_frag�overs.WriteChrom(cID);
#ifdef MY_DEBUG
	bss.Print(cID, _outFName + ".BSS_dump.txt", false);
	//bss.CheckScoreHierarchy();
#endif
	bss.PrintStat();
	_bss.WriteChrom(cID);

	_timer.Stop("Threament: ");	cout << LF;
}
