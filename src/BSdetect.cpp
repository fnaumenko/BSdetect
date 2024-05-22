/************************************************************************************
BSdetect is designed to deconvolve real Binding Sites in NGS alignment

Copyright (C) 2021 Fedor Naumenko (fedor.naumenko@gmail.com)
-------------------------
Last modified: 05/21/2024
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
	{ 'd',"dup-lvl",tOpt::OBLIG,	tINT,	gOTHER,	1, 0, 2, NULL,
	"duplicate reads rejection level:\n-1 - keep all duplicates, 1 - keep one among duplicates, 2 - keep two among duplicates", NULL },
	{ 's',"save-cover",tOpt::NONE,	tENUM,	gOTHER,	FALSE,	NO_VAL,	0, NULL, "save coverage", NULL },
	{ 'i',"save-inter",tOpt::HIDDEN,tENUM,	gOTHER,	FALSE,	NO_VAL,	0, NULL, "save intermediate data", NULL },
	{ 'w', "warn",	tOpt::HIDDEN,tENUM,	gOTHER, FALSE,	NO_VAL, 0, NULL,
	"print each read ambiguity, if they exist" },
	{ 'r',"rank-score",	tOpt::NONE,	tENUM,	gOTHER, TRUE, 0, 2, (char*)Booleans,
	"turn on/off rendering the main result score in greyscale", NULL },
	{ 'O', sOutput,	tOpt::NONE,	tNAME,	gOTHER,	NO_DEF,	0,	0, NULL, "output files common name", NULL },
	{ 't',	sTime,	tOpt::NONE,	tENUM,	gOTHER,	FALSE,	NO_VAL, 0, NULL, sPrTime, NULL },
	{ 'V',"verbose",tOpt::NONE,	tENUM,	gOTHER, Verb::RT, Verb::CRIT, float(Verb::Size()), (char*)Verb::ValTitles, Verb::ValDescr, NULL },
	{ 'v',	sVers,	tOpt::NONE,	tVERS,	gOTHER,	NO_DEF, NO_VAL, 0, NULL, sPrVersion, NULL },
	{ 'h',	sHelp,	tOpt::NONE,	tHELP,	gOTHER,	NO_DEF, NO_VAL, 0, NULL, sPrUsage, NULL },

	{ 'D', "read-cover-dir",	tOpt::NONE,	tNAME,	gOTHER, vUNDEF, 0, 0, NULL, "direct read coverage file", NULL},
	{ 'R', "read-cover-rev",	tOpt::NONE,	tNAME,	gOTHER, vUNDEF, 0, 0, NULL, "reversed read coverage file", NULL},
	{ 'm', "smode",	tOpt::NONE,	tENUM,	gOTHER, 0, 0, ArrCnt(smodes), (char*)smodes,
	"sequencing mode: ? - single end, ? - paired end", NULL },
	{ 'L',"rd-len",	tOpt::NONE,	tINT,	gOTHER, 50, 20, 1000, NULL,
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
		auto ftype = FT::GetType(iName);

		//if (!gName && !FS::HasExt(iName, FT::Ext(FT::BAM)))
		if (!gName && ftype != FT::BAM)
			Err(Options::OptionToStr(oGEN) + " is required while input file is not BAM", iName).Throw();

		BedWriter::SetRankScore(Options::GetBVal(oRANK_SCORE));
		Verb::Set(Options::GetUIVal(oVERB));
		ChromSizes cSizes(gName, true);

		if (ftype == FT::BED) {
			// pre-read first item to check for PE sequence
			RBedReader file(iName, &cSizes, Options::GetIVal(oDUP_LVL), eOInfo::LAC, Verb::Level(Verb::DBG), false, true, true);
			if (Verb::Level(Verb::DBG))	cout << LF;
			file.GetNextItem();		// no need to check for empty sequence
			Detector::IsPEReads = file.IsPaired();

			// detect BS
			Detector bsd(
				file,
				FS::ComposeFileName(Options::GetSVal(oOUTFILE), iName),
				cSizes,
				Options::GetBVal(oSAVE_COVER),
				Options::GetBVal(oSAVE_INTER)
			);
		}
		else {
			Detector::IsPEReads = Options::GetBVal(oSMODE);
			Glob::ReadLen = Options::GetUIVal(oRD_LEN);

			Detector bsd(
				iName,
				FS::CheckedFileName(Options::GetSVal(oPOS_READ_COVER)),
				FS::CheckedFileName(Options::GetSVal(oNEG_READ_COVER)),
				FS::ComposeFileName(Options::GetSVal(oOUTFILE), iName),
				cSizes,
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
bool Detector::IsPEReads = false;

void Detector::CallBS(chrid cID)
{
	const char* noRgnsMsg = "No potential regions found";

	DataSet<TreatedCover>& fragCovers = _frag�overs.ChromData(cID);
	DataSet<TreatedCover>& readCovers = _read�overs.ChromData(cID);
	DataCoverRegions& rgns		 = static_cast<DataCoverRegions&>(_rgns.ChromData(cID));
	DataValuesMap& splines		 = static_cast<DataValuesMap&>(_splines.ChromData(cID));
	DataBoundsValuesMap& derivs = static_cast<DataBoundsValuesMap&>(_derivs.ChromData(cID));
	BS_map& bss = *_bss.ChromData(cID).Data();
	const chrlen cLen = _cSizes[cID];

	if (fragCovers.Empty()) { Verb::PrintMsg(Verb::RT, "Empty fragment coverage"); return; }
	if (readCovers.Empty()) { Verb::PrintMsg(Verb::RT, "Empty read coverage"); return; }
	_timer.Start();
	if (!Glob::ReadLen)	Glob::ReadLen = _file->ReadLength();
	_lineWriter.SetChromID(cID);

	if (IsPEReads) {
		rgns.SetPotentialRegionsPE(fragCovers, cLen, 3);
		if (rgns.Empty()) { Verb::PrintMsg(Verb::RT, noRgnsMsg); return; }

		splines.BuildSplinePE(readCovers, rgns, true);
		_rgns.WriteChrom(cID);	_frag�overs.WriteChrom(cID);	// now we can release both
	}
	else {
		if (IsFragMeanUnset) {
			Verb::PrintMsg(Verb::DBG, "Determine mean fragment length");
			_timer.Start();
			coval maxVal = readCovers.StrandData(POS).GetMaxVal();
			printf("Max cover: %d;  cutoff: %d\n", maxVal, maxVal / 3);
			rgns.SetPotentialRegionsSE(fragCovers, cLen, maxVal / 3, true);	//10
			if (rgns.Empty()) { Verb::PrintMsg(Verb::RT, noRgnsMsg); return; }

			//auto flen = rgns.GetFragMean(fragCovers);
			//printf("\nMass Mean fragment length: %d\n", flen);

			splines.BuildSplineSE(fragCovers, rgns, false, 50);	_frag�overs.Clear();	rgns.Clear();
			Glob::FragLen = splines.GetFragMean();				splines.Clear();
			printf("Mean fragment length: %d\n", Glob::FragLen);

			_frag�overs.Fill(_reads, _saveCover);				_reads.Clear();
			IsFragMeanUnset = false;
			_timer.Stop(0, false, true);
		}

		Verb::PrintMsg(Verb::DBG, "Locate binding sites");
		rgns.SetPotentialRegionsSE(fragCovers, cLen, 3, false);
		if (rgns.Empty()) { Verb::PrintMsg(Verb::RT, noRgnsMsg); return; }

		splines.BuildSplineSE(readCovers, rgns, true);
		_rgns.WriteChrom(cID);	_frag�overs.WriteChrom(cID);	// now we can release both
	}

	splines.Data()->EliminateNonOverlaps();
	splines.Data()->PrintStat(cLen);
	splines.Data()->Numerate();
	//splines.Print(cID);
	derivs.Set(splines);						_splines.WriteChrom(cID);
	//derivs.Print();
	bss.Set(derivs, readCovers, _lineWriter);	_read�overs.WriteChrom(cID); _derivs.WriteChrom(cID);
	//bss.Print(cID, false);
	//return;
	bss.Refine();
	bss.NormalizeScore();
	bss.NormalizeBSwidth();
#ifdef MY_DEBUG
	//bss.Print(cID, false);
	bss.CheckScoreHierarchy();
	//bss.PrintWidthDistrib();
#endif
	bss.PrintStat();
	_bss.WriteChrom(cID);

	_timer.Stop("Threament: ");	cout << LF;
}
