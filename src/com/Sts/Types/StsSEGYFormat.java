package com.Sts.Types;

import com.Sts.Actions.Import.*;
import com.Sts.DB.*;
import com.Sts.IO.*;
import com.Sts.MVC.*;
import com.Sts.UI.*;
import com.Sts.Utilities.*;

import java.io.*;
import java.util.*;

/**
   @(#) StsSEGYFormat.java
   SEGY format definition for associated SEG-Y file, trace, etc..
 */
public final class StsSEGYFormat implements Serializable
{
   public byte dataType = POSTSTACK;
   public byte type = NONE;
   public String name = "SEG-Y";
   public String description = "SEG-Y Variant";
   public boolean isSampleFormatOverride = false;
   public byte overrideSampleFormat = NONE;
   public boolean isLittleEndian = false;
   public int binaryHeaderSize = defaultBinaryHeaderSize;
   public int traceHeaderSize = defaultTraceHeaderSize;
   public int textHeaderSize = defaultTextHeaderSize;
   public int nSamp;
   public float sampleSpacing;
   public String textHeaderFormatString = EBCDIC; // Display name; encoder needs actual name which in this case is "Cp500" (see textHeaderFormatStrings & textHeaderFormats)
   public String EOT = "((SEG:EndText))";
   public byte zDomain = StsParameters.TD_TIME;
   public byte hUnits;
   public byte zUnits;
   public byte tUnits;
   public float startZ = 0.0f;
   public boolean overrideHeader = false;

   transient public String hUnitString;
   transient public String zUnitString;
   transient public String tUnitString;
   transient int overrideNSamples = -1;
   transient float overrideSampleSpacing = -1;

   static public int defaultBinaryHeaderSize = 400;
   static public int defaultTraceHeaderSize = 240;
   static public int defaultTextHeaderSize = 3200;
   static public int defaultHeaderSize = 3600;
   static public final byte defaultSampleFormat = StsSEGYFormat.IBMFLT;
   static public final String defaultTextHeaderFormatString = StsSEGYFormat.EBCDIC;
   static final float nullValue = StsParameters.nullValue;

   /** Initially none (0); in this case sampleFormat will be pulled from binaryHeader using the sampleFormatRec for location and its format.
    *  If user selects a different format, this will be used and saved as desired sampleFormat.  When this new segyFormat is loaded, the
    *  sampleFormatOverride will be used and the binaryHeader value will be ignored.
    */
//   private int sampleFormatOverride = 0;

   // global records
   private StsSEGYFormatRec jobIdNumberRec = new StsSEGYFormatRec("JOBID", 0, INT4, "Job identification number", false);
   private StsSEGYFormatRec lineNumberRec = new StsSEGYFormatRec("LINENUM", 4, INT4, "Line number", false);
   private StsSEGYFormatRec reelNumberRec = new StsSEGYFormatRec("REELNUM", 8, INT4, "Reel number", false);
   private StsSEGYFormatRec tracesPerRec = new StsSEGYFormatRec("NTRACES", 12, INT2, "Trace per record", false);
   private StsSEGYFormatRec auxsPerRec = new StsSEGYFormatRec("NAUXS", 14, INT2, "Auxiliary traces per record", false);
   private StsSEGYFormatRec sampleSpacingRec = new StsSEGYFormatRec("SAMPINT", 16, INT2, "Sample interval in micro-seconds", true);
   private StsSEGYFormatRec origSampleSpacingRec = new StsSEGYFormatRec("OSAMPINT", 18, INT2, "Original sample interval in micro-seconds", true);
   private StsSEGYFormatRec nSamplesRec = new StsSEGYFormatRec("NSAMP", 20, INT2, "Number of samples per trace", true);
   private StsSEGYFormatRec origNSamplesRec = new StsSEGYFormatRec("ONSAMP", 22, INT2, "Original number of samples per trace", false);
   private StsSEGYFormatRec sampleFormatRec = new StsSEGYFormatRec("SAMPFMT", 24, INT2, "Sample format: 1=IBM;2=int;3=short;5=IEEE;6=float8;7=float16;8=byte", true);
   private StsSEGYFormatRec cdpFoldRec = new StsSEGYFormatRec("CDPFOLD", 26, INT2, "CDP Fold (Expected traces per ensemble)", false);
   private StsSEGYFormatRec sortCodeRec = new StsSEGYFormatRec("SORTCODE", 28, INT2, "Trace sorting code: 1=recorded, 2=ensemble, 3=single fold, 4=horizontal stacked", false);
   private StsSEGYFormatRec vertSumCodeRec = new StsSEGYFormatRec("VSUMCODE", 30, INT2, "Vertical sum code: 1=none, 2=two fold...", false);
   private StsSEGYFormatRec startSweepFreqRec = new StsSEGYFormatRec("STARTFREQ", 32, INT2, "Start sweep frequency (hertz)", false);
   private StsSEGYFormatRec endSweepFreqRec = new StsSEGYFormatRec("ENDFREQ", 34, INT2, "End sweep frequency (hertz)", false);
   private StsSEGYFormatRec sweepLengthRec = new StsSEGYFormatRec("SWEEPLEN", 36, INT2, "Sweep length (milliseconds)", false);
   private StsSEGYFormatRec sweepTypeCodeRec = new StsSEGYFormatRec("SWEEPCODE", 38, INT2, "Sweep code: 1-linear, 2-parabolic, 3-exponential, 4-other", false);
   private StsSEGYFormatRec sweepTraceRec = new StsSEGYFormatRec("SWEEPTRC", 40, INT2, "Trace number of sweep channel", false);
   private StsSEGYFormatRec startSweepTaperRec = new StsSEGYFormatRec("STARTTAPER", 42, INT2, "Sweep trace taper length at start (milliseconds)", false);
   private StsSEGYFormatRec endSweepTaperRec = new StsSEGYFormatRec("ENDTAPER", 44, INT2, "Sweep trace taper length at end (milliseconds)", false);
   private StsSEGYFormatRec taperTypeCodeRec = new StsSEGYFormatRec("TAPERCODE", 46, INT2, "Sweep taper code: 1-linear, 2-cosine squared, 3-other", false);
   private StsSEGYFormatRec correlatedRec = new StsSEGYFormatRec("CORRELATED", 48, INT2, "Correlated data traces: 1-no, 2-yes", false);
   private StsSEGYFormatRec binGainRec = new StsSEGYFormatRec("BINGAIN", 50, INT2, "Binary gain recovered: 1-no, 2-yes", false);
   private StsSEGYFormatRec ampRecovMtdRec = new StsSEGYFormatRec("AMPRECOV", 52, INT2, "Amplitude recovery method code: 1-one, 2-spherical divergence, 3-AGC, 4-other", false);
   private StsSEGYFormatRec unitsRec = new StsSEGYFormatRec("UNITS", 54, INT2, "Coordinate units: 1=Meters; 2=Feet", false);
   private StsSEGYFormatRec polarityRec = new StsSEGYFormatRec("POLARITY", 56, INT2, "Impulse signal polarity: 1=Negative; 2=Positive", false);
   private StsSEGYFormatRec vibPolarityRec = new StsSEGYFormatRec("VIBPOL", 58, INT2, "Vibratory polarity code: 1-8", false);

   // Tsunami Multiple PostStack3d SegY Addition
   private StsSEGYFormatRec time01Rec = new StsSEGYFormatRec("TIME1", 200, IBMFLT, "Time for Vertex #1", false);
   private StsSEGYFormatRec value01Rec = new StsSEGYFormatRec("VALUE1", 204, IBMFLT, "Value for Vertex #1", false);
   private StsSEGYFormatRec time02Rec = new StsSEGYFormatRec("TIME2", 208, IBMFLT, "Time for Vertex #2", false);
   private StsSEGYFormatRec value02Rec = new StsSEGYFormatRec("VALUE2", 212, IBMFLT, "Value for Vertex #2", false);
   private StsSEGYFormatRec time03Rec = new StsSEGYFormatRec("TIME3", 216, IBMFLT, "Time for Vertex #3", false);
   private StsSEGYFormatRec value03Rec = new StsSEGYFormatRec("VALUE3", 220, IBMFLT, "Value for Vertex #3", false);
   private StsSEGYFormatRec time04Rec = new StsSEGYFormatRec("TIME4", 224, IBMFLT, "Time for Vertex #4", false);
   private StsSEGYFormatRec value04Rec = new StsSEGYFormatRec("VALUE4", 228, IBMFLT, "Value for Vertex #4", false);
   private StsSEGYFormatRec time05Rec = new StsSEGYFormatRec("TIME5", 232, IBMFLT, "Time for Vertex #5", false);
   private StsSEGYFormatRec value05Rec = new StsSEGYFormatRec("VALUE5", 236, IBMFLT, "Value for Vertex #5", false);
   private StsSEGYFormatRec time06Rec = new StsSEGYFormatRec("TIME6", 240, IBMFLT, "Time for Vertex #6", false);
   private StsSEGYFormatRec value06Rec = new StsSEGYFormatRec("VALUE6", 244, IBMFLT, "Value for Vertex #6", false);
   private StsSEGYFormatRec time07Rec = new StsSEGYFormatRec("TIME7", 248, IBMFLT, "Time for Vertex #7", false);
   private StsSEGYFormatRec value07Rec = new StsSEGYFormatRec("VALUE7", 252, IBMFLT, "Value for Vertex #7", false);
   private StsSEGYFormatRec time08Rec = new StsSEGYFormatRec("TIME8", 256, IBMFLT, "Time for Vertex #8", false);
   private StsSEGYFormatRec value08Rec = new StsSEGYFormatRec("VALUE8", 260, IBMFLT, "Value for Vertex #8", false);
   private StsSEGYFormatRec time09Rec = new StsSEGYFormatRec("TIME9", 264, IBMFLT, "Time for Vertex #9", false);
   private StsSEGYFormatRec value09Rec = new StsSEGYFormatRec("VALUE9", 268, IBMFLT, "Value for Vertex #9", false);
   private StsSEGYFormatRec time10Rec = new StsSEGYFormatRec("TIME10", 272, IBMFLT, "Time for Vertex #10", false);
   private StsSEGYFormatRec value10Rec = new StsSEGYFormatRec("VALUE10", 276, IBMFLT, "Value for Vertex #10", false);
   private StsSEGYFormatRec numVolsRec = new StsSEGYFormatRec("NUMVOLS", 280, IBMFLT, "Number of Velocity Volumes", false);
   private StsSEGYFormatRec velIncRec = new StsSEGYFormatRec("VELINC", 284, IBMFLT, "Velocity Increment between Volumes", false);
   private StsSEGYFormatRec numCDPRec = new StsSEGYFormatRec("NUMCDPS", 288, IBMFLT, "# CDPs Per PostStack3d", false);
   private StsSEGYFormatRec tsuTypeRec = new StsSEGYFormatRec("MVOLTYPE", 292, IBMFLT, "Multi-PostStack3d Attribute Type: 1-pwave, 2-swave 3-eta", false);

   // Remaining text header informaiton
   private StsSEGYFormatRec fmtRev1Rec = new StsSEGYFormatRec("FMTREV1", 300, BYTE, "SEG-Y Format major revision number", false);
   private StsSEGYFormatRec fmtRev2Rec = new StsSEGYFormatRec("FMTREV2", 301, BYTE, "SEG-Y Format minor revision number", false);
   private StsSEGYFormatRec varLenRec = new StsSEGYFormatRec("VARLEN", 302, INT2, "Fixed length/sample interval flag: 0=variable, 1=fixed", false);
   private StsSEGYFormatRec netfhRec = new StsSEGYFormatRec("NETFH", 304, INT2, "Number of external text file headers: -1=variable (scan for ((SEG:EndText)) stanza)", false);

   // trace records - set to definition
   private StsSEGYFormatRec lineTrcRec = new StsSEGYFormatRec(TRACENO, 0, INT4,
      "Trace number in line - unique across files", false);
   private StsSEGYFormatRec reelTrcRec = new StsSEGYFormatRec("REELTRC", 4, INT4, "Trace sequence number within reel", false);
   private StsSEGYFormatRec ffidRec = new StsSEGYFormatRec("FFID", 8, INT4, "Original field record number", false);
   private StsSEGYFormatRec chanRec = new StsSEGYFormatRec("CHAN", 12, INT4,
      "Trace sequence number within original field record", false);
   private StsSEGYFormatRec espNumRec = new StsSEGYFormatRec("ESPNUM", 16, INT4, "Energy source point number", false);
   private StsSEGYFormatRec cdpRec = new StsSEGYFormatRec(CDP, 20, INT4, "CDP ensemble number", true);
   private StsSEGYFormatRec seqNoRec = new StsSEGYFormatRec("SEQNO", 24, INT4,
      "Trace sequence number within CDP ensemble", false);
   private StsSEGYFormatRec traceIDRec = new StsSEGYFormatRec("TRACEID", 28, INT2, "Trace identification code: 1 = seismic data, 2 = dead, 3 = dummy, 4 = time break, 5 = uphole, 6 = sweep, 7 = timing, 8 = water break, 9+ = optional use", false);
   private StsSEGYFormatRec vStackRec = new StsSEGYFormatRec("VSTACK", 30, INT2,
      "Number of vertically summed traces yielding this trace", false);
   private StsSEGYFormatRec foldRec = new StsSEGYFormatRec("FOLD", 32, INT2,
      "Number of horizontally stacked traced yielding this trace", false);
   private StsSEGYFormatRec dataUseRec = new StsSEGYFormatRec("DATAUSE", 34, INT2,
      "Data use: 1 = production, 2 = test", false);
   private StsSEGYFormatRec offsetRec = new StsSEGYFormatRec("OFFSET", 36, INT4,
      "Distance from source point to receiver prefix", false);
   private StsSEGYFormatRec relevRec = new StsSEGYFormatRec("RELEV", 40, INT4, "Receiver prefix elevation", false, "ED-SCAL");
   private StsSEGYFormatRec selevRec = new StsSEGYFormatRec("SELEV", 44, INT4, "Surface elevation at source", false, "ED-SCAL");
   private StsSEGYFormatRec sDepthRec = new StsSEGYFormatRec("SDEPTH", 48, INT4, "Source depth below surface", false, "ED-SCAL");
   private StsSEGYFormatRec rDatumRec = new StsSEGYFormatRec("RDATUM", 52, INT4, "Datum elevation at receiver prefix", false, "ED-SCAL");
   private StsSEGYFormatRec sDatumRec = new StsSEGYFormatRec("SDATUM", 56, INT4, "Datum elevation at source", false, "ED-SCAL");
   private StsSEGYFormatRec wDepthSoRec = new StsSEGYFormatRec("WDEPTHSO", 60, INT4, "Water depth at source", false, "ED-SCAL");
   private StsSEGYFormatRec wDepthRcRec = new StsSEGYFormatRec("WDEPTHRC", 64, INT4, "Water depth at receiver prefix", false, "ED-SCAL");
   private StsSEGYFormatRec edScalRec = new StsSEGYFormatRec("ED-SCAL", 68, INT2,
      "Scalar for elevations and depths: + = multiplier, - = divisor", false);
   private StsSEGYFormatRec coScalRec = new StsSEGYFormatRec("CO-SCAL", 70, INT2,
      "Scalar for coordinates: + = multiplier, - = divisor", false);
   private StsSEGYFormatRec shtXRec = new StsSEGYFormatRec(SHT_X, 72, INT4, "X source coordinate", false, "CO-SCAL");
   private StsSEGYFormatRec shtYRec = new StsSEGYFormatRec(SHT_Y, 76, INT4, "Y source coordinate", false, "CO-SCAL");
   private StsSEGYFormatRec recXRec = new StsSEGYFormatRec(REC_X, 80, INT4, "X receiver prefix coordinate", false, "CO-SCAL");
   private StsSEGYFormatRec recYRec = new StsSEGYFormatRec(REC_Y, 84, INT4, "Y receiver prefix coordinate", false, "CO-SCAL");
   private StsSEGYFormatRec coorUnitRec = new StsSEGYFormatRec("COORUNIT", 88, INT2,
      "Coordinate units: 1 = length in meters or feet, 2 = arc seconds", false);
   private StsSEGYFormatRec wVelRec = new StsSEGYFormatRec("WVEL", 90, INT2, "Weathering velocity", false);
   private StsSEGYFormatRec subVelRec = new StsSEGYFormatRec("SUBVEL", 92, INT2, "Subweathering velocity", false);
   private StsSEGYFormatRec shUpHoleRec = new StsSEGYFormatRec("SHUPHOLE", 94, INT2, "Uphole time at source", false);
   private StsSEGYFormatRec rcUpHoleRec = new StsSEGYFormatRec("RCUPHOLE", 96, INT2, "Uphole time at receiver prefix", false);
   private StsSEGYFormatRec shStatRec = new StsSEGYFormatRec("SHSTAT", 98, INT2, "Source static correction", false);
   private StsSEGYFormatRec rcStatRec = new StsSEGYFormatRec("RCSTAT", 100, INT2, "Receiver prefix static correction", false);
   private StsSEGYFormatRec statApplyRec = new StsSEGYFormatRec("STAPPLY", 102, INT2, "Total static applied", false);
   private StsSEGYFormatRec lagTimeARec = new StsSEGYFormatRec("LAGTIMEA", 104, INT2,
      "Lag time between end of header and time break in milliseconds", false);
   private StsSEGYFormatRec lagTimeBRec = new StsSEGYFormatRec("LAGTIMEB", 106, INT2,
      "Lag time between time break and shot in milliseconds", false);
   private StsSEGYFormatRec delayRec = new StsSEGYFormatRec("DELAY", 108, INT2,
      "Lag time between shot and recording start in milliseconds", false);
   private StsSEGYFormatRec muteStrtRec = new StsSEGYFormatRec("MUTESTRT", 110, INT2, "Start of mute time", false);
   private StsSEGYFormatRec muteEndRec = new StsSEGYFormatRec("MUTEEND", 112, INT2, "End of mute time", false);
   private StsSEGYFormatRec nSampRec = new StsSEGYFormatRec("NSAMPLES", 114, INT2, "Number of samples in this trace", false);
   private StsSEGYFormatRec sRateRec = new StsSEGYFormatRec("SRATE", 116, INT2,
      "Sample interval of this trace in microseconds", false);
   private StsSEGYFormatRec gainTypeRec = new StsSEGYFormatRec("GAINTYPE", 118, INT2,
      "Field instrument gain type code: 1 = fixed, 2 = binary, 3 = floating point, 4+ = optional use", false);
   private StsSEGYFormatRec ingConstRec = new StsSEGYFormatRec("INGCONST", 120, INT2, "Instrument gain constant", false);
   private StsSEGYFormatRec initGainRec = new StsSEGYFormatRec("INITGAIN", 122, INT2,
      "Intrument early gain in decibels", false);
   private StsSEGYFormatRec corrFlagRec = new StsSEGYFormatRec("CORRFLAG", 124, INT2, "Correlated: 1 = no, 2 = yes", false);
   private StsSEGYFormatRec sweepSrtRec = new StsSEGYFormatRec("SWEEPSRT", 126, INT2, "Sweep frequency at start", false);
   private StsSEGYFormatRec sweepEndRec = new StsSEGYFormatRec("SWEEPEND", 128, INT2, "Sweep fequency at end", false);
   private StsSEGYFormatRec sweepLngRec = new StsSEGYFormatRec("SWEEPLNG", 130, INT2, "Sweep length in milliseconds", false);
   private StsSEGYFormatRec sweepTypRec = new StsSEGYFormatRec("SWEEPTYP", 132, INT2,
      "Sweep type code: 1 = linear, 2 = parabolic, 3 = exponential, 4 = other", false);
   private StsSEGYFormatRec sweepStpRec = new StsSEGYFormatRec("SWEEPSTP", 134, INT2,
      "Sweep taper trace length at start in milliseconds", false);
   private StsSEGYFormatRec sweepEtpRec = new StsSEGYFormatRec("SWEEPETP", 136, INT2,
      "Sweep taper trace length at end in milliseconds", false);
   private StsSEGYFormatRec taperTypRec = new StsSEGYFormatRec("TAPERTYP", 138, INT2,
      "Taper type code: 1 = linear, 2 = cosine squared, 3 = other", false);
   private StsSEGYFormatRec aliasFilRec = new StsSEGYFormatRec("ALIASFIL", 140, INT2, "Alias filter frequency", false);
   private StsSEGYFormatRec aliasSlopRec = new StsSEGYFormatRec("ALIASLOP", 142, INT2, "Alias filter slope", false);
   private StsSEGYFormatRec notchFilRec = new StsSEGYFormatRec("NOTCHFIL", 144, INT2, "Notch filter frequency", false);
   private StsSEGYFormatRec notchSlpRec = new StsSEGYFormatRec("NOTCHSLP", 146, INT2, "Notch filter slope", false);
   private StsSEGYFormatRec lowCutRec = new StsSEGYFormatRec("LOWCUT", 148, INT2, "Low cut frequency", false);
   private StsSEGYFormatRec highCutRec = new StsSEGYFormatRec("HIGHCUT", 150, INT2, "High cut frequency", false);
   private StsSEGYFormatRec lowCSlopRec = new StsSEGYFormatRec("LOWCSLOP", 152, INT2, "Low cut slope", false);
   private StsSEGYFormatRec hiCSlopRec = new StsSEGYFormatRec("HICSLOP", 154, INT2, "High cut slope", false);
   private StsSEGYFormatRec yearRec = new StsSEGYFormatRec("YEAR", 156, INT2, "Year data recorded", false);
   private StsSEGYFormatRec dayRec = new StsSEGYFormatRec("DAY", 158, INT2, "Day of year", false);
   private StsSEGYFormatRec hourRec = new StsSEGYFormatRec("HOUR", 160, INT2, "Hour of day: 24-hour clock", false);
   private StsSEGYFormatRec minuteRec = new StsSEGYFormatRec("MINUTE", 162, INT2, "Minute of hour", false);
   private StsSEGYFormatRec secondRec = new StsSEGYFormatRec("SECOND", 164, INT2, "Second of minute", false);
   private StsSEGYFormatRec timeBaseRec = new StsSEGYFormatRec("TIMEBASE", 166, INT2,
      "Time basis: 1 = local, 2 = GMT, 3 = other", false);
   private StsSEGYFormatRec trWeightRec = new StsSEGYFormatRec("TRWEIGHT", 168, INT2,
      "Trace weighting factor for fixed-point format data", false);
   private StsSEGYFormatRec rStaSwp1Rec = new StsSEGYFormatRec("RSTASWP1", 170, INT2,
      "Geophone prefix number of roll switch position one", false);
   private StsSEGYFormatRec rStaTrc1Rec = new StsSEGYFormatRec("RSTATRC1", 172, INT2,
      "Geophone prefix number of first trace of original field record", false);
   private StsSEGYFormatRec rStaTrcNRec = new StsSEGYFormatRec("RSTATRCN", 174, INT2,
      "Geophone prefix number of last trace of original field record", false);
   private StsSEGYFormatRec gapSizeRec = new StsSEGYFormatRec("GAPSIZE", 176, INT2,
      "Gap size: total number of groups dropped", false);
   private StsSEGYFormatRec overTrvlRec = new StsSEGYFormatRec("OVERTRVL", 178, INT2,
      "Overtravel associated with taper: 1 = down/behind, 2 = up/ahead", false);
   private StsSEGYFormatRec op01Rec = new StsSEGYFormatRec("OPTION01", 180, INT4, "Optional 4 byte field", false);
   private StsSEGYFormatRec op02Rec = new StsSEGYFormatRec("OPTION02", 184, INT4, "Optional 4 byte field", false);
   private StsSEGYFormatRec op03Rec = new StsSEGYFormatRec("OPTION03", 188, INT4, "Optional 4 byte field", false);
   private StsSEGYFormatRec op04Rec = new StsSEGYFormatRec("OPTION04", 192, INT4, "Optional 4 byte field", false);
   private StsSEGYFormatRec op05Rec = new StsSEGYFormatRec("OPTION05", 196, INT4, "Optional 4 byte field", false);
   private StsSEGYFormatRec op06Rec = new StsSEGYFormatRec("OPTION06", 200, INT4, "Optional 4 byte field", false);
   private StsSEGYFormatRec op07Rec = new StsSEGYFormatRec("OPTION07", 204, INT4, "Optional 4 byte field", false);
   private StsSEGYFormatRec op08Rec = new StsSEGYFormatRec("OPTION08", 208, INT4, "Optional 4 byte field", false);
   private StsSEGYFormatRec op09Rec = new StsSEGYFormatRec("OPTION09", 212, INT4, "Optional 4 byte field", false);
   private StsSEGYFormatRec op10Rec = new StsSEGYFormatRec("OPTION10", 216, INT4, "Optional 4 byte field", false);
   private StsSEGYFormatRec op11Rec = new StsSEGYFormatRec("OPTION11", 220, INT4, "Optional 4 byte field", false);
   private StsSEGYFormatRec op12Rec = new StsSEGYFormatRec("OPTION12", 224, INT4, "Optional 4 byte field", false);
   private StsSEGYFormatRec op13Rec = new StsSEGYFormatRec("OPTION13", 228, INT4, "Optional 4 byte field", false);
   private StsSEGYFormatRec op14Rec = new StsSEGYFormatRec("OPTION14", 232, INT4, "Optional 4 byte field", false);
   private StsSEGYFormatRec op15Rec = new StsSEGYFormatRec("CDPDATUM", 236, INT4, "Optional 4 byte field", false);

   static public final byte NONE = 0;
   static public final byte PRESTACK_RAW = 1;
   static public final byte PRESTACK_SORTED = 2;
   static public final byte POSTSTACK = 3;
   static public final byte VSP = 4;
   static public final byte POSTSTACK2D = 5;

   // All Binary Records
   private StsSEGYFormatRec[] allBinaryRecs =
      {
      jobIdNumberRec,
      lineNumberRec, reelNumberRec, tracesPerRec, auxsPerRec, sampleSpacingRec,
      origSampleSpacingRec, nSamplesRec, origNSamplesRec, sampleFormatRec,
      cdpFoldRec, sortCodeRec, vertSumCodeRec, startSweepFreqRec, endSweepFreqRec,
      sweepLengthRec, sweepTypeCodeRec, sweepTraceRec, startSweepTaperRec,
      endSweepTaperRec, taperTypeCodeRec, correlatedRec, binGainRec,
      ampRecovMtdRec, unitsRec, polarityRec, vibPolarityRec,
      time01Rec, value01Rec, time02Rec, value02Rec, time03Rec, value03Rec,
      time04Rec, value04Rec, time05Rec, value05Rec, time06Rec, value06Rec,
      time07Rec, value07Rec, time08Rec, value08Rec, time09Rec, value09Rec,
      time10Rec, value10Rec, numVolsRec, velIncRec, numCDPRec, tsuTypeRec,
      fmtRev1Rec, fmtRev2Rec, varLenRec, netfhRec};

   // All Trace Records
   private StsSEGYFormatRec[] allTraceRecs =
      {
      lineTrcRec, reelTrcRec, ffidRec,
      chanRec, espNumRec, cdpRec, seqNoRec, traceIDRec, vStackRec, foldRec,
      dataUseRec, offsetRec, relevRec, selevRec, sDepthRec, rDatumRec, sDatumRec,
      wDepthSoRec, wDepthRcRec, edScalRec, coScalRec, shtXRec, shtYRec, recXRec,
      recYRec, coorUnitRec, wVelRec, subVelRec, shUpHoleRec, rcUpHoleRec, shStatRec,
      rcStatRec, statApplyRec, lagTimeARec, lagTimeBRec, delayRec, muteStrtRec,
      muteEndRec, nSampRec, sRateRec, gainTypeRec, ingConstRec, initGainRec, corrFlagRec,
      sweepSrtRec, sweepEndRec, sweepLngRec, sweepTypRec, sweepStpRec, sweepEtpRec,
      taperTypRec, aliasFilRec, aliasSlopRec, notchFilRec, notchSlpRec, lowCutRec,
      highCutRec, lowCSlopRec, hiCSlopRec, yearRec, dayRec, hourRec, minuteRec,
      secondRec, timeBaseRec, trWeightRec, rStaSwp1Rec, rStaTrc1Rec, rStaTrcNRec,
      gapSizeRec, overTrvlRec, op01Rec, op02Rec, op03Rec, op04Rec, op05Rec,
      op06Rec, op07Rec, op08Rec, op09Rec, op10Rec, op11Rec, op12Rec, op13Rec, op14Rec,
      op15Rec};

   private StsSEGYFormatRec[] coorTraceRecs =
      {shtXRec, shtYRec, recXRec, recYRec};

   public transient byte units = METERS;

   public transient boolean applySCALCO = true;
   // JKF Changed modifier from static to transient on the following line
   private transient boolean formatChanged = false;

   // Word format
   static public final byte FLOAT = -2;
   static public final byte INT = -1;

   // Coordinate units
   static public final byte METERS = 1;
   static public final byte FEET = 2;

   static public final long serialVersionUID = 1L;

   // Unique global attribute IDs.
   // The names and values are arbitrary - I'm using SEG-Y 400 byte binary header locations
   //       static public final int SAMPINT = 16; //Sample interval; over-ridden by D1 if VARLEN is true
   //       static public final int NSAMP = 20;   //# samples/trace; over-ridden by N1 if VARLEN is true
   //       static public final int SAMPFMT = 24; //Sample format; value depends of file format spec.
   //       static public final int UNITS = 54;   //Coordinate units (FEET or METERS)
   static public final int SAMPINT = 0; //Sample interval; over-ridden by D1 if VARLEN is true
   static public final int NSAMP = 1; //# samples/trace; over-ridden by N1 if VARLEN is true
   static public final int SAMPFMT = 1; //Sample format; value depends of file format spec.
   static public final int UNITS = 3; //Coordinate units (FEET or METERS)
   // SEG-Y specific global attribute IDs.  :
   //       static public final int FMTREV1 = 300;
   //       static public final int FMTREV2 = 301;
   //       static public final int VARLEN = 302;  //Variable trace length/sample interval flag
   //       static public final int NETFH = 304;
   static public final int FMTREV1 = 4;
   static public final int FMTREV2 = 5;
   static public final int VARLEN = 6; //Variable trace length/sample interval flag
   static public final int NETFH = 7;

   static public final String XLINE_NO = "XLINE_NO";
   static public final String ILINE_NO = "ILINE_NO";
   static public final String CDP_X = "CDP_X";
   static public final String CDP_Y = "CDP_Y";
   static public final String N1 = "N1";
   static public final String D1 = "D1";
   static public final String TRC_TYPE = "TRC_TYPE";
   static public final String TRACENO = "TRACENO";
   static public final String SCALCO = "SCALCO";
   static public final String SCALEL = "SCALEL";
   static public final String CDP = "CDP";
   static public final String DATASCAL = "DATASCAL";
   static public final String SHT_X = "SHT-X";
   static public final String SHT_Y = "SHT-Y";
   static public final String REC_X = "REC-X";
   static public final String REC_Y = "REC-Y";
   static public final String OFFSET = "OFFSET";
   // SEG-Y sample format codes:
   static public final byte IBMFLT = 1;
   static public final byte INT4 = 2;
   static public final byte INT2 = 3;
   static public final byte FIXED = 4;
   static public final byte IEEEFLT = 5;
   static public final byte FLOAT8 = 6;
   static public final byte FLOAT16 = 7;
   static public final byte BYTE = 8;

   static public final String NONE_STRING = "none/mixed";
   static public final String IBMFLT_STRING = "IBM float";
   static public final String INT4_STRING = "4 byte integer";
   static public final String INT2_STRING = "2 byte integer";
   static public final String FIXED_STRING = "Fixed";
   static public final String IEEEFLT_STRING = "IEEE float";
   static public final String FLOAT8_STRING = "1 byte float";
   static public final String FLOAT16_STRING = "2 byte float";
   static public final String BYTE_STRING = "1 byte integer";

   static public final int maxNSampleFormats = 9;

   static final public String[] sampleFormatStrings = new String[]
   {
      NONE_STRING, IBMFLT_STRING, INT4_STRING, INT2_STRING, FIXED_STRING, IEEEFLT_STRING, FLOAT8_STRING, FLOAT16_STRING, BYTE_STRING
   };

   static final public String[] headerFormatStrings = new String[]
       { IBMFLT_STRING, INT4_STRING, INT2_STRING, IEEEFLT_STRING };

   static public final byte NOTMULTIVOLUME = 0;
   static public final byte PWAVE = 1;
   static public final byte SWAVE = 2;
   static public final byte ETA = 3;
   static public final String[] multiVolumeTypes = new String[] {"None", "P-Wave", "S-Wave", "ETA" };

   static final public String[] distAttributes = new String[] { "OFFSET", "RELEV", "SELEV", "SDEPTH", "CDPDATUM", "RDATUM", "SDATUM", "WDEPTHSO", "WDEPTHRC" };
   static final byte DIST = 0;

   static final public String[] coorAttributes = new String[]
      { SHT_X, SHT_Y, REC_X, REC_Y };
    static final byte COOR = 1;

   static final public String[] timeAttributes = new String[]
      { "SUPHOLE", "RCUPHOLE", "SHSTAT", "RCSTAT", "STAPPLY", "LAGTIMEA", "LAGTIMEB", "MUTESTRT",
      "MUTEEND", "DELAY", "SWEEPSRT", "SWEEPEND", "SWEEPLEN", "FRSTBRK" };
    static final byte TIME = 2;

   static final public String[] datumAttributes = new String[] { "CDPDATUM", "RDATUM", "SDATUM" };
   static final byte DATUM = 3;

   static final public String[] velocityAttributes = new String[] { "WVEL", "SVEL" };
   static final byte VELOCITY = 4;
   static final public String[] scalarAttributes = new String[] { "NONE", "ED-SCAL", "CO-SCAL" };
   static final byte SCALAR = 5;
   static final byte OTHER = 6;

   // Strings displayed in dropdown list
   static public final String ASCII_7BIT = "7 bit ASCII";
   static public final String ASCII_8BIT = "8 bit ASCII";
   static public final String EBCDIC = "EBCDIC";
   static public final String[] textHeaderFormatStrings = new String[] {ASCII_7BIT, ASCII_8BIT, EBCDIC};
   // actual strings required by decoder
   static public final String[] textHeaderFormats = new String[]
      {"US-ASCII", "ISO-8859-1", "Cp500"};

   static public final String group = "segyFormat";
   static public final String format = "obj";

   static private final SegyFormatFilenameFilter segyFormatFilenameFilter = new SegyFormatFilenameFilter(group, format);

   public StsSEGYFormat()
   {
   }

   private StsSEGYFormat(StsModel model, byte type) throws IllegalArgumentException
   {
      this.dataType = type;
      initializeUnits(model);
      switch (type)
      {
         case PRESTACK_RAW:

            // Specific setup for raw shot ordered prestack
            this.getTraceRec(TRACENO).setRequired(true);
            this.getTraceRec("FFID").setUserName(ILINE_NO);
            this.getTraceRec("FFID").setRequired(true);
            this.getTraceRec("CHAN").setUserName(XLINE_NO);
            this.getTraceRec("CHAN").setRequired(true);
            this.getTraceRec("FOLD").setRequired(true);
            this.getTraceRec("TRACEID").setUserName(TRC_TYPE);
            this.getTraceRec("TRACEID").setRequired(true);
            this.getTraceRec("CO-SCAL").setUserName(SCALCO);

//				this.getTraceRec("CO-SCAL").setRequired(true);
            this.getTraceRec("WVEL").setRequired(true);
            this.getTraceRec(SHT_X).setRequired(true);
            this.getTraceRec(SHT_Y).setRequired(true);
            this.getTraceRec(REC_X).setRequired(true);
            this.getTraceRec(REC_Y).setRequired(true);
            this.getTraceRec("NSAMPLES").setUserName(N1);
            this.getTraceRec("NSAMPLES").setRequired(true);
            this.getTraceRec("SRATE").setUserName(D1);
            this.getTraceRec("SRATE").setRequired(true);
            this.getTraceRec("TRWEIGHT").setUserName("B2SCALCO");

//				this.getTraceRec("TRWEIGHT").setRequired(true);
            this.getTraceRec("OPTION05").setUserName(CDP_X);
            this.getTraceRec("OPTION05").setApplyScalar("CO-SCAL");
            this.getTraceRec("OPTION05").setRequired(true);
            this.getTraceRec("OPTION07").setUserName(CDP_Y);
            this.getTraceRec("OPTION07").setApplyScalar("CO-SCAL");
            this.getTraceRec("OPTION07").setRequired(true);
            this.getTraceRec("OFFSET").setRequired(true);
            break;
         case PRESTACK_SORTED:

            // Specific setup for ensemble sorted prestack
            this.getTraceRec(TRACENO).setRequired(true);
            this.getTraceRec("FFID").setUserName(ILINE_NO);
            this.getTraceRec("FFID").setRequired(true);
            this.getTraceRec("CHAN").setUserName(XLINE_NO);
            this.getTraceRec("CHAN").setRequired(true);
            this.getTraceRec("FOLD").setRequired(true);
            this.getTraceRec("TRACEID").setUserName(TRC_TYPE);
            this.getTraceRec("TRACEID").setRequired(true);
            this.getTraceRec("CO-SCAL").setUserName(SCALCO);

//				this.getTraceRec("CO-SCAL").setRequired(true);
            this.getTraceRec("WVEL").setRequired(true);
            this.getTraceRec(SHT_X).setRequired(true);
            this.getTraceRec(SHT_Y).setRequired(true);
            this.getTraceRec(REC_X).setRequired(true);
            this.getTraceRec(REC_Y).setRequired(true);
            this.getTraceRec("NSAMPLES").setUserName(N1);
            this.getTraceRec("NSAMPLES").setRequired(true);
            this.getTraceRec("SRATE").setUserName(D1);
            this.getTraceRec("SRATE").setRequired(true);
            this.getTraceRec("TRWEIGHT").setUserName("B2SCALCO");

//				this.getTraceRec("TRWEIGHT").setRequired(true);
            this.getTraceRec("OPTION05").setUserName(CDP_X);
            this.getTraceRec("OPTION05").setRequired(true);
            this.getTraceRec("OPTION05").setApplyScalar("CO-SCAL");
            this.getTraceRec("OPTION07").setUserName(CDP_Y);
            this.getTraceRec("OPTION07").setRequired(true);
            this.getTraceRec("OPTION07").setApplyScalar("CO-SCAL");
            this.getTraceRec("OFFSET").setRequired(true);
            break;
		 case POSTSTACK2D:
         case POSTSTACK:
            this.getTraceRec(TRACENO).setRequired(true);
            this.getTraceRec("FFID").setUserName(ILINE_NO);
            this.getTraceRec("FFID").setRequired(true);
            this.getTraceRec("CHAN").setUserName(XLINE_NO);
            this.getTraceRec("CHAN").setRequired(true);
            this.getTraceRec("FOLD").setRequired(true);
            this.getTraceRec("STAPPLY").setRequired(true);
            this.getTraceRec("TRACEID").setUserName(TRC_TYPE);
            this.getTraceRec("TRACEID").setRequired(true);
            this.getTraceRec("CO-SCAL").setUserName(SCALCO);

            //			this.getTraceRec("CO-SCAL").setRequired(true);
            this.getTraceRec("WVEL").setRequired(true);
            this.getTraceRec(SHT_X).setUserName(CDP_X);
            this.getTraceRec(SHT_X).setRequired(true);
            this.getTraceRec(SHT_Y).setUserName(CDP_Y);
            this.getTraceRec(SHT_Y).setRequired(true);
            this.getTraceRec("NSAMPLES").setUserName(N1);
            this.getTraceRec("NSAMPLES").setRequired(true);
            this.getTraceRec("SRATE").setUserName(D1);
            this.getTraceRec("SRATE").setRequired(true);
            this.getTraceRec("TRWEIGHT").setUserName("B2SCALCO");

//				this.getTraceRec("TRWEIGHT").setRequired(true);
            break;
        case VSP:
           // Specific setup for vsp post and pre-stack data
            this.getTraceRec(TRACENO).setRequired(true); // 1
            this.getTraceRec("REELTRC").setRequired(true); // 4
            this.getTraceRec("FFID").setRequired(true);    // 9
            this.getTraceRec("CHAN").setRequired(true);    // 13
            this.getTraceRec("ESPNUM").setRequired(true);  // 17
            this.getTraceRec("TRACEID").setRequired(true); // 29
            this.getTraceRec("TRACEID").setUserName(TRC_TYPE);
            this.getTraceRec("VSTACK").setRequired(true);  // 31
            this.getTraceRec("FOLD").setRequired(true);    // 33
            this.getTraceRec("DATAUSE").setRequired(true); // 35
            this.getTraceRec("OFFSET").setRequired(true);  // 37
            this.getTraceRec("RELEV").setRequired(true);   // 41
            this.getTraceRec("SDEPTH").setRequired(true);  // 49
            this.getTraceRec("RDATUM").setRequired(true);  // 53
            this.getTraceRec("WDEPTHSO").setRequired(true);// 61
            this.getTraceRec("WDEPTHRC").setRequired(true);// 65
            this.getTraceRec("ED-SCAL").setRequired(true); // 69
            this.getTraceRec("CO-SCAL").setRequired(true); // 71
            this.getTraceRec("CO-SCAL").setUserName(SCALCO);
            this.getTraceRec(SHT_X).setRequired(true);   // 73
            this.getTraceRec(SHT_Y).setRequired(true);   // 77
            this.getTraceRec(REC_X).setRequired(true);   // 81
            this.getTraceRec(REC_Y).setRequired(true);   // 85
            this.getTraceRec("WVEL").setRequired(true);    // 91
            this.getTraceRec("SUBVEL").setRequired(true);  // 93
            this.getTraceRec("SHUPHOLE").setRequired(true);  // 95
            this.getTraceRec("RCUPHOLE").setRequired(true);  // 97
            this.getTraceRec("OPTION03").setUserName("CDP_X_ARC");
            this.getTraceRec("OPTION03").setRequired(true);  // 189
            this.getTraceRec("OPTION04").setUserName("CDP_Y_ARC");
            this.getTraceRec("OPTION04").setRequired(true);  // 193
            this.getTraceRec("OPTION14").setUserName("FRSTBRK");
            this.getTraceRec("OPTION14").setFormat(INT2);
            this.getTraceRec("OPTION14").setRequired(true);  // 233
            this.getTraceRec("TRWEIGHT").setUserName("B2SCALCO");
            break;
         default:
            new StsMessage(null, StsMessage.WARNING, "StsSEGYFormat.constructTraceAnalyzer called with illegal type: " + type);
            throw new IllegalArgumentException("StsSEGYFormat.constructTraceAnalyzer called with illegal type: " + type);
      }
   }

   static public StsSEGYFormat constructor(StsModel model, byte type)
   {
      try
      {
         return new StsSEGYFormat(model, type);
      }
      catch (Exception e)
      {
         return null;
      }
   }

   static public StsSEGYFormat constructor(StsModel model, String stemname, byte type)
   {
      if (stemname == null)
         return new StsSEGYFormat(model, type);

      String directory = getSegyFormatsDirectory();
      try
      {
         String filename = getFilenameFromStemname(stemname);
         StsSEGYFormat segyFormat = readSegyFormat(directory, filename);
         if (segyFormat == null) segyFormat = new StsSEGYFormat(model, type);
         return segyFormat;
      }
      catch (Exception e)
      {
         StsException.outputException("Failed to read SegyFormat file:", e, StsException.WARNING);
         return constructor(model, type);

      }
   }

   static public byte getFilenameType(String stemname)
   {
      String directory = getSegyFormatsDirectory();
      try
      {
         String filename = getFilenameFromStemname(stemname);
         StsSEGYFormat segyFormat = readSegyFormat(directory, filename);
         if (segyFormat == null) return NONE;
         return segyFormat.dataType;
      }
      catch (Exception e)
      {
         StsException.outputException("Failed to read SegyFormat file:", e, StsException.WARNING);
         return NONE;
      }
   }

   static public StsSEGYFormat readSegyFormat(String directory, String filename)
   {
       String pathname = directory + File.separator + filename;
       StsSEGYFormat segyFormat = new StsSEGYFormat();
       if(!StsDBFileObject.readObjectFile(pathname, segyFormat, null))
       {
           File file = new File(pathname);
           file.delete();
           return null;
       }
       else
           return segyFormat;
   }

   private void initializeUnits(StsModel model)
   {
//	   model.getProject().toggleZDomain(zDomain);
//      zDomain = model.getProject().getZDomain();
      hUnits = model.getProject().getXyUnits();
      hUnitString = StsParameters.DIST_STRINGS[hUnits];
      zUnits = model.getProject().getDepthUnits();
      zUnitString = StsParameters.DIST_STRINGS[zUnits];
      tUnits = model.getProject().getTimeUnits();
      tUnitString = StsParameters.TIME_STRINGS[tUnits];
   }

   static public String getFilenameFromStemname(String stemname)
   {
      return getFilePrefix() + stemname;
   }

   static private StsSEGYFormatRec[] copyRecs(StsSEGYFormatRec[] staticRecs)
   {
      try
      {
         int nRecs = staticRecs.length;
         StsSEGYFormatRec[] outputRecs = new StsSEGYFormatRec[nRecs];
         for (int n = 0; n < nRecs; n++)
            outputRecs[n] = new StsSEGYFormatRec(staticRecs[n]);
         return outputRecs;
      }
      catch (Exception e)
      {
         StsException.systemError("StsSEGYFormat.copyRecs() failed");
         return null;
      }
   }

   public void switchFromPrestackRawToSorted()
   {
      // do what is needed to change PRESTACK_RAW records to PRESTACK_SORTED
   }

   /*
      public void analyzeBinaryHeader(StsSegyVolume segyVolume)
      {
     byte[] hdr = segyVolume.getBinaryHeader();
     boolean getIsLittleEndian = segyVolume.getIsLittleEndian();
     segyVolume.setSampleSpacing((float)getBinaryHdrValue("SAMPINT", hdr, getIsLittleEndian)/1000);
     segyVolume.setNSamples((int)getBinaryHdrValue("NSAMP", hdr, getIsLittleEndian));
//        segyVolume.setSampleFormat((int)getBinaryHdrValue("SAMPFMT", hdr, getIsLittleEndian));
    }
    */
   public float getSampleSpacing(boolean isLittleEndian, byte[] binaryHeader)
   {
      sampleSpacing = (float)getBinaryHdrValue("SAMPINT", binaryHeader, isLittleEndian) / 1000;
      return sampleSpacing;
   }

   public int getNSamples(boolean isLittleEndian, byte[] binaryHeader)
   {
      nSamp = (int)getBinaryHdrValue("NSAMP", binaryHeader, isLittleEndian);
      return nSamp;
   }

   public byte getMultiVolumeType(boolean isLittleEndian, byte[] binaryHeader)
   {
      return (byte)getBinaryHdrValue("MVOLTYPE", binaryHeader, isLittleEndian);
   }

   public boolean isMultiVolume(boolean isLittleEndian, byte[] binaryHeader)
   {
       byte type = (byte)getBinaryHdrValue("MVOLTYPE", binaryHeader, isLittleEndian);
       if((type < PWAVE) || (type > ETA))
           return false;
       else
           return true;
   }

   public void setName(String name)
   {
      this.name = name;
   }

   public String getName()
   {
      return name;
   }

/*
   public void setSampleFormatOverride(int override)
   {
       sampleFormatOverride = override;
       formatChanged = true;
   }

   public int getSampleFormatOverride()
   {
       return sampleFormatOverride;
   }

   public void setSampleFormatString(String sampleFormatString)
   {
      int newSampleFormatOverride = getSampleFormatFromString(sampleFormatString);
      if (newSampleFormatOverride == sampleFormatOverride)
         return;
      sampleFormatOverride = newSampleFormatOverride;
      formatChanged = true;
   }
*/
   public void setZDomainString(String zDomainString)
   {
      setZDomain(StsParameters.getZDomainFromString(zDomainString));
   }

    public void setZDomain(byte zDomain)
    {
        if(this.zDomain == zDomain) return;
        this.zDomain = zDomain;
        setFormatChanged(true);        
    }

   public String getZDomainString()
   {
      return StsParameters.TD_STRINGS[zDomain];
   }

   public StsSEGYFormatRec[] getBinaryRecords()
   {
      return allBinaryRecs;
   }

   public StsSEGYFormatRec[] getAllBinaryRecords()
   {
      return allBinaryRecs;
   }

   public StsSEGYFormatRec[] getAllTraceRecords()
   {
      return allTraceRecs;
   }

	public StsSEGYFormatRec[] getRequiredTraceRecords()
	{
        int nTotalRecords = allTraceRecs.length;
        StsSEGYFormatRec[] requiredRecords = new StsSEGYFormatRec[nTotalRecords];
		int nReqRecords = 0;
		for (int n = 0; n < nTotalRecords; n++)
		{
			if (allTraceRecs[n].isRequired())
				requiredRecords[nReqRecords++] = allTraceRecs[n];
		}
        if(nReqRecords < nTotalRecords)
            requiredRecords = (StsSEGYFormatRec[])StsMath.trimArray(requiredRecords, nReqRecords);
        return requiredRecords;
	}
/*
   public String getSampleFormatString()
   {
      return sampleFormatStrings[sampleFormatOverride];
   }
*/
	public byte getSampleFormat(boolean isLittleEndian, byte[] binaryHeader)
	{
		if (isSampleFormatOverride)
            return overrideSampleFormat;
		byte sampleFormat = (byte)getBinaryHdrValue("SAMPFMT", binaryHeader, isLittleEndian);
        if(sampleFormat >= 1 && sampleFormat <= 8)
            return sampleFormat;
        else
            return defaultSampleFormat;
    }

	public int getBinaryHeaderSampleFormat(boolean isLittleEndian, byte[] binaryHeader)
	{
		return (int)getBinaryHdrValue("SAMPFMT", binaryHeader, isLittleEndian);
	}

   static public int getSampleFormatFromString(String sampleFormatString)
   {
      for (int n = 1; n < sampleFormatStrings.length; n++)
      {
         if (sampleFormatStrings[n].equals(sampleFormatString))
            return n;
      }
      return 0;
   }

   public void setCoordinateFormatString(String sampleFormatString)
   {
      setCoordinteFormat(getSampleFormatFromString(sampleFormatString));
   }

   public void setCoordinteFormat(int format)
   {
      if (format == 0) return;
      StsSEGYFormatRec rec = getTraceRecFromUserName(CDP_X);
      if(rec.getFormat() == format) return;
      rec.setFormat(format);
      setFormatChanged(true);
   }

   public StsSEGYFormatRec getTraceRec(String name)
   {
      for (int n = 0; n < allTraceRecs.length; n++)
         if (allTraceRecs[n].name.equals(name))
            return allTraceRecs[n];
      return null;
   }

   public int getTraceRecIndex(String name)
   {
      for (int n = 0; n < allTraceRecs.length; n++)
         if (allTraceRecs[n].name.equals(name))
            return n;
      return -1;
   }

   public int getNumberOfRequiredRecords()
   {
       int count = 0;
       for (int n = 0; n < allTraceRecs.length; n++)
       {
           if(allTraceRecs[n].required || allTraceRecs[n].userRequired)
               count++;
       }
       return count;
   }

	public StsSEGYFormatRec getTraceRecFromUserName(String name)
	{
		for (int n = 0; n < allTraceRecs.length; n++)
			if (allTraceRecs[n].userName.equals(name))
				return allTraceRecs[n];
		return null;
	}

   public StsSEGYFormatRec getBinaryRec(String name)
   {
      for (int n = 0; n < allBinaryRecs.length; n++)
         if (allBinaryRecs[n].name.equals(name))
            return allBinaryRecs[n];
      return null;
   }

   public StsSEGYFormatRec getBinaryRecFromUserName(String name)
   {
      for (int n = 0; n < allBinaryRecs.length; n++)
         if (allBinaryRecs[n].userName.equals(name))
            return allBinaryRecs[n];
      return null;
   }

   public void setFormatChanged(boolean value)
   {
       formatChanged = value;
   }
    
   public String getCoordinateFormatString()
   {
      int format = getTraceRecFromUserName(CDP_X).getFormat();
      return sampleFormatStrings[format];
   }

   public void setInlineFormatString(String sampleFormatString)
   {
      setInlineFormat(getSampleFormatFromString(sampleFormatString));
   }

   public void setInlineFormat(int format)
   {
      if (format == 0) return;
      StsSEGYFormatRec rec = getTraceRecFromUserName(ILINE_NO);
      if(rec.getFormat() == format) return;
      rec.setFormat(format);
      setFormatChanged(true);
   }

   public String getInlineFormatString()
   {
      int format = getTraceRecFromUserName(ILINE_NO).getFormat();
      return sampleFormatStrings[format];
   }

   public void setXLineFormatString(String sampleFormatString)
   {
      setXLineFormat(getSampleFormatFromString(sampleFormatString));
   }

   public void setXLineFormat(int format)
   {
      if (format == 0) return;
      StsSEGYFormatRec rec = getTraceRecFromUserName(XLINE_NO);
       if(format == rec.format) return;
      rec.setFormat(format);
      setFormatChanged(true);
   }

   public String getXLineFormatString()
   {
      int format = getTraceRecFromUserName(XLINE_NO).getFormat();
      return sampleFormatStrings[format];
   }

   public TraceHeader constructTraceHeader(byte[] hdr, int nTrace, boolean isLittleEndian)
   {
      return new TraceHeader(hdr, nTrace, isLittleEndian);
   }

   public TraceHeader constructTraceHeader(byte[] hdr, int offset, int nTrace, boolean isLittleEndian)
   {
      return new TraceHeader(hdr, offset, nTrace, isLittleEndian);
   }

   final public float getXLine(byte[] hdr, boolean isLittleEndian)
   {
      return (float)getTraceRecFromUserName(XLINE_NO).getHdrValue(hdr, 0, isLittleEndian,
         getCoordinateScale(hdr, isLittleEndian), getElevationScale(hdr, isLittleEndian));
   }

   final public int getDataScale(byte[] hdr, boolean isLittleEndian)
   {
      return (int)getTraceRecFromUserName("B2SCALCO").getHdrValue(hdr, 0, isLittleEndian,
         getCoordinateScale(hdr, isLittleEndian), getElevationScale(hdr, isLittleEndian));
   }

   final public int getCoordinateScale(byte[] hdr, boolean isLittleEndian)
   {
      return (int)getTraceRec("CO-SCAL").getHdrValue(hdr, isLittleEndian);
   }

   final public float getFloatCoordinateScale(byte[] hdr, boolean isLittleEndian)
   {
       return getFloatScale(getCoordinateScale(hdr, isLittleEndian));
   }

   final public float getFloatScale(int scale)
   {
       if(scale < 0) return -1.0f / scale;
	   else if(scale > 0) return (float)scale;
       else return 1.0f;
   }

   final public int getElevationScale(byte[] hdr, boolean isLittleEndian)
   {
      return (int)getTraceRecFromUserName("ED-SCAL").getHdrValue(hdr, isLittleEndian);
   }

   final public float getFloatElevationScale(byte[] hdr, boolean isLittleEndian)
   {
       return getFloatScale(getElevationScale(hdr, isLittleEndian));
   }

   final public int getTraceId(byte[] hdr, boolean isLittleEndian)
   {
      return (int)getTraceRecFromUserName(TRC_TYPE).getHdrValue(hdr, isLittleEndian);
   }

   final public double getBinaryHdrValue(String name, byte[] hdr, boolean isLittleEndian)
   {
      StsSEGYFormatRec record = getBinaryRecFromUserName(name);
      if (record == null)
         return (double)StsParameters.doubleNullValue;
      return record.getHdrValue(hdr, isLittleEndian);
   }

   final public double getTraceHdrValue(String name, byte[] hdr, boolean isLittleEndian)
   {
      return getTraceRecFromUserName(name).getHdrValue(hdr, 0, isLittleEndian, getCoordinateScale(hdr, isLittleEndian), getElevationScale(hdr, isLittleEndian));
   }

   public TraceHeader constructTraceHeader()
   {
      return new TraceHeader();
   }

    public boolean isFormatChanged()
    {
        return formatChanged;
    }

    public class TraceHeader
   {
      public float iLine;
      public float xLine;
      public float x = nullValue;
      public float y = nullValue;
      public int xyScale = -99;
      public int edScale = -99;
      public int dataScale = -99;
      public int nTrace; // index of trace in volume
      public boolean xyFromCDP = true;
	  public int cdp;

      public TraceHeader()
      {
      }

      public TraceHeader(byte[] hdr, int nTrace, boolean littleEndian)
      {
         this(hdr, 0, nTrace, littleEndian);
      }

      public TraceHeader(byte[] hdr, int offset, int nTrace, boolean littleEndian)
      {
         this.nTrace = nTrace;

         xyScale = getCoordinateScale(hdr, littleEndian);
         edScale = getElevationScale(hdr, littleEndian);
/*  this scale check is handled in get__Scale methods above
         if (xyScale == 0)
            xyScale = 1;
        if (edScale == 0)
            edScale = 1;
*/
         if(dataType != VSP)
         {
             iLine = (float)getTraceRecFromUserName(ILINE_NO).getHdrValue(hdr, offset, littleEndian, xyScale, edScale);
             xLine = (float)getTraceRecFromUserName(XLINE_NO).getHdrValue(hdr, offset, littleEndian, xyScale, edScale);
         }
         else
         {
             iLine = (float) getTraceRecFromUserName("FFID").getHdrValue(hdr, offset, littleEndian, xyScale, edScale);
             xLine = (float) getTraceRecFromUserName("CHAN").getHdrValue(hdr, offset, littleEndian, xyScale, edScale);
         }
		 cdp = (int) getTraceRecFromUserName(CDP).getHdrValue(hdr, offset, littleEndian, xyScale, edScale);
//            System.out.println("TraceNum=" + nTrace + "Inline=" + iLine + " Xline=" + xLine);

         // If X & Y values are valid................Else compute X,Y
         StsSEGYFormatRec cdpxRec = getTraceRecFromUserName(CDP_X);
         StsSEGYFormatRec cdpyRec = getTraceRecFromUserName(CDP_Y);
         double xi, yi;
         if (cdpxRec != null && cdpyRec != null)
         {
            x = (float)cdpxRec.getHdrValue(hdr, offset, littleEndian, xyScale, edScale);
            y = (float)cdpyRec.getHdrValue(hdr, offset, littleEndian, xyScale, edScale);
            if((x == 0.0) && (y == 0.0))
                xyFromCDP = false;
            else
                xyFromCDP = true;
         }
         else
             xyFromCDP = false;

         if(!xyFromCDP)
         {
            if (dataType == POSTSTACK)
               return;
            cdpxRec = getTraceRecFromUserName(REC_X);
            if (cdpxRec == null)
               return;
            StsSEGYFormatRec shtxRec = getTraceRecFromUserName(SHT_X);
            if (shtxRec == null)
               return;
            double recx = cdpxRec.getHdrValue(hdr, offset, littleEndian, xyScale, edScale);
            double shtx = shtxRec.getHdrValue(hdr, offset, littleEndian, xyScale, edScale);
//				xi = (recx + shtx) / 2;
            x = (float)((double)(recx + shtx) / (double)2.0);
            cdpyRec = getTraceRecFromUserName(REC_Y);
            if (cdpyRec == null)
               return;
            StsSEGYFormatRec shtyRec = getTraceRecFromUserName(SHT_Y);
            if (shtyRec == null)
               return;
            double recy = cdpyRec.getHdrValue(hdr, offset, littleEndian, xyScale, edScale);
            double shty = shtyRec.getHdrValue(hdr, offset, littleEndian, xyScale, edScale);
//				yi = (recy + shty) / 2;
            y = (float)((double)(recy + shty) / (double)2.0);
         }

         StsSEGYFormatRec scaleRec = getTraceRecFromUserName("B2SCALCO");
         if (scaleRec == null)
            return;
         dataScale = (int)scaleRec.getHdrValue(hdr, offset, littleEndian, xyScale, edScale);
      }

      public boolean isCdpXyOk()
      {
         return x != nullValue;
      }

      public Comparator getXLineComparator()
      {
         return new XLineComparator();
      }

      public class XLineComparator
         implements Comparator
      {
         public XLineComparator()
         {
         }

         public int compare(Object object, Object otherObject)
         {
            float xLine = ((TraceHeader)object).xLine;
            float otherXLine = ((TraceHeader)otherObject).xLine;
            if (xLine > otherXLine)
               return 1;
            else if (xLine < otherXLine)
               return -1;
            else
               return 0;
         }

         public boolean equals(Object otherObject)
         {
            float otherXLine = ((TraceHeader)otherObject).xLine;
            return StsMath.sameAs(xLine, otherXLine);
         }
      }

      public String toString()
      {
         return new String("line: " + iLine + " xline: " + xLine + " x: " + x + " y: " + y + " trace index: " + nTrace);
      }
   }

   public void write(PrintStream out)
   {
      out.println("name: " + name);
      out.println("description: " + description);
      out.println("binaryHeaderSize: " + binaryHeaderSize);
      out.println("traceHeaderSize: " + traceHeaderSize);
      out.println("txtHdrFmt: " + textHeaderFormatString);
      out.println("textHeaderSize: " + textHeaderSize);
//        out.println("netfh: " + netfh);
//        out.println("varlen: " + varlen);
//        out.println("fmtrev1: " + fmtrev1);
//        out.println("fmtrev2: " + fmtrev2);
      out.println("units: " + units);
//        out.println("sampfmt: " + sampleFormat);
      out.println("//Global format records:");
      StsSEGYFormatRec rec = null;
   }

   static public String getFilePrefix()
   {
      return group + "." + format + ".";
   }

   static private String getSegyFormatsDirectory()
   {
      return System.getProperty("user.home") + File.separator + "S2SCache" + File.separator;
   }
   public boolean serialize(StsModel model)
   {
      try
      {
         String directory = getSegyFormatsDirectory();
         File directoryFile = new File(directory);
         if (!directoryFile.exists() && !directoryFile.mkdir())
         {
            new StsMessage(model.win3d, StsMessage.WARNING, "Failed to create directory " + directory);
            return false;
         }

         String[] names = getSegyFormatFilenames();
         StsNameSelectionDialog dialog = new StsNameSelectionDialog(model.win3d,
                                                                    "Save SEGY Format",
                                                                    "Format Name",
                                                                    names,
                                                                    name,
                                                                    StsNameSelectionDialog.QUESTION_OVERWRITE,
                                                                    StsNameSelectionDialog.ALLOW_CANCEL);


         String filePrefix = getFilePrefix();

         dialog.setVisible(true);
         if (dialog.wasCanceled())
         {
            return false;
         }

         name = dialog.getSelectedName();
		 String filePathname = directory + File.separator + filePrefix + name;
		 File f = new File(filePathname);
		 if (f.exists())
		 {
			f.delete();
		 }
		 serializeFormat(filePathname);
		 return true;
	  }
	  catch (Exception e)
	  {
		 StsException.outputException("StsSEGYFormat.serialize() failed.", e, StsException.WARNING);
		 return false;
	  }
   }

   public boolean serializeSilently(String name)
   {
	  try
	  {
		  this.name = name;
		 String directory = getSegyFormatsDirectory();
		 File directoryFile = new File(directory);
		 if (!directoryFile.exists() && !directoryFile.mkdir())
		 {
			throw new Exception("Failed to create directory " + directory + ".");
		 }

		 String filePrefix = getFilePrefix();
         String filePathname = directory + File.separator + filePrefix + name;
         File f = new File(filePathname);
         if (f.exists())
         {
            f.delete();
         }
         serializeFormat(filePathname);
         return true;
      }
      catch (Exception e)
      {
         StsException.outputException("StsSEGYFormat.serialize() failed.", e, StsException.WARNING);
         return false;
      }
   }

   public void serializeFormat(String pathname)
   {
      try
      {
 			File f = new File(pathname);
			if(f.exists())f.delete();
			StsDBFileObject.writeObjectFile(pathname, this, null);
//			StsParameterFile.writeSerializedFields(this, parameterFileName);
			StsMessageFiles.logMessage("Views settings saved in " + pathname);
      }
      catch (Exception e)
      {
         e.printStackTrace();
      }
   }
/*
   public void serializeFormat(String pathname)
   {
      try
      {
         FileOutputStream out = new FileOutputStream(pathname);
         ObjectOutputStream os = new ObjectOutputStream(out);
         os.writeObject(this);
         os.flush();
         out.close();
      }
      catch (Exception e)
      {
         e.printStackTrace();
      }
   }
*/
   static public String[] getSegyFormatFilenames()
   {
      String[] segyFormatFilenames;

      // Create list of templates from Project directory
      String directory = getSegyFormatsDirectory();
      StsFileSet availableSegyFormatFileSet = StsFileSet.constructor(directory, segyFormatFilenameFilter);
      StsAbstractFile[] availableSegyFormatFiles = availableSegyFormatFileSet.getFiles();
      int nFiles = availableSegyFormatFiles.length;
      if (nFiles > 0)
      {
         segyFormatFilenames = new String[nFiles];
         for (int n = 0; n < nFiles; n++)
            segyFormatFilenames[n] = segyFormatFilenameFilter.getFilenameEnding(availableSegyFormatFiles[n].
               getFilename());
         return segyFormatFilenames;
      }
      else
         return new String[0];
   }

   static final class SegyFormatFilenameFilter
      implements FilenameFilter
   {
      String group;
      String format;
      int length;

      public SegyFormatFilenameFilter(String group, String format)
      {
         this.group = group;
         this.format = format;
         length = group.length() + format.length() + 2; // string length of "prefix.format."

      }

      public boolean accept(File dir, String name)
      {
         StsKeywordIO.parseBinaryFilename(name);
         String keywordGroup = StsKeywordIO.group;
         String keywordFormat = StsKeywordIO.format;
         boolean groupOK = keywordGroup.equals(group);
         boolean formatOK = keywordFormat.equals(format);
         return StsKeywordIO.group.equals(group) && StsKeywordIO.format.equals(format);
      }

      public String getFilenameEnding(String filename)
      {
         return filename.substring(length);
      }
   }

   public int getTraceHeaderSize()
   {
      return this.traceHeaderSize;
   }

   public void setTraceHeaderSize(int size)
   {
       if(this.traceHeaderSize == size) return;
      this.traceHeaderSize = size;
      setFormatChanged(true);
   }

   public int getTextHeaderSize()
   {
      return this.textHeaderSize;
   }

   public void setTextHeaderSize(int size)
   {
      if(this.textHeaderSize == size) return;
      this.textHeaderSize = size;
      setFormatChanged(true);
   }

   public int getBinaryHeaderSize()
   {
      return this.binaryHeaderSize;
   }

   public void setBinaryHeaderSize(int size)
   {
      if(this.binaryHeaderSize == size) return;
      this.binaryHeaderSize = size;
      setFormatChanged(true);
   }

   public String getLineNumberFormatString()
   {
      int format = getTraceRecFromUserName(XLINE_NO).getFormat();
      return sampleFormatStrings[format];
   }

   public void setLineNumberFormat(String formatString)
   {
      int format = getSampleFormatFromString(formatString);
      getTraceRecFromUserName(ILINE_NO).setFormat(format);
      getTraceRecFromUserName(XLINE_NO).setFormat(format);
   }

   public int getInlineLoc()
   {
      StsSEGYFormatRec rec = getTraceRecFromUserName(ILINE_NO);
      if (rec == null)
         return -1;
      return rec.getLoc();
   }

   public void setType(int objectType)
   {
       type = (byte)objectType;
       setFormatChanged(true);
   }

   public int getType() { return type; }

   public void setInlineLoc(int loc)
   {
      StsSEGYFormatRec rec = getTraceRecFromUserName(ILINE_NO);
      if (rec == null) return;
      if(rec.loc == loc) return;
      setFormatChanged(true);
      rec.setLoc(loc);
   }

   public int getXlineLoc()
   {
      StsSEGYFormatRec rec = getTraceRecFromUserName(XLINE_NO);
      if (rec == null) return -1;
      return rec.getLoc();
   }

   public void setXlineLoc(int loc)
   {
      StsSEGYFormatRec rec = getTraceRecFromUserName(XLINE_NO);
      if (rec == null) return;
      if(rec.loc == loc) return;
      setFormatChanged(true);
      rec.setLoc(loc);
   }

   public boolean isCoordinateAttribute(StsSEGYFormatRec rec)
   {
      for (int i = 0; i < coorTraceRecs.length; i++)
      {
         if (rec == coorTraceRecs[i])
            return true;
      }
      return false;
   }

   public int getXLoc()
   {
      StsSEGYFormatRec rec = getTraceRecFromUserName(CDP_X);
      if (rec == null)
         return -1;
      return rec.getLoc();
   }

   public void setXLoc(int loc)
   {
      StsSEGYFormatRec rec = getTraceRecFromUserName(CDP_X);
      if (rec == null) return;
      if(rec.loc == loc) return;
      setFormatChanged(true);
      rec.setLoc(loc);
   }

   public int getYLoc()
   {
      StsSEGYFormatRec rec = getTraceRecFromUserName(CDP_Y);
      if (rec == null) return -1;
      return rec.getLoc();
   }

   public void setYLoc(int loc)
   {
      StsSEGYFormatRec rec = getTraceRecFromUserName(CDP_Y);
      if (rec == null) return;
      if(rec.loc == loc) return;
      setFormatChanged(true);
      rec.setLoc(loc);
   }

   public float getSRate()
   {
      return sampleSpacing;
   }

   public void setSRate(float value)
   {
       if(this.sampleSpacing == value) return;
      this.sampleSpacing = value;
      setFormatChanged(true);
   }

    public void setOverrideHeader(boolean override)
    {
        overrideHeader = override;
    }
    
    public boolean getOverrideHeader() { return overrideHeader; }

   public int getNSamp(boolean isLittleEndian, byte[] binaryHeader)
   {
      return nSamp;
   }

    public void setStartZ(float z) { this.startZ = z; }
    public float getStartZ() { return startZ; }

   public void setNSamp(int value)
   {
      if(this.nSamp == value) return;
      this.nSamp = value;
      setFormatChanged(true);
   }

   public void setIsLittleEndian(boolean isLittleEndian)
   {
      if(this.isLittleEndian == isLittleEndian) return;
      this.isLittleEndian = isLittleEndian;
      setFormatChanged(true);
   }

   public boolean getIsLittleEndian()
   {
      return isLittleEndian;
   }

    public void setIsSampleFormatOverride(boolean b)
    {
        isSampleFormatOverride = b;
        if(!isSampleFormatOverride)
            setOverrideSampleFormat(NONE);
    }

    public boolean getIsSampleFormatOverride()
    {
        return isSampleFormatOverride;
    }

//	public void setStartZ(float startZ) { this.startZ = startZ; }
//	public float getStartZ() { return startZ; }

   public void setHUnitString(String hUnitString)
   {
      setHUnits(StsParameters.getDistanceUnitsFromString(hUnitString));
   }

   public void setHUnits(byte hUnits)
   {
       if(this.hUnits == hUnits) return;
       this.hUnits = hUnits;
       this.hUnitString = StsParameters.DIST_STRINGS[hUnits];
       setFormatChanged(true);
   }

   public void setTUnitString(String vUnitString)
   {
      setTUnits(StsParameters.getTimeUnitsFromString(vUnitString));

   }

   public void setTUnits(byte tUnits)
   {
       if(this.tUnits == tUnits) return;
       this.tUnits = tUnits;
       this.tUnitString = StsParameters.TIME_STRINGS[tUnits];
       setFormatChanged(true);
   }

   public void setZUnitString(String zUnitString)
   {
      setZUnits(StsParameters.getDistanceUnitsFromString(zUnitString));
   }

   public void setZUnits(byte zUnits)
   {
       if(this.zUnits == zUnits) return;
       this.zUnits = zUnits;
       this.zUnitString = StsParameters.DIST_STRINGS[zUnits];
       setFormatChanged(true);
   }

   public byte getHUnits()
   {
      return hUnits;
   }

   public byte getZUnits()
   {
      return zUnits;
   }

   public byte getTUnits()
   {
      return tUnits;
   }

   public int getOverrideNSamples()
   {
      return overrideNSamples;
   }

   public float getOverrideSampleSpacing()
   {
      return overrideSampleSpacing;
   }

   public String getOverrideSampleFormatString()
   {
       return sampleFormatStrings[overrideSampleFormat];
   }

   public void setOverrideNSamples(int ns)
   {
      overrideNSamples = ns;
//      formatChanged = true;
   }

   public void setOverrideSampleSpacing(float ss)
   {
      overrideSampleSpacing = ss;
//      formatChanged = true;
   }

    public void setOverrideSampleFormat(byte sampleFormat)
    {
        if(sampleFormat == overrideSampleFormat) return;
        setFormatChanged(true);
        overrideSampleFormat = sampleFormat;
    }

   public String getTextHeaderFormatString()
   {
      return textHeaderFormatString;
   }

   public void setTextHeaderFormatString(String format)
   {
      if(textHeaderFormatString == format) return;
      textHeaderFormatString = format;
      setFormatChanged(true);
   }

   public String getTextHeaderFormat()
   {
      for (int n = 0; n < textHeaderFormatStrings.length; n++)
         if (textHeaderFormatStrings[n].equals(textHeaderFormatString))
            return textHeaderFormats[n];
      return null;
   }
   static public boolean isOtherAttribute(String attribute)
   {
       if(!isTimeAttribute(attribute) && !isDistanceAttribute(attribute) && !isVelocityAttribute(attribute)
          && !isCoordinateAttribute(attribute) && !isDatumAttribute(attribute))
          return true;
       else
          return false;
   }
   static public boolean isTimeAttribute(String attribute)
   {
       return isInList(timeAttributes, attribute);
   }
   static public boolean isDistanceAttribute(String attribute)
   {
       return isInList(distAttributes, attribute);
   }
   static public boolean isVelocityAttribute(String attribute)
   {
       return isInList(velocityAttributes, attribute);
   }
   static public boolean isCoordinateAttribute(String attribute)
   {
       return isInList(coorAttributes, attribute);
   }
   static public boolean isDatumAttribute(String attribute)
   {
       return isInList(datumAttributes, attribute);
   }
   static public boolean isInList(String[] list, String attribute)
   {
       int idx = java.util.Arrays.binarySearch(list, attribute);
       if(idx < 0) return false;
       else return true;
   }
   static public String[] getOtherAttributes(String[] attributes)
   {
       return getAttributes(attributes, OTHER);
   }
   static public String[] getTimeAttributes(String[] attributes)
   {
       return getAttributes(attributes, TIME);
   }
   static public String[] getDistanceAttributes(String[] attributes)
   {
       return getAttributes(attributes, DIST);
   }
   static public String[] getVelocityAttributes(String[] attributes)
   {
       return getAttributes(attributes, VELOCITY);
   }
   static public String[] getScalarAttributes(String[] attributes)
   {
       return getAttributes(attributes, SCALAR);
   }
   static public String[] getCoordinateAttributes(String[] attributes)
   {
       return getAttributes(attributes, COOR);
   }
   static public String[] getDatumAttributes(String[] attributes)
   {
       return getAttributes(attributes, DATUM);
   }
   static public String[] getAttributes(String[] attributes, byte type)
   {
       if(attributes == null) return new String[0];
       String[] atts = new String[attributes.length];
       String[] list = null;
       switch(type)
       {
           case TIME:
               list  = timeAttributes;
               break;
            case DIST:
                list  = distAttributes;
                break;
            case VELOCITY:
                list  = velocityAttributes;
                break;
            case COOR:
               list  = coorAttributes;
               break;
            case DATUM:
                list  = datumAttributes;
                break;
            case SCALAR:
                list  = scalarAttributes;
                break;
            case OTHER:
                list = null;
                break;
            default:
                return null;
       }

       int cnt = 0;
       // Find all attributes in list that are not classified.
       if(list == null)
       {
          for(int i=0; i<attributes.length; i++)
          {
              if(!isTimeAttribute(attributes[i]) && !isDistanceAttribute(attributes[i]) &&
                 !isVelocityAttribute(attributes[i]) && !isCoordinateAttribute(attributes[i]) &&
                 !isDatumAttribute(attributes[i]))
                  atts[cnt++] = attributes[i];
          }
       }
       else
       {
           for(int i=0; i<attributes.length; i++)
           {
               for(int j=0; j<list.length; j++)
               {
                   if(attributes[i].equals(list[j]))
                   {
                       atts[cnt++] = attributes[i];
                       break;
                   }
               }
           }
       }
       return (String[])StsMath.trimArray(atts,cnt);
   }

   static public String toStringStatic(TraceHeader traceHeader)
   {
       if(traceHeader == null) return "nulltraceHeader";
       else return traceHeader.toString();
   }
}