package com.Sts.Types;

import com.Sts.IO.*;
import com.Sts.Utilities.*;

import java.awt.*;

/**
 *
 * <p>Title: S2S Development</p>
 * <p>Description: </p>
 * <p>Copyright: Copyright (c) 2002-2003</p>
 * <p>Company: S2S Systems LLC</p>
 * @author TJLasseter
 * @version beta 1.0
 */

public class StsCoherence
{
    static final public byte DIP_BIN = 0;
    static final public byte DIP_FX = 1;
    static final public byte DIP_USER_DEFINED = 2;

    static final public byte LOW = 0;
    static final public byte MEDIUM = 1;
    static final public byte HIGH = 2;
    static final public byte USER_DEFINED = 3;

    static final public byte TA_SAMPLES = 0;
    static final public byte TA_USER_DEFINED = 1;
    static final public byte SA_RECTANGLE = 0;
    static final public byte SA_CIRCLE = 1;
    static final public byte SA_USER_DEFINED = 2;

    transient public StsSegyVolume segyVolume = null;
    transient public byte dipMethod = DIP_BIN;
    transient public float maxInlineDip = 1;
    transient public float maxXlineDip = 1;
    transient public byte dipSampling = LOW;

    transient public byte temporalAperture = TA_SAMPLES;
    transient public int temporalSamples = 5;
    transient public byte spatialAperature = SA_RECTANGLE;
    transient public boolean excludeDiagnol = false;
    transient public int inlineTracePerSide = 1;
    transient public int xlineTracePerSide = 1;
    transient public float radius = 1;

    transient public boolean adaptiveEigen = false;
    transient public byte adaptiveEigenSmoothing = LOW;
    transient public byte adaptionRate = LOW;
    transient public byte highResSharpening = LOW;
    transient public String cmdLineArgs = null;

    transient public int outputFormat = StsSEGYFormat.BYTE;
    transient public float outputMin = -20.0f;
    transient public float outputMax = 127.0f;

    transient private String executable = "cohv41.exe";

    transient public boolean[] outputVolumes = new boolean[] {false,true,false,false,false,false,false,false,false,
        false,false,false,false, false, false,false, false};

    transient public static String[] volumeTypes = new String[] {"Semblance","Eigen","Eigen High Resolution","Dip","Azimuth","Instantaneous Amplitude Envelope",
        "Instantaneous Frequency","Instantaneous Phase","Carrier (Cosine of Phase)","Response Amplitude",
        "Response Frequency","Response Phase","Response Zero Phase", "Response All Non-Zero Phase",
        "Instantaneous Bandwidth", "Response Bandwidth","Omega"};
    transient public static String[] volumeTips = new String[] {"Semblance","Eigen","Eigen High Resolution","Dip","Azimuth","Instantaneous Amplitude Envelope",
        "Instantaneous Frequency","Instantaneous Phase","Carrier (Cosine of Phase)","Response Amplitude",
        "Response Frequency","Response Phase", "Response Zero Phase", "Response All Non-Zero Phase",
        "Instantaneous Bandwidth", "Response Bandwidth","Omega"};
    transient public static String[] volumeFilename = new String[] {"Semblance","Eigen","EigenHighRes","Dip",
        "Azimuth","InstAmpEnv","InstFreq","InstPhase","Carrier","RespAmp","RespFreq","RespPhase",
        "RespZeroPhase", "RespAllPhase","InstBw", "RespBw","Omega"};

    transient public static boolean[] volumeState = new boolean[] {false,true,false,false,false,false,false,false,false,
        false,false,false,false, false, false,false, false};

    transient public static boolean[] volumeValid = new boolean[] {true,true,true,true,true,true,true,true,true,
        true,true,true,true,true,true, true, true};

    transient public static String[] volumeFlag = new String[] {"-S","-E","-HRE","-D","-A","-IE","-IF","-IP","-C",
        "-RE","-RF","-RP","-P0","-P90","-IBW","-RBW","-O"};

    public StsCoherence() {}

    private StsCoherence(Frame frame, StsSegyVolume segyVolume)
    {
        try
        {
            this.segyVolume = segyVolume;
        }
        catch (Exception e)
        {
            StsException.outputException("StsCoherence(frame, segyVol) failed.", e, StsException.WARNING);
        }
    }

    static public StsCoherence constructor(StsSegyVolume segyVol, Frame frame)
    {
        try
        {
            return new StsCoherence(frame, segyVol);
        }
        catch (Exception e)
        {
            StsException.outputException("StsCoherence.constructTraceAnalyzer(segyVol, frame) failed.", e, StsException.WARNING);
            return null;
        }
    }

    public void prepareDryRunFile()
    {
        StsFile file = StsFile.constructor(segyVolume.stsDirectory + segyVolume.stemname + ".asc");
        StsAsciiFile asciiFile = new StsAsciiFile(file);
        if(!asciiFile.openWriteWithErrorMessage()) return;

        buildFile(true, true, asciiFile);

        asciiFile.close();
    }

    public void prepareExecutionFile(boolean isForeground)
    {
        StsFile file = StsFile.constructor(segyVolume.stsDirectory + segyVolume.stemname + ".asc");
        StsAsciiFile asciiFile = new StsAsciiFile(file);
        if(!asciiFile.openWriteWithErrorMessage()) return;

        buildFile(false , isForeground, asciiFile);

        asciiFile.close();
    }

    public boolean verifyExecutableExist()
    {
        try
        {
            Runtime rt = Runtime.getRuntime();
            Process prcs = rt.exec(executable);
            prcs.destroy();
            return true;
        }
        catch(Exception e)
        {
            return false;
        }
    }

    public String getDryRunCmd()
    {
        String cmd = null;

        if(verifyExecutableExist())
            cmd = new String(executable + " " + segyVolume.stsDirectory + segyVolume.stemname + ".asc");
        else
            cmd = "ERROR: Unable to locate executable (" + executable +")";
        return cmd;
    }

    public String getExecutionCmd()
    {
        String cmd = null;

        if(verifyExecutableExist())
            cmd = new String(executable + " " + segyVolume.stsDirectory + segyVolume.stemname + ".asc");
        else
            cmd = "ERROR: Unable to locate executable (" + executable +")";

        return cmd;
    }

    public void buildFile(boolean isDryRun, boolean isForeground, StsAsciiFile asciiFile)
    {
        try
        {
            asciiFile.writeLine("-fileType Segy");
            asciiFile.writeLine("-N " + segyVolume.segyDirectory + segyVolume.segyFilename);
            for(int i=0; i<volumeState.length; i++)
            {
                if(volumeState[i] == true)
                    asciiFile.writeLine(volumeFlag[i] + " " + segyVolume.segyDirectory + segyVolume.getName() + volumeFilename[i]  + ".sgy");
            }

            if(isDryRun)
                asciiFile.writeLine("-dryrun");

            if(!isForeground)
                asciiFile.writeLine("-pfile " + segyVolume.segyDirectory + segyVolume.stemname + ".log");

            int fmt = segyVolume.segyData.getSampleFormat();
			String fmtString = StsSEGYFormat.sampleFormatStrings[fmt];
            asciiFile.writeLine("-sgy_informat " + fmtString);
            fmtString = StsSEGYFormat.sampleFormatStrings[outputFormat];
            asciiFile.writeLine("-outformat " + fmtString);

            float num = segyVolume.getRowNumMin() + (segyVolume.cropBox.rowMin - segyVolume.getRowMin()) * segyVolume.getRowNumInc();
            asciiFile.writeLine("-first_lineout " + num);
            num = segyVolume.getRowNumMin() + (segyVolume.cropBox.rowMax - segyVolume.getRowMin()) * segyVolume.getRowNumInc();
            asciiFile.writeLine("-last_lineout " + num);
            num = segyVolume.getColNumMin() + (segyVolume.cropBox.colMin - segyVolume.getColMin()) * segyVolume.getColNumInc();
            asciiFile.writeLine("-first_traceout " + num);
            num = segyVolume.getColNumMin() + (segyVolume.cropBox.colMax - segyVolume.getColMin()) * segyVolume.getColNumInc();
            asciiFile.writeLine("-last_traceout " + num);
            asciiFile.writeLine("-lineinc " + segyVolume.rowNumInc);
            asciiFile.writeLine("-traceinc " + segyVolume.colNumInc);
            num = segyVolume.cropBox.sliceMin * segyVolume.getZInc();
            asciiFile.writeLine("-tstart " + num);
            num = segyVolume.cropBox.sliceMax * segyVolume.getZInc();
            asciiFile.writeLine("-tend " + num);
            asciiFile.writeLine("-tskip " + (segyVolume.cropBox.sliceInc * segyVolume.getZInc()));
            asciiFile.writeLine("-ilazim " + segyVolume.getAngle());

            if((segyVolume.getAngle()+90) > 360.0f)
                asciiFile.writeLine("-clazim " + (segyVolume.getAngle()+90-360));
            else
                asciiFile.writeLine("-clazim " + (segyVolume.getAngle()+90));

            asciiFile.writeLine("-ildm " + segyVolume.getXInc());
            asciiFile.writeLine("-cldm " + segyVolume.getYInc());
            asciiFile.writeLine("-linebyte " + (segyVolume.getSegyFormat().getTraceRecFromUserName("ILINE_NO").getLoc()+1));
            asciiFile.writeLine("-tracebyte " + (segyVolume.getSegyFormat().getTraceRecFromUserName("XLINE_NO").getLoc()+1));
            if(dipMethod == DIP_BIN)
            {
                asciiFile.writeLine("-dipmethod bins");
                asciiFile.writeLine("-ildip_max " + maxInlineDip);
                asciiFile.writeLine("-cldip_max " + maxXlineDip);
            }
            else if(dipMethod == DIP_FX)
                asciiFile.writeLine("-dipmethod fx");
            else
                asciiFile.writeLine("-dipmethod userdefined");

            if(dipSampling == LOW)
                asciiFile.writeLine("-dip_sample low");
            else if(dipMethod == MEDIUM)
                asciiFile.writeLine("-dip_sample med");
            else if(dipMethod == HIGH)
                asciiFile.writeLine("-dip_sample high");

            if(temporalAperture == TA_SAMPLES)
            {
                asciiFile.writeLine("-w " + temporalSamples);
                asciiFile.writeLine("-lentgate " + temporalSamples);
            }
            else if(temporalAperture == TA_USER_DEFINED)
            {
                ;
            }

            if(spatialAperature == SA_RECTANGLE)
            {
                asciiFile.writeLine("-ilhap " + inlineTracePerSide);
                asciiFile.writeLine("-clhap " + xlineTracePerSide);
            }
            else if(spatialAperature == SA_CIRCLE)
            {
                asciiFile.writeLine("-rad " + radius);
            }
            else
            {
                ;
            }

            if(excludeDiagnol)
                asciiFile.writeLine("-excludediag");

            asciiFile.writeLine("-att_outmin " + outputMin);
            asciiFile.writeLine("-att_outmax " + outputMax);
            asciiFile.writeLine("-outmin -128");
            asciiFile.writeLine("-outmax 127");

            if(adaptiveEigen)
            {
                asciiFile.writeLine("-adaptive_yes");
                if(adaptiveEigenSmoothing == LOW)
                    asciiFile.writeLine("-adapt_smooth low");
                else if(adaptiveEigenSmoothing == MEDIUM)
                    asciiFile.writeLine("-adapt_smooth med");
                else if(adaptiveEigenSmoothing == HIGH)
                    asciiFile.writeLine("-adapt_smooth high");

                if(adaptionRate == LOW)
                    asciiFile.writeLine("-adapt_rate low");
                else if(adaptionRate == MEDIUM)
                    asciiFile.writeLine("-adapt_rate med");
                else if(adaptionRate == HIGH)
                    asciiFile.writeLine("-adapt_rate high");
            }
            if(highResSharpening == LOW)
                asciiFile.writeLine("-hrlevel low");
            else if(highResSharpening == MEDIUM)
                asciiFile.writeLine("-hrlevel med");
            else if(highResSharpening == HIGH)
                asciiFile.writeLine("-hrlevel high");

            asciiFile.writeLine("-input_file_tstart " + segyVolume.getZMin());
            asciiFile.writeLine("-input_file_tend " + segyVolume.getZMax());

            asciiFile.writeLine(cmdLineArgs);

            asciiFile.close();
        }
        catch (Exception e)
        {
            StsException.outputException("StsCoherence.buildFile() failed.",
                e, StsException.WARNING);
            return;
        }
    }

    public void validateOutputTypes()
    {
        if(adaptiveEigen)
        {
            for(int i=0; i<volumeValid.length; i++)
            {
                if((i == 1) || (i == 2))
                {
                    volumeState[i] = true;
                    continue; // Eigen is only valid attribute with adaptive enabled.
                }
                volumeValid[i] = false;
                volumeState[i] = false;
            }
        }
    }
/*
    private String getFormat(int val)
    {
        try
        {
            switch (val)
            {
                case StsSEGYFormat.BYTE:
                    return "int_8";
                case StsSEGYFormat.INT2:
                    return "int_16";
                case StsSEGYFormat.INT4:
                    return "int_32";
                case StsSEGYFormat.IBMFLT:
                    return "float_ibm";
                case StsSEGYFormat.IEEEFLT:
                    return "float_32";
                case StsSEGYFormat.FLOAT8:
                    return "float_8";
                case StsSEGYFormat.FLOAT16:
                    return "float_16";
                default:
                    System.err.println("Format " +  segyVolume.segyData.getSampleFormat() + " unrecognized");
                    return null;
            }
        }
        catch (Exception e)
        {
            StsException.outputException("StsSegyVolume.getTraceHeaderBinary() failed.",
                e, StsException.WARNING);
            return null;
        }
    }
*/
}
