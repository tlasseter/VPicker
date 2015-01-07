package com.Sts.DBTypes;

import com.Sts.Actions.Import.*;
import com.Sts.IO.*;
import com.Sts.MVC.*;
import com.Sts.MVC.View3d.*;
import com.Sts.Types.*;
import com.Sts.UI.Beans.*;
import com.Sts.UI.ObjectPanel.*;
import com.Sts.Utilities.DateTime.*;
import com.Sts.Utilities.*;

import javax.media.opengl.*;
import java.text.*;
import java.util.*;

/**
 * Created by IntelliJ IDEA.
 * User: Tom Lasseter
 * Date: Feb 13, 2009
 * Time: 12:00:21 AM
 * To change this template use File | Settings | File Templates.
 */
public class StsTimeLogCurve extends StsLogCurve // implements StsTreeObjectI
{
    public String directory;
    public String filename;
    public StsWell well;
    protected StsLongVector timesVector;
    protected StsLogVector[] surveyVectors;
    private float minValue;
    private float maxValue;
    protected boolean asCurtain = true;

    transient byte[] data;
    transient public int nSelectedSurvey;
    transient long selectedSurveyTime;
    transient long projectTime;
    transient StsTimeLogCurtain logCurtain;

    static final String GROUP = StsLogVector.WELL_TIME_LOG_PREFIX;
    static final String FORMAT = StsLogVector.FORMAT_BIN;
    static final String MDEPTH = StsWellIO.MDEPTH;

    static public SimpleDateFormat dateTimeFormat;
    static public StsFieldBean[] displayFields = null;
    static protected StsObjectPanel objectPanel = null;

    public StsTimeLogCurve()
    {
    }

    public StsTimeLogCurve(String directory, String filename, String curvename, String curveTypeName, StsSpectrum spectrum, float dataMin, float dataMax)
    {
        super(false);
        String stemname = StsStringUtils.trimSuffix(filename);
        setName(stemname);
        this.directory = directory;
        this.filename = filename;
        this.curvename = curvename;
        logCurveType = (StsTimeLogCurveType)getStsClassObjectWithName(StsTimeLogCurveType.class, curveTypeName);
        if(logCurveType == null)
            logCurveType = new StsTimeLogCurveType(curveTypeName, dataMin, dataMax, spectrum);
        else
            logCurveType.adjustRange(dataMin, dataMax);
        dateTimeFormat = new SimpleDateFormat(currentModel.getProject().getTimeDateFormatString());
        addToModel();
        addToLogCurveType();
    }

    public StsTimeLogCurve(StsFile file, String curvename, String curveTypeName, StsSpectrum spectrum, float dataMin, float dataMax)
    {
        this(file.getDirectory(), file.getFilename(), curvename, curveTypeName, spectrum, dataMin, dataMax);
    }

    private void addToLogCurveType()
    {
        logCurveType.addLogCurve(this);
    }

    public boolean initialize(StsModel model)
    {
        initialize();
        dateTimeFormat = new SimpleDateFormat(model.getProject().getTimeDateFormatString());
        return true;
    }

    public long getDeltaT()
    {
        projectTime = currentModel.getProject().getProjectTime();
        return projectTime - bornDate;
    }

    public float getDeltaTFraction()
    {
        projectTime = currentModel.getProject().getProjectTime();
        float f = (float)((double)(projectTime - bornDate)/(deathDate - bornDate));
        return StsMath.minMax(f, 0.0f, 1.0f);
    }

    public void setDeltaT(float deltaT)
    {
        projectTime = (long)deltaT + bornDate;
        currentModel.getProject().setProjectTime(projectTime, true);
    }

    public long getTimeFromFraction(double fraction)
    {
        return (long)(bornDate + fraction*(deathDate - bornDate));
    }

    public boolean initialize(float[] mdepths, String[] timeStrings, float[][] vectorData)
    {
        long[] times = constructTimes(timeStrings);
        int nTimes = times.length;
        bornDate = times[0];
        dbFieldChanged("bornDate", bornDate);
        deathDate = times[nTimes - 1];
        dbFieldChanged("deathDate", deathDate);
        // currentModel.project.setProjectTime(bornDate);

        // mdepthVector is nontransient, so if this object is reloaded from db, we just need to fill in the values
        String wellCurveName = wellname + "." + logCurveType.name;
        String binFilename = StsKeywordIO.constructFilename(GROUP, FORMAT, wellCurveName, MDEPTH, version);
        mdepthVector = new StsLogVector(filename, binFilename, MDEPTH);
        dbFieldChanged("mdepthVector", mdepthVector);
        mdepthVector.setValues(mdepths);
        mdepthVector.deleteWriteBinaryFile();
        // depthVector = well.getDepthsFromMDepths(mdepthVector);
        // depthVector.deleteWriteBinaryFile();
        // timesVector and valueVectors are always transient, so construct and fill with values
        timesVector = new StsLongVector(times);
        dbFieldChanged("timesVector", timesVector);
        int nSurveys = vectorData.length;
        surveyVectors = new StsLogVector[nSurveys];
        dbFieldChanged("surveyVectors", surveyVectors);
        for(int n = 0; n < nSurveys; n++)
        {
            binFilename = StsKeywordIO.constructFilename(GROUP, FORMAT, wellCurveName, timeStrings[n], version);
            surveyVectors[n] = new StsLogVector(filename, binFilename, curvename);
            surveyVectors[n].setValues(vectorData[n]);
            surveyVectors[n].deleteWriteBinaryFile();
        }
        return true;
    }

    public int checkSetNSelectedSurvey()
    {
        if(surveyVectors == null) return -1;
        if(projectTime == selectedSurveyTime) return nSelectedSurvey;
        
        if(selectedSurveyTime == 0 || selectedSurveyTime != projectTime)
        {
            selectedSurveyTime = projectTime;
            long[] times = timesVector.getValues();
            int nTimes = times.length;
            if(selectedSurveyTime <= bornDate)
            {
                selectedSurveyTime = bornDate;
                nSelectedSurvey = 0;
            }
            else if(selectedSurveyTime >= deathDate)
            {
                selectedSurveyTime = deathDate;
                nSelectedSurvey = nTimes -1;
            }
            else
            {
                nSelectedSurvey = StsMath.binarySearch(times, selectedSurveyTime);
                long timeAbove = times[nSelectedSurvey];
                long timeBelow = times[nSelectedSurvey - 1];
                float f = (selectedSurveyTime - timeBelow) / (timeAbove - timeBelow);
                if(f <= 0.5f)
                    nSelectedSurvey--;
            }
        }
        return nSelectedSurvey;
    }

    public StsLogVector getValueVector()
    {
        checkSetNSelectedSurvey();
        if(nSelectedSurvey == -1) return null;
        StsLogVector selectedSurveyVector = surveyVectors[nSelectedSurvey];
        selectedSurveyVector.checkLoadVector(stsDirectory);
        return selectedSurveyVector;
    }

    public float getValueAtTimeAndMdepth(long time, float mdepth)
    {
        float indexF = mdepthVector.getIndexF(mdepth);
        if(indexF < 0.0f || indexF > mdepthVector.getSize()-1)
            return StsParameters.nullValue;
        else
        {
            long[] times = timesVector.getValues();
            int nSelectedSurvey = StsMath.binarySearch(times, time);
            StsLogVector surveyVector = surveyVectors[nSelectedSurvey];
            if(surveyVector == null) return StsParameters.nullValue;
            return surveyVector.getValue(indexF);
        }
    }

    public StsFloatVector getValuesFloatVector()
    {
        StsLogVector surveyVector = getValueVector();
        if(surveyVector == null) return null;
        return surveyVector.getValues();
    }

    public long[] getTimes() { return timesVector.getValues(); }

    public boolean initializeData()
    {
        float[] range = getColorscaleDisplayRange();
        float scale = StsMath.floatToUnsignedByteScale(range[0], range[1]);
        float offset = StsMath.floatToUnsignedByteScaleOffset(scale, range[0]);
        int nVectors = surveyVectors.length;
        int nVectorValues = mdepthVector.getSize();
        data = new byte[nVectors * nVectorValues];
        int i = 0;
        for(int n = 0; n < nVectors; n++)
        {
            float[] values = surveyVectors[n].getFloats();
            for(int v = 0; v < nVectorValues; v++)
                data[i++] = StsMath.floatToUnsignedByte254WithScale(values[v], scale, offset);
        }
        return true;
    }

    private float[] getColorscaleDisplayRange()
    {
        return ((StsTimeLogCurveType)logCurveType).getColorscale().getEditRange();
    }

    public byte[] getData()
    {
        if(data != null) return data;
        if(!initializeData()) return null;
        return data;
    }

    public boolean readFile()
    {
        // Process the files.
        StsFile file = StsFile.constructor(directory, filename);
        String[][] fileStrings = StsAsciiTokensFile.getStringTable(directory, filename, ",");
        if(fileStrings == null) return false;
        int nRows = fileStrings.length;
        int nCols = fileStrings[0].length;
        int nTimes = nCols - 1;
        int nDepths = nRows - 3;

        // mdepthVector = new StsLogVector();
        float[] mdepths = new float[nDepths];
        for(int row = 3, nDepth = 0; row < nRows; row++, nDepth++)
        {
            try
            {
                mdepths[nDepth] = Float.parseFloat(fileStrings[row][0]);
            }
            catch(Exception e)
            {
                StsException.systemError(this, "run", "Couldn't parse to float: " + fileStrings[row][0]);
            }
        }
        // mdepthVector.setValues(mdepths);
        String[] timeStrings = getDateAndTimeStrings(fileStrings[0], fileStrings[1], 1);
        float[][] vectorData = new float[nTimes][nDepths];
        int nGoodTimes = 0;
        minValue = StsParameters.largeFloat;
        maxValue = -StsParameters.largeFloat;
        for(int col = 1, nTime = 0; col < nCols; col++, nTime++)
        {
            if(timeStrings[nTime] != null)
            {
                for(int row = 3, nDepth = 0; row < nRows; row++, nDepth++)
                {
                    try
                    {
                        float value = Float.parseFloat(fileStrings[row][col]);
                        minValue = Math.min(minValue, value);
                        maxValue = Math.max(maxValue, value);
                        vectorData[nGoodTimes][nDepth] = Float.parseFloat(fileStrings[row][col]);
                    }
                    catch(Exception e)
                    {
                        StsException.systemError(this, "run", "Couldn't parse to float: " + fileStrings[row][col]);
                    }
                }
                nGoodTimes++;
            }
        }
        if(nGoodTimes < nTimes)
            vectorData = (float[][])StsMath.trimArray(vectorData, nGoodTimes);

        initialize(mdepths, timeStrings, vectorData);
        return true;
    }

    private String[] getDateAndTimeStrings(String[] dates, String[] clockTimes, int nSkip)
    {
        int nStrings = dates.length;
        int nTimes = nStrings - nSkip;
        String[] dateTimes = new String[nTimes];
        for(int n = 0, i = nSkip; n < nTimes; n++, i++)
        {
            dateTimes[n] = dates[i] + " " + clockTimes[i];
            try
            {
                Calendar cal = CalendarParser.parse(dateTimes[n], CalendarParser.YY_MM_DD, true);
                dateTimes[n] = CalendarParser.prettyString(cal);
            }
            catch(Exception e)
            {
                StsException.systemError(this, "getDateAndTimeString", "failed to parse date/time: " + dateTimes[n]);
                dateTimes[n] = null;                          
            }
        }
        return dateTimes;
    }

    private long[] constructTimes(String[] dateAndTimes)
    {
        int nTimes = dateAndTimes.length;
        long[] times = new long[nTimes];
        for(int n = 0; n < nTimes; n++)
        {
            try
            {
                if(dateAndTimes[n] == null)
                    times[n] = 0;
                else
                {
                    Calendar cal = CalendarParser.parse(dateAndTimes[n], CalendarParser.YY_MM_DD, true);
                    times[n] = cal.getTimeInMillis();
                }
            }
            catch(Exception e)
            {
                times[n] = 0;
            }
        }
        return times;
    }

    public int getNSurveys()
    {
        return surveyVectors.length;
    }

    public int getNSamples()
    {
        return mdepthVector.getSize();
    }

    public StsWell getWell()
    {
        return well;
    }

    public void setWell(StsWell well)
    {
        this.well = well;
        if(well == null)
            wellname = null;
        else
            wellname = well.getName();
    }

    public String getName()
    {
        return name;
    }

    public StsTimeLogCurtain getTimeLogCurtain()
    {
        if(logCurtain == null) constructLogCurtain();
        return logCurtain;
    }

    private boolean constructLogCurtain()
    {
        logCurtain = StsTimeLogCurtain.constructor(currentModel, this);
        return logCurtain != null;
    }

    public void display3d(StsGLPanel3d glPanel3d, StsWell well, float origin)
    {
        display3d(glPanel3d, well, origin, well.getMaxMDepth());
    }

    public void display3d(StsGLPanel3d glPanel3d, StsWell well, float origin, float mdLimit)
    {
        //if(!isAlive(currentModel.getProject().getProjectTime()))
        //    return;
        String displayTypeString = wellClass.getDisplayTypeString();
        if(displayTypeString == StsWellClass.DISPLAY_CURTAIN)
        {
            int logCurveWidth = wellClass.getLogCurveDisplayWidth();
            if(logCurtain == null && !constructLogCurtain()) return;
            logCurtain.displayTextureCurtain(glPanel3d, logCurveWidth, origin);   // ToDo: Need to be able to md limit these
        }
        else if(displayTypeString == StsWellClass.DISPLAY_CURVE)
            displayLog3d(glPanel3d, well, origin, mdLimit);
        else if(displayTypeString == StsWellClass.DISPLAY_LINE)
        {
            int logLineWidth = wellClass.getLogLineDisplayWidth();
            if(logCurtain == null && !constructLogCurtain()) return;
            logCurtain.displayTextureLine(glPanel3d, logLineWidth, origin);   // ToDo: Need to be able to md limit these
        }
    }

    public float getMinValue()
    {
        return minValue;
    }

    public float getMaxValue()
    {
        return maxValue;
    }

    public String getBornTimeString() { return dateTimeFormat.format(bornDate); }
    public String getDeathTimeString() { return dateTimeFormat.format(deathDate); }
    public String getProjectTimeString() { return dateTimeFormat.format(projectTime); }

    public StsPoint[] getWellPoints()
    {
        float[] mdepths = mdepthVector.getFloats();
        int nPoints = mdepths.length;
        StsPoint[] points = new StsPoint[nPoints];
        for(int n = 0; n < nPoints; n++)
            points[n] = well.getPointAtMDepth(mdepths[n], true);
        return points;
    }

    public StsPoint[] getWellSlopes()
    {
        float[] mdepths = mdepthVector.getFloats();
        return well.getSlopesAtMDepths(mdepths);
    }

    /** override in subclass if this feature is available */
    public boolean getIsPixelMode() { return false; }
/*
    public void initializeColorscale()
    {
        try
        {
            if(colorscale == null)
            {
                constructColorscale();
                dbFieldChanged("colorscale", colorscale);
            }
            colorList = new StsColorList(colorscale);
            colorscale.addActionListener(this);
        }
        catch (Exception e)
        {
            StsException.outputException("StsTimeLogCurve.initializeColorscale() failed.", e, StsException.WARNING);
        }
    }

    private void constructColorscale()
    {
        String name = logCurveType.getName();
        colorscale = new StsColorscale(name, spectrum, minValue, maxValue);
    }
*/
    public boolean setGLColorList(GL gl, boolean nullsFilled, int shader)
	{
        return ((StsTimeLogCurveType)logCurveType).setGLColorList(gl, nullsFilled, shader);
	}

    public void colorscaleChanged()
    {
        data = null;
        if(logCurtain != null) logCurtain.setTextureChanged();
    }

    public String getDateTimeString(long dateTime)
    {
        return dateTimeFormat.format(dateTime);
    }

    public boolean hasColorscale() { return true; }

    public StsColorscale getColorscale() { return ((StsTimeLogCurveType)logCurveType).getColorscale(); }

    // not currently displaying instance of timeLogCurve on object panel
/*
    public StsFieldBean[] getDisplayFields()
    {
       try
       {
           if (displayFields == null)
           {
               displayFields = new StsFieldBean[]
               {
                    new StsBooleanFieldBean(StsTimeLogCurve.class, "isVisible", "Enable"),
                    new StsDateFieldBean(StsTimeLogCurve.class, "bornDate", false, "Minimum Time:"),
                    new StsDateFieldBean(StsTimeLogCurve.class, "deathDate", false, "Maximum Time:"),
                    new StsFloatFieldBean(StsTimeLogCurve.class, "minValue", "Min Temp:"),
                    new StsFloatFieldBean(StsTimeLogCurve.class, "maxValue", "Max Temp:")
                    //new StsEditableColorscaleFieldBean(StsTimeLogCurve.class, "colorscale")
               };
           }
            return displayFields;
        }
        catch (Exception e)
        {
            StsException.outputException("StsTimeLogCurve.getDisplayFields() failed.", e, StsException.WARNING);
            return null;
        }
    }
*/
    // not currently displaying instance of timeLogCurve on objectPanel
    /*
    public StsFieldBean[] getPropertyFields() { return null; }

    public StsFieldBean[] getDefaultFields() { return null; }

    public Object[] getChildren() { return new Object[0];}

    public StsObjectPanel getObjectPanel()
    {
        if (objectPanel == null)
        {
            objectPanel = StsObjectPanel.constructor(this, true);
        }
        return objectPanel;
    }

    public boolean anyDependencies() { return false; }

    public boolean canExport() { return false; }

    public void treeObjectSelected() { }
    */
}
