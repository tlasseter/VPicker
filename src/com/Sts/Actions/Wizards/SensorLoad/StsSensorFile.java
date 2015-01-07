

//Title:        S2S: Seismic-to-simulation
//Version:
//Copyright:    Copyright (c) 2001
//Author:       Tom Lasseter
//Company:      4D Systems LLC
//Description:  Web-enabled Integrated Interpretation System

package com.Sts.Actions.Wizards.SensorLoad;

import com.Sts.Actions.Import.*;
import com.Sts.DBTypes.*;
import com.Sts.IO.*;
import com.Sts.MVC.*;
import com.Sts.UI.Beans.StsDateFieldBean;
import com.Sts.Utilities.DateTime.*;
import com.Sts.Utilities.*;

import java.io.*;
import java.util.*;

public class StsSensorFile
{
	public boolean relativeCoordinates = false;
    public double relativeX = 0.0f;
    public double relativeY = 0.0f;
    public double relativeZ = 0.0f;

	public int numStages = 1;
	public int[] stageOffsets = null;
    public boolean overrideStages = false;

	public byte positionType = StsSensor.STATIC;
    public double staticX = 0.0f;
    public double staticY = 0.0f;
    public double staticZ = 0.0f;

	public byte timeType = StsSensorKeywordIO.NO_TIME;
	public long currentTime = 0L;
    public long startTime = 0L;
    private int dateOrder = CalendarParser.DD_MM_YY;

	public StsFile file = null;
	public int[] attIndices = null;
	public String[] curveNames = null;
	public String[] validCurves = null;
	public double[] currentValues = null;

    public boolean hasHeader = true;
    public int numHeaderRows = 0;
    public int numCols = 0;

    public byte type = VARIABLE;
    public static byte FIXED = 0;
    public static byte VARIABLE = 1;
    public String delimiters = ";,";

    public static byte X = 0;
    public static byte Y = 1;
    public static byte Z = 2;
    public static byte TIME = 3;
    public static byte DATE = 4;
    int[] colOrder = new int[] { -1, -1, -1, -1, -1 };
    int numColsPerRow = 5;

    static public char[] invalidChars = new char[] {'#',' ',':',';','"','\'',',', ')', '(', '[', ']', '{', '}', '=', '*', '@', '?','/', '\\', '$', '%', '&', '*'};

    public StsSensorFile(StsFile file)
    {
    	this.file = file;
    }

    public void analyzeFile()
    {
    	analyzeFile(null, type);
    }
    public void analyzeFile(StsModel model)
    {
        analyzeFile(model, type);
    }
    public boolean getAttributeValues(StsModel model, String line, int dateOrder)
    {
        this.dateOrder = dateOrder;
        return getAttributeValues(model, line, null);
    }
    public boolean getAttributeValues(StsModel model, String line)
    {
        return getAttributeValues(model, line, null);
    }
    public boolean getAttributeValues(StsModel model, String line, String dateToken)
    {
        StringTokenizer stok = null;
        String token = null;
        double value = 0.0f;

        String timeToken = null;
        try
        {
            line = line.trim();
            line = StsStringUtils.detabString(line);
            line = StsStringUtils.deQuoteString(line);
            line = StsStringUtils.cleanLine(line, delimiters);

            stok = new StringTokenizer(line, delimiters);
            int nTokens = stok.countTokens();
            if(nTokens < curveNames.length)
                return false;
            currentValues = new double[nTokens];
            validCurves = new String[nTokens];
            //dateToken = null;
            int cnt = 0;
            for (int i = 0; i < nTokens; i++)
            {
                token = stok.nextToken();
                if(colOrder[TIME] == i)
                	timeToken = StsStringUtils.deWhiteSpaceString(token);
                else if(colOrder[DATE] == i)
                	dateToken = StsStringUtils.deWhiteSpaceString(token);
                else
                {
                    try
                    {
                    	value = Double.parseDouble(token);
                    }
                    catch (Exception e)
                    {
                    	if(model != null)
                            value = model.getProject().getLogNull();      // Not be a number
                        else
                            value = 0.0f;
                    }
                    if(Double.valueOf(value).isNaN())
                    {
                    	if(model != null)
                            value = model.getProject().getLogNull();      // Not be a number
                        else
                            value = 0.0f;
                    }
                    currentValues[cnt] = value;
                    //System.out.println("curveName=" + curveNames[cnt] + " Value=" + value);
                    validCurves[cnt] = curveNames[cnt];
                    cnt++;
                }
            }
            if(model != null)
                dateOrder = model.getProject().getDateOrder();

            String dateTime = dateToken;
            if(timeType == StsSensorKeywordIO.ELAPSED_TIME)
            {
                currentTime = startTime + (long)(Double.valueOf(timeToken) * 60000); // Assuming minutes
            }
            else if((timeType != StsSensorKeywordIO.TIME_AND_DATE) && (timeType != StsSensorKeywordIO.TIME_OR_DATE))
            {
                if(dateToken == null)
                {
                    dateToken = "01-01-71 00:00:00.0";
                    if(model != null)
                        dateToken = model.getProject().getDateFromLong(startTime);
                }
                dateTime = timeToken;
                if(dateToken != null)
                    dateTime = dateToken + " " + timeToken;
                Calendar cal = CalendarParser.parse(dateTime, dateOrder , true);
                currentTime = cal.getTimeInMillis();
            }
            else
            {
                dateTime = timeToken;
                if(dateToken != null)
                    dateTime = dateToken + " " + timeToken;
                Calendar cal = CalendarParser.parse(dateTime, dateOrder , true);
                currentTime = cal.getTimeInMillis();
            }

        	currentValues = (double[])StsMath.trimArray(currentValues, cnt);
        	validCurves = (String[])StsMath.trimArray(validCurves, cnt);
        }
        catch(Exception e)
        {
            StsMessageFiles.logMessage("Failed to find or parse time/date column(s) in file - check the DATE ORDER paramter: " + file.getFilename());
            return false;
        }
        return true;
    }

    public boolean analyzeFile(StsModel model, byte type)
    {
        if(type == FIXED)
            return analyzeFixedColumnFile(model);
        else
            return analyzeVariableColumnFile(model);
    }

    public boolean analyzeFixedColumnFile(StsModel model)
    {
        BufferedReader bufRdr = null;
        String line = null;
        StringTokenizer stok = null;
        String token = null;
        delimiters = " ";
        try
        {
            bufRdr = new BufferedReader(new FileReader(file.getPathname()));
            line = StsSensorFile.readLine(bufRdr);
            line = line.trim();
            line = StsStringUtils.detabString(line);
            line = StsStringUtils.deQuoteString(line);
            line = StsStringUtils.cleanLine(line, delimiters);

            if(line == null) return false;

            stok = new StringTokenizer(line, delimiters);
            int nTokens = stok.countTokens();
            int hdrTokens = nTokens;
            numCols = nTokens - 1; // Have time and date columns but only time header
            int subtract = 0;
            for (int i = 0; i < hdrTokens; i++)
            {
            	token = stok.nextToken().trim();
                if(token.equalsIgnoreCase("SRC"))   // This format always has date and time in column 0 and 1 but only one header ....time for both
                {
                    nTokens = nTokens - 1;
                    subtract = 1;
                    colOrder[DATE] = 0;
                    colOrder[TIME] = 1;
                    break;
                }
                token = token.replace(".", " ").trim();
            	colOrder[TIME] = checkIndex(token, i - subtract, colOrder[TIME], StsSensorKeywordIO.TIME_KEYWORDS, null);
            	colOrder[DATE] = checkIndex(token, i - subtract, colOrder[DATE], StsSensorKeywordIO.DATE_KEYWORDS, null);
            }

            // Parse the header line to find key indices and attribute names
            parseHdrAttributes(line, subtract, nTokens);

            // Determine if time and date columns exist and what type
            line = determineTimeType(model, bufRdr);

            // If file has both a time and a date column need to trim arrays
            if(timeType == StsSensorKeywordIO.TIME_OR_DATE)
            {
            	curveNames = (String[])StsMath.arrayDeleteLastElement(curveNames);
            	attIndices = (int[])StsMath.arrayDeleteLastElement(attIndices);
            }
            if(model != null)
            {
                if(line == null)
                    return false;
                line = line.trim();
                line = StsStringUtils.detabString(line);
                line = StsStringUtils.deQuoteString(line);
                line = StsStringUtils.cleanLine(line, delimiters);

                StsSensorClass sensorClass = getSensorClass(model);
                relativeCoordinates = relativeCoordinates(line, sensorClass.getXRelativeCriteria());
            }
            else
                relativeCoordinates = false;
        }
        catch(Exception e)
        {
            StsMessageFiles.logMessage("Failed to find time column in file: " + file.getFilename());
            return false;
        }
        return true;
    }

    public boolean analyzeVariableColumnFile(StsModel model)
    {
        BufferedReader bufRdr = null;
        String line = null;
        StringTokenizer stok = null;
        String token = null;
        delimiters = ";,";
        try
        {
            boolean colHeaderFound = false;
            bufRdr = new BufferedReader(new FileReader(file.getPathname()));
            while((line = StsSensorFile.readLine(bufRdr)) != null)
            {
                for(int i=0; i<StsSensorKeywordIO.TIME_KEYWORDS.length; i++)
                {
                    if(line.toLowerCase().indexOf(StsSensorKeywordIO.TIME_KEYWORDS[i].toLowerCase()) != -1)
                    {
                        colHeaderFound = true;
                        break;
                    }
                }
                if(colHeaderFound) break;
            }
            if(!colHeaderFound)
            {
                hasHeader = false;
                return false;
            }
            line = line.trim();
            line = StsStringUtils.detabString(line);
            line = StsStringUtils.deQuoteString(line);
            line = StsStringUtils.cleanLine(line, delimiters);

            if(line == null) return false;


            stok = new StringTokenizer(line, delimiters, false);
            int nTokens = stok.countTokens();
            curveNames = new String[nTokens-1];
            attIndices = new int[nTokens-1];
            numCols = nTokens;
            int cnt = 0;
            for (int i = 0; i < nTokens; i++)
            {
            	token = stok.nextToken().trim();
                colOrder[TIME] = checkIndex(token, i, colOrder[TIME], StsSensorKeywordIO.TIME_KEYWORDS, null);
                colOrder[DATE] = checkIndex(token, i, colOrder[DATE], StsSensorKeywordIO.DATE_KEYWORDS, null);
            }

            if(colOrder[TIME] == -1 && colOrder[DATE] == -1)
                hasHeader = false;

            // Parse the header line to find key indices and attribute names
            parseHdrAttributes(line, 0, nTokens);

            // Determine if time and date columns exist and what type
            line = determineTimeType(model, bufRdr);

            // If file has both a time and a date column need to trim arrays
            if(timeType == StsSensorKeywordIO.TIME_OR_DATE)
            {
            	curveNames = (String[])StsMath.arrayDeleteLastElement(curveNames);
            	attIndices = (int[])StsMath.arrayDeleteLastElement(attIndices);
            }
            if(model != null)
            {
                StsSensorClass sensorClass = getSensorClass(model);
                if(line == null)       // Nothing but header in file
                    return false;
                line = line.trim();
                line = StsStringUtils.detabString(line);
                line = StsStringUtils.deQuoteString(line);
                line = StsStringUtils.cleanLine(line, delimiters);

                relativeCoordinates = relativeCoordinates(line, sensorClass.getXRelativeCriteria());
            }
            else
                relativeCoordinates = false;
        }
        catch(Exception e)
        {
            StsMessageFiles.logMessage("Failed to find time column in file: " + file.getFilename());
            return false;
        }
        return true;
    }
    
    public byte getType() { return type; }
    public void setType(byte val) { type = val; }
    public int getNumCols()
    {
        return numCols;
    }
    public void setNumCols(int num)    
    {
        numCols = num;
    }
    public int getColLocation(byte type)
    {
        return colOrder[type];
    }
    public void setColLocation(byte type, int location)
    {
        colOrder[type] = location;
    }
    public String getName() { return file.filename; }
    public void parseHdrAttributes(String line, int subtract, int nTokens)
    {
        String token = null;
        StringTokenizer stok = null;
        int cnt = 0;

        curveNames = new String[nTokens-1];
        attIndices = new int[nTokens-1];

        if(colOrder[TIME] == -1 && colOrder[DATE] == -1)
            hasHeader = false;

        stok = new StringTokenizer(line, delimiters);
        nTokens = stok.countTokens();
        boolean foundSrcTkn = false;
        for (int i = 0; i < nTokens; i++)
        {
          	token = stok.nextToken().trim();
            if(foundSrcTkn)
            {
                foundSrcTkn = false;
                continue;
            }
            if(token.equalsIgnoreCase("SRC"))
            {
                foundSrcTkn = true;
                continue;
            }
            token = token.replace(".", " ").trim();
          	colOrder[X] = checkIndex(token, i - subtract, colOrder[X], StsSensorKeywordIO.X_KEYWORDS, null);
           	colOrder[Y] = checkIndex(token, i - subtract, colOrder[Y], StsSensorKeywordIO.Y_KEYWORDS, null);
           	colOrder[Z] = checkIndex(token, i - subtract, colOrder[Z], StsSensorKeywordIO.Z_KEYWORDS, null);
            if(!hasHeader)
            {
                colOrder[TIME] = checkIndex(token, i - subtract, colOrder[TIME], StsSensorKeywordIO.TIME_KEYWORDS, new String[] {":"});
                colOrder[DATE] = checkIndex(token, i - subtract, colOrder[DATE], StsSensorKeywordIO.DATE_KEYWORDS, new String[] {"/"});
            }
           	if(((i - subtract )!= colOrder[TIME]) && ((i - subtract) != colOrder[DATE]))
           	{
           		attIndices[cnt] =  i - subtract;
           		if(colOrder[Z] ==  (i - subtract))
           			curveNames[cnt++] = StsSensorKeywordIO.DEPTH;
           		else if(colOrder[X] ==  (i - subtract))
           			curveNames[cnt++] = StsSensorKeywordIO.X;
           		else if(colOrder[Y] ==  (i - subtract))
           			curveNames[cnt++] = StsSensorKeywordIO.Y;
           		else
                {
                   if(hasHeader)
            		    curveNames[cnt++] = StsStringUtils.cleanString(token, invalidChars);
                   else
                        curveNames[cnt++] = "Att" + cnt;
                }
            }
        }

        if((colOrder[X] != -1) && (colOrder[Y] != -1) && (colOrder[Z] != -1))
            positionType = StsSensor.DYNAMIC;
        else
        	positionType = StsSensor.STATIC;
    }

    public boolean determineTimeType(StsModel model)
    {
        String line = null;
        StringTokenizer stok = null;
        try
        {
            BufferedReader bufRdr = new BufferedReader(new FileReader(file.getPathname()));
            for(int i=0; i<getNumHeaderRows(); i++)
            {
                bufRdr.readLine();
            }
            line = StsSensorFile.readLine(bufRdr);
            line = line.trim();
            line = StsStringUtils.detabString(line);
            line = StsStringUtils.deQuoteString(line);
            line = StsStringUtils.cleanLine(line, delimiters);
            stok = new StringTokenizer(line, delimiters);
            dateOrder = model.getProject().getDateOrder();
            if(timeTypeAndFormat(model, line, dateOrder, delimiters) == StsSensorKeywordIO.NO_TIME)
                return false;
            else
                return true;
        }
        catch(Exception ex)
        {
            StsException.outputException(ex, StsException.WARNING);
            return false;
        }
    }

    public String determineTimeType(StsModel model, BufferedReader bufRdr)
    {
        String line = null;
        StringTokenizer stok = null;
        String token = null;
        timeType = StsSensorKeywordIO.NO_TIME;

        while(timeType == StsSensorKeywordIO.NO_TIME)
        {
            line = StsSensorFile.readLine(bufRdr);
            line = line.trim();
            line = StsStringUtils.detabString(line);
            line = StsStringUtils.deQuoteString(line);
            line = StsStringUtils.cleanLine(line, delimiters);

            stok = new StringTokenizer(line, delimiters);
            int nTokens = stok.countTokens();
            if(model != null)
            {
                dateOrder = model.getProject().getDateOrder();
                timeType = timeTypeAndFormat(model, line, dateOrder, delimiters);
            }
            else
            {
                if((colOrder[DATE] != -1) && (colOrder[TIME] != -1))
                    timeType = StsSensorKeywordIO.TIME_OR_DATE;
                else
                    timeType = StsSensorKeywordIO.TIME_ONLY;
            }
        }
        return line;
    }

    static public StsSensorClass getSensorClass(StsModel model)
    {
        return (StsSensorClass)model.getCreateStsClass(StsSensor.class);
    }

    static public String readLine(BufferedReader bufRdr)
    {
    	try
    	{
    		String line = bufRdr.readLine();
    		while(line.startsWith("#") || line.trim().length() == 0)
    			line = bufRdr.readLine();
    		return line;
    	}
    	catch(Exception ex)
    	{
    		return null;
    	}

    }

	public int checkIndex(String token, int currentIndex, int index, String[] keywords, String[] keychars)
	{
        if(index == -1)
        {
            for(int j=0; j<keywords.length; j++)
            {
                if(token.equalsIgnoreCase(keywords[j]))
                    return currentIndex;
            }
            // Unable to find header word....possibly no header so see if token has key characters
            if(keychars != null)
            {
                for(int j=0; j<keychars.length; j++)
                {
                    if(token.contains(keychars[j]))
                        return currentIndex;
                }
            }
            return index;
        }
        else
        	return index;
	}

	public boolean positionToEndOfStage(BufferedReader bufRdr, int stage)
	{
		int num = 0;
    	try
    	{
    		while(num < stageOffsets[stage]-1)
    		{
    			StsSensorFile.readLine(bufRdr);
    			num++;
    		}
    		return true;
    	}
    	catch(Exception ex)
    	{
    		return false;
    	}
	}

	public boolean positionToStage(BufferedReader bufRdr, int stage)
	{
		int num = 0;
    	try
    	{
    		while(num < stageOffsets[stage-1])
    		{
    			StsSensorFile.readLine(bufRdr);
    			num++;
    		}
    		return true;
    	}
    	catch(Exception ex)
    	{
    		return false;
    	}
	}

	public int getLinesInAllStages()
	{
 		if(numStages > 1)
			return stageOffsets[numStages] - stageOffsets[0];
		else
			return stageOffsets[1];
	}

	public int getLinesInStage(int stage)
	{
		if(stage >= 1)
			return stageOffsets[stage] - stageOffsets[stage-1];
		else
			return stageOffsets[0];
	}

    public int getNumLinesInFile()
    {
        int lineNum = 0;
    	try
    	{
    		BufferedReader bufRdr = new BufferedReader(new FileReader(file.getPathname()));
    		String line = StsSensorFile.readLine(bufRdr); // Header
            if(!hasHeader) lineNum++;
    		while(line != null)
    		{
    			line = StsSensorFile.readLine(bufRdr);
    			if(line != null) lineNum++;
            }
    		return lineNum;
        }
        catch(Exception ex)
        {
        	return 0;
        }
    }

    public int[] getStageOffsets()
    {
        return stageOffsets;
    }
    public void setStageOffsets(int[] offsets)
    {
        stageOffsets = offsets;
    }

    public void setOverrideStages(boolean val)
    {
        overrideStages = val;
    }
    public int numStages(StsModel model)
    {
		long previousTime = -1, currentTime = -1;
        double previousStage = -1.f;
		String line = null;
		int nTokens = 0;
		String time = null, date = null;
		StringTokenizer stok = null;
		int lineNum = 0, stageIdx = -1, stage = -1;
        
        if(overrideStages)
            return numStages;

		stageOffsets = new int[5000];

    	try
    	{
    		BufferedReader bufRdr = new BufferedReader(new FileReader(file.getPathname()));
    		StsProject project = model.getProject();
            line = this.readLine(bufRdr);
    		if(hasHeader)
            {
    		    lineNum++;
                line = StsSensorFile.readLine(bufRdr);

                // Determine if stage column exists and use it
                for(int i=0; i<curveNames.length; i++)
                {
                    if(curveNames[i].toLowerCase().startsWith("stage"))
                    {
                        stageIdx = i;
                        break;
                    }
                }
            }
            // Use the stage column to determine stage limits
            if(stageIdx != -1)
            {
    		    while(line != null)
    		    {
                    lineNum++;
                    if(lineNum < numHeaderRows)
                    {
                        line = StsSensorFile.readLine(bufRdr);
                        continue;                        
                    }
                    if(!getAttributeValues(model, line))
                    {
                        line = StsSensorFile.readLine(bufRdr);
                        continue;
                    }

    			    if(currentValues[stageIdx] != previousStage)
    			    {
    					stageOffsets[numStages-1] = lineNum-1;
                        numStages++;
    				    previousStage = currentValues[stageIdx];
    			    }
                    line = StsSensorFile.readLine(bufRdr);
                }
                numStages--;
            }
            // No stage column so try to determine from time breaks
            else
            {
                stageOffsets[numStages-1] = lineNum;
    		    date = project.getDateFromLong(startTime);
    		    while(line != null)
    		    {
    			    lineNum++;
    			    stok = new StringTokenizer(line, delimiters);
    			    nTokens = stok.countTokens();
                    if(nTokens < 2)
                    {
                        line = StsSensorFile.readLine(bufRdr);
                        continue;
                    }
    			    for(int i=0; i<nTokens; i++)
    			    {
    				    String token = stok.nextToken();
    				    if(i == colOrder[TIME])
    					    time = token;

    				    if(i == colOrder[DATE])
    					    date = token;
    			    }

    			    if(previousTime != -1)
    			    {
    				    currentTime = StsSensorKeywordIO.getTime(timeType, time, date, startTime, project.getDateOrder());
    				    StsSensorClass sensorClass = getSensorClass(model);
                        float multiStageCriteria = sensorClass.getMultiStageCriteria();
                        if((currentTime - previousTime) > multiStageCriteria*1000)
    				    {
    					    numStages++;
    					    stageOffsets[numStages-1] = lineNum-1;
    				    }
    				    previousTime = currentTime;
    			    }
    			    else
    				    previousTime = StsSensorKeywordIO.getTime(timeType, time, date, startTime, project.getDateOrder());
                    line = StsSensorFile.readLine(bufRdr);
                }
            }
            StsMessageFiles.infoMessage("Found " + numStages + " stages in the file.");
    		return numStages;
        }
        catch(Exception ex)
        {
        	return numStages;
        }
        finally
        {
        	stageOffsets[numStages] = lineNum;   // lines in file
        	stageOffsets = (int[])StsMath.trimArray(stageOffsets, numStages+1);
        }
    }

    private boolean relativeCoordinates(String line, double criteria)
    {
        String xStg = null;
        StringTokenizer stok = new StringTokenizer(line, delimiters);
        int nTokens = stok.countTokens();

        for(int i=0; i<nTokens; i++)
        {
        	xStg = stok.nextToken();
        	if(i != colOrder[X])
        		continue;

        	// Add logic to determine if relative
        	double x = Double.parseDouble(xStg);
        	if(x < criteria)
        		return true;
        	else
        		return false;
        }
    	return false;
    }
    public void setNumStages(int num)
    {
        if(num != numStages)
            numStages = num;
        // ToDo: Need to reconfigure the stageOffsets to reasonable defaults based on increase or decrease of number of stages.
        

    }
    public long getStartTimeForStage(StsModel model, int stageNum)
    {
        try
        {
            BufferedReader bufRdr = new BufferedReader(new FileReader(file.getPathname()));
            positionToStage(bufRdr, stageNum);
            String line = bufRdr.readLine();
            getAttributeValues(model, line);
            return currentTime;
        }
        catch(Exception ex)
        {
        	return 0l;
        }
    }

    public long getEndTimeForStage(StsModel model, int stageNum)
    {
        try
        {
            BufferedReader bufRdr = new BufferedReader(new FileReader(file.getPathname()));
            positionToEndOfStage(bufRdr, stageNum);
            String line = bufRdr.readLine();
            getAttributeValues(model, line);
            return currentTime;
        }
        catch(Exception ex)
        {
        	return 0l;
        }
    }
    private byte timeTypeAndFormat(StsModel model, String line, int dateOrder, String delimiters)
    {
        String time = null, date = null;
        StringTokenizer stok = new StringTokenizer(line, delimiters);
        int nTokens = stok.countTokens();

        try
        {
            for(int i=0; i<nTokens; i++)
            {
        	    String token = stok.nextToken();
        	    if(i == colOrder[TIME])
        		    time = token.trim();

        	    if(i == colOrder[DATE])
        		    date = token.trim();
            }

            if(!isDateInString(time))
            {
                if((time.indexOf(":") < 0) && (time.indexOf(".") > 0))
                {
                    //StsMessageFiles.infoMessage("File contains unsupported time format (" + time + ").");
                    return StsSensorKeywordIO.ELAPSED_TIME;
                }
                else
                {
                    if(colOrder[DATE] == -1)
            	        return StsSensorKeywordIO.TIME_ONLY;
                    else
                    {
                        Calendar cal = CalendarParser.parse(date + " " + time, dateOrder , true);
                        startTime = cal.getTimeInMillis();
            	        return StsSensorKeywordIO.TIME_OR_DATE;
                    }
                }
            }
            else
            {
                Calendar cal = CalendarParser.parse(time, dateOrder , true);
                startTime = cal.getTimeInMillis();
                return StsSensorKeywordIO.TIME_AND_DATE;
            }
        }
        catch(Exception ex)
        {
            StsMessageFiles.errorMessage("Error analyzing time type in sensor file. Processing will abort.");
            //StsException.outputException("Error determining time type in sensor file: ", ex, StsException.WARNING);
            return StsSensorKeywordIO.NO_TIME;
        }
    }

    public boolean isDateInString(String dt)
    {
    	String[] dateStrings = new String[] {"/","-","Jan","Feb","Mar","Apr",
    			                     "May","Jun","Jul","Aug","Sep","Oct","Nov",
    			                     "Dec","Mon","Tue","Wed","Thu","Fri","Sat","Sun"};
    	if((dt.indexOf(".") >= 0) && (dt.indexOf(".") < dt.indexOf(":")))  // May be . in date and time so non-standard check required
    		return true;
    	for(int i=0; i<dateStrings.length; i++)
    		if(dt.indexOf(dateStrings[i]) >= 0)
    			return true;
    	return false;
    }

    public void setNumHeaderRows(int n)
    {
        numHeaderRows = n;
    }
    public int getNumHeaderRows() { return numHeaderRows; }
    public String toString() { return file.filename; }
    public void setRelativeX(double x) { this.relativeX = x; }
    public void setRelativeY(double y) { this.relativeY = y; }
    public void setRelativeZ(double z) { this.relativeZ = z; }
    public double getRelativeX() { return relativeX; }
    public double getRelativeY() { return relativeY; }
    public double getRelativeZ() { return relativeZ; }

    public void setStaticX(double x) { this.staticX = x; }
    public void setStaticY(double y) { this.staticY = y; }
    public void setStaticZ(double z) { this.staticZ = z; }
    public double getStaticX() { return staticX; }
    public double getStaticY() { return staticY; }
    public double getStaticZ() { return staticZ; }

    public long getStartTime() { return startTime; }
    public void setStartTime(long time) { this.startTime = time; }
}
