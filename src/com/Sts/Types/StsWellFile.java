

//Title:        S2S: Seismic-to-simulation
//Version:
//Copyright:    Copyright (c) 2001
//Author:       Tom Lasseter
//Company:      4D Systems LLC
//Description:  Web-enabled Integrated Interpretation System

package com.Sts.Types;

import com.Sts.Actions.Import.StsSensorKeywordIO;
import com.Sts.Actions.Import.StsWellKeywordIO;
import com.Sts.DBTypes.StsLogVector;
import com.Sts.DBTypes.StsSensor;
import com.Sts.DBTypes.StsSensorClass;
import com.Sts.IO.StsFile;
import com.Sts.MVC.StsMessageFiles;
import com.Sts.MVC.StsModel;
import com.Sts.MVC.StsProject;
import com.Sts.Utilities.DateTime.CalendarParser;
import com.Sts.Utilities.StsException;
import com.Sts.Utilities.StsMath;
import com.Sts.Utilities.StsStringUtils;

import java.io.BufferedReader;
import java.io.FileReader;
import java.util.Calendar;
import java.util.StringTokenizer;

public class StsWellFile
{
	public byte timeType = StsSensorKeywordIO.NO_TIME;
	public long currentTime = 0L;
    public long startTime = 0L;
    private int dateOrder = CalendarParser.DD_MM_YY;

    public double xOrigin = 0.0f;
    public double yOrigin = 0.0f;
    public float nullValue = -999.0f;

	public StsFile file = null;
	public int[] attIndices = null;
	public String[] curveNames = null;
	public String[] validCurves = null;
	public double[] currentValues = null;

    public int numCols = 0;
    public int numLines = 0;
    public String wellName = "None";
    public String delimiters = " ";

    transient BufferedReader bufRdr = null;

    public static byte X = 0;
    public static byte Y = 1;
    public static byte Z = 2;
    public static byte TIME = 3;
    public static byte DATE = 4;
    int[] colOrder = new int[] { -1, -1, -1, -1, -1 };

    static public char[] invalidChars = new char[] {'#',' ',':',';','"','\'',',', ')', '(', '[', ']', '{', '}', '=', '*', '@', '?','/', '\\', '$', '%', '&', '*'};

    public StsWellFile(StsFile file)
    {
    	this.file = file;
    }

    public void analyzeFile()
    {
    	analyzeFile(null);
    }

    public boolean getAttributeValues(StsModel model, String line, int dateOrder)
    {
        this.dateOrder = dateOrder;
        return getAttributeValues(model, line);
    }

    public boolean getAttributeValues(StsModel model, String line)
    {
        StringTokenizer stok = null;
        String token = null;
        double value = 0.0f;

        String timeToken = null;
        String dateToken = null;
        try
        {
            line = line.trim();
            line = StsStringUtils.detabString(line);
            stok = new StringTokenizer(line, " ");
            int nTokens = stok.countTokens();
            if(nTokens < curveNames.length)
                return false;
            currentValues = new double[nTokens];
            validCurves = new String[nTokens];
            int cnt = 0;
            for (int i = 0; i < nTokens; i++)
            {
                token = stok.nextToken();
                if(colOrder[TIME] == i)
                	timeToken = token;
                else if(colOrder[DATE] == i)
                	dateToken = token;
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
                    validCurves[cnt] = curveNames[cnt];
                    cnt++;
                }
            }
            if(model != null)
                dateOrder = model.getProject().getDateOrder();
            if(timeToken != null)
            {
                String dateTime = dateToken;

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
            StsMessageFiles.logMessage("Failed to find or parse time/date column(s) in file - check the DATE ORDER parameter");
            return false;
        }
        return true;
    }

    public String readLine()
    {
        try
        {
            return bufRdr.readLine().trim();
        }
        catch(Exception ex)
        {
            StsException.outputException("StsWellFile.readLine():", ex, StsException.WARNING);
            return null;
        }
    }
    /* read the well name */
    public boolean analyzeFile(StsModel model)
    {
        String line = null;
        StringTokenizer stok = null;
        String token = null;
        delimiters = ";,";
        try
        {
            bufRdr = new BufferedReader(new FileReader(file.getPathname()));
            while(true)
			{
			    line = bufRdr.readLine().trim();
                line = StsStringUtils.detabString(line);
		        if(line.endsWith(StsWellKeywordIO.WELL_NAME))
			        wellName =  new String(bufRdr.readLine().trim());
			    else if(line.indexOf(StsWellKeywordIO.ORIGIN) >= 0)  // is origin keyword there?
			    {
				    boolean xyOrder = true;
				    if (line.indexOf(StsWellKeywordIO.YX) >= 0)
                        xyOrder = false;  // determine x-y order
			        line = bufRdr.readLine().trim();  // get the next line

				    // tokenize the x-y values and convert to a point object
				    stok = new StringTokenizer(line);
				    xOrigin = Double.valueOf(stok.nextToken()).doubleValue();
				    yOrigin = Double.valueOf(stok.nextToken()).doubleValue();

				    if(!xyOrder)
				    {
					    double temp = xOrigin;
					    xOrigin = yOrigin;
					    yOrigin = temp;
				    }
                }
			    else if(line.endsWith(StsWellKeywordIO.NULL_VALUE))
			    {
				    line = bufRdr.readLine().trim();  // get the next line
				    stok = new StringTokenizer(line);
				    nullValue = Float.valueOf(stok.nextToken()).floatValue();
			    }
			    else if(line.endsWith(StsWellKeywordIO.CURVE))
			    {
				    String[] curveNames = new String[0];
				    line = bufRdr.readLine().trim();
                    int currentIdx = 0;
				    while(!StsWellKeywordIO.lineHasKeyword(line, StsWellKeywordIO.VALUE))
				    {
                        colOrder[X] = StsWellKeywordIO.checkIndex(line, currentIdx, colOrder[X], StsWellKeywordIO.X_KEYWORDS);
                        colOrder[Y] = StsWellKeywordIO.checkIndex(line, currentIdx, colOrder[Y], StsWellKeywordIO.Y_KEYWORDS);
                        colOrder[Z] = StsWellKeywordIO.checkIndex(line, currentIdx, colOrder[Z], StsWellKeywordIO.Z_KEYWORDS);
                        colOrder[TIME] = StsWellKeywordIO.checkIndex(line, currentIdx, colOrder[TIME], StsWellKeywordIO.TIME_KEYWORDS);
                        colOrder[DATE] = StsWellKeywordIO.checkIndex(line, currentIdx, colOrder[DATE], StsWellKeywordIO.DATE_KEYWORDS);
			    	    if((currentIdx != colOrder[TIME]) && (currentIdx != colOrder[DATE]))
                            curveNames = (String[])StsMath.arrayAddElement(curveNames, line);
                        currentIdx++;
					    line = bufRdr.readLine().trim();
				    }
                    this.curveNames = curveNames;
                    return true;
			    }
			    else if(line.endsWith(StsWellKeywordIO.VALUE))
                    return false;
			    else
				    return false;
            }
        }
        catch(Exception e)
        {
            StsException.outputException("StsWellKeywordIO.readFileHeader() failed.",
				e, StsException.WARNING);
			return false;
        }
    }

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

    public int getNumLinesInFile()
    {
        int lineNum = 0;
    	try
    	{
    		BufferedReader bufRdr = new BufferedReader(new FileReader(file.getPathname()));
    		String line = bufRdr.readLine();
    		while(line != null)
    		{
                if(line.endsWith(StsWellKeywordIO.VALUE))
                {
                    line = bufRdr.readLine();
                    while(line != null)
                    {
    			        lineNum++;
                        line = bufRdr.readLine();
                    }
                }
                else
                    line = bufRdr.readLine();
            }
    		return lineNum;
        }
        catch(Exception ex)
        {
        	return 0;
        }
    }

    public String toString() { return file.filename; }

    public long getStartTime() { return startTime; }
    public void setStartTime(long time) { this.startTime = time; }
}