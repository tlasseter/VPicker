

//Title:        S2S: Seismic-to-simulation
//Version:
//Copyright:    Copyright (c) 2001
//Author:       Tom Lasseter
//Company:      4D Systems LLC
//Description:  Web-enabled Integrated Interpretation System

package com.Sts.Actions.Wizards.FracSetMVLoad;

import com.Sts.Actions.Import.*;
import com.Sts.IO.*;
import com.Sts.MVC.*;
import com.Sts.UI.*;
import com.Sts.Utilities.*;

import java.io.*;
import java.util.*;

public class StsMVFracSetFile
{
	public boolean relativeCoordinates = false;
	public StsFile file = null;
	public String[] fieldNames = null;
	public String[] validCurves = null;
	public double[] currentValues = null;
	public String name = null;
	public int numHeaderLines = 0;
	public int numAttributes = 0;
    public long currentTime = 0L;
	public int timeIdx = -1, dateIdx = -1;;
    public byte timeType = StsSensorKeywordIO.NO_TIME;

	static public final String[] TYPE_KEYWORDS = { "Real", "Integer" };
	static public final String[] UNITS_KEYWORDS = { "Metre", "Meter", "Degree", "Dimensionless", "Feet" };
    static public final String[] PROPERTY_KEYWORDS = { "East","North","Z","dip","azimuth","strike",
    	"aperture","aspect_ratio","length", "X", "Y", "Depth"};
    static public final String[] PROP_TYPE_KEYWORDS = { "xy","elevation","dimension","dimensionless"};

    public String[] currentTypes = null;
    public String[] currentUnits = null;
    BufferedReader bufRdr = null;
	
    public StsMVFracSetFile(StsFile file)
    {
    	this.file = file;
    	String fn = file.getFilename();
    	StringTokenizer sTokens = new StringTokenizer(fn, ".");
        if(sTokens.hasMoreTokens())
           name = sTokens.nextToken();
    }
    
    public boolean getPropertyValues(StsModel model, int dateOrder)
    {
    	if (bufRdr == null) analyzeFile();
        StringTokenizer stok = null;
        String token = null;
        double value = 0.0f;
        String timeToken = null;
        String dateToken = null;
        String line;
        try
        {   
        	line = StsMVFracSetFile.readLine(bufRdr);
        	if (line == null) return false;
            timeType = timeTypeAndFormat(model, line);
            if((timeType == StsSensorKeywordIO.NO_TIME) && (timeIdx != -1))
                return false;

            stok = new StringTokenizer(line, ",");
            int nTokens = stok.countTokens();
            currentValues = new double[nTokens];
            validCurves = new String[nTokens];
            int cnt = 0;
            for (int i = 0; i < nTokens; i++)
            {
                token = stok.nextToken();
                try 
                {
                    if(i == timeIdx)
                    {
                        timeToken = token;
                        continue;
                    }
                    if(i == dateIdx)
                    {
                        dateToken = token;
                        continue;
                    }
                    else
                	    value = Double.parseDouble(token);
                }
                catch (Exception e) 
                { 
                	continue;      // Not be a number so skip it
                }
                currentValues[cnt] = value;
                validCurves[cnt] = fieldNames[cnt];
                cnt++;           	
            }
            if(timeType != StsSensorKeywordIO.NO_TIME)
                currentTime = StsSensorKeywordIO.getTime(timeType, timeToken, dateToken, 0, dateOrder);
            currentValues = (double[])StsMath.trimArray(currentValues, cnt);
            validCurves = (String[])StsMath.trimArray(validCurves, cnt);

        }
        catch(Exception e)
        {
            StsMessageFiles.logMessage("Failed to find time column in file: " + file.getFilename());
            return false;
        }
        return true; 		    	
    } 
    
    /*
     * read the header lines
     */
    public void analyzeFile()
    {
        String line = null;
        StringTokenizer stok = null;
        String token = null;
        try
        {
            bufRdr = new BufferedReader(new FileReader(file.getPathname()));
            boolean readingHeader = true;
            boolean hasAttributes = false;
            while (readingHeader)
            {
	            line = StsMVFracSetFile.readLine(bufRdr);
	            numHeaderLines++;
	            
	            stok = new StringTokenizer(line, ",");
	            int nTokens = stok.countTokens();
	            
	            String[] currentKeywords = new String[nTokens];
	            for (int i = 0; i < nTokens; i++)
	            {
	            	token = stok.nextToken().trim();
	            	currentKeywords[i] = token.toString();
	            }
	            if (currentKeywords[0].equalsIgnoreCase("Midland Valley")) continue;
	            if (currentKeywords[0].equalsIgnoreCase("Vertex Data"))
	            {
	            	if (nTokens > 9) hasAttributes = true;
	            	numAttributes = nTokens - 9;
	            	continue;
	            }
	            if (checkKeyword(currentKeywords[0], TYPE_KEYWORDS))
	            {
	            	currentTypes = new String[nTokens];
	            	for (int i=0; i<nTokens; i++)
	            		currentTypes[i] = currentKeywords[i];
	            	continue;
	            }
	            if (checkKeyword(currentKeywords[0], PROPERTY_KEYWORDS))
	            {
                    timeIdx = locateColumn(currentKeywords, StsSensorKeywordIO.TIME_KEYWORDS);
            	    dateIdx = locateColumn(currentKeywords, StsSensorKeywordIO.DATE_KEYWORDS);
	            	if((timeIdx != -1) && (dateIdx != -1))
                        fieldNames = new String[nTokens-2];
                    else if((timeIdx != -1) || (dateIdx != -1))
                        fieldNames = new String[nTokens-1];
                    else
                        fieldNames = new String[nTokens];

                    int aCnt = 0;
	            	for (int i=0; i<nTokens; i++)
                    {
                        if((i != timeIdx) && (i != dateIdx))
                        {
	            		    fieldNames[aCnt] = currentKeywords[i];
                            aCnt++;
                        }
                    }
	            	continue;
	            }
	            if (currentKeywords[0].startsWith("#")) continue;	
	            if ( checkKeyword(currentKeywords[0], UNITS_KEYWORDS))
	            {
	            	currentUnits = new String[nTokens];
	            	for (int i=0; i<nTokens; i++)
	            		currentUnits[i] = currentKeywords[i];
	            	continue;
	            }
	            if ( checkKeyword(currentKeywords[0], PROP_TYPE_KEYWORDS))
	            {
	            	readingHeader = false;
	            	continue;
	            }
            }
            if(timeIdx != -1) numAttributes--;
            if(dateIdx != -1) numAttributes--;
        }
         catch(Exception e)
        {
            StsMessageFiles.logMessage("Failed to find time column in file: " + file.getFilename());
            return;
        }
        return; 		    	
    }

//    static public StsMVFractureSetClass getFracSetClass(StsModel model)
//    {
//        return (StsMVFractureSetClass)model.getCreateStsClass(StsMVFractureSet.class);
//    }

    public int locateColumn(String[] tokens, String[] keywords)
    {
        for(int i=0; i<tokens.length; i++)
        {
            for(int j=0; j<keywords.length; j++)
            {
                if(tokens[i].equalsIgnoreCase(keywords[j]))
                    return i;
            }
        }
        return -1;
    }
    static public String readLine(BufferedReader bufRdr)
    {   	
    	try
    	{
    		String line = bufRdr.readLine();
//    		while(line.startsWith("#") || line.trim().length() == 0)
//    			line = bufRdr.readLine();
    		return line;
    	}
    	catch(Exception ex)
    	{
    		return null;
    	}
    	
    }
    
    public boolean checkKeyword(String token, String[] keywords)
	{
        for(int j=0; j<keywords.length; j++)
        {
            if(token.equalsIgnoreCase(keywords[j]))
                return true;
        }
        return false;
	}
    
	public int checkIndex(String token, int currentIndex, int index, String[] keywords)
	{
        if(index == -1)
        {
            for(int j=0; j<keywords.length; j++)
            {
                if(token.equalsIgnoreCase(keywords[j]))
                    return currentIndex;
            }
            return index;
        } 
        else
        	return index;
	}
	

    private byte timeTypeAndFormat(StsModel model, String line)
    {
        String time = null, date = null;
        StringTokenizer stok = new StringTokenizer(line, ",");
        int nTokens = stok.countTokens();

        if((timeIdx == -1) && (dateIdx == -1))
           return StsSensorKeywordIO.NO_TIME;

        for(int i=0; i<nTokens; i++)
        {
        	String token = stok.nextToken();
        	if(i == timeIdx)
        		time = token;

        	if(i == dateIdx)
        		date = token;
        }

        if(!isDateInString(time))
        {
           if(time.indexOf(":") < 0)
           {
               new StsMessage(model.win3d, StsMessage.ERROR, "File contains unsupport time format (" + time + ").");
               return StsSensorKeywordIO.NO_TIME;
           }
           else
           {
               if(dateIdx == -1)
                  return StsSensorKeywordIO.TIME_ONLY;
               else
            	  return StsSensorKeywordIO.TIME_OR_DATE;
           }
        }
        else
        {
           return StsSensorKeywordIO.TIME_AND_DATE;
        }
    }
     
    private boolean relativeCoordinates(String line, double criteria)
    {  
//        String xStg = null;
//        StringTokenizer stok = new StringTokenizer(line, ",");
//        int nTokens = stok.countTokens();
//        
//        for(int i=0; i<nTokens; i++)
//        {
//        	xStg = stok.nextToken();
//        	if(i != xIdx)
//        		continue;  
//        	
//        	// Add logic to determine if relative
//        	double x = Double.parseDouble(xStg);
//        	if(x < criteria)
//        		return true;
//        	else
//        		return false;
//        }
    	return false;
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
    
}
