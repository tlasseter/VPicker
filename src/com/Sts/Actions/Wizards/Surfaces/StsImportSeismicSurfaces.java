
//Title:        S2S: Seismic-to-simulation
//Version:
//Copyright:    Copyright (c) 2001
//Author:       Tom Lasseter
//Company:      4D Systems LLC
//Description:  Web-enabled Integrated Interpretation System

package com.Sts.Actions.Wizards.Surfaces;

import com.Sts.DBTypes.*;
import com.Sts.IO.*;
import com.Sts.MVC.*;
import com.Sts.UI.Progress.*;
import com.Sts.UI.*;
import com.Sts.Utilities.*;

import java.util.*;

public class StsImportSeismicSurfaces extends StsImportSurfaces
{
    static float minLine, maxLine;
    static float minXline, maxXline;
    static float lineInc, xlineInc;
    static int nCols;
    static int nRows;
    static boolean isTrueDepth;
    static double xOrigin, yOrigin;
    static double xInc, yInc;
    static double angle;
    static boolean isCrosslineCCW;
    static double modelXOrigin, modelYOrigin;
    static boolean hasAttributes = false;
    static int horzIndex = -1;

    static public final String[] Z_KEYWORDS = { "DEPTH", "Z" };
    static public final String[] T_KEYWORDS = { "TIME", "T", "TM", "MSEC", "SEC", "SECONDS" };
    static public final String[] X_KEYWORDS = { "X", "EASTING", "DX" };
    static public final String[] Y_KEYWORDS = { "Y", "NORTHING", "DY" };
    static public final String[] INLINE_KEYWORDS = { "INLINE", "IN-LINE", "IL", "ILINE", "LINE" };
    static public final String[] XLINE_KEYWORDS = { "CROSSLINE", "CROSS-LINE", "XL", "XLINE", "TRACE" };
    static public final String[] HORIZON_KEYWORDS = { "HRZ", "SURFACE", "HORZ", "HORIZON", "SRF" };

    // Column Order
    public static byte X = 0;
    public static byte Y = 1;
    public static byte T = 2;
    public static byte ILINE = 3;
    public static byte XLINE = 4;
    static int[] colOrder = new int[] { X, Y, T, ILINE, XLINE };
    static int numColsPerRow = 5;
    static int skippedLines = 0;
    static int linesInFile = 0;
    static int filePosition = 0;
    static boolean scanned = false;

    static byte XYZIX = 0;
    static byte XYZ = 1;
    static byte gridType = XYZIX;
    static int numInlines = -1;
    static int numXlines = -1;

    static final float largeFloat = StsParameters.largeFloat;

    static public void setColLocation(int val, int col)
    {
            colOrder[val] = col;
    }
    static public void setNumCols(int numCols)
    {
        numColsPerRow = numCols;
    }
    static public void setSkippedLines(int headerLines)
    {
        skippedLines = headerLines;
    }
    
    static public void determineColumnOrder(String hdrLine, boolean isDepth)
    {
    	int col = 0;
        String[] TD_KEYWORDS = Z_KEYWORDS;
        if(!isDepth)
            TD_KEYWORDS = T_KEYWORDS;
    	String[][] typeKeywords =  { X_KEYWORDS, Y_KEYWORDS, TD_KEYWORDS, INLINE_KEYWORDS, XLINE_KEYWORDS };

    	hdrLine = StsStringUtils.detabString(hdrLine);
    	StringTokenizer stok = new StringTokenizer(hdrLine, ", ");

        // Gridded seismic may have all 5 columns or just XYZ.
        if(stok.countTokens() == 2) // NumInline and numXlines
        {
            String token = stok.nextToken();
        	token = StsStringUtils.deQuoteString(token);
            numInlines = Integer.valueOf(token);
            token = stok.nextToken();
        	token = StsStringUtils.deQuoteString(token);
            numXlines = Integer.valueOf(token);
            gridType = XYZ;
            return;
        }
        else
            gridType = XYZIX;

        // Initial all columns to -1 so user will now which columns were not found
        for(int i=0; i<colOrder.length; i++)
            colOrder[i] = -1;

        // Determine column locations for key values.
    	while(stok.hasMoreTokens())
    	{
    		String token = stok.nextToken();
        	token = StsStringUtils.deQuoteString(token);    		
    		for(int j=0; j<typeKeywords.length; j++)
    		{
    			if(isType(token, typeKeywords[j]))
                {
    				colOrder[j] = col;
                    break;
                }
    		}
			col++;
    	}
    }
    
    static public boolean isType(String token, String[] strings)
    {
    	for(int i=0; i<strings.length; i++)
    	{
    		if(token.equalsIgnoreCase(strings[i]))
    			return true;
    	}
    	return false;
    }

    static public void initializeFile()
    {
        scanned = false;
        filePosition = 0;
    }
    /** Scan the first line in the file to determine if the units are reasonable in x,y,t,line,crossline
     * The line and crossline correspond to row and column of grid
     * If not every line and crossline picked, discern average pick interval and
     * build a new grid with these picked lines/crosslines as new rows/columns
     */
    static public boolean scanSurface(StsModel model, StsAbstractFile file, String name, float datum,
                                      byte hUnits, byte vUnits, byte tUnits, boolean concatLines, boolean isTrueDepth_)
	{  
        double x, y, t;
        float nLine = 0, nCrossline = 0;
        String[] tokens;

        if(scanned == true)
            return true;
    	StsAsciiFile asciiFile = new StsAsciiFile(file);
        if(!asciiFile.openReadWithErrorMessage()) return false;

        datumShift = datum;
        scanned = true;
        isTrueDepth = isTrueDepth_;
        linesInFile = 0;
        
        horizontalUnits = hUnits;
        if(StsImportSeismicSurfaces.isTrueDepth)
            verticalUnits = vUnits;
        else
            verticalUnits = tUnits;
        
        try
        {
            int nFileLinesRead = 0;
        	boolean unitsChecked = false;
            nLine = 1;
            nCrossline = 0;
            while((tokens = asciiFile.getTokens(new String[] {" ",","})) != null)
            {
                //Skip over user defined header
                if(nFileLinesRead < skippedLines)
                {
                    nFileLinesRead++;
                    if((nFileLinesRead == skippedLines) && (gridType != XYZ)) // Assume last header line is column names
                    {
                        // Concatenate lines to reach num of columns specified.
                        tokens = concatenateLines(tokens, asciiFile, concatLines);

                        // Determine if a Horizon Column exists
                        for(int i=0; i<tokens.length; i++)
                        {
                            String aname = StsStringUtils.detabString(StsStringUtils.deQuoteString(tokens[i]));
                            // Determine if a horizon column exists
                            for(int ii=0; ii<HORIZON_KEYWORDS.length; ii++)
                            {
                                if(aname.equalsIgnoreCase(HORIZON_KEYWORDS[ii]))
                                {
                                    horzIndex = i;
                                    break;
                                }
                            }
                        }
                    }
                    linesInFile++;
                    continue;
                }
                
                // Concatenate lines to reach num of columns specified.
                tokens = concatenateLines(tokens, asciiFile, concatLines);
                if(tokens == null)
                   continue;

                int nTokens = tokens.length;
                if(nTokens > numColsPerRow)
                	nTokens = numColsPerRow;
                
                x = ((Double.parseDouble(tokens[colOrder[X]])) - modelXOrigin);
                y = ((Double.parseDouble(tokens[colOrder[Y]])) - modelYOrigin);
                t = Double.parseDouble(tokens[colOrder[T]]) + datumShift;
                if(gridType == XYZIX)
                {
                    nLine = (int)Float.parseFloat(tokens[colOrder[ILINE]]);
                    nCrossline = (int)Float.parseFloat(tokens[colOrder[XLINE]]);
                }
                else
                {
                    if(nCrossline == numXlines)
                    {
                        nCrossline = 1;
                        nLine++;
                    }
                    else
                        nCrossline++;
                }
                // Validate the coordinates and grid against the project.
                if(StsImportSeismicSurfaces.isTrueDepth)
                {
                	;  // Not sure how to run an all inclusive check on ft versus meters.
                }
                else if(!unitsChecked)
                {
                	if((verticalUnits == StsParameters.TIME_MSECOND) && (t < 20))
                    {
                        if(StsYesNoDialog.questionValue(model.win3d, "Selected file " + asciiFile.getFilename() + " appears to have time in seconds, you selected milliseconds.\n\n Do you want to change it to seconds?\n"))
                			verticalUnits = StsParameters.TIME_SECOND;
                        unitsChecked = true;
                    }
                    if((verticalUnits == StsParameters.TIME_SECOND) && (t > 20))
                    {
                        if(StsYesNoDialog.questionValue(model.win3d, "Selected file " + asciiFile.getFilename() + " appears to have time in milliseconds, you selected seconds.\n\n Do you want to change it to milliseconds?\n"))
                			verticalUnits = StsParameters.TIME_MSECOND;
                        unitsChecked = true;
                    }
                }
                linesInFile++;
            }
            asciiFile.close();
        }
        catch(Exception e)
       	{
        	new StsMessage(model.win3d, StsMessage.WARNING, "File read error during scan of " + asciiFile.getFilename() + "\n" + "Data: " + asciiFile.getLine());
            if(asciiFile != null) asciiFile.close();
            return false;
        }    	    	
    	return true;
	}

    static private String[] concatenateLines(String[] toks, StsAsciiFile file, boolean concat)
    {
        try
        {
            if(concat)
            {
                while(toks.length < numColsPerRow)
                {
                   System.out.println("Found " + toks.length + " tokens....Needed " + numColsPerRow);
                   String[] newToks = file.getTokens(new String[] {" ",","});
                   if(newToks == null)
                       return null;
                   toks = (String[]) StsMath.arrayAddArray(toks, newToks);
                   if(toks.length > numColsPerRow)
                       StsMath.trimArray(toks,numColsPerRow);
                }
            }
            return toks;
        }
        catch (Exception ex)
        {
             StsException.outputException("Failed to concat lines together to meet user specified number of columns.", ex, StsException.WARNING);
             return null;
        }
    }
    /** Read in x,y,t,line,crossline
     * The line and crossline correspond to row and column of grid
     * If not every line and crossline picked, discern average pick interval and
     * build a new grid with these picked lines/crosslines as new rows/columns
     */
    static public StsSurface createSurface(StsModel model, StsAbstractFile file, String name, boolean concatLines, StsProgressPanel panel)
	{
        StsAsciiFile asciiFile = null;
        StsSurface newSurface;
        float[] lineNumbers, crosslineNumbers;
        SeismicPoint[] points;
        int arraySizeInc = 100000;
        float[][] zValues;
        float[][][] attributes;
        SeismicPoint[][] pointArray;
        boolean hasNulls = false;
        double x, y, t;
        float nLine = 0, nCrossline = 0;
        String hzName = "none";

        asciiFile = new StsAsciiFile(file);
        if(debug) asciiFile.debug = true;

        String fullPathname = file.getURLPathname();

        // see if we already have this surface read in
        StringBuffer surfaceName = new StringBuffer();  // object to receive surface name
        String filename = file.getFilename();
        parseFilename(filename, surfaceName);  // get surface name

      	StsMessageFiles.logMessage("Loading grid from: " + fullPathname + " ...");

        byte zDomain = StsParameters.TD_TIME;
		if(isTrueDepth)
            zDomain = StsParameters.TD_DEPTH;

        modelXOrigin = model.getXOrigin();
        modelYOrigin = model.getYOrigin();

        int nLinesRead = 0;
        String[] tokens;
        TreeSet sortedLines, sortedCrosslines;

        points = new SeismicPoint[arraySizeInc];

        int nAttributes = 0;
        String[] attributeNames = new String[0];
        float[] attMaxs = null;
        float[] attMins = null;
        boolean hasAttributes = false;

        minLine = StsParameters.largeFloat;
        maxLine = -StsParameters.largeFloat;
        minXline = StsParameters.largeFloat;
        maxXline = -StsParameters.largeFloat;

        sortedLines = new TreeSet();
        sortedCrosslines = new TreeSet();
        long fileLength = asciiFile.length();
        int nFileLinesRead = 0;
        
        if(filePosition == linesInFile)
            return null;

        if(horzIndex != -1)
            name = "none";

        // Get a color from the basic spectrum
        StsColor color = model.getCurrentSpectrumColor("Basic");
        int nPoints = 0;

        if(!asciiFile.openReadWithErrorMessage()) return null;

        try
        {
            // progressBar has an interval for each surface
            // allocate 60% here for reading ascii file
            long bytesRead = 0;
            panel.progressBar.setSubInterval(0.0f, 0.6f, fileLength);
            panel.appendLine("Reading ASCII seismic file....");
            boolean checked = false;
            float[] xy = new float[2];
        	double[] lattributes = null;
            nLine = 1.0f;
            nCrossline = 0.0f;
            while((tokens = asciiFile.getTokens(new String[] {" ",","})) != null)
            {
                if(nFileLinesRead < skippedLines)
                {
                    nFileLinesRead++;
                    if((nFileLinesRead == skippedLines) && (gridType != XYZ)) // Assume last header line is column names
                    {
                         // Concatenate lines to reach num of columns specified.
                        tokens = concatenateLines(tokens, asciiFile, concatLines);

                    	if(numColsPerRow != tokens.length) continue;  // Obviously not a column header
                    	{
                    		for(int i=0; i<tokens.length; i++)
                    		{
                    			if(StsMath.arrayContains(colOrder, i))
                    				continue;
                                String aname = StsStringUtils.detabString(StsStringUtils.deQuoteString(tokens[i]));
                    			attributeNames = (String[])StsMath.arrayAddElement(attributeNames, aname);
                    			nAttributes++;
                    		}
                    	}
                    }
                    continue;
                }
                nFileLinesRead++;
                tokens = concatenateLines(tokens, asciiFile, concatLines);
                if(tokens == null)
                   continue;

                // Position to the correct location in the file.
                if(nFileLinesRead <= filePosition)
                   continue;

                // Allocate the attribute arrays.
                lattributes = new double[nAttributes];
                if(nLinesRead == 0)
                {
                	nAttributes = tokens.length - 5;
                	if(nAttributes > 0)
                	{
                		hasAttributes = true;
                		lattributes = new double[nAttributes];
                		if(attributeNames.length == 0)
                		{
                			attributeNames = new String[nAttributes];
                			for(int i=0;i<nAttributes; i++)
                				attributeNames[i] = new String("Att00" + i);
                		}
                	}
                    else
                        nAttributes = 0;
                }
                nLinesRead++;
                int nTokens = tokens.length;
                if(nTokens > numColsPerRow)
                	nTokens = numColsPerRow;
                
                if(debug)
                {
                    System.out.print("    tokens:");
                    for(int n = 0; n < nTokens; n++)
                        System.out.print(" " + tokens[n]);
                    System.out.println();
                }
                x = ((Double.parseDouble(tokens[colOrder[X]])) - modelXOrigin);
                y = ((Double.parseDouble(tokens[colOrder[Y]])) - modelYOrigin);
                t = Double.parseDouble(tokens[colOrder[T]]);

                if(gridType == XYZ)
                {
                    // Determine the inline and crossline
                    StsObject[] objects = model.getObjectList(StsSeismicVolume.class);
                    if(objects.length != 0)
                    {
                        nCrossline = (int)((StsSeismicVolume)objects[0]).getColCoor((float)x);
                        nLine = (int) ((StsSeismicVolume)objects[0]).getRowCoor((float)y);
                    }
                    else if(numXlines != -1)
                    {
                        if(nCrossline == numXlines)
                        {
                            nCrossline = 1;
                            nLine++;
                        }
                        else
                            nCrossline++;
                     }
                    else
                    {
                        StsMessageFiles.errorMessage("No seismic grid defined in Project and no inline, crossline columns in file. Unable to load.");
                        return null;
                    }
                }
                else
                {
                        nLine = (int)Float.parseFloat(tokens[colOrder[ILINE]]);
                        nCrossline = (int)Float.parseFloat(tokens[colOrder[XLINE]]);
                }

                if((float)t == model.getProject().getMapGenericNull())
                    continue;
                t = t + datumShift;
                
                //System.out.println("Inline= " + nLine + " Crossline= " + nCrossline);
                // Parse out the horizon name if column exists

                if(horzIndex != -1)
                {
                    String token = StsStringUtils.detabString(StsStringUtils.deQuoteString(tokens[horzIndex]));
                    if(!name.equalsIgnoreCase("none") && !name.equalsIgnoreCase(token))
                    {
                        panel.appendLine("Found start of new surface (" + token + "), will complete processing of " + name + " first.");
                        nFileLinesRead--;
                        break;
                    }
                    name = token;
                }

                // Parse out any other columns as attributes.
                int cnt = 0;
                if(nAttributes > 0) 		    // Do attribute columns exist?
                {
                	for(int i=0; i<tokens.length; i++)
                	{
                		if(StsMath.arrayContains(colOrder, i))
                			continue;
                		try 
                		{ 
                			lattributes[cnt] = Double.parseDouble(tokens[i]);
                		}
                		catch(Exception ex) 
                		{ 
                			lattributes[cnt] = 0.0f;
                		}
                		cnt++;
                	}
                	if(cnt != nAttributes)
                	{
                		for(int i=cnt; i<nAttributes; i++)
                			lattributes[i] = 0.0f;
                	}
                }
                
                if(nPoints < 1000 && !checked)
                {
                    if (x < -largeFloat || x > largeFloat)
                    {
                        if(!valueOK(model, "File line " + nLinesRead + " X: " + tokens[colOrder[X]])) return null;
                        checked = true;
                    }
                    if (y < -largeFloat || y > largeFloat)
                    {
                        if(!valueOK(model, "File line " + nLinesRead + " Y: " + tokens[colOrder[Y]])) return null;
                    }
                    if (t < -largeFloat || t > largeFloat)
                    {
                        if (!valueOK(model, "File line " + nLinesRead + " T: " + tokens[colOrder[T]])) return null;
                        checked = true;
                    }
                    if(nLine < -100000 || nLine > 100000)
                    {
                        if (!valueOK(model, "File line " + nLinesRead + " ILINE: " + tokens[colOrder[ILINE]])) return null;
                        checked = true;
                   }
                    if(nCrossline < -100000 || nCrossline > 100000)
                    {
                        if (!valueOK(model, "File line " + nLinesRead + " XLINE: " + tokens[colOrder[XLINE]]))  return null;
                         checked = true;
                   }
                 }
                SeismicPoint point = new SeismicPoint(x, y, t, nLine, nCrossline, lattributes);
                points = (SeismicPoint[])StsMath.arrayAddElement(points, point, nPoints, arraySizeInc);
                if(points == null)
                {
                    StsException.systemError("Cannot create seismic surface for: " + name +
                                             " due to StsMath.arrayAddElement error.");
                    return null;
                }
                nPoints++;

                sortedLines.add(nLine);
                sortedCrosslines.add(nCrossline);
                panel.progressBar.setCount(asciiFile.getNBytes());
            }
            panel.appendLine("Processed " + nPoints + " for surface " + name);            
            asciiFile.close();
        }
        catch(Exception e)
       	{
        	new StsMessage(model.win3d, StsMessage.WARNING, "File read error for file " + fullPathname + " line " +  nLinesRead + "\n" + "Data: " + asciiFile.getLine());
            if(asciiFile != null) asciiFile.close();
            return null;
        }

        try
        {
            points = (SeismicPoint[])StsMath.trimArray(points, nPoints);

            lineNumbers = getLineNumbers(sortedLines);
            crosslineNumbers = getLineNumbers(sortedCrosslines);

            minLine = (Float)sortedLines.first();
            maxLine = (Float)sortedLines.last();
            minXline = (Float)sortedCrosslines.first();
            maxXline = (Float)sortedCrosslines.last();

            lineInc = getMinimumLineIncrement(lineNumbers);
            xlineInc = getMinimumLineIncrement(crosslineNumbers);

            nRows = Math.round((maxLine - minLine)/lineInc + 1);
            nCols =  Math.round((maxXline - minXline)/xlineInc + 1);

            // allocate 20 progressBar units for this
            pointArray = new SeismicPoint[nRows][nCols];
            attributes = new float[nAttributes][nRows][nCols];
            attMins = new float[nAttributes];
            attMaxs = new float[nAttributes];
            zValues = new float[nRows][nCols];

            panel.progressBar.setSubInterval(0.6f, 0.8f, nPoints);
            panel.appendLine("Assigning seismic points to array....");
            for(int n = 0; n < nPoints; n++)
            {
                int row = Math.round((points[n].nLine - minLine)/lineInc);
                int col = Math.round((points[n].nCrossline - minXline)/xlineInc);
                pointArray[row][col] = points[n];
                panel.progressBar.incrementCount();
            }

            // compute angle and location of origin
            if(!computeGeometry(model, pointArray, nRows, nCols))
            {
                new StsMessage(model.win3d, StsMessage.WARNING, "Failed to construct grid from seismic pick data.");
                return null;
            }

            // set zValues
            // allocate 10 progressBar units for this
            panel.progressBar.setSubInterval(0.8f, 0.9f, nRows);
            panel.appendLine("Completing surface construction for " + name + "....");
            for(int nn=0; nn<nAttributes; nn++)
            {
            	attMaxs[nn] = -StsParameters.largeFloat;
            	attMins[nn] = StsParameters.largeFloat;
            }
            for(int row = 0; row < nRows; row++)
            {
                for(int col = 0; col < nCols; col++)
                {
                    if(pointArray[row][col] != null)
                    {
                        zValues[row][col] = (float)pointArray[row][col].t;
                        for(int nn=0; nn<nAttributes; nn++)
                        {
                        	attributes[nn][row][col] = (float)pointArray[row][col].attributes[nn]; 
                        	if(attributes[nn][row][col] > attMaxs[nn]) attMaxs[nn] = attributes[nn][row][col];
                        	if(attributes[nn][row][col] < attMins[nn]) attMins[nn] = attributes[nn][row][col];
                        }
                    }
                    else
                    {
                        hasNulls = true;
                        zValues[row][col] = StsParameters.nullValue;
                        for(int nn=0; nn<nAttributes; nn++)
                        	attributes[nn][row][col] = StsParameters.nullValue;                        
                    }
                }
                panel.progressBar.incrementCount();
            }

            if (debug)
            {
                float xSize = (float)((nCols-1)*xInc);
                float ySize = (float)((nRows-1)*yInc);

                System.out.println("xOrigin = " + xOrigin + ", xSize = " + xSize
                    + ", xInc = " + xInc + ", yOrigin = " + yOrigin + ", ySize = " + ySize
                    + ", yInc = " + yInc + ", angle = " + angle  + ", nCols = " + nCols + ", nRows = " + nRows);
            }

            // create the surface
            StsMessageFiles.logMessage("Adding surface to model ...");
    /*
			String stringName = name.toString();

			boolean isDepth;
			if(StsStringUtils.stringContainsString(stringName, "depth"))
				isDepth = true;
			else if(StsStringUtils.stringContainsString(stringName, "time"))
				isDepth = false;
			else
				isDepth = StsMessage.questionValue(model.win3d, "Is " + stringName + " a depth surface (otherwise we will assume time)?");

			byte zDomain = StsParameters.TD_TIME;
			if(isDepth) zDomain = StsParameters.TD_DEPTH;
    */
            panel.progressBar.setSubInterval(0.9f, 1.0f, nRows);
            String stringName = name.toString();
            newSurface = StsSurface.constructSurface(stringName, color, StsSurface.IMPORTED, nCols, nRows, xOrigin, yOrigin, (float)xInc, (float)yInc, 0.0f, 0.0f,
            			(float)angle, zValues, hasNulls, StsParameters.nullValue, zDomain, verticalUnits, horizontalUnits, panel);
            newSurface.setRowColNumRange(minLine, maxLine, minXline, maxXline, lineInc, xlineInc);

            StsSurfaceAttribute[] surfaceAttributes = new StsSurfaceAttribute[attributes.length];
            if(hasAttributes)
            {
            	for(int i=0; i<attributes.length; i++)
            		surfaceAttributes[i] = StsSurfaceAttribute.constructor(newSurface, attributeNames[i], attMins[i], attMaxs[i], attributes[i]);
            }
            if(newSurface == null) return null;

            newSurface.setNativeUnits(verticalUnits, horizontalUnits);
            newSurface.addSurfaceAttributes(surfaceAttributes);

            writeBinaryFile(model, newSurface, StsSurface.seismicGrp);
            
            model.incrementSpectrumColor("Basic");

            // Translate Surface to Model Units
            newSurface.toProjectUnits();

            if(!model.getProject().addToProject(newSurface))
            {
                newSurface.delete();
                new StsMessage(model.win3d, StsMessage.WARNING, "Failed to load seismic surface: " + stringName);
                return null;
            }
            model.getProject().setCursorDisplayXYAndGridCheckbox(false);
        }
       	catch (Exception e)
        {
            StsException.outputException("StsImportSeismicSurfaces.createSurface() failed.", e, StsException.WARNING);
            return null;
        }

        System.out.println("Loaded file: " + fullPathname);
        StsMessageFiles.logMessage("Loaded map: " + name.toString());
        filePosition = nFileLinesRead;
        return newSurface;
    }

    static public void initializeColOrder()
    {
        colOrder = new int[] { X, Y, T, ILINE, XLINE };
    }

    static private boolean valueOK(StsModel model, String string)
    {
        return StsYesNoDialog.questionValue(model.win3d, "Data format seems wrong. Questionable value: " + string + "\nDo you wish to continue?");
    }

    static class SeismicPoint
    {
        public double x, y, t;
        public float nLine, nCrossline;
        public double[] attributes = null;

        SeismicPoint(double x, double y, double t, int nLine, int nCrossline)
        {
            this.x = x;
            this.y = y;
            this.t = t;
            this.nLine = nLine;
            this.nCrossline = nCrossline;
        }
        
        SeismicPoint(double x, double y, double t, float nLine, float nCrossline, double[] attributes)
        {
            this.x = x;
            this.y = y;
            this.t = t;
            this.nLine = nLine;
            this.nCrossline = nCrossline;
            this.attributes = attributes;
        } 
        
        public int getNumberAttributes() { return attributes.length; }
    }

    static float[] getLineNumbers(TreeSet sortedLines)
    {
        int nLines = sortedLines.size();
        float[] lineNumbers = new float[nLines];
        Iterator iter = sortedLines.iterator();
        int n = 0;
        while(iter.hasNext())
            lineNumbers[n++] = ((Float)iter.next());
        return lineNumbers;
    }

    static float getMinimumLineIncrement(float[] lineNumbers)
    {
        // get min increment from ordered list of lineNumbers
        float minInc = StsParameters.largeFloat;
        int nLines = lineNumbers.length;
        for(int n = 1; n < nLines; n++)
        {
            minInc = Math.min(minInc, Math.abs( (lineNumbers[n] - lineNumbers[n-1]) ) );
            if(minInc == 1) return 1;
        }

        // check that increment is divisible into each line spacing
        boolean incOK = false;
        while(!incOK)
        {
            incOK = true;
            for(int n = 1; n < nLines; n++)
            {
                float dif = Math.abs( (lineNumbers[n] - lineNumbers[n-1]) );
                if(!StsMath.isIntegralRatio(dif, minInc, 1.0e-3f))
                {
                    new StsMessage(StsModel.getCurrentModel().win3d,  StsMessage.ERROR, "Adjacent line numbers " + lineNumbers[n-1] + " and " + lineNumbers[n] +
                                    " are not divisible by minimum line number spacing " + minInc + ".\nWill use min spacing.");
                    return minInc;
                }
            }
        }
        return minInc;
    }

    static private boolean computeGeometry(StsModel model, SeismicPoint[][] pointArray, int nRows, int nCols)
    {
        double angleX = 0.0, angleY = 0.0;
        int firstColFirstRow, firstRowFirstCol;
        int firstColLastRow, firstRowLastCol;
        int firstNonNullCol;
        double dx, dy;
        SeismicPoint point0, point1;
        double distance;
        double cosA, sinA;
        double angDif;
        int row, col;

        try
        {
            /* Search x dimension for two valid points */

            firstRowFirstCol = -1;
            firstRowLastCol = -1;
            for(row = 0; row < nRows; row++)
            {
                for(col = 0; col < nCols; col++)
                {
                    if(pointArray[row][col] != null)
                    {
                        firstRowFirstCol = col;
                        break;
                    }
                }

                for(col = nCols-1; col >= 0; col--)
                {
                    if(pointArray[row][col] != null)
                    {
                        firstRowLastCol = col;
                        break;
                    }
                }

            /* Compute rotation angle and increment */
                if((firstRowFirstCol >= 0 && firstRowLastCol >= 0) && ((firstRowLastCol - firstRowFirstCol) > (nCols/4)))
                {
                    if (firstRowFirstCol == firstRowLastCol)
                    {
                        firstRowFirstCol = -1;
                        firstRowLastCol = -1;
                    }
                    else
                    {
                        point0 = pointArray[row][firstRowFirstCol];
                        point1 = pointArray[row][firstRowLastCol];

                        dx = point1.x - point0.x;
                        dy = point1.y - point0.y;

                    /* angle_rad_x is ccw angle between row direction and positive x axis	*/

                        angleX = StsMath.atan2(dx, dy);
                        distance = Math.sqrt(dx*dx + dy*dy);
                        xInc = distance/(firstRowLastCol - firstRowFirstCol);

                        break;
                    }
                }
            }

            /* Search in y dimension for two valid points */

            firstColFirstRow = -1;
            firstColLastRow = -1;
            int firstCol = -1;
            for(col = 0; col < nCols; col++)
            {
                for(row = 0; row < nRows; row++)
                 {
                    if(pointArray[row][col] != null)
                    {
                        firstColFirstRow = row;
                        firstCol = col;
                        break;
                    }
                }

                for(row = nRows-1; row >= 0; row--)
                {
                    if(pointArray[row][col] != null)
                    {
                        firstColLastRow = row;
                        break;
                    }
                }

                if(firstColFirstRow >= 0 && firstColLastRow >= 0)
                {
                    if (firstColFirstRow == firstColLastRow)
                    {
                        firstColFirstRow = -1;
                        firstColLastRow = -1;
                    }
                    else
                    {
                        point0 = pointArray[firstColFirstRow][col];
                        point1 = pointArray[firstColLastRow][col];

                        dx = point1.x - point0.x;
                        dy = point1.y - point0.y;

                    /* angle_radY is ccw angle between col direction and positive x axis	*/

                        angleY = StsMath.atan2(dx, dy);
                        distance = Math.sqrt(dx*dx + dy*dy);
                        yInc = distance/(firstColLastRow - firstColFirstRow);

                        angDif = angleY - angleX;
                        if(angDif < 0.0f) angDif += 360.0f;

                        if(Math.abs(angDif  - 90.0f) < 1.0f)
                            isCrosslineCCW = true;
                        else if(Math.abs(angDif - 270.0f) < 1.0f)
                             isCrosslineCCW = false;
                        else
                        {
                            boolean doIt = StsYesNoDialog.questionValue(model.win3d, "Horizon line and trace angles" +
                                    " nonorthogonal by more than 1 degree.\n Might need to adjust column ordering. Load anyways?");
                            if(!doIt) return false;

                            if(angDif <= 180.0f)
                                isCrosslineCCW = true;
                            else
                                isCrosslineCCW = false;
                        }
                        break;
                    }
                }
            }

            if (xInc == 0.0 || yInc == 0.0)
            {
                StsMessage.printMessage("Horizon pick increment can not be computed from pick file.");
                return false;
            }

            angle = angleX;
            cosA = StsMath.cosd(angle);
            sinA = StsMath.sind(angle);

            // if crossline direction is not 90 degrees CCW from inline direction, reverse crosslines (columns)
            if(!isCrosslineCCW)
            {
                for(col = 0; col < nCols; col++)
                {
                    int r0 = 0;
                    int r1 = nRows-1;
                    while(r0 < r1)
                    {
                        SeismicPoint p0 = pointArray[r0][col];
                        SeismicPoint p1 = pointArray[r1][col];
                        pointArray[r0++][col] = p1;
                        pointArray[r1--][col] = p0;
                    }
                }
                firstColFirstRow = nRows-1 - firstColLastRow;
                float minLineTemp = minLine;
                minLine = maxLine;
                maxLine = minLineTemp;
                lineInc = -lineInc;
            }

            // adjust the origin so it is at row=0, col=0
            // firstRowFirstCol is the first col with a non-null in it at row firstColFirstRow
            xOrigin = modelXOrigin + pointArray[firstColFirstRow][firstCol].x - (firstCol *cosA*xInc - firstColFirstRow *sinA*yInc);
            yOrigin = modelYOrigin + pointArray[firstColFirstRow][firstCol].y - (firstCol *sinA*xInc + firstColFirstRow *cosA*yInc);


            return true;
        }
        catch(Exception e)
        {
            StsException.outputException("StsImportSeismicSurfaces.computeGeometry() failed.",
                e, StsException.WARNING);
            return false;
        }
    }
}