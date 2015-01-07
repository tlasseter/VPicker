package com.Sts.Actions.Wizards.FaultSticks;

import com.Sts.Actions.Wizards.*;
import com.Sts.DBTypes.*;
import com.Sts.IO.*;
import com.Sts.MVC.*;
import com.Sts.Types.*;
import com.Sts.UI.Progress.*;
import com.Sts.UI.*;
import com.Sts.Utilities.*;

import java.awt.*;
import java.util.*;

/**
 * <p>Title: </p>
 *
 * <p>Description: </p>
 *
 * <p>Copyright: Copyright (c) 2002</p>
 *
 * <p>Company: </p>
 *
 * @author not attributable
 * @version 1.0
 */
public class StsFaultSticksWizard extends StsWizard
{
    StsAbstractFileSet fileSet;
    StsFaultStickSet faultStickSet;
    String faultStickSetName = null;

    byte vUnits = StsParameters.DIST_NONE;
    byte hUnits = StsParameters.DIST_NONE;
    byte tUnits = StsParameters.TIME_NONE;
    float datumShift = 0.0f;
    byte zDomain = StsProject.TD_TIME;

    static final byte COOR = 0;
    static final byte GRID = 1;
    byte coorOrGrid = COOR;

    public StsFaultSticksSelect faultStickSelect = new StsFaultSticksSelect(this);
    public StsFaultSticksAttributes faultStickAttributes = new StsFaultSticksAttributes(this);
    public StsFaultSticksLoad faultStickLoad = new StsFaultSticksLoad(this);
    private StsWizardStep[] mySteps =
    {
    		faultStickSelect, faultStickAttributes, faultStickLoad
    };

    static public final double rmsNull = 999.0f;

    public StsFaultSticksWizard(StsActionManager actionManager)
    {
        super(actionManager,500,500);
        addSteps(mySteps);
    }

    public boolean start()
    {
        System.runFinalization();
        System.gc();
        dialog.setTitle("Process & Load Fault Stick Sets");
        disableFinish();
        return super.start();
    }

    public boolean end()
    {
        if (success) model.setActionStatus(getClass().getName(), StsModel.STARTED);
        return super.end();
    }

    public void previous()
    {
        gotoPreviousStep();
    }

    public void next()
    {
        if(currentStep == faultStickAttributes)
            faultStickLoad.constructPanel();
        gotoNextStep();
    }

    public void finish()
    {
        super.finish();
    }

    public StsAbstractFile[] getSelectedFiles()
    {
        return faultStickSelect.getSelectedFiles();
    }

    public boolean createFaultStickSet(StsProgressPanel panel)
    {
       if( model == null ) return false;

       try
       {
           StsAbstractFile[] selectedFiles = getSelectedFiles();
           model.disableDisplay();

           disablePrevious();

           // turn on the wait cursor
           StsCursor cursor = new StsCursor(model.win3d, Cursor.WAIT_CURSOR);
           int nSelected = selectedFiles.length;
           int nLoaded = 0;

            String timeOrDepth = faultStickSelect.panel.verticalUnitsPanel.timeDepthString;
            boolean isTrueDepth = timeOrDepth.equals(StsParameters.TD_DEPTH_STRING);
            if(isTrueDepth) zDomain = StsProject.TD_DEPTH;

            String tvdOrSsString = faultStickSelect.panel.verticalUnitsPanel.tvdOrSsString;
            boolean isTvd = tvdOrSsString.equals(StsParameters.TVD_STRING);
            vUnits = faultStickSelect.panel.unitsPanel.vUnits;
            hUnits = faultStickSelect.panel.unitsPanel.hUnits;
            tUnits = faultStickSelect.panel.unitsPanel.tUnits;
            datumShift = faultStickSelect.panel.verticalUnitsPanel.datumShift;

           StsProject project = model.getProject();
           panel.initialize(nSelected);
           for(int n = 0; n < nSelected; n++)
           {
               StsAbstractFile file = selectedFiles[n];
               if (selectedFiles[n] == null) continue;
               setFileType(file);
               panel.appendLine("Scanning file named: " + selectedFiles[n].getFilename());
               computeRanges(file, panel);
            	   
               panel.appendLine("Loading file named: " + selectedFiles[n].getFilename());
               if(!loadFaultSticks(file, panel))
               {
            	   panel.appendLine("Failed to load file named: " + selectedFiles[n].getFilename());
            	   panel.setLevel(StsProgressBar.WARNING);
               }
               else
            	   nLoaded++;

               panel.setValue(n+1);
               panel.setDescription("Fault stick set " + nLoaded + " of " + nSelected + " loaded.");
           }
           panel.setDescription("Loaded " + nLoaded + " fault stick files of " + nSelected + " selected.");
           
           // turn off the wait cursor
           if( cursor!=null ) cursor.restoreCursor();

           disableCancel();
           enableFinish();

           if(nLoaded > 0)
           {
               StsSeismicVelocityModel velocityModel = project.getSeismicVelocityModel();
               if(velocityModel != null) checkSetTimeOrDepth(velocityModel);
               project.adjustBoundingBoxes(true, false);
        	   project.checkAddUnrotatedClass(StsFaultStickSet.class);
        	   project.rangeChanged();
        	   model.win3d.cursor3d.initialize();
        	   model.win3d.cursor3dPanel.setSliderValues();
        	   model.getGlPanel3d().setDefaultView();          
        	   model.enableDisplay();
        	   model.win3dDisplay(); // display the surfaces
           }

           return true;
       }
       catch (Exception e)
       {
           panel.appendLine("StsFaultStickWizard.createFaultStickSets() failed." + e.getMessage());
           panel.setDescriptionAndLevel("Exception thrown.", StsProgressBar.ERROR);
           StsException.outputException("StsFaultStickWizard.createFaultStickSets() failed.", e, StsException.WARNING);
           return false;
       }
    }

    private void setFileType(StsAbstractFile file)
    {
        if(file.format != StsAbstractFile.asciiFormat) return;
        file.fileType = file.group;
    }

    private void checkSetTimeOrDepth(StsSeismicVelocityModel velocityModel)
    {
        StsObjectRefList faultSticksList = faultStickSet.getFaultStickList();
        int nFaultSticks = faultSticksList.getSize();
        for (int n = 0; n < nFaultSticks; n++)
        {
            StsFaultLine faultStick = (StsFaultLine)faultSticksList.getElement(n);
            faultStick.adjustFromVelocityModel(velocityModel);
        }
    }

    private boolean computeRanges(StsAbstractFile file, StsProgressPanel panel)
    {
        String group = file.group;
        if(group.equals(StsFaultStickSet.lmkGrp))
            return computeRangesLmk(file, panel);
        else if(group.equals(StsFaultStickSet.rmsGrp))
            return computeRangesRms(file, panel);
        else
            return false;
    }

    private boolean computeRangesLmk(StsAbstractFile file, StsProgressPanel panel)
    {
    	int nLoaded = 0;
    	float x, y, zt;
    	float zMin = StsParameters.largeFloat, zMax = -StsParameters.largeFloat;
    	float xMin = StsParameters.largeFloat, xMax = -StsParameters.largeFloat;
    	float yMin = StsParameters.largeFloat, yMax = -StsParameters.largeFloat;
    	
        StringTokenizer st;
        String strToken;
        
        try 
        {
            StsAsciiFile asciiFile = new StsAsciiFile(file);
            if(!asciiFile.openReadWithErrorMessage()) return false;

            // Read and parse the selected file.
            String aline = asciiFile.readLine();
            while(aline != null)
            {
                st = new StringTokenizer(aline, ", ");
                x = Float.parseFloat(st.nextToken());
                if(x < xMin) xMin = x;
                if(x > xMax) xMax = x;                
                y = Float.parseFloat(st.nextToken());
                if(y < yMin) yMin = y;
                if(y > yMax) yMax = y;                
                
                zt = Float.parseFloat(st.nextToken());
                if(zt < zMin) zMin = zt;
                if(zt > zMax) zMax = zt;
                aline = asciiFile.readLine();
            } 
            asciiFile.close();
            panel.appendLine("   XRange --> " + xMin + " ," + xMax);
            panel.appendLine("   YRange --> " + yMin + " ," + yMax);
            panel.appendLine("   ZTRange --> " + zMin + " ," + zMax);
            
            if(!model.getProject().isInProjectBounds(xMax, yMax) && !model.getProject().isInProjectBounds(xMin, yMin))
            {
            	panel.appendLine("Positioning data appears to be inline and crosslines.");
            	coorOrGrid = GRID;            
            }
            else
            {
            	panel.appendLine("Positioning data appears to be X & Y.");
            	coorOrGrid = COOR;            
            }           	
        } 
        catch (Exception e) 
        {
        	StsException.outputException("Unable to scan selected file", e, StsException.FATAL);
        	e.printStackTrace();
        	return false;
        }
        if(nLoaded == 0) 
        	return false;
        else
        {
        	panel.appendLine("Scanned " + nLoaded + " fault sticks from file: " + file.getFilename());
        	return true;
        }
    }
    private boolean computeRangesRms(StsAbstractFile file, StsProgressPanel panel)
    {
    	int nLoaded = 0;
    	float x, y, zt;
    	float zMin = StsParameters.largeFloat, zMax = -StsParameters.largeFloat;
    	float xMin = StsParameters.largeFloat, xMax = -StsParameters.largeFloat;
    	float yMin = StsParameters.largeFloat, yMax = -StsParameters.largeFloat;

        StringTokenizer st;
        String strToken;

        try
        {
            StsAsciiFile asciiFile = new StsAsciiFile(file);
            if(!asciiFile.openReadWithErrorMessage()) return false;

            // Read and parse the selected file.
            String aline = asciiFile.readLine();
            while(aline != null)
            {
                st = new StringTokenizer(aline, ", ");
                x = Float.parseFloat(st.nextToken());
                if(x < xMin) xMin = x;
                if(x > xMax) xMax = x;
                y = Float.parseFloat(st.nextToken());
                if(y < yMin) yMin = y;
                if(y > yMax) yMax = y;

                zt = Float.parseFloat(st.nextToken());
                if(zt < zMin) zMin = zt;
                if(zt > zMax) zMax = zt;
                aline = asciiFile.readLine();
            }
            asciiFile.close();
            panel.appendLine("   XRange --> " + xMin + " ," + xMax);
            panel.appendLine("   YRange --> " + yMin + " ," + yMax);
            panel.appendLine("   ZTRange --> " + zMin + " ," + zMax);
        }
        catch (Exception e)
        {
        	StsException.outputException("Unable to scan selected file", e, StsException.FATAL);
        	e.printStackTrace();
        	return false;
        }
        if(nLoaded == 0)
        	return false;
        else
        {
        	panel.appendLine("Scanned " + nLoaded + " fault sticks from file: " + file.getFilename());
        	return true;
        }
    }

    private boolean loadFaultSticks(StsAbstractFile file, StsProgressPanel panel)
    {
        String fileType = file.fileType;
        if(fileType.equals(StsFaultStickSet.lmkGrp))
            return loadFaultSticksLmk(file, panel);
        else if(fileType.equals(StsFaultStickSet.rmsGrp))
            return loadFaultSticksRms(file, panel);
        else
            return false;
    }

    private boolean loadFaultSticksLmk(StsAbstractFile file, StsProgressPanel panel)
    {
    	int nLoaded = 0;
    	float c1, c2, zt;
    	float[] xyzmt;
    	int pos;
    	StsPoint point;
    	StsFaultLine line = null;
        StsObjectRefList lineVertices = null;
        StringTokenizer st;
        StsFaultStickSet faultStickSet = StsFaultStickSet.constructor(faultStickSetName);
        StsProject project = model.getProject();
        StsColor faultColor = faultStickAttributes.panel.getStsColor();
        try 
        {
            StsAsciiFile asciiFile = new StsAsciiFile(file);
            if(!asciiFile.openReadWithErrorMessage()) return false;

            // Read and parse the selected file.
            String aline = asciiFile.readLine();
            while(aline != null)
            {
                st = new StringTokenizer(aline, ", ");

                c1 = Float.parseFloat(st.nextToken());
                c2 = Float.parseFloat(st.nextToken());
                zt = Float.parseFloat(st.nextToken());
                
                xyzmt = new float[] { c1, c2, zt, 0.0f, zt };
                if(coorOrGrid == GRID)
                {
                	float rotatedX = project.getRotatedBoundingBox().getXFromColNum(c2);
                	float rotatedY = project.getRotatedBoundingBox().getYFromRowNum(c1);
                    float[] xy = project.getUnrotatedRelativeXYFromRotatedXY(rotatedX, rotatedY);
                    xyzmt[0] = xy[0];
                    xyzmt[1] = xy[1];
                }
                point = new StsPoint(xyzmt);

                int id = Integer.parseInt(st.nextToken());  // not used
                int id2 = Integer.parseInt(st.nextToken()); // not used
                pos = Integer.parseInt(st.nextToken());
                switch(pos)
                {
                	case 1: // First point in new line
                		line = StsFaultLine.buildImportedFault(faultStickSet);
                        line.setXOrigin(project.getXOrigin());
                        line.setYOrigin(project.getYOrigin());
                        lineVertices = line.getLineVertices();
                        lineVertices.add(new StsSurfaceVertex(point, line));
                		break;
                	case 2: // Middle point
                		lineVertices.add(new StsSurfaceVertex(point, line));
                		break;
                	case 3: // Last point in line
                		lineVertices.add(new StsSurfaceVertex(point, line));
                        line.initialize(model);
                        line.addToSet(faultStickSet);
                        line.computeMDepths();
                        line.addToProject(zDomain);
                		nLoaded++;
                        break;
                }
                aline = asciiFile.readLine();
            }            
        } 
        catch (Exception e) 
        {
        	StsException.outputException("Unable to load selected file", e, StsException.FATAL);
        	e.printStackTrace();
        	return false;
        }
        if(nLoaded == 0) 
        	return false;
        else
        {
        	panel.appendLine("Loaded " + nLoaded + " fault sticks from file: " + file.getFilename());
        	faultStickSet.setStsColor(faultColor);
        	return true;
        }	
    }

    private boolean loadFaultSticksRms(StsAbstractFile file, StsProgressPanel panel)
    {
    	int nLoaded = 0;
    	double x, y;
        float zt;
    	int pos;
    	StsPoint point = null;
    	StsFaultLine line = null;
        StsObjectRefList lineVertices = null;
        StringTokenizer st;
        StsFaultStickSet faultStickSet = StsFaultStickSet.constructor(faultStickSetName);
        StsProject project = model.getProject();
        StsColor faultColor = faultStickAttributes.panel.getStsColor();
        pos = 1;
        try
        {
            StsAsciiFile asciiFile = new StsAsciiFile(file);
            if(!asciiFile.openReadWithErrorMessage()) return false;

            // Read and parse the selected file.
            String aline = asciiFile.readLine();
            while(aline != null)
            {
                st = new StringTokenizer(aline, ", ");

                x = Double.parseDouble(st.nextToken());
                y = Float.parseFloat(st.nextToken());
                zt = Float.parseFloat(st.nextToken());
                float[] xy = project.getRelativeXY(x, y); // returns xy in relative unrotated coor system with respect to origin
                if(isLastPoint(x, y, zt))
                    pos = 3;
                else
                    point = new StsPoint(new float[] { xy[0], xy[1], zt, 0.0f, zt });

                switch(pos)
                {
                	case 1: // First point in new line
                		line = StsFaultLine.buildImportedFault();
                        line.setXOrigin(project.getXOrigin());
                        line.setYOrigin(project.getYOrigin());
                        lineVertices = line.getLineVertices();
                        lineVertices.add(new StsSurfaceVertex(point, line));
                        pos = 2;
                        break;
                	case 2: // Middle point
                        lineVertices.add(new StsSurfaceVertex(point, line));
                		break;
                	case 3: // Last point in line
                        line.checkSortLineVertices();
                        line.removeClosePoints();
                        // line.checkIsVertical();
                        line.computeMDepths();
                        line.initialize(model);
                        line.addToSet(faultStickSet);
                        line.addToProject(zDomain);
                		nLoaded++;
                        pos = 1;
                        break;
                }
                aline = asciiFile.readLine();
            }
        }
        catch (Exception e)
        {
        	StsException.outputException("Unable to load selected file", e, StsException.FATAL);
        	e.printStackTrace();
        	return false;
        }
        if(nLoaded == 0)
        	return false;
        else
        {
        	panel.appendLine("Loaded " + nLoaded + " fault sticks from file: " + file.getFilename());
        	faultStickSet.setStsColor(faultColor);
        	return true;
        }
    }

    private boolean isLastPoint(double x, double y, float z)
    {
        return x == rmsNull && y == rmsNull;
    }

    public boolean checkPrerequisites() 
    { 
    	Object[] objects = (Object[])model.getCastObjectList(StsSurface.class);
    	if(objects.length > 0) return true;
    	
    	objects = (Object[])model.getCastObjectList(StsSeismicVolume.class);
    	if(objects.length > 0) return true; 
    	
    	new StsMessage(frame,StsMessage.WARNING, "Must load a surface or seismic volume prior to loading fault sticks.");
    	return false; 
    }
    
    public String getFaultStickSetName() { return faultStickAttributes.panel.getName(); }
    public void setFaultStickSetName(String name)
    {
        faultStickSetName = name;
    }
}
