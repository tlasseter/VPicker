
//Title:        S2S: Seismic-to-simulation
//Version:
//Copyright:    Copyright (c) 2001
//Author:       Tom Lasseter
//Company:      4D Systems LLC
//Description:  Web-enabled Integrated Interpretation System

package com.Sts.Actions.Build;

import com.Sts.Actions.*;
import com.Sts.DBTypes.*;
import com.Sts.MVC.*;
import com.Sts.RFUI.*;

public class StsBuildSurfacesFromMarkers extends StsAction implements Runnable
{
//    private SDSFile sdsFile = null;
//    private SDSFile rfHdfSdsFile = null;
    private int nSurfacesBuilt = 0;
    private boolean success = false;


 	public StsBuildSurfacesFromMarkers(StsActionManager actionManager)
    {
        super(actionManager);
    }

    public void run()
    {
        StsObject[] markers = null;
        StsModelSurface xySurface = null;
        try
        {
            statusArea.setTitle("Build surfaces:");
            logMessage("Choose markers to derive surfaces from.");

            // get markers
            StsSelectStsObjects selector = StsSelectStsObjects.constructor(model,
                    StsMarker.class, model.getName(), "Select markers:", false, false);
			if(selector == null) return;
            markers = selector.selectObjects();
            if (markers==null)
            {
                statusArea.textOnly();
                logMessage("No markers selected.");
    	        actionManager.endCurrentAction();
                return;
            }

            // get a surface
            StsObject[] surfaces = model.getObjectListOfType(StsModelSurface.class, StsModelSurface.IMPORTED);
            if (surfaces != null)  // pick an existing grid
            {
                selector = StsSelectStsObjects.constructor(model, surfaces, model.getName(), "Select a surface grid:", true);
				if(selector == null) return;
                StsModelSurface surface = (StsModelSurface)selector.selectObject();
                if(surface == null)
                {
                    statusArea.textOnly();
                    logMessage("No surface grid selected.");
                    return;
                }
				actionManager.endCurrentAction();
           }
        }
        catch(Exception e) { return; }
/*
        // look in surfaces.set files for existing surfaces
        HDFSurfaceFile hdfSurfaceFile = null;
        try
        {
            File f = new File(model.getProject().getDataFullDirString()
                    + File.separator + StsImportHDFSurfaces.HDF_SURFACE_SET);
 	        sdsFile = new SDSFile(f, SDSFile.ACCESS_READ);
            hdfSurfaceFile = new HDFSurfaceFile(sdsFile);
        }
        catch (Exception e)  { if (sdsFile!=null) sdsFile.end(); }
        HDFSurfaceFile rfHdfSurfaceFile = null;
        try
        {
            File f = new File(model.getProject().getBinaryFullDirString()
                    + File.separator + StsImportHDFSurfaces.HDF_SURFACE_SET);
       	    // open the HDF file
            rfHdfSdsFile = new SDSFile(f, SDSFile.ACCESS_READ);
            rfHdfSurfaceFile = new HDFSurfaceFile(rfHdfSdsFile);
        }
        catch (Exception e)  { if (rfHdfSdsFile!=null) sdsFile.end(); }

        try
        {
            // build the surfaces
            status.setMaximum(markers.length);
            StsStatusArea.staticSetProgress(0);
            status.addProgress();
            status.sleep(10);
            for (int i=0; i<markers.length; i++)
            {
                StsMarker marker = (StsMarker)markers[i];
                String name = marker.getName();
                if (model.getSurface(name) != null)
                {
                    logMessage("Already have surface: " + name +
                            ".  Skipping to next marker.");
                    StsStatusArea.staticSetProgress(i+1); // update progress
                    continue;
                }
                if ((hdfSurfaceFile != null &&  hdfSurfaceFile.getSurface(name) != null) ||
                    (rfHdfSurfaceFile != null && rfHdfSurfaceFile.getSurface(name) != null))
                {
                    logMessage("The surface: " + name + " is already saved in "
                                + "the Sts data store.  Skipping to next marker.");
                    StsStatusArea.staticSetProgress(i+1); // update progress
                    continue;
                }

                StsSurface newSurface = new StsSurface(model, name,
                        marker.getStsColor(), StsSurface.IMPORTED);
                model.refreshModelState();
                if (xyGrid == null) // build hardwired grid from model x-y range
                {
                    int nX = 100;
                    int nY = 100;
                    StsProject p = model.getProject();
                    float xMin = p.getXMin();
                    float xMax = p.getXMax();
                    float yMin = p.getYMin();
                    float yMax = p.getYMax();
                    float xInc = (xMax - xMin) / nX;
                    float yInc = (yMax - yMin) / nY;
                    xyGrid = new StsGrid(model, newSurface, nX, nY,
                            xMin, yMin, xInc, yInc, null,
                            true, StsParameters.nullValue);
                }
                StsGrid newGrid = xyGrid.cloneParameters(newSurface);
                newSurface.setGrid(newGrid);
                newGrid.setPrefixString("Surface: " + name);
                if (!newGrid.setGridPoints(marker, xyGridCopied ? xyGrid : null,
                        StsParameters.largeFloat, 4.0f))
                {
                    logMessage("Grid creation for surface:  " +
                            name + " failed.");
                    model.delete(newSurface);
                }
                newGrid.persistPoints();
                newGrid.constructGrid();
                //newGrid.saveValuesToDataStore();
                newSurface.setMarker(marker);
                marker.setParentSurface(newSurface);
                nSurfacesBuilt++;
                logMessage("Built new surface: " + name);
                StsStatusArea.staticSetProgress(i+1); // update progress
            }

            success = true;
            closeSdsFiles();
        	glPanel.getActionManager().endCurrentAction();
        }
        catch(Exception e)
        {
            StsStatusArea.getStatusArea().removeProgress();
            StsException.outputException("StsBuildSurfacesFromMarkers failed.",
                    e, StsException.WARNING);
            closeSdsFiles();
        }
*/
    }

    public boolean end()
    {
        statusArea.textOnly();
        if (success)
        {
            logMessage("Built " + nSurfacesBuilt +
                    " new surfaces from markers.");
        }
        else
        {
            logMessage("Unable to build surfaces from markers.");
        }
        return success;
    }
/*
    private void closeSdsFiles()
    {
        if (sdsFile != null) sdsFile.end();
        if (rfHdfSdsFile != null) rfHdfSdsFile.end();
    }
*/
}
