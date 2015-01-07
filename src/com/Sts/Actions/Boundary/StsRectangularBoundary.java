

//Title:        S2S: Seismic-to-simulation
//Version:
//Copyright:    Copyright (c) 2001
//Author:       Tom Lasseter
//Company:      4D Systems LLC
//Description:  Web-enabled Integrated Interpretation System

package com.Sts.Actions.Boundary;

import com.Sts.DBTypes.*;
import com.Sts.MVC.*;
import com.Sts.MVC.View3d.*;
import com.Sts.Types.*;
import com.Sts.Utilities.*;


/** build a model boundary with row-col alignment */
public class StsRectangularBoundary extends StsPolygonBoundary
{
    int rowMin = 0, rowMax = 0, colMin = 0, colMax = 0;

 	public StsRectangularBoundary(StsActionManager actionManager)
    {
        super(actionManager);
    }

	public boolean start()
    {
        try
        {
            if (!actionOK())  return false;
            statusArea.setTitle("Build row-col boundary: ");
            gridPoints = new StsGridPoint[4];
            logMessage("Select one rectangle corner by pressing mouse; drag to opposite corner.");
            return true;
        }
        catch(Exception e){ return false; }
    }

    /** mouse action for 3d window */
   	public boolean performMouseAction(StsMouse mouse, StsGLPanel glPanel)
	{
        StsGridSectionPoint gridSectionPoint;
        StsPoint pickPoint;

    	int leftButtonState = mouse.getButtonStateCheckClear(StsMouse.LEFT);

        if(leftButtonState == StsMouse.PRESSED || leftButtonState == StsMouse.DRAGGED)
        {
            StsGridPoint gridPoint = surface.getSurfacePosition(mouse, true, (StsGLPanel3d)glPanel);
            if(gridPoint == null) return true;
            gridPoint.adjustToNearestRowCol();

            if(gridPoints[0] == null)
                gridPoints[0] = gridPoint;
            else
            {
                gridPoints[1] = gridPoint;
                computeRowColRange();
            }
        }
        else if (StsMouse.isButtonStateReleasedOrClicked(leftButtonState))  // add a point or finish picking
       	{
            isLoop = true;
            constructCornerGridPoints();
            glPanel.actionManager.endCurrentAction();
        }

        if(gridPoints[0] == null || gridPoints[1] == null) return true;

        if(boundaryEdge == null)
        {
            boundaryEdge = new StsSectionEdge(StsParameters.BOUNDARY, surface, false);
            model.addDisplayableInstance(boundaryEdge);
        }
        StsPoint[] points = surface.getRectanglePoints(rowMin, rowMax, colMin, colMax);
        boundaryEdge.setPoints(points);
        model.win3d.win3dDisplay();
        return true;
	}

    private void computeRowColRange()
    {
        if(gridPoints[0] == null || gridPoints[1] == null) return;

        int row0 = gridPoints[0].row;
        int col0 = gridPoints[0].col;
        int row1 = gridPoints[1].row;
        int col1 = gridPoints[1].col;

        if(row0 <= row1)
        {
            rowMin = row0;
            rowMax = row1;
        }
        else
        {
            rowMin = row1;
            rowMax = row0;
        }

        if(col0 <= col1)
        {
            colMin = col0;
            colMax = col1;
        }
        else
        {
            colMin = col1;
            colMax = col0;
        }
    }

    private void constructCornerGridPoints()
    {
        StsGridDefinition gridDef = model.getGridDefinition();
        gridPoints[0] = new StsGridPoint(rowMin, colMin, gridDef);
        gridPoints[1] = new StsGridPoint(rowMax, colMin, gridDef);
        gridPoints[2] = new StsGridPoint(rowMax, colMax, gridDef);
        gridPoints[3] = new StsGridPoint(rowMin, colMax, gridDef);
    }
}


