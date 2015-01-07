
//Title:        S2S: Seismic-to-simulation
//Version:
//Copyright:    Copyright (c) 2001
//Author:       Tom Lasseter
//Company:      4D Systems LLC
//Description:  Web-enabled Integrated Interpretation System

package com.Sts.Actions.Boundary;

import com.Sts.Actions.*;
import com.Sts.DBTypes.*;
import com.Sts.MVC.*;
import com.Sts.MVC.View3d.*;
import com.Sts.Types.*;
import com.Sts.UI.Toolbars.*;
import com.Sts.Utilities.*;

/** build a model boundary with row-col alignment */
public class StsPolygonBoundary extends StsAction
{
    protected StsModelSurface surface = null;
    protected boolean visible = true;
    private boolean success = false;

    protected StsGridPoint[] gridPoints;
    protected boolean isLoop = false;
    protected StsEdge boundaryEdge;

    // convenience copies of flags
/*
    static protected final int VERTICAL_ROW_MINUS = StsParameters.VERTICAL_ROW_MINUS;
    static protected final int VERTICAL_COL_MINUS = StsParameters.VERTICAL_COL_MINUS;
    static protected final int VERTICAL_ROW_PLUS = StsParameters.VERTICAL_ROW_PLUS;
    static protected final int VERTICAL_COL_PLUS = StsParameters.VERTICAL_COL_PLUS;
*/
    static protected final int NONE = StsParameters.NONE;

    static private float gridNearSq = 1.0f; // min distance between points; otherwise the same

 	public StsPolygonBoundary(StsActionManager actionManager)
    {
        super(actionManager);
    }

	public boolean actionOK()
    {
        StsModelSurfaceClass surfaceClass = (StsModelSurfaceClass)model.getStsClass(StsModelSurface.class);
        surface = (StsModelSurface)surfaceClass.getTopVisibleSurface();
        if(surface == null)
        {
            statusArea.textOnly();
            logMessage("No horizons available!");
            return false;
        }
        return true;
    }

	public boolean start()
    {
		if( !actionOK() )  return false;

        try
        {
            statusArea.setTitle("Build arbitrary boundary: ");

            logMessage("Pick points in clockwise sequence; finish by selecting first again" + " or with Done button.");
            addEndButton();
            gridPoints = new StsGridPoint[0];

            model.win3d.win3dDisplay();
            return true;
        }
        catch(Exception e) { return false;}
    }

    /** mouse action for 3d window */
   	public boolean performMouseAction(StsMouse mouse, StsGLPanel glPanel)
	{
    	int leftButtonState = mouse.getButtonStateCheckClear(StsMouse.LEFT);

        if (leftButtonState == StsMouse.PRESSED || leftButtonState == StsMouse.DRAGGED)
        {
            surface.getSurfacePosition(mouse, true, (StsGLPanel3d)glPanel);
        }
       	else if (StsMouse.isButtonStateReleasedOrClicked(leftButtonState))  // add a point or finish picking
       	{
            StsGridPoint gridPoint = surface.getSurfacePosition(mouse, true, (StsGLPanel3d)glPanel);
            if(gridPoint == null) return true;

            int nPicks = StsMath.arrayLength(gridPoints);
            if(nPicks == 0)
            {
                boundaryEdge = new StsSectionEdge(StsParameters.BOUNDARY, surface, false); // not persistent
                boundaryEdge.setDrawPoints(true);
				model.addDisplayableInstance(boundaryEdge);
                gridPoints = (StsGridPoint[])StsMath.arrayAddElement(gridPoints, gridPoint);
            }
            else if(nPicks > 0 && nearLastPick(gridPoint))
            {
                logMessage("Picked point identical to previous one. Please reselect.");
                return true;
            }
            else if(nPicks >= 3 && nearFirstPick(gridPoint))  // finished a closed loop
            {
                isLoop = true;
                glPanel.actionManager.endCurrentAction();
            }
            else
                gridPoints = (StsGridPoint[])StsMath.arrayAddElement(gridPoints, gridPoint);

            StsPoint point = surface.getPointZorT(gridPoint.row, gridPoint.col);
            boundaryEdge.addGridEdgePoint(new StsGridSectionPoint(point, surface, false));
			boundaryEdge.getPointsFromGridEdgePoints();
			model.win3d.win3dDisplay();
        }
        return true;
	}

    public boolean end()
    {
        boolean success = constructBoundarySections();

        surface.setIsVisible(visible);
        statusArea.textOnly();
        if (success)
        {
            if (isLoop) logMessage("Closed-boundary successfully created.");
            else logMessage("Boundary successfully created.  You must close "
                    + "the boundary loop with picked fault sections.");
            model.setActionStatus(StsBuildBoundary.class.getName(), StsModel.ENDED);
        }
        else logMessage("Boundary picking unsuccessful.");
        model.win3d.removeToolbar(StsBuildBoundaryToolbar.NAME);
		model.removeDisplayableInstance(boundaryEdge);
        model.win3d.win3dDisplay();
        return success;
    }

    protected boolean nearFirstPick(StsGridPoint gridPoint)
    {
        return nearPick(gridPoint, gridPoints[0]);
    }

    protected boolean nearLastPick(StsGridPoint gridPoint)
    {
        return nearPick(gridPoint, gridPoints[gridPoints.length-1]);
    }

    protected boolean nearPick(StsGridPoint gridPoint0, StsGridPoint gridPoint1)
    {
        float dRowF = gridPoint0.rowF - gridPoint1.rowF;
        float dColF = gridPoint0.colF - gridPoint1.colF;
        return dRowF*dRowF + dColF*dColF < gridNearSq;
    }

    // Assume that all picks are aligned on rows and columns
    protected boolean constructBoundarySections()
    {
        int i, row, col, nSections;
        try
        {
            int nPoints = (gridPoints != null ? gridPoints.length : 0);
            if(nPoints == 0) return false;

            // calculate alignments for each grid point
            StsSectionGeometry[] alignments = new StsSectionGeometry[nPoints];
            for (i=0; i<nPoints; i++)
                alignments[i] = new StsSectionGeometry(gridPoints[i], gridPoints[(i+1)%nPoints], StsSectionGeometry.RIBBON);

            // construct vertical pseudo-wells at each grid point
            StsLine[] pseudos = new StsLine[nPoints];
            for (i=0; i<nPoints; i++)
                pseudos[i] = StsLine.buildVertical(gridPoints[i], StsParameters.BOUNDARY);

            // Build sections from vertical pseudos
            nSections = isLoop ? nPoints : nPoints - 1;

			boolean constructionOK = true;
            for (i=0; i<nSections; i++)
            {
				StsSection section = StsSection.constructor(null, StsParameters.BOUNDARY,
									    pseudos[i], pseudos[(i+1)%nPoints], alignments[i]);
				if(section == null)
				{
					logMessage("PolygonBoundary section constructor failed for side " +
								    i + ". : see error log for details.");
					constructionOK = false;
				}

            }
		    if(!constructionOK) return false;
            return true;
        }
        catch(Exception e)
        {
            StsException.outputException("StsPolygonBoundary.constructBoundarySections() failed.",
                                            e, StsException.WARNING);
            return false;
        }
    }
}



