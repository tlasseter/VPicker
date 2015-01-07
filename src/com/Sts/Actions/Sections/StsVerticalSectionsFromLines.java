

//Title:        S2S: Seismic-to-simulation
//Version:
//Copyright:    Copyright (c) 2001
//Author:       Tom Lasseter
//Company:      4D Systems LLC
//Description:  Web-enabled Integrated Interpretation System

package com.Sts.Actions.Sections;

import com.Sts.Actions.*;
import com.Sts.DBTypes.*;
import com.Sts.MVC.*;
import com.Sts.MVC.View3d.*;
import com.Sts.Types.*;
import com.Sts.Utilities.*;

public class StsVerticalSectionsFromLines extends StsAction
{
    private StsModelSurface surface;
    private boolean visible = true;
    private StsSectionEdge edge;
    private StsSection section;
    private StsLine[] pickedLines = new StsLine[2];
    private int nLinesPicked = 0;
    private String surfaceName = null;

    static int nPointsPerVertex = 1;
    static final byte FAULT = StsParameters.FAULT;
	static final byte BOUNDARY = StsParameters.BOUNDARY;

       /** create a new section and attach a mouse listener to it */
     public StsVerticalSectionsFromLines(StsActionManager actionManager)
    {
        super(actionManager);
    }

    public void initializeRepeatAction(StsAction lastAction)
    {
        surface = ((StsVerticalSectionsFromLines)lastAction).getSurface();
    }

    public StsModelSurface getSurface() { return surface; }

    private boolean actionOK()
    {
        try
        {
            StsSurfaceClass surfaceClass = (StsSurfaceClass)model.getCreateStsClass(StsSurface.class);
            if (surfaceClass.getSize() <= 0)
            {
                logMessage("No surfaces have been loaded: terminating action.");
                return false;
            }
            return true;
        }
        catch(Exception e)
        {
            StsException.outputException("StsVerticalSectionsFromLines.actionOK() failed.",
                e, StsException.WARNING);
            return false;
        }
    }

    public boolean start()
    {
        if( !actionOK() )  return false;

        try
        {

            // select a horizon
            if(surface == null)
            {
                StsModelSurfaceClass surfaceClass = (StsModelSurfaceClass)model.getCreateStsClass(StsModelSurface.class);
                surface = (StsModelSurface)surfaceClass.selectSurface(true, StsModelSurface.class); // try visible 1st
            }
            if (surface == null)
            {
                statusArea.textOnly();
                logMessage("No horizon selected.");
                return false;
            }
            visible = surface.getIsVisible();
            surface.setIsVisible(true);

            surfaceName = surface.getName();

            if(!surface.checkIsLoaded())
            {
                statusArea.textOnly();
                logMessage("Unable to load " + surfaceName + ".");
                return false;
            }
            visible = surface.getIsVisible();
            surface.setIsVisible(true);


            logMessage("Select a fault line, digitize fault cut on " + surfaceName
                    + ", then finish by selecting another fault line.");

            addAbortButton();
            addEndAllButton();

            /** Allocate a section; this may be temporary if we find
             * that we are adding an edge to an existing section */

//            model.viewProperties.save();
//            model.viewProperties.set("Sections", true);
//            model.viewProperties.set("Patches", false);

//            model.viewProperties.set("EditingSections", true);

            for(int n = 0; n < 2; n++) pickedLines[n] = null;

//            win3d.glPanel3d.set3dOverlay(true);
            model.win3d.win3dDisplay();
            return true;
        }
        catch(Exception e)
        {
            StsException.outputException("StsVerticalSectionsFromLines.start() failed.",
                e, StsException.WARNING);
            return false;
        }
    }

    /** mouse action for 3d window */
       public boolean performMouseAction(StsMouse mouse, StsGLPanel glPanel)
    {
        StsLine line;
        StsLine firstLinePicked = null;
        try
        {
            int leftButtonState = mouse.getButtonStateCheckClear(StsMouse.LEFT);
            if(leftButtonState == StsMouse.PRESSED || leftButtonState == StsMouse.DRAGGED)
            {
                if(nLinesPicked == 1)
                {
                    StsMousePoint mousePoint = mouse.getMousePoint();
                    surface.getSurfacePosition(mouse, true, (StsGLPanel3d)glPanel);
                }
            }
            else if(StsMouse.isButtonStateReleasedOrClicked(leftButtonState))  // add a point
            {
                /** first check if line has been picked; if so, add it to lineEnties list.
                  * If no line pick and one line has been picked,  compute intersection on surface
                  */
                StsGLPanel3d glPanel3d = (StsGLPanel3d)glPanel;
                StsObject[] visibleLines = getVisibleLines();
                line = (StsLine)StsJOGLPick.pickClass3d(glPanel3d, visibleLines, StsJOGLPick.PICKSIZE_MEDIUM, StsJOGLPick.PICK_CLOSEST);
                boolean lineVisible = ( line != null && line.getIsVisible() );

                if(lineVisible && line != pickedLines[0])
                {
                    pickedLines[nLinesPicked] = line;

                    // find intersection of line line with grid, then snap to line line
                    StsGridPoint gridPoint = line.computeGridIntersect(surface);
                    StsPoint linePoint = line.getXYZPointAtZorT(gridPoint.getZorT(), true); // snap to line line

                    nLinesPicked++;

                    StsSurfaceVertex vertex = new StsSurfaceVertex(linePoint, line, surface, null, true);
//                    line.insertSurfaceVertex(vertex);

                    if(nLinesPicked == 1)
                    {
                        firstLinePicked = line;
                        // section.addEndLine(line);
                        edge = new StsSectionEdge(FAULT, null, surface, nPointsPerVertex);
                        edge.setPointsFromEdgePoints();
                        model.addDisplayableInstance(edge);
                        edge.addPrevVertex(vertex);
                        edge.setDrawPoints(true);
                    }
                    else
                    {
                        if(section == null)
                            section = new StsSection(FAULT, firstLinePicked, line);
                        edge.addNextVertex(vertex);
                        section.addSectionEdge(edge);

                        section.completeSection();
                        // section.constructSection();
                        glPanel.actionManager.endCurrentAction();
                    }
                }
                else
                    StsMessageFiles.logMessage("No line selected to make section.");
            }
            return true;
        }
        catch(Exception e)
        {
            StsException.outputException("StsVerticalSectionsFromLines.performMouseAction() failed.",
                e, StsException.WARNING);
            return false;
        }
    }

    private StsObject[] getVisibleLines()
    {
        StsObject[] visibleLines = model.getVisibleObjectList(StsLine.class);
        StsObject[] visibleFaultLines = model.getVisibleObjectList(StsFaultLine.class);
        return (StsObject[])StsMath.arrayAddArray(visibleLines, visibleFaultLines);
    }

    public boolean end()
    {
        boolean ok = true;

        if(nLinesPicked == 0) return true; // nothing built; everthing ok
        surface.setIsVisible(visible);
        model.removeDisplayableInstance(edge);
//        model.viewProperties.restore();
        StsLine.clearHighlightedLines();
        statusArea.textOnly();

        if(nLinesPicked != 2)
        {
            section.delete();
            ok = false;
            logMessage("Fault section not built.");
        }
        else
        {
            logMessage("Section:  " + section.getLabel() + " built successfully.");
            ok = true;
        }
//        win3d.glPanel3d.set3dOverlay(false);
        model.win3d.win3dDisplay();
        return ok;
    }
}


