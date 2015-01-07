

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

public class StsVerticalSectionsOnSurface extends StsAction
{
    private StsSurface surface;
    private boolean visible = true;
    private StsSectionEdge edge;
    private StsSection section;
    private StsLine prevLine = null;
    private StsLine nextLine = null;
    private int nLinesPicked = 0;
    private String surfaceName = null;

    static private int nPointsPerVertex = 1;

       /** create a new section and attach a mouse listener to it */
     public StsVerticalSectionsOnSurface(StsActionManager actionManager)
    {
        super(actionManager);
    }

    public void initializeRepeatAction(StsAction lastAction)
    {
        surface = ((StsVerticalSectionsOnSurface)lastAction).getSurface();
    }

    public StsSurface getSurface() { return surface; }

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
            StsException.outputException("StsVerticalSectionsOnSurface.actionOK() failed.",
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
                StsSurfaceClass surfaceClass = (StsSurfaceClass)model.getCreateStsClass(StsSurface.class);
                surface = (StsSurface)surfaceClass.selectSurface(true, StsSurface.class); // try visible 1st
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

            prevLine = null;
            nextLine = null;

//            win3d.glPanel3d.set3dOverlay(true);
            model.win3d.win3dDisplay();
            return true;
        }
        catch(Exception e)
        {
            StsException.outputException("StsVerticalSectionsOnSurface.start() failed.",
                e, StsException.WARNING);
            return false;
        }
    }

    /** mouse action for 3d window */
    public boolean performMouseAction(StsMouse mouse, StsGLPanel glPanel)
    {
        try
        {
            StsSurfaceClass surfaceClass = (StsSurfaceClass) model.getStsClass(StsSurface.class);
            StsSurface surface = surfaceClass.getTopVisibleSurface();
            if (surface == null)return false;

            int leftButtonState = mouse.getButtonStateCheckClear(StsMouse.LEFT);

            if (leftButtonState == StsMouse.PRESSED || leftButtonState == StsMouse.DRAGGED)
            {
                surface.getSurfacePosition(mouse, true, (StsGLPanel3d)glPanel);
            }
            else if (StsMouse.isButtonStateReleasedOrClicked(leftButtonState)) // add a point
            {
                StsGridPoint gridPoint = surface.getSurfacePosition(mouse, true, (StsGLPanel3d)glPanel);
                if (gridPoint == null)return true;
                StsLine line = StsLine.buildVertical(gridPoint, StsParameters.FAULT);

                if (line != null)
                {
                    prevLine = nextLine;
                    nextLine = line;
                    if (prevLine != null)
                    {
                        StsSection section = StsSection.constructor(null, StsParameters.AUXILIARY,
                                                prevLine, nextLine, StsSectionGeometry.RIBBON);
                        section.completeSection();
                        // section.constructSection();
                    }
                }
                model.win3d.win3dDisplay();
            }
            return true;
        }
        catch(Exception e)
        {
            StsException.outputException("StsVerticalSectionsOnSurface.performMouseAction() failed.",
                e, StsException.WARNING);
            return false;
        }
    }

    public boolean end()
    {
        return true;
    }
}


