package com.Sts.Actions.Wizards.EditWellMarkers;

/**
 * <p>Title: jS2S development</p>
 * <p>Description: </p>
 * <p>Copyright: Copyright (c) 2004</p>
 * <p>Company: S2S Systems LLC</p>
 * @author Tom Lasseter
 * @version 1.0
 */

import com.Sts.Actions.Wizards.*;
import com.Sts.Actions.Wizards.WizardHeaders.*;
import com.Sts.DBTypes.*;
import com.Sts.MVC.*;
import com.Sts.MVC.View3d.*;
import com.Sts.Types.*;
import com.Sts.Utilities.*;

public class StsEditWellMarkers extends StsWizardStep
{
    StsEditWellMarkersPanel panel;
    StsHeaderPanel header, info;
    StsWell selectedWell = null;
    StsPoint selectedPoint = null;

    int mode = SELECT_WELL;

    static int SELECT_WELL = 0;
    static int SELECT_POINT = 1;

    static final boolean debug = false;

    public StsEditWellMarkers(StsEditWellMarkersWizard wizard)
    {
        super(wizard);
        panel = new StsEditWellMarkersPanel(this);
        header = new StsHeaderPanel();
        setPanels(panel, header);
        header.setTitle("Well Marker Editing");
        header.setSubtitle("Selecting Well");
        header.setLink("http://www.s2ssystems.com/Protected/Marketing/AppLinks.html#EditWellMarkers");
        header.setInfoText(wizardDialog,
                           "(1) Select the well from list or in 3D window to edit.\n" +
                           "(2) Select operation (Add or Edit).\n" +
                           "------------------- Adding a Marker -----------------\n" +
                           "(3) Select the pick type\n" +
                           "(4) Using the mouse in the well track set the measured depth.\n" +
                           "(5) Define remaining type specific properties\n" +
                           "(6) Press the Accept Button\n" +
                           "------------------ Editing a Marker -----------------\n" +
                           "(3) Select the marker to edit.\n" +
                           "(4) Adjust the measured depth using the mouse in the well track.\n" +
                           "(5) Adjust remaining type specific properties.\n" +
                           "(6) Press Accept Button.\n\n" +
                           "(7) Press Finish when all edits and additions are complete.");
    }

    public boolean start()
    {
        panel.initialize(getVisibleWellList());
        wizard.enableFinish();
        ((StsEditWellMarkersWizard)wizard).disableCancelBtn();
        StsPoint.setCompareIndex(StsPoint.mIndex);
        return true;
    }

    public boolean end()
    {
        return true;
    }

    public StsObject[] getVisibleWellList()
    {
        StsObject[] wells = model.getVisibleObjectList(StsWell.class);
        int nWells = wells.length;
        StsObject[] wellList = new StsObject[nWells + 1];
        StsWell nullWell = StsWell.nullWellConstructor("none");
        wellList[0] = nullWell;
        for (int n = 1; n <= nWells; n++)
        {
            wellList[n] = wells[n-1];
        }
		selectedWell = (StsWell) wellList[0];

        StsObject[] lwells = model.getVisibleObjectList(StsLiveWell.class);
        if(lwells.length == 0)
            return wellList;

        int nLiveWells = lwells.length;
        StsObject[] lwellList = new StsObject[nLiveWells];
        for (int n = 0; n < nLiveWells; n++)
        {
            lwellList[n] = lwells[n];
        }
        wellList = (StsObject[])StsMath.arrayAddArray(wellList, lwellList);
        return wellList;
    }

    public void setSelectedWell(Object object)
    {
        selectedWell = (StsWell) object;
        model.win3dDisplay();
    }

    public Object getSelectedWell()
    {
        return selectedWell;
    }

    /** mouse action for 3d window */

    public boolean performMouseAction(StsMouse mouse, StsGLPanel glPanel)
    {
        try
        {
            int leftButtonState = mouse.getButtonStateCheckClear(StsMouse.LEFT);

            if (leftButtonState == StsMouse.PRESSED)
            {
                StsGLPanel3d glPanel3d = (StsGLPanel3d)glPanel;
                getPickedWell(glPanel3d, mouse);
            }
            else if (leftButtonState == StsMouse.DRAGGED)
            {
//                    if(mode == SELECT_WELL) return true;

                if (selectedWell == null)
                {
                    StsMessageFiles.errorMessage("No well has been selected.  Try again.");
                    return true;
                }
                StsGLPanel3d glPanel3d = (StsGLPanel3d)glPanel;
                StsPoint pickedPoint = getPointOnLineNearestMouse(glPanel3d, mouse);
                if (pickedPoint != null)
                {
                    model.win3dDisplay();
                }
            }
            else if (StsMouse.isButtonStateReleasedOrClicked(leftButtonState))
            {
                if (selectedWell == null)
                {
                    return true;
                }
                selectedPoint = null;
            }
            return true;
        }
        catch (Exception e)
        {
            StsException.outputException("Failed during pick.", e,
                                         StsException.WARNING);
            return false;
        }
    }

    private StsPoint getPointOnLineNearestMouse(StsGLPanel3d glPanel3d, StsMouse mouse)
    {
        if (selectedPoint == null)
        {
            return null;
        }
        StsPoint point0 = new StsPoint(selectedPoint);
        StsPoint point1 = new StsPoint(selectedPoint);
        float t = selectedPoint.getT();
        point0.setT(t - 10.0f);
        point1.setT(t + 10.0f);
        return glPanel3d.getPointOnLineNearestMouse(mouse, point0, point1);
    }

    private boolean getPickedWell(StsGLPanel3d glPanel3d, StsMouse mouse)
    {
        try
        {
            StsObject[] visibleWells = (StsObject[])getVisibleWellList();          // ToDo: Debug this...
            StsWell pickedWell = (StsWell)StsJOGLPick.pickClass3d(glPanel3d, visibleWells, StsJOGLPick.PICKSIZE_MEDIUM, StsJOGLPick.PICK_CLOSEST);
            if (pickedWell == null)
            {
                StsMessageFiles.logMessage("No well picked.");
                return true;

            }
            if (pickedWell != selectedWell)
            {
                if (debug)
                {
                    System.out.println("New well selected: " + pickedWell.getName());
                }
                mode = SELECT_WELL;
                panel.setSelectedWell(pickedWell);
                //               setSelectedWell(pickedWell);
                return true;
            }
            mode = SELECT_POINT;

            // save this pick info in case we need to add a new selected point
            StsPickItem pickItem = StsJOGLPick.pickItems[0];
            int nSegment = pickItem.names[1];
/*
            if (!glPanel.glc.pick(picker))
            {
                selectedPoint = selectedWell.getPointOnLineNearestMouse(nSegment, mouse, (StsGLPanel3d) glPanel);
                if (selectedPoint == null)
                {
                    return true;
                }
                if (mainDebug)
                {
                    System.out.println("New point being added.");
                }
                model.win3dDisplay();
            }
            else
            {
*/
                pickItem = StsJOGLPick.pickItems[0];
                int pointIndex = pickItem.names[1];
                if (debug)
                {
                    System.out.println("Picked point " + pointIndex);
                }
//            }
            return true;
        }
        catch (Exception e)
        {
            StsException.outputException("StsEditTDCurve.getPickedWell() failed.",
                                         e, StsException.WARNING);
            return false;
        }
    }

    public void addMarker()
    {
    }

    public void editMarker()
    {
    }

    public void deleteMarker()
    {
    }
}
