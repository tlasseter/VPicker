package com.Sts.Actions.Wizards.EditTd;

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

 public class StsEditTdCurve extends StsWizardStep
{
	public StsEditTdCurve()
	{
		try
		{
			jbInit();
		}
		catch(Exception ex)
		{
			ex.printStackTrace();
		}
	}

	StsEditTdCurvePanel panel;
    StsHeaderPanel header, info;

    StsWell selectedWell = null;
    StsEditTdSet editTdSet = null;
    StsPoint selectedPoint = null;
//    StsLine lineThruPoint = null;

    int mode = SELECT_WELL;

    static int SELECT_WELL = 0;
    static int SELECT_POINT = 1;

    static final boolean debug = false;

    public StsEditTdCurve(StsWizard wizard)
    {
        super(wizard);
        panel = new StsEditTdCurvePanel(this);
        header = new StsHeaderPanel();
        setPanels(panel, header);
//        int height = header.displayableHeight + panel.displayableHeight;
//        panel.setPreferredSize(new Dimension(300, height));
        header.setTitle("TD Editing");
        header.setSubtitle("Selecting Well");
        header.setLink("http://www.s2ssystems.com/Protected/Marketing/AppLinks.html#EditTD");
        header.setInfoText(wizardDialog,"(1) Select the well from list or in 3D window to edit.\n" +
                           "(2) Display or add depth markers on well.\n" +
                           " (3) Move marker in 3D window vertically to match correct seismic event.\n");
    }

    public boolean start()
     {
         panel.initialize(getVisibleWellList());
         wizard.enableFinish();
         StsPoint.setCompareIndex(StsPoint.mIndex);
         return true;
     }

     public boolean end()
     {
         if(editTdSet == null) return true;
         editTdSet.endTransaction();
         editTdSet.deleteSeismicCurtain();
         model.removeDisplayableInstance(editTdSet);
         return true;
     }

    public Object[] getVisibleWellList()
    {
        Object[] wells = model.getVisibleObjectList(StsWell.class);
        int nWells = wells.length;
        Object[] wellList = new Object[nWells+1];
        StsWell nullWell = StsWell.nullWellConstructor("none");
        wellList[0] = nullWell;
        for(int n = 0; n < nWells; n++)
        {
            wellList[n+1] = wells[n];
        }
		setSelectedWell(wellList[0]);
        return wellList;
    }

    public void setSelectedWell(Object object)
    {
        StsWell newSelectedWell = (StsWell)object;
        if(editTdSet != null) editTdSet.deleteSeismicCurtain();
        if(newSelectedWell.getName().equals("none"))
        {
            selectedWell = null;
        }
        else
        {
            selectedWell = newSelectedWell;
            editTdSet = getTdEdits();
            editTdSet.createSeismicCurtain(selectedWell);
            model.addDisplayableInstance(editTdSet);
        }
        model.win3dDisplay();
    }

    StsEditTdSet getEditTdSet(StsWell well)
    {
        if(!well.isPersistent()) return null;
        if(editTdSet != null) editTdSet.endTransaction();
        editTdSet = StsEditTdSet.getCreateEditTdSet(well);
        return editTdSet;
    }

    public StsEditTdSet getTdEdits()
    {
        if(selectedWell == null) return null;
        return getEditTdSet(selectedWell);
    }

    public StsEditTd createTdEdit()
    {
        StsEditTdSet editTdSet = getTdEdits();
        if (editTdSet == null)return null;
        StsEditTd newEditTd = editTdSet.createTdEditFromCurrent();
        if(newEditTd == null) return null;
        panel.initializeTdEditListBean();
        return newEditTd;
    }

    public Object getSelectedWell()
    {
            return selectedWell;
    }

    public void setTdEdit(Object tdEdit)
    {
        if(editTdSet == null) return;
        editTdSet.setTdEdit(tdEdit);
    }

    public Object getTdEdit()
    {
        if(editTdSet == null) return null;
        return editTdSet.getTdEdit();
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
				if(debug) System.out.println(mouse.toString() + " " + mouse.mousePoint.toString());
            }
            else if (leftButtonState == StsMouse.DRAGGED)
            {
                if(mode == SELECT_WELL) return true;

                if (selectedWell == null)
                {
                    StsMessageFiles.errorMessage("No well has been selected.  Try again.");
                    return true;
                }
                StsGLPanel3d glPanel3d = (StsGLPanel3d)glPanel;
                StsPoint pickedPoint = getPointOnLineNearestMouse(glPanel3d, mouse);
                if(pickedPoint != null)
                {
					if(debug) System.out.println(mouse.toString() + " " + mouse.mousePoint.toString());

                    if(debug) System.out.println("Moved point from " + selectedPoint.getT() + " to " + pickedPoint.getT());
                    editTdSet.movePoint(selectedPoint, pickedPoint);
                    model.win3dDisplay();
					panel.wellViewModel.wellWindowPanel.repaint();
                }

//				selectedWell.editTDCurve(selectedPoint, pickedPoint);
            }
            else if (StsMouse.isButtonStateReleasedOrClicked(leftButtonState))
            {
                if (selectedWell == null)
                {
                    return true;
                }
//				StsPoint pickedPoint = selectedWell.getPointOnWellCurtain(mouse, selectedPoint);
//				selectedWell.editTDCurve(selectedPoint, pickedPoint);
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
        if(selectedPoint == null) return null;
        StsPoint point0 = new StsPoint(selectedPoint);
        StsPoint point1 = new StsPoint(selectedPoint);
        float t = selectedPoint.getT();
        point0.setT(t-10.0f);
        point1.setT(t+10.0f);
        return glPanel3d.getPointOnLineNearestMouse(mouse, point0, point1);
    }

    private boolean getPickedWell(StsGLPanel3d glPanel3d, StsMouse mouse)
    {
        try
        {
            StsMainObject[] visibleWells = model.getVisibleObjectList(StsWell.class);
            StsWell pickedWell = (StsWell)StsJOGLPick.pickClass3d(glPanel3d, visibleWells, StsJOGLPick.PICKSIZE_MEDIUM, StsJOGLPick.PICK_CLOSEST);
            if(pickedWell == null)
            {
                StsMessageFiles.logMessage("No well picked.");
                return true;

            }
            if (pickedWell != selectedWell)
            {
                if(debug) System.out.println("New well selected: " + pickedWell.getName());
                mode = SELECT_WELL;
                panel.setSelectedWell(pickedWell);
 //               setSelectedWell(pickedWell);
                return true;
            }
            mode = SELECT_POINT;

            // save this pick info in case we need to add a new selected point
            StsPickItem pickItem = StsJOGLPick.pickItems[0];
            int nSegment = pickItem.names[1];

            StsMainObject[] tdEditsList = new StsMainObject[] { this.editTdSet };
            StsObject pickedObject = StsJOGLPick.pickClass3d(glPanel3d, tdEditsList, StsJOGLPick.PICKSIZE_MEDIUM, StsJOGLPick.PICK_CLOSEST);
            if (pickedObject == null)
            {
                selectedPoint = selectedWell.getPointOnLineNearestMouse(nSegment, mouse, glPanel3d);
                if (selectedPoint == null) return true;
                if(debug) System.out.println("New point being added.");
                editTdSet.addPoint(selectedPoint);
//                lineThruPoint = StsLine.buildVertical(selectedPoint, false);
                model.win3dDisplay();
           }
           else
           {
               pickItem = StsJOGLPick.pickItems[0];
               int pointIndex = pickItem.names[2];
               selectedPoint = editTdSet.setSelectedPoint(pointIndex);
               if(debug) System.out.println("Picked point " + pointIndex);
           }
		   if(debug) System.out.println("selectPoint set. T:" + selectedPoint.getT());
           return true;
       }
        catch (Exception e)
        {
            StsException.outputException("StsEditTDCurve.getPickedWell() failed.",
                                         e, StsException.WARNING);
            return false;
        }
    }

	public void jbInit() throws Exception
	{
	}
}
