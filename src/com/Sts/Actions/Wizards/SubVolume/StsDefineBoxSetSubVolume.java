package com.Sts.Actions.Wizards.SubVolume;

import com.Sts.Actions.Wizards.CombinationVolume.*;
import com.Sts.Actions.Wizards.*;
import com.Sts.Actions.Wizards.WizardHeaders.*;
import com.Sts.DBTypes.*;
import com.Sts.MVC.*;
import com.Sts.MVC.View3d.*;
import com.Sts.Types.*;
import com.Sts.UI.*;
import com.Sts.Utilities.*;

import java.awt.*;

/**
 * <p>Title: S2S Development</p>
 * <p>Description: </p>
 * <p>Copyright: Copyright (c) 2002</p>
 * <p>Company: S2S Systems LLC</p>
 * @author unascribed
 * @version 1.1
 */

public class StsDefineBoxSetSubVolume extends StsWizardStep
{
    StsDefineBoxSetSubVolumePanel panel;
    StsHeaderPanel header;

    StsBoxSetSubVolume boxSetSubVolume;
    StsColor stsColor;
    private boolean editCube = false;
    private StsCursorPoint cubePoint = null;
    private int state = ACTION_NONE;

    public static final byte ACTION_NONE = 0;
    public static final byte ACTION_DEFINE = 1;
    public static final byte ACTION_MOVE = 2;
    public static final byte ACTION_EDIT = 3;
    public static final byte ACTION_DELETE = 4;

    StsBoxSubVolume currentBox = null;
    boolean visible = false;

    static int suffix = 1;
    static StsBoxSetSubVolumeClass boxSetClass = null;

    public StsDefineBoxSetSubVolume(StsWizard wizard)
    {
        super(wizard);
        panel = new StsDefineBoxSetSubVolumePanel(wizard, this);
        header = new StsHeaderPanel();
        setPanels(panel, header);
        panel.setPreferredSize(new Dimension(500, 500));
        header.setTitle("Sub-Volume Definition");
        header.setSubtitle("Defining a subvolume consisting of a set of boxes");
        header.setLink("http://www.s2ssystems.com/Protected/Marketing/AppLinks.html#SubVolume");
        header.setInfoText(wizardDialog, "(1) Specify the name of the box set sub-volume.\n" +
                           "(2) Select the desired color for all the hexahedron in the set.\n" +
                           "(3) Specify the default size of the hexahedron when they are created.\n" +
                           "    **** Once created, the size and location can be graphically edited.\n" +
                           "(4a) Press the Create Button to begin creating hexahedron by selecting a center point on the seismic.\n" +
                           "(4b) Press the Delete Button to enter delete mode, select existing center point to delete hexahedron.\n" +
                           "(4c) Press the Move Button to enter move mode, select existing center point and move mouse to move hexahedron.\n" +
                           "(4d) Press the Edit Button to enter edit mode, select white handles on any face to stretch hexahedron in that direction.\n" +
                           "    **** This set can be edited later by re-entering this wizard and selecting this set ****\n" +
                           "(5) Press the Next> Button when all the hexahedron in this set have been created, to proceed to the setup screen\n");
    }

    public boolean start()
    {
        panel.initializeButtonGroup();
        state = ACTION_DEFINE;
        boxSetClass = (StsBoxSetSubVolumeClass) model.getStsClass(StsBoxSetSubVolume.class);

        StsSubVolumeClass subVolumeClass = (StsSubVolumeClass)model.getStsClass(StsSubVolume.class);
        subVolumeClass.setIsVisible(false);
        subVolumeClass.setIsApplied(false);

        boxSetSubVolume = (StsBoxSetSubVolume)boxSetClass.getCurrentObject();
        if (boxSetSubVolume == null)
        {
            boxSetSubVolume = new StsBoxSetSubVolume(false, ((StsSubVolumeWizard)wizard).getSubVolumeName());
            boxSetSubVolume.setStsColor(panel.getStsColor());
            boxSetClass.setCurrentObject(boxSetSubVolume);
			boxSetSubVolume.addToModel();
        }
        boxSetSubVolume.setIsVisible(true);
        boxSetSubVolume.setIsApplied(false);

        updateAll();
        return true;
    }

    private void createBoxSet(String name)
    {
    	StsSubVolumeWizard svWizard = (StsSubVolumeWizard)wizard;
        boxSetSubVolume = new StsBoxSetSubVolume(svWizard.getSubVolumeName());
        byte zDomain = model.getProject().getZDomain();
        boxSetSubVolume.setZDomain(zDomain);
        if(wizard instanceof StsSubVolumeWizard)
            ((StsSubVolumeWizard)wizard).setSubVolume(boxSetSubVolume);
        else
            ((StsCombinationVolumeWizard)wizard).setSubVolume(boxSetSubVolume);
    }

    public boolean end()
    {
        if(wizard instanceof StsSubVolumeWizard)
        {
            if(((StsSubVolumeWizard)wizard).isCanceled())
            {
                if(boxSetSubVolume != null)
                {
                    boxSetSubVolume.deleteBoxes();
                    boxSetSubVolume.delete();
                }
                boxSetSubVolume = null;
                model.refreshObjectPanel();
                return false;
            }
        }
        StsSubVolumeClass subVolumeClass = (StsSubVolumeClass)model.getStsClass(StsSubVolume.class);
        subVolumeClass.setIsVisible(false);
        subVolumeClass.setIsApplied(true);

        if (boxSetSubVolume == null)
            return false;

        if(boxSetSubVolume.getNBoxes() <= 0)
            boxSetSubVolume.delete();
        else
            boxSetSubVolume.setName(((StsSubVolumeWizard)wizard).getSubVolumeName());

        boxSetSubVolume = null;
		model.subVolumeChanged();
        return true;
    }

    public void setBoxSet(StsBoxSetSubVolume boxSet)
    {
		if(boxSetSubVolume != null) boxSetSubVolume.setIsVisible(false);
		if(boxSet == null) return;
        boxSetSubVolume = boxSet;
		boxSetSubVolume.setIsVisible(true);
    }

    public void setDefineCube()
    {
        state = ACTION_DEFINE;
        wizard.enableFinish();
    }

    public void setMoveCube()
    {
        state = ACTION_MOVE;
        boxSetSubVolume.setAction(ACTION_NONE);
    }

    public void setEditCube()
    {
        state = ACTION_EDIT;
        boxSetSubVolume.setAction(ACTION_NONE);
    }

    public void setDeleteCube()
    {
        state = ACTION_DELETE;
    }

    public void resetCubeState()
    {
        state = ACTION_NONE;
    }

    public boolean performMouseAction(StsMouse mouse, StsGLPanel glPanel)
    {
        if(boxSetSubVolume == null) return true;
        StsGLPanel3d glPanel3d = (StsGLPanel3d)glPanel;
        if(state == ACTION_DEFINE) return defineCube(mouse, glPanel3d);
        else if(state == ACTION_MOVE) return moveCube(mouse, glPanel3d);
        else if(state == ACTION_EDIT) return editCube(mouse, glPanel3d);
        else if(state == ACTION_DELETE)
        {
            boolean ok = deleteCube(mouse, glPanel3d);
            if(ok)
            {
                int nBoxes = boxSetSubVolume.getBoxes().getSize();
                if(nBoxes == 0)
                {
                    boxSetSubVolume.delete();
                    boxSetSubVolume = (StsBoxSetSubVolume)boxSetClass.getCurrentObject();
                    updateAll();
                }
                return ok;
            }
        }
        return true;
    }

    private void updateAll()
    {
        panel.setPanelObject(boxSetSubVolume);
        panel.updateBoxProperties();
    }

    private boolean defineCube(StsMouse mouse, StsGLPanel3d glPanel3d)
    {
        if( mouse.getCurrentButton() != StsMouse.LEFT ) return true;
        int buttonState = mouse.getButtonStateCheckClear(StsMouse.LEFT);

        // set temporary cursorPoint on mouse PRESS or DRAG
        if (buttonState == StsMouse.PRESSED || buttonState == StsMouse.DRAGGED)
        {
            StsView currentView = glPanel3d.getView();
            if (currentView instanceof StsView3d)
            {
                StsCursor3d cursor3d = glPanel3d.getCursor3d();
                cubePoint = cursor3d.getCursorPoint(glPanel3d, mouse);
                if(cubePoint == null) return true;
                cursor3d.setCurrentDirNo(cubePoint.dirNo);
            }
            else if (currentView instanceof StsViewCursor)
            {
                cubePoint = ( (StsViewCursor) currentView).getCursorPoint(mouse);
                if(cubePoint == null) return true;
            }
            else
                return true;
        }
        // permanently add this point when mouse.RELEASED
        else if (buttonState == StsMouse.RELEASED)
        {
            StsSeismicVolume currentVolume = (StsSeismicVolume)model.getStsClass(StsSeismicVolume.class).getCurrentObject();
            if(currentVolume == null)
            {
                new StsMessage(model.win3d, StsMessage.WARNING, "No seismic volumes loaded.");
                return true;
            }
            if(cubePoint == null) return true;
            currentBox = new StsBoxSubVolume(cubePoint, currentVolume, panel.getNBoxRows(), panel.getNBoxCols(), panel.getNBoxSlices());
            boxSetSubVolume.add(currentBox);
            updateAll();
            model.subVolumeChanged();
        }
        model.win3dDisplay();
        return true;
    }

    private boolean moveCube(StsMouse mouse, StsGLPanel3d glPanel3d)
    {
        boolean moved = boxSetSubVolume.move(mouse, glPanel3d);
        if(moved)
        {
            updateAll();
            model.subVolumeChanged();
        }
        return moved;
    }

    private boolean editCube(StsMouse mouse, StsGLPanel3d glPanel3d)
    {
        boolean edited = boxSetSubVolume.edit(mouse, glPanel3d);
        currentBox = boxSetSubVolume.getCurrentBox();
        if(edited)
        {
            updateAll();
            model.subVolumeChanged();
        }
        return edited;
    }

    private boolean deleteCube(StsMouse mouse, StsGLPanel3d glPanel3d)
    {
        boolean deleted = boxSetSubVolume.delete(mouse, glPanel3d);
        if(deleted)
        {
            updateAll();
            model.subVolumeChanged();
        }
        return deleted;
    }


    public void setStsColor(StsColor color) { if(boxSetSubVolume != null) boxSetSubVolume.setStsColor(color); }
    public StsColor getStsColor()
    {
        if(boxSetSubVolume == null) return null;
        return boxSetSubVolume.getStsColor();
    }

    public void setBoxName(String name) { if(boxSetSubVolume != null) boxSetSubVolume.setName(name); }
    public void setBoxXMin(float xMin) { if(boxSetSubVolume != null) currentBox.adjustXMin(xMin); model.win3dDisplay(); }
    public void setBoxXMax(float xMax) { if(boxSetSubVolume != null) currentBox.adjustXMax(xMax); model.win3dDisplay(); }
    public void setBoxYMin(float yMin) { if(boxSetSubVolume != null) currentBox.adjustYMin(yMin); model.win3dDisplay(); }
    public void setBoxYMax(float yMax) { if(boxSetSubVolume != null) currentBox.adjustYMax(yMax); model.win3dDisplay(); }
    public void setBoxZMin(float zMin) { if(boxSetSubVolume != null) currentBox.adjustZMin(zMin); model.win3dDisplay(); }
    public void setBoxZMax(float zMax) { if(boxSetSubVolume != null) currentBox.adjustZMax(zMax); model.win3dDisplay(); }

    public String getBoxName() 
    { 
    	StsBoxSetSubVolumeClass subVolumeClass = (StsBoxSetSubVolumeClass)model.getStsClass(StsBoxSetSubVolume.class);    	
    	if(boxSetSubVolume == null) 
    		return "boxes" + subVolumeClass.getSize(); 
    	else 
    		return boxSetSubVolume.getName(); 
    }
    public float getBoxXMin() { if(currentBox == null) return StsParameters.nullValue; else return currentBox.getXMin(); }
    public float getBoxXMax() { if(currentBox == null) return StsParameters.nullValue; else return currentBox.getXMax(); }
    public float getBoxYMin() { if(currentBox == null) return StsParameters.nullValue; else return currentBox.getYMin(); }
    public float getBoxYMax() { if(currentBox == null) return StsParameters.nullValue; else return currentBox.getYMax(); }
    public float getBoxZMin() { if(currentBox == null) return StsParameters.nullValue; else return currentBox.getZMin(); }
    public float getBoxZMax() { if(currentBox == null) return StsParameters.nullValue; else return currentBox.getZMax(); }

}
