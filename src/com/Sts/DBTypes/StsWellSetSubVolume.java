package com.Sts.DBTypes;

import com.Sts.Actions.Wizards.SubVolume.*;
import com.Sts.DB.DBCommand.*;
import com.Sts.Interfaces.*;
import com.Sts.MVC.*;
import com.Sts.MVC.View3d.*;
import com.Sts.Reflect.*;
import com.Sts.Types.*;
import com.Sts.UI.Beans.*;
import com.Sts.UI.ObjectPanel.*;
import com.Sts.UI.*;
import com.Sts.Utilities.*;

import javax.media.opengl.*;
import java.awt.*;

/**
 * <p>Title: </p>
 * <p>Description: </p>
 * <p>Copyright: Copyright (c) 2002</p>
 * <p>Company: </p>
 * @author not attributable
 * @version 1.0
 */

public class StsWellSetSubVolume extends StsSubVolume implements StsTreeObjectI
{
    private StsObjectRefList wellSvs;
    protected StsColor stsColor = new StsColor(Color.RED);
    transient boolean isEditing = false;
    transient StsWellSubVolume currentSubVol = null;
    transient int pickedPointIndex = -1;
    transient byte action = ACTION_NONE;
    transient StsSeismicVolume currentVolume = null;
    transient boolean honorVolume = false;

    static byte ACTION_NONE = StsDefineBoxSetSubVolume.ACTION_NONE;
    static byte ACTION_MOVE_CYL = StsDefineBoxSetSubVolume.ACTION_MOVE;
    static byte ACTION_MOVE_POINT = StsDefineBoxSetSubVolume.ACTION_EDIT;
    static byte ACTION_DELETE = StsDefineBoxSetSubVolume.ACTION_DELETE;

    static protected StsObjectPanel objectPanel = null;

    static public StsFieldBean[] displayFields = null;

    static StsJOGLPick pickCurrentWellSetCenterPoint = null;
    static StsJOGLPick pickCurrentWellSetFacePoint = null;

    public StsWellSetSubVolume()
    {
		this(null, true);
    }

    public StsWellSetSubVolume(String name)
    {
		this(name, true);
    }

    public StsWellSetSubVolume(String name, boolean persistent)
    {
        super(persistent);
		setName(name);
		initializeVisibleFlags();
    }

    public void add(StsWellSubVolume wellsv)
    {
        if(wellSvs == null)
        {
            wellSvs = StsObjectRefList.constructor(4, 4, "wellSvs", this);
        }
        wellSvs.add(wellsv);
        currentSubVol = wellsv;
    }
       
    public boolean delete()
    {
        wellSvs.deleteAll();
        boolean success = super.delete();
        currentModel.viewObjectChangedAndRepaint(this, this);
        return success;        
    }

    public boolean anyDependencies()
    {
        return false;
    }

    public StsFieldBean[] getDisplayFields()
    {
        if(displayFields == null)
        {
            displayFields = new StsFieldBean[]
            {
                new StsStringFieldBean(StsWellSetSubVolume.class, "name", "Name"),
                new StsBooleanFieldBean(StsWellSetSubVolume.class, "isVisible", "Visible"),
                new StsBooleanFieldBean(StsWellSetSubVolume.class, "isApplied", "Applied"),
                new StsBooleanFieldBean(StsWellSetSubVolume.class, "isInclusive", "Inclusive"),
                new StsColorComboBoxFieldBean(StsWellSetSubVolume.class, "stsColor", "Color", currentModel.getSpectrum("Basic").getStsColors())
            };
        }
        return displayFields;
    }

    public StsFieldBean[] getPropertyFields()
    {
		return null;
    }

    public Object[] getChildren()
    {return new Object[0];
    }

    public StsObjectPanel getObjectPanel()
    {
        if (objectPanel == null)
        {
            objectPanel = StsObjectPanel.constructor(this, true);
        }
        return objectPanel;
    }

    public void display(StsGLPanel glPanel, boolean isCurrentObject)
    {
        int nWellSubVols = wellSvs.getSize();
        for (int n = 0; n < nWellSubVols; n++)
        {
            StsWellSubVolume wellsv = (StsWellSubVolume) wellSvs.getElement(n);
            boolean editing = currentSubVol == wellsv;
            wellsv.display(glPanel, stsColor, action, editing);
        }
    }

    public void setStsColor(StsColor stsColor)
    {
        this.stsColor = stsColor;
        currentModel.addTransactionCmd("SubVolume color change", new StsChangeCmd(this, stsColor, "stsColor", false));
        currentModel.win3dDisplayAll();
    }

    public StsColor getStsColor() { return stsColor; }

    public StsGridBoundingBox getBoundingBox()
    {
        StsGridBoundingBox gridBoundingBox = new StsGridBoundingBox(false);
        int nSubVols = wellSvs.getSize();
        for (int n = 0; n < nSubVols; n++)
        {
            StsBoxSubVolume wellSv = (StsBoxSubVolume) wellSvs.getElement(n);
            gridBoundingBox.addUnrotatedBoundingBox(wellSv);
        }
        return gridBoundingBox;
    }

    public void setCurrentVolume(StsSeismicVolume vol)
    {
        currentVolume = vol;
    }

    public void setHonorVolumes(boolean honor)
    {
        honorVolume = honor;
    }

    public void addUnion(byte[] subVolumePlane, int dir, float dirCoordinate, StsSeismicVolume vol, boolean honorVols, byte zDomainData)
    {
        currentVolume = vol;
        honorVolume = honorVols;
        addUnion(subVolumePlane, dir, dirCoordinate, vol, zDomainData);
    }

    public void addUnion(byte[] subVolumePlane, int dir, float dirCoordinate, StsRotatedGridBoundingBox cursor3dBoundingBox, byte zDomainData)
    {

        if (!isApplied) return;

        int nWellSubVols = wellSvs.getSize();
        switch (dir)
        {
            case StsCursor3d.XDIR:
                for (int n = 0; n < nWellSubVols; n++)
                {
                    StsWellSubVolume wellSv = (StsWellSubVolume) wellSvs.getElement(n);
                    if((wellSv.getVolume() != currentVolume) && (honorVolume))
                        continue;

                    // Add code to determine X plane area that bisects current cylinder

                }
                break;
            case StsCursor3d.YDIR:
                for (int n = 0; n < nWellSubVols; n++)
                {
                    StsWellSubVolume wellSv = (StsWellSubVolume) wellSvs.getElement(n);
                    if((wellSv.getVolume() != currentVolume) && (honorVolume))
                        continue;

                    // Add code to determine Y plane area that bisects current cylinder

                }
                break;
            case StsCursor3d.ZDIR:
                for (int n = 0; n < nWellSubVols; n++)
                {
                    StsWellSubVolume wellSv = (StsWellSubVolume) wellSvs.getElement(n);
                    if((wellSv.getVolume() != currentVolume) && (honorVolume))
                        continue;

                    // Add code to determine Z plane area that bisects current cylinder

                }
                break;
        }
    }

    public void setIsEditing(boolean isEditing)
    {this.isEditing = isEditing;
    }

    public int getNSubVols()
    {return wellSvs.getSize();
    }

    public void setAction(byte action)
    {
        this.action = action;
    }

    public boolean move(StsMouse mouse, StsGLPanel3d glPanel3d)
    {
        if (mouse.getCurrentButton() != StsMouse.LEFT)
        {
            return false;
        }
        int buttonState = mouse.getButtonStateCheckClear(StsMouse.LEFT);

        if (buttonState == StsMouse.CLEARED)
        {
            return true;
        }

        checkInitializePickMethods(glPanel3d);

        if (buttonState == StsMouse.PRESSED)
        {
            currentSubVol = pickCenter(glPanel3d);
            if (currentSubVol != null)
            {
                action = ACTION_MOVE_CYL;
            }
            else
            {
                action = ACTION_NONE;
            }
            currentModel.win3dDisplay();
        }
        else if (buttonState == StsMouse.DRAGGED)
        {
            if (currentSubVol == null)
            {
                return true;
            }
            if (action != ACTION_MOVE_CYL)
            {
                StsMessageFiles.errorMessage("StsWellSetSubVolume.action not set to ACTION_MOVE_CYL.");
                return true;
            }
            currentSubVol.moveCylinder(mouse, glPanel3d);
            currentModel.win3dDisplay();
        }
        else if (buttonState == StsMouse.RELEASED)
        {
            action = ACTION_NONE;
        }
        return true;
    }

    private void checkInitializePickMethods(StsGLPanel3d glPanel3d)
    {
        if (pickCurrentWellSetCenterPoint == null)
        {
            GL gl = glPanel3d.getGL();
            StsMethod method = new StsMethod(StsWellSetSubVolume.class, "pickCurrentWellSetCenterPoint", gl, GL.class);
            pickCurrentWellSetCenterPoint = new StsJOGLPick(glPanel3d, method, 10, StsMethodPick.PICK_CLOSEST);
        }
        if (pickCurrentWellSetFacePoint == null)
        {
            GL gl = glPanel3d.getGL();
            StsMethod method = new StsMethod(StsWellSetSubVolume.class, "pickCurrentWellSetFacePoint", gl, GL.class);
            pickCurrentWellSetFacePoint = new StsJOGLPick(glPanel3d, method, 10, StsMethodPick.PICK_CLOSEST);
        }
    }


    public boolean edit(StsMouse mouse, StsGLPanel3d glPanel3d)
    {
        if (mouse.getCurrentButton() != StsMouse.LEFT)
        {
            return false;
        }
        int buttonState = mouse.getButtonStateCheckClear(StsMouse.LEFT);

        if (buttonState == StsMouse.CLEARED)
        {
            return true;
        }

        checkInitializePickMethods(glPanel3d);

        if (buttonState == StsMouse.PRESSED)
        {
            if (action != ACTION_MOVE_POINT)
            {
                StsWellSubVolume pickedBox = pickCenter(glPanel3d);
                if (pickedBox == null)
                {
                    currentSubVol = null;
                    action = ACTION_NONE;
                }
                else
                {
                    currentSubVol = pickedBox;
                    action = ACTION_MOVE_POINT;
                }
            }
            else
            {
                boolean facePicked = pickFace(glPanel3d);
                if (!facePicked)
                {
                    StsWellSubVolume pickedWellSv = pickCenter(glPanel3d);
                    if (pickedWellSv != null && currentSubVol != pickedWellSv)
                    {
                        currentSubVol = pickedWellSv;
                    }
                    else if(pickedWellSv == null)
                    {
                        currentSubVol = null;
                        action = ACTION_NONE;
                    }
                }
            }
            currentModel.win3dDisplay();
            return true;
        }
        else if (buttonState == StsMouse.DRAGGED)
        {
            if (currentSubVol == null)
            {
                return true;
            }
            if (action != ACTION_MOVE_POINT)
            {
                StsMessageFiles.errorMessage("StsWellSetSubVolume.action not set to ACTION_MOVE_POINT.");
                return true;
            }
            currentSubVol.movePoint(mouse, glPanel3d, pickedPointIndex);
            currentModel.win3dDisplay();
        }
        else if (buttonState == StsMouse.RELEASED)
        {
            action = ACTION_NONE;
        }
        return true;
    }

    public boolean delete(StsMouse mouse, StsGLPanel3d glPanel3d)
    {
        if (mouse.getCurrentButton() != StsMouse.LEFT) return false;
        int buttonState = mouse.getButtonStateCheckClear(StsMouse.LEFT);

        if (buttonState == StsMouse.CLEARED) return false;

        checkInitializePickMethods(glPanel3d);

        if (buttonState == StsMouse.RELEASED)
        {
            currentSubVol = pickCenter(glPanel3d);
            if (currentSubVol == null)
            {
                new StsMessage(currentModel.win3d, StsMessage.WARNING, "No well cylinder selected.  Try again.");
                return false;
            }
            else
            {
                wellSvs.delete(currentSubVol);
                currentSubVol.delete();
                currentModel.win3dDisplay();
                return true;
            }
        }
        else
            return false;
    }

    private StsWellSubVolume pickCenter(StsGLPanel3d glPanel3d)
    {
        try
        {
            if (!pickCurrentWellSetCenterPoint.methodPick3d()) return null;
            StsPickItem pickItem = pickCurrentWellSetCenterPoint.pickItems[0];
            int nWellSv = pickItem.names[0];
            return (StsWellSubVolume) wellSvs.getElement(nWellSv);
        }
        catch(Exception e) { return null; }
    }

    private boolean pickFace(StsGLPanel3d glPanel3d)
    {
        try
        {
            if(!pickCurrentWellSetFacePoint.methodPick3d()) return false;
            StsPickItem pickItem = pickCurrentWellSetFacePoint.pickItems[0];
            int nWellSv = pickItem.names[0];
            StsWellSubVolume pickedSubVol = (StsWellSubVolume) wellSvs.getElement(nWellSv);
            if (pickedSubVol != currentSubVol)
            {
                currentSubVol = pickedSubVol;
                return true;
            }
            pickedPointIndex = pickItem.names[1];
            currentModel.win3dDisplay();
            return true;
        }
        catch(Exception e) { return false; }
    }
    static public void pickCurrentWellSetCenterPoint(GL gl)
    {
        StsWellSetSubVolumeClass wellSetClass = (StsWellSetSubVolumeClass)StsObject.getCurrentModel().getStsClass(StsWellSetSubVolume.class);
        StsWellSetSubVolume currentWellSet = (StsWellSetSubVolume)wellSetClass.getCurrentObject();
        if(currentWellSet == null) return;
        currentWellSet.pickCenterPoint(gl);
    }

    public void pickCenterPoint(GL gl)
    {
        int nWellSvs = wellSvs.getSize();
        for (int n = 0; n < nWellSvs; n++)
        {
            StsWellSubVolume wellSv = (StsWellSubVolume) wellSvs.getElement(n);
            boolean editing = currentSubVol == wellSv;
            wellSv.pickCenterPoint(gl, editing, n);
        }
    }

    static public void pickCurrentWellSetFacePoint(GL gl)
    {
        StsWellSetSubVolumeClass wellSetClass = (StsWellSetSubVolumeClass)StsObject.getCurrentModel().getStsClass(StsWellSetSubVolume.class);
        StsWellSetSubVolume currentWellSet = (StsWellSetSubVolume)wellSetClass.getCurrentObject();
        if(currentWellSet == null) return;
        currentWellSet.pickFacePoint(gl);
    }
    
    public void pickFacePoint(GL gl)
    {
        int nWelSvs = wellSvs.getSize();
        for (int n = 0; n < nWelSvs; n++)
        {
            StsWellSubVolume wellSv = (StsWellSubVolume) wellSvs.getElement(n);
            boolean editing = currentSubVol == wellSv;
            wellSv.pickFacePoint(gl, n);
        }
    }
    public StsObjectRefList getWellSubVols() { return wellSvs; }

    public void setCylinderTop(float top) { if(currentSubVol != null) currentSubVol.setTop(top); currentModel.win3dDisplay(); }
    public void setCylinderBtm(float btm) { if(currentSubVol != null) currentSubVol.setBtm(btm); currentModel.win3dDisplay(); }
    public void setCylinderRadius(float radius) { if(currentSubVol != null) currentSubVol.setRadius(radius); currentModel.win3dDisplay(); }

    public float getCylinderTop() { if(currentSubVol == null) return StsParameters.nullValue; else return currentSubVol.getTop(); }
    public float getCylinderBtm() { if(currentSubVol == null) return StsParameters.nullValue; else return currentSubVol.getBtm(); }
    public float getCylinderRadius() { if(currentSubVol == null) return StsParameters.nullValue; else return currentSubVol.getRadius(); }
}
