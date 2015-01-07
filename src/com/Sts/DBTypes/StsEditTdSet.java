package com.Sts.DBTypes;

/**
 * <p>Title: jS2S development</p>
 * <p>Description: </p>
 * <p>Copyright: Copyright (c) 2004</p>
 * <p>Company: S2S Systems LLC</p>
 * @author Tom Lasseter
 * @version 1.0
 */

import com.Sts.Interfaces.*;
import com.Sts.MVC.*;
import com.Sts.MVC.View3d.*;
import com.Sts.Types.*;
import com.Sts.UI.Beans.*;
import com.Sts.UI.ObjectPanel.*;
import com.Sts.Utilities.*;

import javax.media.opengl.*;

public class StsEditTdSet extends StsMainObject implements StsTreeObjectI, StsInstance3dDisplayable
{
    StsWell well;
    StsObjectRefList tdEdits;
    transient StsEditTd tdEdit = null;

    /** a curtain of vertical seismic traces passing down thru this line path */
    transient protected StsSeismicCurtain seismicCurtain = null;

    static StsObjectPanel objectPanel = null;

    static StsComboBoxFieldBean tdEditsListBean;
    static public StsFieldBean[] displayFields = null;

    public StsEditTdSet()
    {
    }

    /** We have added a tdEdit point to this well.  Created a "no edit" tdEdit
     *  curve and our first tdEdit curve with just a tdEdit point at the top.
     */
    private StsEditTdSet(StsWell well)
    {
        this.well = well;
        setName(well.getName());
        StsEditTd noTdEdit = new StsEditTd(well, "no edits");
        StsEditTd firstTdEdit = new StsEditTd(well);
        tdEdits = StsObjectRefList.constructor(2, 2, "tdEdits", this);
        tdEdits.add(noTdEdit);
        tdEdits.add(firstTdEdit);
        Object[] tdEditsList = getTdEdits();
//        setLogCurve(firstTdEdit);
        StsPoint topPoint = well.getTopPoint();
        firstTdEdit.addPoint(topPoint);
        firstTdEdit.resetTimes();
        this.tdEdit = firstTdEdit;
        if(well.wellViewModel != null) well.wellViewModel.display();
        if(tdEditsListBean == null) return;
        tdEditsListBean.setListItems(tdEditsList);
        tdEditsListBean.setSelectedItem(tdEdit);
    }


    static public StsEditTdSet getCreateEditTdSet(StsWell well)
    {
        try
        {
            StsEditTdSet tdEditSet = (StsEditTdSet) currentModel.getObjectWithName(StsEditTdSet.class, well.getName());
            if (tdEditSet != null)return tdEditSet;
            return new StsEditTdSet(well);
        }
        catch (Exception e)
        {
            StsException.outputWarningException(StsEditTdSet.class, "getCreateEditTdSet", e);
            return null;
        }
    }

    public boolean initialize(StsModel model)
    {
        return true;
    }

    public void createSeismicCurtain(StsWell well)
    {
        isVisible = true;
        well.isDrawingCurtain = true;
        StsSeismicVolume seismicVolume = (StsSeismicVolume) currentModel.getCurrentObject(StsSeismicVolume.class);
        StsPoint[] rotatedPoints = well.getRotatedPoints();
		if (seismicVolume != null)
        seismicCurtain = new StsSeismicCurtain(currentModel, rotatedPoints, seismicVolume);
    }

    public boolean delete()
    {
        if (seismicCurtain != null)
        {
            deleteSeismicCurtain();
         }
        if (!super.delete())
        {
            return false;
        }
        StsObjectRefList.deleteAll(tdEdits);
        return true;
    }

    public void deleteSeismicCurtain()
    {
        if (seismicCurtain != null)
        {
            well.isDrawingCurtain = false;
            seismicCurtain.textureChanged();
            currentModel.win3dDisplay();
        }
        seismicCurtain = null;
    }

    public void setTdEdits(Object[] tdEditObjects)
    {
        if(tdEdits == null) return;

        if(tdEdits != null) tdEdits.deleteAll();
        else tdEdits = new StsObjectRefList();

        for(int n = 0; n < tdEditObjects.length; n++)
            tdEdits.add((StsEditTd)tdEditObjects[n]);
    }

    public Object[] getTdEdits()
    {
        return tdEdits.getList().getTrimmedList();
    }

    public StsEditTd createTdEditFromCurrent()
    {
        if(tdEdit == null) return null;
        StsEditTd newTdEdit = new StsEditTd(tdEdit);
        tdEdits.add(newTdEdit);
        return newTdEdit;
    }

    public StsFieldBean[] getDisplayFields()
    {
        if(displayFields == null)
        {
            tdEditsListBean = new StsComboBoxFieldBean(StsEditTdSet.class, "tdEdit", "TD Edits");
            displayFields = new StsFieldBean[]
            {
                new StsBooleanFieldBean(StsEditTdSet.class, "isVisible", "Well Seismic"),
                tdEditsListBean
            };
        }
        tdEditsListBean.setListItems(getTdEdits());
        tdEditsListBean.setValueObject(tdEdit);
        return displayFields;
    }

    public StsFieldBean[] getPropertyFields()
    {
        return null;
    }

    public Object[] getChildren()
    {
        return new Object[0];
    }

    public boolean anyDependencies()
    {
        return false;
    }

    public StsObjectPanel getObjectPanel()
    {
        if (objectPanel == null)
        {
            objectPanel = StsObjectPanel.constructor(this, true);
        }
        return objectPanel;
    }

    public void treeObjectSelected()
    {
        currentModel.getStsClass(StsEditTdSet.class).selected(this);
    }

    public StsPoint getTDPointOnWell(StsGLPanel glPanel)
    {
        try
        {
            StsEditTdSet pickedTdEdits = (StsEditTdSet)StsJOGLPick.pickClass3d(glPanel, new StsObject[] { this }, StsJOGLPick.PICKSIZE_MEDIUM, StsJOGLPick.PICK_CLOSEST);
            if(pickedTdEdits == null) return null;
            StsPickItem[] pickItems = StsJOGLPick.pickItems;
            int names[] = pickItems[0].names;
            int pointIndex = names[2];
            return tdEdit.tdPoints[pointIndex];
        }
        catch(Exception e)
        {
            return null;
        }
    }
    public void pick(GL gl, StsGLPanel glPanel)
    {
        if(tdEdit == null) return;
        gl.glPushName(getIndex());
        tdEdit.pick(gl, glPanel);
        gl.glPopName();
    }

    public void display(StsGLPanel3d glPanel3d)
    {
        if(tdEdit == null) return;
        if(!isVisible) return;
        if (seismicCurtain != null)
        {
            seismicCurtain.displayTexture3d(glPanel3d);
            seismicCurtain.display(glPanel3d);
        }
        tdEdit.display(glPanel3d);
    }

    public Object getTdEdit() { return tdEdit; }

    public void setTdEdit(Object tdEdit)
    {
        if(this.tdEdit == tdEdit) return;
//        tdEditsListBean.setValueObject(tdEdit);
//        tdEditsListBean.getComboBox().getModel().setSelectedItem(tdEdit);
        this.tdEdit = (StsEditTd)tdEdit;
        this.tdEdit.resetTimes();
        if(well.wellViewModel != null) well.wellViewModel.display();
        if(tdEditsListBean != null)
            tdEditsListBean.setSelectedItem(tdEdit);
    }

    public void addPoint(StsPoint point)
    {
        if(tdEdit == null) return;
        tdEdit.addPoint(point);
    }

    public void movePoint(StsPoint selectedPoint, StsPoint pickedPoint)
    {
        if(tdEdit == null) return;
        tdEdit.adjustWellPath(well, selectedPoint, pickedPoint, well.wellViewModel);
    }

    public void endTransaction()
    {
        tdEdits.forEach("endTransaction");
    }

    public StsPoint setSelectedPoint(int pointIndex)
    {
        return tdEdit.setSelectedPoint(pointIndex);
    }
}
