package com.Sts.DBTypes;

/**
 * <p>Title: S2S Development</p>
 * <p>Description: </p>
 * <p>Copyright: Copyright (c) 2002</p>
 * <p>Company: S2S Systems LLC</p>
 * @author TJLasseter
 * @version beta 1.0
 */

import com.Sts.Interfaces.*;
import com.Sts.MVC.*;
import com.Sts.MVC.View3d.*;
import com.Sts.Types.*;
import com.Sts.UI.Beans.*;
import com.Sts.UI.ObjectPanel.*;
import com.Sts.Utilities.*;

public class StsCropVolume extends StsRotatedGridBoundingBox implements StsTreeObjectI, Cloneable
{
    public boolean applyCrop = false;

    transient StsRotatedGridBoundingBox initialBoundingBox;
    transient boolean valueChanged = false;

    /** List of display beans for StsProject */
    static public StsBooleanFieldBean applyCropBean;
    static public StsFieldBean[] displayFields = null;

    static public StsFloatFieldBean cropXMinBean;
    static public StsFloatFieldBean cropXMaxBean;
    static public StsFloatFieldBean cropYMinBean;
    static public StsFloatFieldBean cropYMaxBean;
    static public StsFloatFieldBean cropZMinBean;
    static public StsFloatFieldBean cropZMaxBean;
    static public StsFloatFieldBean cropRowNumMinBean;
    static public StsFloatFieldBean cropRowNumMaxBean;
    static public StsFloatFieldBean cropColNumMinBean;
    static public StsFloatFieldBean cropColNumMaxBean;
    static public StsFieldBean[] propertyFields;

    static StsObjectPanel objectPanel = null;

    static final long serialVersionUID = Main.getTime(103, 7, 3, 0, 0, 0);

	public StsCropVolume()
	{
	}

	public StsCropVolume(boolean persistent)
	{
		super(persistent);
        initialize();
    }

    public boolean initialize(StsModel model)
    {
        objectPanel = null;
        initialBoundingBox = new StsRotatedGridBoundingBox(false);
        setInitialToProject();
        StsFieldBean[] fieldBeans = getPropertyFields();
        StsFieldBean.setBeanObject(fieldBeans, this);
        return true;
    }

    public void initialize()
    {
        objectPanel = null;
        initialBoundingBox = new StsRotatedGridBoundingBox(false);
        setInitialToProject();
        setCropToInitial();
        StsFieldBean[] fieldBeans = getPropertyFields();
        StsFieldBean.setBeanObject(fieldBeans, this);
    }

    public void setInitialToProject()
    {
        if(!StsToolkit.copySubToSuperclass(getProjectBoundingBox(), initialBoundingBox, StsBoundingBox.class)) return;
//        setCropToInitial();
    }

    private StsRotatedGridBoundingBox getProjectBoundingBox()
    {
        if(currentModel == null || currentModel.getProject() == null) return null;
        return currentModel.getProject().getRotatedBoundingBox();
    }

    private void setCropToInitial()
    {
        StsToolkit.copySubToSuperclass(initialBoundingBox, this, StsBoundingBox.class);
    }

    public void setMinMaxRanges()
    {
        if(xMin != largeFloat) cropXMinBean.setValueAndRange(xMin, initialBoundingBox.xMin, initialBoundingBox.xMax);
        if(xMax != -largeFloat) cropXMaxBean.setValueAndRange(xMax, initialBoundingBox.xMin, initialBoundingBox.xMax);

        if(yMin != largeFloat) cropYMinBean.setValueAndRange(yMin, initialBoundingBox.yMin, initialBoundingBox.yMax);
        if(yMax != -largeFloat) cropYMaxBean.setValueAndRange(yMax, initialBoundingBox.yMin, initialBoundingBox.yMax);

        if(zMin != largeFloat) cropZMinBean.setValueAndRange(zMin, initialBoundingBox.zMin, initialBoundingBox.zMax);
        if(zMax != -largeFloat) cropZMaxBean.setValueAndRange(zMax, initialBoundingBox.zMin, initialBoundingBox.zMax);

		if(rowNumMin != nullValue && rowNumMax != nullValue)
		{
			if (initialBoundingBox.rowNumInc > 0.0f)
			{
				cropRowNumMinBean.setValueAndRange(rowNumMin, initialBoundingBox.rowNumMin, initialBoundingBox.rowNumMax);
				cropRowNumMaxBean.setValueAndRange(rowNumMax, initialBoundingBox.rowNumMin, initialBoundingBox.rowNumMax);
			}
			else
			{
				cropRowNumMinBean.setValueAndRange(rowNumMin, initialBoundingBox.rowNumMax, initialBoundingBox.rowNumMin);
				cropRowNumMaxBean.setValueAndRange(rowNumMax, initialBoundingBox.rowNumMax, initialBoundingBox.rowNumMin);
			}
		}
		if(colNumMin != nullValue && colNumMax != nullValue)
		{
			if (initialBoundingBox.colNumInc > 0.0f)
			{
				cropColNumMinBean.setValueAndRange(colNumMin, initialBoundingBox.colNumMin, initialBoundingBox.colNumMax);
				cropColNumMaxBean.setValueAndRange(colNumMax, initialBoundingBox.colNumMin, initialBoundingBox.colNumMax);
			}
			else
			{
				cropColNumMinBean.setValueAndRange(colNumMin, initialBoundingBox.colNumMax, initialBoundingBox.colNumMin);
				cropColNumMaxBean.setValueAndRange(colNumMax, initialBoundingBox.colNumMax, initialBoundingBox.colNumMin);
			}
		}
    }

    public float getCropXMin() { return xMin; }
    public float getCropXMax() { return xMax; }
    public float getCropYMin() { return yMin; }
    public float getCropYMax() { return yMax; }
    public float getCropZMin() { return zMin; }
    public float getCropZMax() { return zMax; }

    public int getCropRowMin() { return initialBoundingBox.getNearestBoundedRowCoor(yMin); }
    public int getCropRowMax() { return initialBoundingBox.getNearestBoundedRowCoor(yMax); }
    public int getCropColMin() { return initialBoundingBox.getNearestBoundedColCoor(xMin); }
    public int getCropColMax() { return initialBoundingBox.getNearestBoundedColCoor(xMax); }
    public int getCropSliceMin() { return initialBoundingBox.getNearestSliceCoor(zMin); }
    public int getCropSliceMax() { return initialBoundingBox.getNearestSliceCoor(zMax); }

    public int getNCropRows() { return getCropRowMax() - getCropRowMin() + 1; }
    public int getNCropCols() { return getCropColMax() - getCropColMin() + 1; }
    public int getNCropSlices() { return getCropSliceMax() - getCropSliceMin() + 1; }

    public float getCropRowNumMin() { return (float)rowNumMin; }
    public float getCropRowNumMax() { return (float)rowNumMax; }
    public float getCropColNumMin() { return (float)colNumMin; }
    public float getCropColNumMax() { return (float)colNumMax; }

    public void setCropRowNumMin(float cropRowNumMin)
    {
        if(this.rowNumMin == cropRowNumMin) return;
        setRowNumMin(cropRowNumMin);
        setYMin(initialBoundingBox.getYFromRowNum(cropRowNumMin));
//        setRowMin(Math.round(projectBoundingBox.getRowCoor(yMin)));
        cropYMinBean.setValue(yMin);
        valueChanged();
    }

    public void setRowNumMin(float rowNumMin)
    {
        dbFieldChanged("rowNumMin", rowNumMin);
        super.setRowNumMin(rowNumMin);
    }

    public void setRowNumMax(float rowNumMax)
    {
        dbFieldChanged("rowNumMax", rowNumMax);
        super.setRowNumMax(rowNumMax);
    }

    public void setColNumMin(float colNumMin)
    {
        dbFieldChanged("colNumMin", colNumMin);
        super.setColNumMin(colNumMin);
    }

    public void setColNumMax(float colNumMax)
    {
        dbFieldChanged("colNumMax", colNumMax);
        super.setColNumMax(colNumMax);
    }

    public void setXMin(float xMin)
    {
        dbFieldChanged("xMin", xMin);
        super.setXMin(xMin);
    }

    public void setXMax(float xMax)
    {
        dbFieldChanged("xMax", xMax);
        super.setXMax(xMax);
    }

    public void setYMin(float yMin)
    {
        dbFieldChanged("yMin", yMin);
        super.setYMin(yMin);
    }

    public void setYMax(float yMax)
    {
        dbFieldChanged("yMax", yMax);
        super.setYMax(yMax);
    }

    public void setZMin(float zMin)
    {
        dbFieldChanged("zMin", zMin);
        super.setZMin(zMin);
    }

    public void setZMax(float zMax)
    {
        dbFieldChanged("zMax", zMax);
        super.setZMax(zMax);
    }

    public void setCropRowNumMax(float cropRowNumMax)
    {
        if(this.rowNumMax == cropRowNumMax) return;
        setRowNumMax(cropRowNumMax);
        setYMax(initialBoundingBox.getYFromRowNum(cropRowNumMax));
//        setRowMax(Math.round(projectBoundingBox.getRowCoor(yMax)));
        cropYMaxBean.setValue(yMax);
        valueChanged();
    }

    public void setCropColNumMin(float cropColNumMin)
    {
        if(this.colNumMin == cropColNumMin) return;
        setColNumMin(cropColNumMin);
        setXMin(initialBoundingBox.getXFromColNum(cropColNumMin));
//        setColMin(Math.round(projectBoundingBox.getColCoor(xMin)));
        cropXMinBean.setValue(xMin);
        valueChanged();
    }

    public void setCropColNumMax(float cropColNumMax)
    {
        if(this.colNumMax == cropColNumMax) return;
        setColNumMax(cropColNumMax);
        setXMax(initialBoundingBox.getXFromColNum(cropColNumMax));
//        setColMax(Math.round(projectBoundingBox.getColCoor(xMax)));
        cropXMaxBean.setValue(xMax);
        valueChanged();
    }

    public void setCropXMin(float cropXMin)
    {
        if(this.xMin == cropXMin) return;
        setXMin(cropXMin);
//        setColMin(Math.round(projectBoundingBox.getColCoor(cropXMin)));
        setColNumMin(initialBoundingBox.getNearestBoundedColNumFromX(cropXMin));
        cropColNumMinBean.setValue(colNumMin);
        valueChanged();
    }

    public void setCropXMax(float cropXMax)
    {
        if(this.xMax == cropXMax) return;
        setXMax(cropXMax);
//        setColMax(Math.round(projectBoundingBox.getColCoor(cropXMax)));
        setColNumMax(initialBoundingBox.getNearestBoundedColNumFromX(cropXMax));
        cropColNumMaxBean.setValue(colNumMax);
        valueChanged();
    }

    public void setCropYMin(float cropYMin)
    {
        if(this.yMin == cropYMin) return;
        setYMin(cropYMin);
//        setRowMin(Math.round(projectBoundingBox.getRowCoor(cropYMin)));
        setRowNumMin(initialBoundingBox.getNearestBoundedRowNumFromY(cropYMin));
        cropRowNumMinBean.setValue(rowNumMin);
        valueChanged();
    }

    public void setCropYMax(float cropYMax)
    {
        if(this.yMax == cropYMax) return;
        setYMax(cropYMax);
//        setRowMax(Math.round(projectBoundingBox.getRowCoor(cropYMax)));
        setRowNumMax(initialBoundingBox.getNearestBoundedRowNumFromY(cropYMax));
        cropRowNumMaxBean.setValue(rowNumMax);
        valueChanged();
    }

    public void setCropZMin(float cropZMin)
    {
        setZMin(cropZMin);
//        setSliceMin(Math.round(projectBoundingBox.getSliceCoor(cropZMin)));
        valueChanged();
    }

    public void setCropZMax(float cropZMax)
    {
        setZMax(cropZMax);
//        setSliceMax(Math.round(projectBoundingBox.getSliceCoor(cropZMax)));
        valueChanged();
    }

    public int getCursorRowMin(int dir)
    {
        if     (dir == XDIR) return initialBoundingBox.getNearestBoundedRowCoor(yMin);
        else if(dir == YDIR) return initialBoundingBox.getNearestBoundedColCoor(xMin);
        else if(dir == ZDIR) return initialBoundingBox.getNearestBoundedRowCoor(yMin);
        else                 return 0;
    }

    public int getCursorColMin(int dir)
    {
        if     (dir == XDIR) return initialBoundingBox.getNearestSliceCoor(zMin);
        else if(dir == YDIR) return initialBoundingBox.getNearestSliceCoor(zMin);
        else if(dir == ZDIR) return initialBoundingBox.getNearestBoundedColCoor(xMin);
        else                 return 0;
    }

    public int getCursorRowMax(int dir)
    {
        if     (dir == XDIR) return initialBoundingBox.getNearestBoundedRowCoor(yMax);
        else if(dir == YDIR) return initialBoundingBox.getNearestBoundedColCoor(xMax);
        else if(dir == ZDIR) return initialBoundingBox.getNearestBoundedRowCoor(yMax);
        else                 return 0;
    }

    public int getCursorColMax(int dir)
    {
        if     (dir == XDIR) return initialBoundingBox.getNearestSliceCoor(zMax);
        else if(dir == YDIR) return initialBoundingBox.getNearestSliceCoor(zMax);
        else if(dir == ZDIR) return initialBoundingBox.getNearestBoundedColCoor(xMax);
        else                 return 0;
    }

    public StsFieldBean[] getDisplayFields()
    {
        if(displayFields == null)
        {
            applyCropBean = new StsBooleanFieldBean(StsCropVolume.class, "applyCrop", "Crop Volumes:");
            displayFields = new StsFieldBean[]
            {
                applyCropBean
            };
        }
        return displayFields;
    }

    public StsFieldBean[] getPropertyFields()
    {
        if(propertyFields == null)
        {
            cropXMinBean = new StsFloatFieldBean(StsCropVolume.class, "cropXMin", "Cropped Min X");
            cropXMaxBean = new StsFloatFieldBean(StsCropVolume.class, "cropXMax", "Cropped Max X");
            cropYMinBean = new StsFloatFieldBean(StsCropVolume.class, "cropYMin", "Cropped Min Y");
            cropYMaxBean = new StsFloatFieldBean(StsCropVolume.class, "cropYMax", "Cropped Max Y");
            cropZMinBean = new StsFloatFieldBean(StsCropVolume.class, "cropZMin", "Cropped Min Z or T");
            cropZMaxBean = new StsFloatFieldBean(StsCropVolume.class, "cropZMax", "Cropped Max Z or T");
            cropRowNumMinBean = new StsFloatFieldBean(StsCropVolume.class, "cropRowNumMin", "Min Cropped Row/Line");
            cropRowNumMaxBean = new StsFloatFieldBean(StsCropVolume.class, "cropRowNumMax", "Max Cropped Row/Line");
            cropColNumMinBean = new StsFloatFieldBean(StsCropVolume.class, "cropColNumMin", "Min Cropped Col/XLine");
            cropColNumMaxBean = new StsFloatFieldBean(StsCropVolume.class, "cropColNumMax", "Max Cropped Col/XLine");
            propertyFields = new StsFieldBean[] { cropXMinBean, cropXMaxBean, cropYMinBean, cropYMaxBean,
                                                  cropZMinBean, cropZMaxBean, cropRowNumMinBean, cropRowNumMaxBean,
                                                  cropColNumMinBean, cropColNumMaxBean };
            setMinMaxRanges();
        }
        return propertyFields;
    }

    public Object[] getChildren() { return new Object[0]; }
    public StsObjectPanel getObjectPanel()
    {
        if(objectPanel == null) objectPanel = StsObjectPanel.constructor(this, true);
        return objectPanel;
    }
    public void treeObjectSelected() { }

    public boolean anyDependencies()
    {
        return true;
    }

    public boolean isDirCoordinateCropped(int dir, float coor)
    {
        switch (dir)
        {
            case StsCursor3d.XDIR:
                return coor < xMin || coor > xMax;
            case StsCursor3d.YDIR:
                return coor < yMin || coor > yMax;
            case StsCursor3d.ZDIR:
                return coor < zMin || coor > zMax;
            default:
                return false;
        }
    }

    public void valueChanged()
    {
        valueChanged = true;
        checkApplyCrop();
 //       currentModel.cropChanged();
 //       currentModel.win3dDisplayAll();
//        redisplay();
    }

    private void checkApplyCrop()
    {
       if(applyCrop) return;
       applyCropBean.setValue(true);
       setApplyCrop(true);
    }

    public boolean isCropped()
    {
        if(!applyCrop) return false;
        if(xMin > initialBoundingBox.xMin) return true;
        if(xMax < initialBoundingBox.xMax) return true;
        if(yMin > initialBoundingBox.yMin) return true;
        if(yMax < initialBoundingBox.yMax) return true;
        if(zMin > initialBoundingBox.zMin) return true;
        if(zMax < initialBoundingBox.zMax) return true;
        return false;
    }

	public boolean isZCropped()
	{
		if(!applyCrop) return false;
		if(zMin > initialBoundingBox.zMin) return true;
		if(zMax < initialBoundingBox.zMax) return true;
		return false;
	}

    public boolean getValueChanged() { return valueChanged; }


    /** @param crop indicates whether we want cropping isVisible */
    public void setApplyCrop(boolean crop)
    {
        if(crop == applyCrop) return;
        applyCrop = crop;
        dbFieldChanged("applyCrop", applyCrop);
        currentModel.cropChanged();
        currentModel.win3dDisplayAll();
    }

    public boolean getApplyCrop() { return applyCrop; }
}
