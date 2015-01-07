package com.Sts.DBTypes;

/**
 * <p>Title: S2S Development</p>
 * <p>Description: </p>
 * <p>Copyright: Copyright (c) 2002</p>
 * <p>Company: S2S Systems LLC</p>
 * @author unascribed
 * @version 1.1
 */

import com.Sts.DB.*;
import com.Sts.Interfaces.*;
import com.Sts.MVC.*;
import com.Sts.MVC.View3d.*;
import com.Sts.UI.Beans.*;
import com.Sts.UI.*;
import com.Sts.UI.Toolbars.*;

import java.awt.*;
import java.beans.*;

public class StsCrossplotClass extends StsClass implements StsSerializable, StsClassCursor3dTextureDisplayable
{
    protected boolean verticalLinesOnCrossplot = false;
    protected boolean horizontalLinesOnCrossplot = false;

    protected boolean isVisibleOnCursor = true;
    protected boolean displayOnSubVolumes = true;
    protected boolean displayEntireVolume = false;
    protected boolean displayGridLinesOnCrossplot = true;
    protected boolean displayCrossplotAxis = true;
    protected boolean isPixelMode = true;

    private String xplotSpectrumName = StsSpectrumClass.SPECTRUM_RWB;

    public StsCrossplotClass()
    {
        userName = "Volume Crossplot";
        addPropertyChangeListener();
    }

    private void addPropertyChangeListener()
    {
        StsSubVolumeClass subVolumeClass = (StsSubVolumeClass)currentModel.getCreateStsClass(StsSubVolume.class);
        subVolumeClass.addPropertyChangeListener
        (
            new PropertyChangeListener()
            {
                public void propertyChange(PropertyChangeEvent e)
                {
                    propertyChanged(e);
                }
            }
        );
    }

    private void propertyChanged(PropertyChangeEvent e)
    {
        String propertyName = e.getPropertyName();
        if(propertyName.equals("StsSubVolumeClass.isApplied"))
        {
             forEach("subVolumeVisibilityChanged", e.getNewValue());
        }
    }

    public void initializeDisplayFields()
    {
        displayFields = new StsFieldBean[]
        {
            new StsBooleanFieldBean(this, "isVisibleOnCursor", "Visible on Cursor"),
            new StsBooleanFieldBean(this, "displayOnSubVolumes", "Filter by SubVolumes"),
            new StsBooleanFieldBean(this, "displayEntireVolume", "Crossplot Entire PostStack3d"),
            new StsBooleanFieldBean(this, "displayGridLinesOnCrossplot", "Grid Lines:"),
            new StsBooleanFieldBean(this, "displayCrossplotAxis", "Plot Axis:"),
            new StsBooleanFieldBean(this, "isPixelMode", "Display Pixel Mode:"),
        };
    }

    public void initializeDefaultFields()
    {
       defaultFields = new StsFieldBean[]
       {
            new StsComboBoxFieldBean(this, "xplotSpectrumName", "Spectrum:", StsSpectrumClass.cannedSpectrums),
       };

    }
    
    public boolean setCurrentObject(StsObject object)
    {
        if(currentModel.win3d == null) 
        	return super.setCurrentObject(object);
        super.setCurrentObject(object);
        StsCrossplotToolbar crossplotToolbar = (StsCrossplotToolbar)currentModel.win3d.getToolbarNamed(StsCrossplotToolbar.NAME);
        if(crossplotToolbar != null)
        	crossplotToolbar.enableDensityToggle(((StsCrossplot)object).hasAttributeData());        
        return true;
    }
/*    
    public void selected(StsCrossplot crossplot)
    {
        super.selected(crossplot);
        setCurrentObject(crossplot);    
    }
*/
    public boolean getIsVisibleOnCursor() { return isVisible && isVisibleOnCursor; }
    public void setIsVisibleOnCursor(boolean isVisible)
    {
        if(this.isVisibleOnCursor == isVisible) return;
        this.isVisibleOnCursor = isVisible;
//        setDisplayField("isVisibleOnCursor", isVisibleOnCursor);
		currentModel.win3dDisplayAll();
    }

    public boolean getDisplayGridLinesOnCrossplot() { return displayGridLinesOnCrossplot; }
    public void setDisplayGridLinesOnCrossplot(boolean b)
    {
        if(this.displayGridLinesOnCrossplot == b) return;
        this.displayGridLinesOnCrossplot = b;
//        setDisplayField("displayGridLinesOnCrossplot", displayGridLinesOnCrossplot);
		currentModel.win3dDisplayAll();
    }

    public boolean getDisplayCrossplotAxis() { return displayCrossplotAxis; }
    public void setDisplayCrossplotAxis(boolean b)
    {
        if(this.displayCrossplotAxis == b) return;
        this.displayCrossplotAxis = b;
//        setDisplayField("displayCrossplotAxis", displayCrossplotAxis);
		currentModel.win3dDisplayAll();
    }

    public boolean getDisplayOnSubVolumes() { return displayOnSubVolumes; }
    public void setDisplayOnSubVolumes(boolean b)
    {
        if(this.displayOnSubVolumes == b) return;
        this.displayOnSubVolumes = b;
//        setDisplayField("displayOnSubVolumes", displayOnSubVolumes);
        currentModel.subVolumeChanged();
        currentModel.win3dDisplayAll();
    }

    public boolean getDisplayEntireVolume() { return displayEntireVolume; }
    public void setDisplayEntireVolume(boolean b)
    {
        if(this.displayEntireVolume == b) return;
        this.displayEntireVolume = b;
//        setDisplayField("displayEntireVolume", displayEntireVolume);
        StsCursor cursor = new StsCursor(currentModel.win3d, Cursor.WAIT_CURSOR);
        currentModel.clearTextureClassDisplays(StsCrossplot.class);
        currentModel.win3dDisplayAll();
        cursor.restoreCursor();
    }

/*
    public boolean getCropMode() { return cropMode; }
    public void setCropMode(boolean mode)
    {
        cropMode = mode;
        int size = list.getSize();
        for (int i = 0; i < size; i++)
        {
            StsCrossplot crossplot = (StsCrossplot)list.getElement(i);
            crossplot.clearColors();
        }
        currentModel.win3dDisplayAll();
    }
*/
    public boolean getIsPixelMode() { return isPixelMode; }
    public void setIsPixelMode(boolean mode)
    {
        if(this.isPixelMode == mode) return;
        this.isPixelMode = mode;
//        setDisplayField("isPixelMode", isPixelMode);

        int size = list.getSize();
        for (int i = 0; i < size; i++)
        {
            StsCrossplot crossplot = (StsCrossplot)list.getElement(i);
            crossplot.clearColors();
        }
        currentModel.win3dDisplayAll();
    /*
        for(int i=0; i<currentModel.getNumberFamilies(); i++)
        {
            for (int n = 0; n < currentModel.getWindows(i).size(); n++)
            {
                StsView view = (StsView) ((StsWin3dBase)currentModel.getWindows(i).get(n)).glPanel3d.getCurrentView();
                if (view instanceof StsViewXP)
                {
                    StsViewXP viewXP = (StsViewXP)view;
                    viewXP.checkClearColors();
                }
            }
        }
        currentModel.win3dDisplayAll();
    */
    }

    public void previousCrossplot()
    {
        StsCrossplot cp;

        int size = list.getSize();
        if (size < 1) return;
        String currentName = currentObject.getName();

        cp = (StsCrossplot)list.getElement(0);
        if(currentName.equals(cp.getName()))
        {
            setCurrentObject((StsCrossplot)list.getElement(size - 1));
//            setCurrentCrossplot((StsCrossplot)list.getElement(size - 1));
            return;
        }
        for (int i = 1; i < size; i++)
        {
            cp = (StsCrossplot)list.getElement(i);
            if(currentName.equals(cp.getName()))
                setCurrentObject((StsCrossplot)list.getElement(i - 1));
//                setCurrentCrossplot((StsCrossplot)list.getElement(i - 1));
        }
    }

    public void nextCrossplot()
    {
        StsCrossplot cp;

        int size = list.getSize();
        if (size < 1) return;
        String currentName = currentObject.getName();

        cp = (StsCrossplot)list.getElement(size-1);
        if(currentName.equals(cp.getName()))
        {
            setCurrentObject((StsCrossplot)list.getElement(0));
//            setCurrentCrossplot((StsCrossplot)list.getElement(0));
            return;
        }
        for (int i = 0; i < size-1; i++)
        {
            cp = (StsCrossplot)list.getElement(i);
            if(currentName.equals(cp.getName()))
                setCurrentObject((StsCrossplot)list.getElement(i + 1));
//                setCurrentCrossplot((StsCrossplot)list.getElement(i + 1));
        }
    }

    public String getCurrentCrossplotName()
    {
        if(currentObject == null) return new String("");
        else return currentObject.getName();
    }

    public void setCurrentCrossplotName(String name)
    {
        StsCrossplot newCrossplot = (StsCrossplot)getObjectWithName(name);
        setCurrentObject(newCrossplot);
//        setCurrentCrossplot(newCrossplot);
    }

    public StsCrossplot getCurrentCrossplot()
    {
    /*
        if(currentCrossplot == null) return null;
        if(!currentCrossplot.indexOK())
        {
            currentCrossplot = (StsCrossplot)getElement(0);
            if(currentCrossplot != null) currentCrossplot.processPolygons();
        }
    */
        return (StsCrossplot)currentObject;
    }
/*
    public void setCurrentCrossplot(StsCrossplot crossplot)
    {
        if(currentObject == crossplot) return;
        currentObject = crossplot;
        StsClassChangeCmd cmd = new StsClassChangeCmd(this, currentObject, "currentObject", false);
        currentModel.addTransactionCmd("setCurrentCrossplot", cmd);
    }
*/
    public boolean getVerticalLinesOnCrossplot() { return verticalLinesOnCrossplot; }
    public boolean getHorizontalLinesOnCrossplot() { return horizontalLinesOnCrossplot; }

    public void setVerticalLinesOnCrossplot(boolean b){ verticalLinesOnCrossplot = b; }
    public void setHorizontalLinesOnCrossplot(boolean b) { horizontalLinesOnCrossplot = b; }

    public StsCursor3dTexture constructDisplayableSection(StsModel model, StsCursor3d cursor3d, int dir)
    {
        return new StsCrossplotCursorSection(model, (StsCrossplot)currentObject, cursor3d, dir);
    }

    public boolean drawLast() { return true; }
/*
    public void displayOnCursor(StsCursor3d.StsCursor3dTexture cursorSection, StsGLPanel3d glPanel3d, int dirNo, float dirCoordinate, boolean is3d)
      {
          if(currentObject == null || !((StsCrossplot)currentObject).getIsVisibleOnCursor()) return;
          StsCrossplotCursorSection crossplotCursorSection = cursorSection.hasCursor3dDisplayable(StsCrossplotCursorSection.class);
          if(crossplotCursorSection == null)
          {
              crossplotCursorSection = new StsCrossplotCursorSection((StsCrossplot)currentObject, glPanel3d, dirNo, dirCoordinate);
              cursorSection.addDisplayable(crossplotCursorSection);
          }
          crossplotCursorSection.display(this, glPanel3d, is3d);
      }
  */
    public String getXplotSpectrumName()  { return xplotSpectrumName; }
    public void setXplotSpectrumName(String value)
    {
        if(this.xplotSpectrumName.equals(value)) return;
        this.xplotSpectrumName = value;
//        setDisplayField("xplotSpectrumName", xplotSpectrumName);
    }
}
