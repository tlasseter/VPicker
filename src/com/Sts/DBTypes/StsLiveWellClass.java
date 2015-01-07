package com.Sts.DBTypes;

import com.Sts.DB.StsSerializable;
import com.Sts.Interfaces.*;
import com.Sts.MVC.View3d.StsGLPanel3d;
import com.Sts.UI.Beans.*;
import com.Sts.Utilities.StsSumTimer;

import java.util.Iterator;

/**
 * <p>Title: S2S Development</p>
 * <p>Description: </p>
 * <p>Copyright: Copyright (c) 2002</p>
 * <p>Company: S2S Systems LLC</p>
 * @author unascribed
 * @version 1.1
 */

public class StsLiveWellClass extends StsWellClass implements StsSerializable, StsClassTimeDisplayable, StsClassCursorDisplayable, StsClassTimeSeriesDisplayable, StsClassViewSelectable
{
    public boolean realtimeTrack = true;
    public StsLiveWellClass()
    {
        super();
        userName = "Time-Dependent Wells";
    }

    public void initializeDisplayFields()
    {
        StsFloatFieldBean perfScaleBean = new StsFloatFieldBean(this,"perfScale", 2.0f, 50.0f, "Perforation Scale:", true);
        perfScaleBean.setRangeFixStep(2.0f, 50.0f, 1.0f);
        StsFloatFieldBean equipScaleBean = new StsFloatFieldBean(this,"equipmentScale", 2.0f, 50.0f, "Equipment Scale:", true);
        equipScaleBean.setRangeFixStep(2.0f, 50.0f, 1.0f);
        displayFields = new StsFieldBean[]
        {
            new StsBooleanFieldBean(this, "displayNames", "Names"),
            new StsBooleanFieldBean(this, "realtimeTrack", "Realtime Tracking"),
            new StsBooleanFieldBean(this, "displayLines", "Show Path"),
            new StsBooleanFieldBean(this, "enableTime", "Enable Time"),
            new StsBooleanFieldBean(this, "enableSound", "Enable Sound"),
            new StsBooleanFieldBean(this, "displayIn2d", "Show Markers in 2D"),
            new StsBooleanFieldBean(this, "displayMarkers", "Show Geologic Markers"),
            new StsBooleanFieldBean(this, "displayPerfMarkers", "Show Perforations"),
            new StsBooleanFieldBean(this, "displayEquipMarkers", "Show Equipment"),
            new StsBooleanFieldBean(this, "displayFmiMarkers", "Show FMI Markers"),
            new StsFloatFieldBean(this,"fmiScale", true, "FMI Scale:"),
            perfScaleBean, equipScaleBean,
            new StsFloatFieldBean(this,"curtainStep", true, "Curtain Step:"),
            new StsStringFieldBean(this, "labelFormatAsString", "Label Format:"),
            new StsComboBoxFieldBean(this, "logTypeDisplay3dLeft", "Left Log Type:", "logTypeDisplayList"),
            new StsComboBoxFieldBean(this, "logTypeDisplay3dRight", "Right Log Type:", "logTypeDisplayList"),
            new StsComboBoxFieldBean(this, "displayTypeString", "Select log type display", displayTypeStrings),
            new StsIntFieldBean(this, "logCurveDisplayWidth", true, "Log display width, pixels"),
            new StsIntFieldBean(this, "LogCurveLineWidth", 1, 4, "Log line width, pixels", true),
            new StsIntFieldBean(this, "logLineDisplayWidth", true, "Width of colored log line, pixels")
        };
    }

	public void displayClass(StsGLPanel3d glPanel3d)
	{
        super.displayClass(glPanel3d);
    }

    public void displayTimeClass(StsGLPanel3d glPanel3d, long time)
    {
        StsObject[] wells = getVisibleObjectList();
        for(int i=0; i<wells.length; i++)
        {
            StsLiveWell well = (StsLiveWell)wells[i];
            if(well.timeVector == null)    // Nothing in well yet.
                continue;
            if(!well.setTimeIndex(time))
                continue;
            if((enableTime && well.isAlive(time)) || (!enableTime))
                well.display(glPanel3d, displayNames, displayMarkers, displayPerfMarkers, displayFmiMarkers, this);
        }
    }

    public void drawOnCursor2d(StsGLPanel3d glPanel3d, int dirNo, float dirCoordinate, boolean axesFlipped, boolean xAxisReversed, boolean yAxisReversed)
    {
        Iterator iter = getVisibleObjectIterator();
        long time = currentModel.getProject().getProjectTime();
        while(iter.hasNext())
        {
            StsLiveWell well = (StsLiveWell)iter.next();
            well.setTimeIndex(time);
            if(displayLines)
            {
                if(!currentModel.getProject().isDepth())
                {
                    if(well.isInTime())
                        well.display2d(glPanel3d, displayNames, dirNo, dirCoordinate, axesFlipped,
                                xAxisReversed, yAxisReversed, displayIn2d,
                                displayMarkers, displayPerfMarkers, displayFmiMarkers);
                }
                else
                    well.display2d(glPanel3d, displayNames, dirNo, dirCoordinate, axesFlipped,
                            xAxisReversed, yAxisReversed, displayIn2d,
                            displayMarkers, displayPerfMarkers, displayFmiMarkers);
            }
        }
    }

    public boolean viewObjectChanged(Object source, Object object)
    {
        for(int n = 0; n < getSize(); n++)
        {
            StsWellFrameViewModel wellViewModel = ((StsWell)getElement(n)).wellViewModel;
            if(wellViewModel == null) continue;
            if(wellViewModel.isVisible)
                wellViewModel.viewObjectChanged(source, object);
        }
        return false;
    }

    public boolean getRealtimeTrack() { return realtimeTrack; }
    public void setRealtimeTrack(boolean track)
    {
        if(realtimeTrack == track) return;
        realtimeTrack = track;
    }
}