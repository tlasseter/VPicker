package com.Sts.UI;

import com.Sts.DB.*;
import com.Sts.DBTypes.*;
import com.Sts.UI.Beans.*;
import com.Sts.Utilities.*;

/**
 * <p>Title: S2S development</p>
 *
 * <p>Description: Integrated seismic to simulation software</p>
 *
 * <p>Copyright: Copyright (c) 2004</p>
 *
 * <p>Company: S2S Systems LLC</p>
 *
 * @author not attributable
 * @version c51c
 */

public class StsSuperGatherProperties extends StsPanelProperties implements StsSerializable
{
    boolean ignoreSuperGather = false;
    int nSuperGatherCols = 3; //default should use supergather - SWC 6/19/09
    int nSuperGatherRows = 1;
    byte gatherType = StsPreStackLineSetClass.SUPER_INLINE; //default should use supergather - SWC 6/19/09
    boolean is3d;
    transient StsComboBoxFieldBean gatherTypeBean;
    transient StsIntFieldBean nRowsBean;
    transient StsIntFieldBean nColsBean = null;
    transient StsBooleanFieldBean ignoreSuperGatherBean;
    transient StsSuperGatherConfigPanel gatherGraphicPanel = new StsSuperGatherConfigPanel();

	static private final String title = "Super Gather Properties";
	public StsSuperGatherProperties()
	{
	}

	public StsSuperGatherProperties(StsPreStackLineSetClass prestackClass, String fieldName)
	{
		super(title, fieldName);
        is3d = prestackClass instanceof StsPreStackLineSet3dClass;
    }

	public StsSuperGatherProperties(StsObject parentObject, StsSuperGatherProperties defaultProperties, String fieldName)
	{
        super(parentObject, title, fieldName);
        initializeDefaultProperties(defaultProperties);
    }
    public void initializeDefaultProperties(Object defaultProperties)
    {
        if(defaultProperties != null)
            StsToolkit.copyObjectNonTransientFields(defaultProperties, this);
    }

	public void initializeBeans()
	{
       gatherTypeBean = new StsComboBoxFieldBean(this, "gatherTypeString", "Gather Type:", StsPreStackLineSet3dClass.GATHER_TYPE_STRINGS);
       ignoreSuperGatherBean = new StsBooleanFieldBean(this, "ignoreSuperGather", "Ignore Super Gather:");
       if(is3d)
       {
            nRowsBean = new StsIntFieldBean(this,"nSuperGatherRows", 1, 999, "Number of inlines:",true);
            nColsBean = new StsIntFieldBean(this,"nSuperGatherCols", 1, 999, "Number of crosslines:",true);
            nRowsBean.setStep(2);
            nColsBean.setStep(2);
            propertyBeans = new StsFieldBean[] { ignoreSuperGatherBean, gatherTypeBean, nRowsBean, nColsBean };
       }
       else
       {
            nColsBean = new StsIntFieldBean(this,"nSuperGatherCols", 1, 999, "Number of crosslines:",true);
            nColsBean.setStep(2);
            propertyBeans = new StsFieldBean[] { ignoreSuperGatherBean, gatherTypeBean, nColsBean };
       }
       configureGatherTypeStatus();
       configureIgnoreStatus();
	}

    public StsSuperGatherConfigPanel getGraphicsPanel()
    {
        return gatherGraphicPanel;
    }

    public byte getGatherType() { return gatherType; }   
    public void setGatherType(byte type)
    { 
        this.gatherType = type;
    }

    public String getGatherTypeString()
    {
        return StsPreStackLineSetClass.GATHER_TYPE_STRINGS[gatherType];
    }

    public void setGatherTypeString(String option)
    {
        for (byte i = 0; i < StsPreStackLineSetClass.GATHER_TYPE_STRINGS.length; i++)
        {
            if (StsPreStackLineSetClass.GATHER_TYPE_STRINGS[i] == option)
            {
                if(gatherType == i) return;
                gatherType = i;
                configureGatherTypeStatus();
//                changed = true;
                return;
            }
        }
    }

    public void configureGatherTypeStatus()
    {
        if(nRowsBean == null) return;
        switch(gatherType)
        {
            case StsPreStackLineSetClass.SUPER_SINGLE:
                if(is3d)
                {
                    nRowsBean.setEditable(false);
                    nRowsBean.setValue(1);
                    setNSuperGatherRows(1);
                }
                nColsBean.setEditable(false);
                nColsBean.setValue(1);
                setNSuperGatherCols(1);
                break;
             case StsPreStackLineSetClass.SUPER_CROSS:
             case StsPreStackLineSetClass.SUPER_RECT:
             if(is3d)
                 nRowsBean.setEditable(true);
             nColsBean.setEditable(true);
             break;
             case StsPreStackLineSetClass.SUPER_INLINE:
                if(is3d) nRowsBean.setEditable(false);
                break;
             case StsPreStackLineSetClass.SUPER_XLINE:
                 if(is3d) nRowsBean.setEditable(true);
                 nColsBean.setEditable(false);
                 break;
        }
    }

    public int getNSuperGatherCols()
    {
        return nSuperGatherCols;
    }
    public void setNSuperGatherCols(int val)
    {
        this.nSuperGatherCols = val;
    }

    public int getNSuperGatherRows()
    {
        return nSuperGatherRows;
    }
    public void setNSuperGatherRows(int val)
    {
        this.nSuperGatherRows = val;
    }

    public void setIgnoreSuperGather(boolean use)
    {
        ignoreSuperGather = use;
        configureIgnoreStatus();
    }

    public void configureIgnoreStatus()
    {
        if(nRowsBean == null) return;
        if(ignoreSuperGather)
        {
            ignoreSuperGatherBean.setSelected(true);
            nRowsBean.setEditable(false);
            if(nColsBean != null)
                nColsBean.setEditable(false);
            gatherTypeBean.setEditable(false);
        }
        else
        {
            ignoreSuperGatherBean.setSelected(false);
            nRowsBean.setEditable(true);
            if(nColsBean != null)
                nColsBean.setEditable(true);
            gatherTypeBean.setEditable(true);
        }
    }

    public boolean getIgnoreSuperGather() { return ignoreSuperGather; }

}
