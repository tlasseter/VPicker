package com.Sts.Actions.Wizards.SuperGather;

import com.Sts.Actions.Wizards.*;
import com.Sts.DBTypes.*;
import com.Sts.UI.Beans.*;
import com.Sts.UI.*;

import java.awt.*;

/**
 * <p>Title: </p>
 * <p>Description: </p>
 * <p>Copyright: Copyright (c) 2002</p>
 * <p>Company: </p>
 * @author unascribed
 * @version 1.0
 */

public class StsStackParameterPanel extends StsFieldBeanPanel
{
    private StsSuperGatherWizard wizard;
    private StsStackParameter wizardStep;

	StsPreStackLineSet volume;
    StsPreStackLineSetClass volumeClass;

	StsPreStackVelocityModel velocityModel;

    boolean changed = false;

    private StsGroupBox defineGroupBox = new StsGroupBox("Traces and Samples to Stack");
    private StsFloatFieldBean neighborRadiusBean;
    private StsFloatFieldBean maxZBean;

    StsGroupBox gatherDefineBox = new StsGroupBox("Gather Configuration");
    StsGroupBox gatherGraphicsBox = new StsGroupBox();

    StsComboBoxFieldBean gatherTypeBean;
    StsIntFieldBean xExtentBean;
	StsIntFieldBean yExtentBean;
    StsSuperGatherConfigPanel gatherGraphicPanel = new StsSuperGatherConfigPanel();

    public StsStackParameterPanel(StsWizard wizard, StsWizardStep wizardStep)
    {
        this.wizard = (StsSuperGatherWizard)wizard;
        this.wizardStep = (StsStackParameter)wizardStep;
        try
        {
            buildBeans();
            jbInit();
        }
        catch(Exception ex)
        {
            ex.printStackTrace();
        }
    }

    private void buildBeans()
    {
        neighborRadiusBean = new StsFloatFieldBean(StsPreStackLineSetClass.class, "stackNeighborRadius", true, "Stack neighbor radius:",true);
        maxZBean = new StsFloatFieldBean(StsPreStackLineSetClass.class, "maxZToStack", true, "Maximum Z:",true);

        gatherTypeBean = new StsComboBoxFieldBean(this, "gatherTypeString", "Gather Type:", StsPreStackLineSetClass.GATHER_TYPE_STRINGS);
        yExtentBean = new StsIntFieldBean(this, "yGatherExtent", 1, 9, "Number CrossLine Gathers:", true);
        xExtentBean = new StsIntFieldBean(this, "xGatherExtent", 1, 9, "Number Inline Gathers:", true);
    }

    public void initialize()
    {
        volume = this.wizard.getPreStackVolume();
        volumeClass = volume.lineSetClass;
		velocityModel = this.wizard.getVelocityModel();

		defineGroupBox.setPanelObject(volumeClass);
        maxZBean.setValueAndRange(volumeClass.getMaxZToStack(), volume.getZMin(), volume.getZMax());
        maxZBean.fixStep(volume.getZInc());

        gatherGraphicPanel.setPreStackSeismic(volume);

        gatherTypeBean.setSelectedItem(volume.superGatherProperties.getGatherTypeString());
        xExtentBean.setValue(volume.superGatherProperties.getNSuperGatherCols());
        xExtentBean.setStep(2);
        yExtentBean.setValue(volume.superGatherProperties.getNSuperGatherRows());
        yExtentBean.setStep(2);

        if(volume.superGatherProperties.getGatherType() == volumeClass.SUPER_SINGLE)
        {
            xExtentBean.setEditable(false);
            yExtentBean.setEditable(false);
        }
        if(volume instanceof StsPreStackLineSet2d)
        {
            yExtentBean.setEditable(false);
        }
    }

    void jbInit() throws Exception
    {
        defineGroupBox.gbc.fill = gbc.HORIZONTAL;
        gbc.anchor = gbc.NORTH;
        defineGroupBox.addEndRow(neighborRadiusBean);
        defineGroupBox.addEndRow(maxZBean);

        gatherDefineBox.gbc.fill = gbc.HORIZONTAL;
        gatherDefineBox.gbc.anchor = gbc.NORTH;
        gatherDefineBox.gbc.weighty = 0.0;
        gatherDefineBox.addEndRow(gatherTypeBean);
        gatherDefineBox.addEndRow(xExtentBean);
        gatherDefineBox.addEndRow(yExtentBean);

        gatherGraphicPanel.setBackground(Color.WHITE);
        gatherGraphicPanel.setSize(300,300);
        gatherGraphicsBox.gbc.fill = gbc.BOTH;
        gatherGraphicsBox.addEndRow(gatherGraphicPanel);
        gatherDefineBox.gbc.fill = gbc.BOTH;
        gatherDefineBox.gbc.gridwidth = 2;
        gatherDefineBox.gbc.weighty = 1.0;
        gatherDefineBox.add(gatherGraphicsBox);

        gbc.weighty = 0.0;
        gbc.fill = gbc.HORIZONTAL;
        add(defineGroupBox);

        gbc.fill = gbc.BOTH;
        gbc.weighty = 1.0;
        add(gatherDefineBox);
    }

    public String getGatherTypeString()
    {
        if(volume == null)
            return StsPreStackLineSet3dClass.SINGLE_GATHER_STRING;
        return volume.superGatherProperties.getGatherTypeString();
    }
    public void setGatherTypeString(String option)
    {
        if(volume == null) return;
        volume.superGatherProperties.setGatherTypeString(option);
        changed = true;
        if(volume.superGatherProperties.getGatherType() == volumeClass.SUPER_SINGLE)
        {
            xExtentBean.setEditable(false);
            yExtentBean.setEditable(false);
            xExtentBean.setValue(1);
            volume.superGatherProperties.setNSuperGatherCols(1);
            yExtentBean.setValue(1);
            volume.superGatherProperties.setNSuperGatherRows(1);
        }
        else
        {
           xExtentBean.setEditable(true);
           yExtentBean.setEditable(true);
       }
       if(volume instanceof StsPreStackLineSet2d)
            yExtentBean.setEditable(false);
    }
    public int getXGatherExtent()
    {
        if(volume == null) return 1;
        return volume.superGatherProperties.getNSuperGatherCols();
    }
    public void setXGatherExtent(int val)
    {
        if(volume == null) return;
        volume.superGatherProperties.setNSuperGatherCols(val);
        changed = true;
        gatherGraphicPanel.repaint();
    }

    public int getYGatherExtent()
    {
        if(volume == null) return 1;
        return volume.superGatherProperties.getNSuperGatherRows();
    }

    public void setYGatherExtent(int val)
    {
        if(volume == null) return;
        volume.superGatherProperties.setNSuperGatherRows(val);
        changed = true;
        gatherGraphicPanel.repaint();
	}

    public void checkSaveSuperGatherProperties()
    {
        if(changed) volume.superGatherProperties.commitChanges();
    }
}
