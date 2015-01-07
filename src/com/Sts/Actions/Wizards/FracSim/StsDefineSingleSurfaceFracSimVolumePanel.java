package com.Sts.Actions.Wizards.FracSim;

import com.Sts.Actions.Wizards.*;
import com.Sts.DBTypes.*;
import com.Sts.MVC.*;
import com.Sts.UI.Beans.*;

import javax.swing.*;

/**
 * <p>Title: </p>
 * <p>Description: </p>
 * <p>Copyright: Copyright (c) 2002</p>
 * <p>Company: </p>
 * @author unascribed
 * @version 1.0
 */
public class StsDefineSingleSurfaceFracSimVolumePanel extends JPanel
{
    private StsWizard wizard;
    private StsDefineSingleSurfaceFracSimVolume defineSurfaceSubVolume;

    private StsModelSurface[] availableSurfaces = null;
    private StsModel model = null;

    StsFieldBeanPanel beanPanel = new StsFieldBeanPanel();
    StsStringFieldBean nameBean = new StsStringFieldBean();
    StsComboBoxFieldBean surfaceListBean = new StsComboBoxFieldBean();
    StsComboBoxFieldBean botSurfaceListBean = new StsComboBoxFieldBean();
    StsFloatFieldBean topOffsetBean = new StsFloatFieldBean();
    StsFloatFieldBean botOffsetBean = new StsFloatFieldBean();

    public StsDefineSingleSurfaceFracSimVolumePanel(StsWizard wizard, StsDefineSingleSurfaceFracSimVolume defineSurfaceSubVolume)
    {
        this.wizard = wizard;
        this.defineSurfaceSubVolume = defineSurfaceSubVolume;
        try
        {
            jbInit();
        }
        catch(Exception ex)
        {
            ex.printStackTrace();
        }
    }

    public void initialize()
    {
        model = wizard.getModel();
        availableSurfaces = (StsModelSurface[])model.getCastObjectList(StsModelSurface.class);
		if(availableSurfaces != null) defineSurfaceSubVolume.setSurface(availableSurfaces[0]);
        surfaceListBean.setListItems(availableSurfaces);
        topOffsetBean.setValue(0.0f);
        botOffsetBean.setValue(0.0f);
    }

    void jbInit() throws Exception
    {
        nameBean.initialize(defineSurfaceSubVolume, "subVolumeName", true, "Name:");
        surfaceListBean.initialize(defineSurfaceSubVolume, "surface", "Surface: ");
        topOffsetBean.initialize(defineSurfaceSubVolume, "topOffset", -10000.0f, 10000.0f, "Top Offset: ");
        botOffsetBean.initialize(defineSurfaceSubVolume, "botOffset", -10000.0f, 10000.0f, "Bottom Offset: ");

        beanPanel.add(nameBean);
        beanPanel.add(surfaceListBean);
        beanPanel.add(topOffsetBean);
        beanPanel.add(botOffsetBean);

        beanPanel.setBorder(BorderFactory.createEtchedBorder());
        this.add(beanPanel);
    }
}