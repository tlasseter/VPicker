package com.Sts.Actions.Wizards.SubVolume;

import com.Sts.Actions.Wizards.*;
import com.Sts.DBTypes.*;
import com.Sts.MVC.*;
import com.Sts.UI.Beans.*;
import com.Sts.Utilities.*;

import javax.swing.*;

/**
 * <p>Title: </p>
 * <p>Description: </p>
 * <p>Copyright: Copyright (c) 2002</p>
 * <p>Company: </p>
 * @author unascribed
 * @version 1.0
 */
public class StsDefineSurfaceSubVolumePanel extends JPanel
{
    private StsWizard wizard;
    private StsDefineSurfaceSubVolume defineSurfaceSubVolume;

    private StsModelSurface[] availableSurfaces = null;
    private StsModel model = null;

    StsFieldBeanPanel beanPanel = new StsFieldBeanPanel();
    StsComboBoxFieldBean surfaceListBean = new StsComboBoxFieldBean();
    StsComboBoxFieldBean botSurfaceListBean = new StsComboBoxFieldBean();
    StsFloatFieldBean topOffsetBean = new StsFloatFieldBean();
    StsFloatFieldBean botOffsetBean = new StsFloatFieldBean();
    StsComboBoxFieldBean offsetDomainBean = new StsComboBoxFieldBean();

    public StsDefineSurfaceSubVolumePanel(StsWizard wizard, StsDefineSurfaceSubVolume defineSurfaceSubVolume)
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
        byte offsetDomain = model.getProject().getZDomain();
        offsetDomainBean.setSelectedIndex(offsetDomain);
        if(model.getProject().getZDomainSupported() != StsParameters.TD_TIME_DEPTH)
        	offsetDomainBean.setEditable(false);
    }

    void jbInit() throws Exception
    {
        surfaceListBean.initialize(defineSurfaceSubVolume, "surface", "Surface: ");
        topOffsetBean.initialize(defineSurfaceSubVolume, "topOffset", -10000.0f, 10000.0f, "Top Offset: ");
        botOffsetBean.initialize(defineSurfaceSubVolume, "botOffset", -10000.0f, 10000.0f, "Bottom Offset: ");
        offsetDomainBean.initialize(defineSurfaceSubVolume, "offsetDomainString", "Offset domain:", StsParameters.TD_STRINGS);
        offsetDomainBean.setToolTipText("Select the domain of the offsets if used.");
        beanPanel.add(surfaceListBean);
        beanPanel.add(topOffsetBean);
        beanPanel.add(botOffsetBean);
         beanPanel.add(offsetDomainBean);

        beanPanel.setBorder(BorderFactory.createEtchedBorder());
        this.add(beanPanel);
    }
}
