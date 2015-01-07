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

public class StsDefineDualSurfaceSubVolumePanel extends JPanel
{
    private StsWizard wizard;
    private StsDefineDualSurfaceSubVolume defineDualSurfaceSubVolume;

    private StsModelSurface[] availableSurfaces = null;
    private StsModel model = null;

    StsFieldBeanPanel beanPanel = new StsFieldBeanPanel();
    StsComboBoxFieldBean topSurfaceListBean = new StsComboBoxFieldBean();
    StsComboBoxFieldBean botSurfaceListBean = new StsComboBoxFieldBean();
    StsFloatFieldBean topOffsetBean = new StsFloatFieldBean();
    StsFloatFieldBean botOffsetBean = new StsFloatFieldBean();
    StsComboBoxFieldBean offsetDomainBean = new StsComboBoxFieldBean();

    public StsDefineDualSurfaceSubVolumePanel(StsWizard wizard, StsDefineDualSurfaceSubVolume defineDualSurfaceSubVolume)
    {
        this.wizard = wizard;
        this.defineDualSurfaceSubVolume = defineDualSurfaceSubVolume;
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
        topSurfaceListBean.setListItems(availableSurfaces);
        botSurfaceListBean.setListItems(availableSurfaces);
        topSurfaceListBean.setSelectedIndex(0);
        botSurfaceListBean.setSelectedIndex(1);
        topOffsetBean.setValue(0.0f);
        botOffsetBean.setValue(0.0f);
        byte offsetDomain = model.getProject().getZDomain();
        offsetDomainBean.setSelectedIndex(offsetDomain);
        if(model.getProject().getZDomainSupported() != StsParameters.TD_TIME_DEPTH)
        	offsetDomainBean.setEditable(false);
    }

    void jbInit() throws Exception
    {
        topSurfaceListBean.initialize(defineDualSurfaceSubVolume, "topSurface", "Top Surface: ");
        botSurfaceListBean.initialize(defineDualSurfaceSubVolume, "botSurface", "Bottom Surface: ");
        topOffsetBean.initialize(defineDualSurfaceSubVolume, "topOffset", -10000.0f, 10000.0f, "Top Offset: ");
        botOffsetBean.initialize(defineDualSurfaceSubVolume, "botOffset", -10000.0f, 10000.0f, "Bottom Offset: ");
        offsetDomainBean.initialize(defineDualSurfaceSubVolume, "offsetDomainString", "Offset domain:", StsParameters.TD_STRINGS);
        offsetDomainBean.setToolTipText("Select the domain of the offsets if used.");
        beanPanel.add(topSurfaceListBean);
        beanPanel.add(botSurfaceListBean);
        beanPanel.add(topOffsetBean);
        beanPanel.add(botOffsetBean);
        beanPanel.add(offsetDomainBean);
        beanPanel.setBorder(BorderFactory.createEtchedBorder());
        this.add(beanPanel);
    }

    public void checkSetTopSurfaceSelected()
    {
        checkSetSurfaceSelected(topSurfaceListBean, botSurfaceListBean);
    }

    public void checkSetBotSurfaceSelected()
    {
        checkSetSurfaceSelected(botSurfaceListBean, topSurfaceListBean);
    }

    private void checkSetSurfaceSelected(StsComboBoxFieldBean selectedListBean, StsComboBoxFieldBean otherListBean)
    {
        int selectedIndex = selectedListBean.getSelectedIndex();
        int otherIndex = otherListBean.getSelectedIndex();

        if(selectedIndex != otherIndex) return;

        for(int n = 0; n < availableSurfaces.length; n++)
        if(n != selectedIndex) otherListBean.setSelectedIndex(n);
    }
}
