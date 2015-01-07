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

public class StsDefineDualSurfaceFracSimVolumePanel extends JPanel
{
    private StsWizard wizard;
    private StsDefineDualSurfaceFracSimVolume defineDualSurfaceFracSimVolume;

    private StsModelSurface[] availableSurfaces = null;
    private StsModel model = null;
    StsFieldBeanPanel beanPanel = new StsFieldBeanPanel();
    StsStringFieldBean nameBean = new StsStringFieldBean();
    StsComboBoxFieldBean topSurfaceListBean = new StsComboBoxFieldBean();
    StsComboBoxFieldBean botSurfaceListBean = new StsComboBoxFieldBean();
    StsFloatFieldBean topOffsetBean = new StsFloatFieldBean();
    StsFloatFieldBean botOffsetBean = new StsFloatFieldBean();

    public StsDefineDualSurfaceFracSimVolumePanel(StsWizard wizard, StsDefineDualSurfaceFracSimVolume defineDualSurfaceSubVolume)
    {
        this.wizard = wizard;
        this.defineDualSurfaceFracSimVolume = defineDualSurfaceSubVolume;
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
    }

    void jbInit() throws Exception
    {
//        this.setLayout(new GridBagLayout());
        nameBean.initialize(defineDualSurfaceFracSimVolume, "subVolumeName", true, "Name:");
        topSurfaceListBean.initialize(defineDualSurfaceFracSimVolume, "topSurface", "Top Surface: ");
        botSurfaceListBean.initialize(defineDualSurfaceFracSimVolume, "botSurface", "Bottom Surface: ");
        topOffsetBean.initialize(defineDualSurfaceFracSimVolume, "topOffset", -10000.0f, 10000.0f, "Top Offset: ");
        botOffsetBean.initialize(defineDualSurfaceFracSimVolume, "botOffset", -10000.0f, 10000.0f, "Bottom Offset: ");
/*
        beanPanel.add(nameBean, new GridBagConstraints(0, 0, 1, 1, 1.0, 0.0
            , GridBagConstraints.WEST, GridBagConstraints.HORIZONTAL, new Insets(10, 5, 0, 0), 0, 0));
        beanPanel.add(topSurfaceListBean, new GridBagConstraints(0, 1, 1, 1, 1.0, 0.0
            , GridBagConstraints.WEST, GridBagConstraints.HORIZONTAL, new Insets(5, 5, 0, 0), 0, 0));
        beanPanel.add(botSurfaceListBean, new GridBagConstraints(0, 2, 1, 1, 1.0, 0.0
            , GridBagConstraints.WEST, GridBagConstraints.HORIZONTAL, new Insets(5, 5, 0, 0), 0, 0));
        beanPanel.add(topOffsetBean, new GridBagConstraints(0, 3, 1, 1, 1.0, 0.0
            , GridBagConstraints.WEST, GridBagConstraints.HORIZONTAL, new Insets(5, 5, 0, 0), 0, 0));
        beanPanel.add(botOffsetBean, new GridBagConstraints(0, 4, 1, 1, 1.0, 0.0
            , GridBagConstraints.WEST, GridBagConstraints.HORIZONTAL, new Insets(5, 5, 0, 0), 0, 0));
        */
        beanPanel.add(nameBean);
        beanPanel.add(topSurfaceListBean);
        beanPanel.add(botSurfaceListBean);
        beanPanel.add(topOffsetBean);
        beanPanel.add(botOffsetBean);

        beanPanel.setBorder(BorderFactory.createEtchedBorder());
        this.add(beanPanel);
//        this.add(beanPanel, new GridBagConstraints(0, 0, 1, 1, 1.0, 0.0
//            , GridBagConstraints.WEST, GridBagConstraints.HORIZONTAL, new Insets(2, 5, 0, 0), 200, 0));
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